#!/usr/bin/env python
#
#  Module to deal with running STAR/StarSolo
#

import gzip
import logging
import os
import requests
import subprocess

from scqc.utils  import gzip_decompress, download_ftpurl

# inputs should be runs  identified as 'some10x'
# fastq files should already downloaded.
# srrid and species needs to be passed in from the dataframe
# can pass star parameters from config
class Align10xSTAR(object):
    '''
        - Identifies 10x version
        - gets the path to fastq files
        - 
        Simple wrapper for STAR - 10x input
    '''

    def __init__(self, config, srrid, species, outlist):
        self.log = logging.getLogger('star')
        self.config = config

        self.tempdir = os.path.expanduser(
            self.config.get('analysis', 'tempdir'))
        self.srrid = srrid
        self.log.debug(f'aligning id {srrid}')
        self.staroutdir = os.path.expanduser(
            self.config.get('analysis', 'staroutdir'))
        self.resourcedir = os.path.expanduser(
            self.config.get('analysis', 'resourcedir'))
        self.species = species
        self.outlist = outlist
        self.num_streams = self.config.get('analysis', 'num_streams')




    def _get_10x_STAR_parameters(self, tech):
        d = {
            "10xv1":
                {
                    "solo_type": "CB_UMI_Simple",
                    "white_list_path": f'{self.resourcedir}/whitelists/whitelist_10xv1.txt',
                    "CB_length": "14",
                    "UMI_start": "15",
                    "UMI_length": "10"
                },
            "10xv2":
                {
                    "solo_type": "CB_UMI_Simple",
                    "white_list_path": f'{self.resourcedir}/whitelists/whitelist_10xv2.txt',
                    "CB_length": "16",
                    "UMI_start": "17",
                    "UMI_length": "10"
                },
            "10xv3":
                {
                    "solo_type": "CB_UMI_Simple",
                    "white_list_path": f'{self.resourcedir}/whitelists/whitelist_10xv3.txt',
                    "CB_length": "16",
                    "UMI_start": "17",
                    "UMI_length": "12"
                }
        }
        return(d[tech])

    # tested on SRR14633482 - did not get a Solo.out directory? Ran with 10xv3 params (though umi+cb=30)
    def execute(self):
        read_bio, read_tech, tech = self._impute_10x_version()
        star_param = self._get_10x_STAR_parameters(tech)  # as dictionary

        if tech != 'unknown':
            cmd = ['STAR',
                   '--runMode', 'alignReads',
                   '--runThreadN', f'{self.num_streams}',
                   '--genomeDir', f'{self.resourcedir}/genomes/{self.species}/STAR_index',
                   '--outFileNamePrefix', f'{self.staroutdir}/{self.srrid}_{tech}_',
                   '--soloType', star_param["solo_type"],
                   '--soloCBwhitelist', star_param["white_list_path"],
                   '--soloCBlen', star_param["CB_length"],
                   '--soloUMIstart', f'{int(star_param["CB_length"]) + 1}',
                   '--soloUMIlen', star_param["UMI_length"],
                   '--soloFeatures', 'Gene',
                   '--readFilesIn', read_bio, read_tech,
                   '--outSAMtype', 'None']

            cp = subprocess.run(cmd)

            cmdstr = " ".join(cmd)
            logging.debug(f"Fasterq-dump command: {cmdstr} running...")
            cp = subprocess.run(cmd)
            logging.debug(
                f"Ran cmd='{cmdstr}' returncode={cp.returncode} {type(cp.returncode)} ")
            # successful runs - append to outlist.
            if str(cp.returncode) == "0":
                self.outlist.append(self.srrid)

        else:
            pass


# currently, if solo.out dir exsts, save results to the temp directory
# then merge results with previous star runs (to do)
# input should be a project accession
# no real way to verify that the reads are indeed smartseq...
# can pass star parameters from config
class AlignSmartSeqSTAR(object):

    def __init__(self, config, species, srpid, outlist):
        self.log = logging.getLogger('sra')
        self.config = config

        self.tempdir = os.path.expanduser(
            self.config.get('analysis', 'tempdir'))
        self.srpid = srpid
        self.log.debug(f'aligning smartseq run from {srpid}')
        self.staroutdir = os.path.expanduser(
            self.config.get('analysis', 'staroutdir'))
        self.resourcedir = os.path.expanduser(
            self.config.get('analysis', 'resourcedir'))
        self.species = species
        self.outlist = outlist
        self.num_streams = self.config.get('analysis', 'num_streams')

    def _make_manifest(self):
        # search for all fastq files with <run>_[0-9].fastq
        manipath = f"{self.metadir}/{self.srpid}_smartseq_manifest.tsv"

        # metadata here should only contain smart seq data.
        df = pd.read_csv(
            f'{self.metadir}/{self.srpid}_smartseq_metadata.tsv', sep="\t")

        df.runs = df.runs.apply(ast.literal_eval)
        df = df.explode('runs')
        # runlist = df.runs.values

        # runid = runlist[1]
        allRows = []
        for runid in df.runs:
            fqs = glob.glob(f'{self.cachedir}/{runid}*.fastq').sort()

            if len(fqs) > 0 and len(fqs) < 3:  # fastq files found
                if len(fqs) == 1:
                    fqs.append(['-', runid])
                elif len(fqs) == 2:
                    fqs.append(runid)

                allRows.append(fqs)

        manifest = pd.DataFrame(allRows, columns=['read1', 'read2', 'run'])

        # overwrite
        manifest.to_csv(manipath, sep="\t", header=None, index=False, mode="w")

        return(manipath, manifest)

    # to do;
    def _merge_solo_out_results(self, final_path, temp_path):
        '''
        Merge the Solo.out results for two star runs on different parts of the data
        '''
        # merge data,

        # delete from temp

        pass

    def execute(self):

        ss_params = {"solo_type": "SmartSeq",
                     "soloUMIdedup": "Exact",
                     "soloStrand": "Unstranded"}

        manipath, manifest = self._make_manifest()

        # do we already have star output for this project.
        # If so, save starout results in a temp directory, merge data/stats
        # then delete tmp
        runs_done_file = f'{self.staroutdir}/{self.srpid}_smartseq_Solo.out/Gene/raw/barcodes.tsv'
        if os.path.isfile(runs_done_file):
            runs_done = open(runs_done_file).read().strip().split('\n')

            # of the runs that i find in the manifest, which have already been aligned?
            new_runs = listdiff(manifest.run.values, runs_done)
            #list(set(manifest.run.values) - set(runs_done)).sort()
            tmp_mani = manifest.loc[new_runs == manifest.run, :]

            # save tmp manifest to temp directory - may be empty
            tmp_manipath = manipath.replace(
                f'{self.metadir}', f'{self.tempdir}')
            tmp_mani.to_csv(tmp_manipath, sep="\t")
            out_file_prefix = f'{self.tempdir}/{self.srpid}_smartseq_'
        else:
            tmp_manipath = manipath
            out_file_prefix = f'{self.staroutdir}/{self.srpid}_smartseq_'

        cmd = ['STAR',
               '--runMode', 'alignReads',
               '--runThreadN', f'{self.num_streams}',
               '--genomeDir', f'{self.resourcedir}/genomes/{self.species}/STAR_index',
               '--outFileNamePrefix', out_file_prefix,
               '--soloType', ss_params["solo_type"],
               '--soloFeatures', 'Gene',
               '--readFilesManifest', f'{tmp_manipath}',
               '--soloUMIdedup', ss_params["soloUMIdedup"],
               '--soloStrand', ss_params["soloStrand"],
               '--outSAMtype', 'None']

        cmdstr = " ".join(cmd)

        logging.debug(f"STAR command: {cmdstr} running...")
        cp = subprocess.run(cmd)

        logging.debug(
            f"Ran cmd='{cmdstr}' returncode={cp.returncode} {type(cp.returncode)} ")
        # successful runs - append to outlist.
        if str(cp.returncode) == "0":
            self.outlist.append(self.srrid)

        # did we write to a temp directory?
        if out_file_prefix.startswith(f'{self.tempdir}'):
            self._merge_solo_out_results(
                f'{self.staroutdir}/{self.srpid}_smartseq_',    # starout direc
                out_file_prefix)                                # temp direc



def setup(config, force=False):
    '''
    Builds directories in config file 
    Download the appropriate supplement data.
    Only needs to be done once.
    '''
    log = logging.getLogger('star')
    metadir = os.path.expanduser(config.get('star', 'metadir'))
    cachedir = os.path.expanduser(config.get('star', 'cachedir'))
    tempdir = os.path.expanduser(config.get('star', 'tempdir'))
    resourcedir = os.path.expanduser(config.get('star', 'resourcedir'))
    #n_core = config.get('star', 'n_core')
    
    for d in [metadir, cachedir, tempdir, resourcedir]:
        try:
            log.debug(f"making directory: {d}")
            os.makedirs(d)
        except FileExistsError:
            pass
  
    get_whitelists(config, force)
    get_genome_data(config, force)
    build_genome_indices(config, force)


def get_whitelists(config):
    """
    Get cellranger tag whitelists. Assumes resourcedir exists. 
    """
    log = logging.getLogger('star')
    wls = [ '10x_v1_whitelist','10x_v2_whitelist', '10x_v3_whitelist']
    outdir = os.path.expanduser(config.get('star','resourcedir'))
    
    for key in wls:
        url = config.get('star', key)
        log.debug(f'getting cellranger whitelist from {url}...')
        r = requests.get(url)
        if r.status_code == 200:
            with open("/".join([outdir, key]), "w") as f:
                if url.endswith(".gz"):
                    f.write(gzip.decompress(r.content).decode())
                else:
                    f.write(r.text)
        else:
            self.log.warning(f"Retrieving {url} failed with status_code: {str(r.status_code)}")


def get_genome_data(config):
    """
    Download required genome data
    """
    log = logging.getLogger('star')
    resourcedir = os.path.expanduser(config.get('star','resourcedir'))
    speciesnames = config.get('star','species')
    specieslist = [i.strip() for i in speciesnames.split(',')]
    log.debug(f'got species {specieslist}')
    
    for species in specieslist:
        outdir = "/".join([resourcedir, species])
        try:
            os.makedirs(outdir)
        except FileExistsError:
            pass

        fa_url=config.get('star',f'{species}_fa')
        log.debug(f'{species} fa_url={fa_url}')
        download_ftpurl(fa_url, outdir, 'genome.fa')
        
        gtf_url=config.get('star',f'{species}_gtf')
        log.debug(f'{species} gtf_url={gtf_url}')
        download_ftpurl(gtf_url, outdir, 'annotation.gtf')

      
def build_genome_indices(config, force=False):
    log = logging.getLogger('star')
    
    n_core = int(config.get('star','n_core'))
    resourcedir = os.path.expanduser(config.get('star','resourcedir'))
    speciesnames = config.get('star','species')
    specieslist = [i.strip() for i in speciesnames.split(',')]
    log.debug(f'got species {specieslist}')
    
    for species in specieslist:
        outdir = "/".join([resourcedir, species])
        cmd = ["STAR",
           "--runMode", "genomeGenerate",
           "--genomeSAsparseD", "3" ,   # for low memory
           "--genomeSAindexNbases", "12" # for low memory
           "--runThreadN", f'{n_core}',
           "--genomeDir", f'{outdir}',
           "--genomeFastaFiles", f'{outdir}/genome.fa',
           "--sjdbGTFfile", f'{outdir}/annotation.gtf']

        cmdstr = " ".join(cmd)
        log.info(f'running {cmdstr} ...')
        cp = subprocess.run(cmdstr, 
                            shell=True, 
                            universal_newlines=True, 
                            stdout=subprocess.PIPE, 
                            stderr=subprocess.PIPE)
        log.debug(f'Ran cmd={cmdstr}  returncode={cp.returncode} ')

    
    logging.debug(
        f"Ran cmd='{cmdstr}' returncode={cp.returncode} {type(cp.returncode)} ")
    if str(cp.returncode) == "0":
        log.warning(f"non-zero return code from STAR. See Log.out...")



# to be done - get marker sets for mouse brain - see bin/getMarkers.R
def get_meta_marker_sets(config):
    pass








# to do: include drivers for 10x and ss alignments
if __name__ == "__main__":

    gitpath = os.path.expanduser("~/git/scqc")
    sys.path.append(gitpath)

    FORMAT = '%(asctime)s (UTC) [ %(levelname)s ] %(filename)s:%(lineno)d %(name)s.%(funcName)s(): %(message)s'
    logging.basicConfig(format=FORMAT)
    logging.getLogger().setLevel(logging.DEBUG)

    parser = argparse.ArgumentParser()

    parser.add_argument('-d', '--debug',
                        action="store_true",
                        dest='debug',
                        help='debug logging')

    parser.add_argument('-v', '--verbose',
                        action="store_true",
                        dest='verbose',
                        help='verbose logging')

    parser.add_argument('-s', '--setup',
                        action='store_true',
                        dest='setup',
                        help='Set up directories and downloads supplemental data')

    parser.add_argument('-F', '--force',
                        action='store_true',
                        dest='force',
                        default=False, 
                        help='Re-do, overwriting what is done.')

    parser.add_argument('-f', '--fasterq',
                        metavar='fasterq',
                        type=str,
                        nargs='+',
                        required=False,
                        default=None,
                        help='Download args with fasterq-dump. e.g. SRR14584407')


    parser.add_argument('-tx', '--tenx',
                        metavar='tenx_align',
                        type=str,
                        nargs='+',
                        required=False,
                        default=None,
                        help='Align 10x args with STAR. e.g. SRR14584407')

    parser.add_argument('-ss', '--smartseq',
                        metavar='ss_align',
                        type=str,
                        nargs='+',
                        required=False,
                        default=None,
                        help='Align SmartSeq args with STAR. e.g. SRP308826')


    args = parser.parse_args()

    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)

    cp = get_default_config()
    cs = get_configstr(cp)

    logging.debug(f"got config: {cs}")

    if args.setup:
        s = SetUp(cp, force=args.force)
        s.execute()





