
# old
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


    # tested on SRR14633482 - did not get a Solo.out directory? Ran with 10xv3 params (though umi+cb=30)
    def execute(self):
        read_bio, read_tech, tech = self._impute_10x_version()
        star_param = self._get_10x_STAR_parameters(tech)  # as dictionary

        if tech != 'unknown':
            cmd = ['STAR',
                   '--runMode', 'alignReads',
                   '--runThreadN', f'{self.num_streams}',
                   '--genomeDir', f'{self.resourcedir}/genomes/{self.species}',
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
            # if str(cp.returncode) == "0":
            #     self.outlist.append(self.srrid)

        else:
            pass


# currently, if solo.out dir exsts, save results to the temp directory
# then merge results with previous star runs (to do)
# input should be a project accession
# no real way to verify that the reads are indeed smartseq...
# can pass star parameters from config


# XXX moved to utils
# def get_default_config():
#     cp = ConfigParser()
#     cp.read(os.path.expanduser("~/git/scqc/etc/scqc.conf"))
#     return cp


# # john lee is satisfied with this class 6/3/2021
# def get_configstr(cp):
#     with io.StringIO() as ss:
#         cp.write(ss)
#         ss.seek(0)  # rewind
#         return ss.read()


# old
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
               '--genomeDir', f'{self.resourcedir}/genomes/{self.species}',
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




        # TODO filter to the runs that we havne't aligned yet. 
        #   adjust manifest accordingly.

        # runs_done_file = f'{self.outputdir}/{self.srpid}_smartseq_Solo.out/Gene/raw/barcodes.tsv'
        # if os.path.isfile(runs_done_file):
        #     runs_done = open(runs_done_file).read().strip().split('\n')

        #     # of the runs that i find in the manifest, which have already been aligned?
        #     new_runs = listdiff(manifest.run.values, runs_done)
        #     #list(set(manifest.run.values) - set(runs_done)).sort()
        #     tmp_mani = manifest.loc[new_runs == manifest.run, :]

        #     # save tmp manifest to temp directory - may be empty
        #     tmp_manipath = manipath.replace(
        #         f'{self.metadir}', f'{self.tempdir}')
        #     tmp_mani.to_csv(tmp_manipath, sep="\t")
        #     out_file_prefix = f'{self.tempdir}/{self.srpid}_smartseq_'
        # else:
        #     tmp_manipath = manipath
        #     out_file_prefix = f'{self.staroutdir}/{self.srpid}_smartseq_'
        