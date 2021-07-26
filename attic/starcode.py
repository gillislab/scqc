
if flag < 2:  # at least one file was updated. Let's re-generate the genome
    self.log.info(f"... Generating the genome for {self.species}")

    cmd = ["STAR",
           "--runMode", "genomeGenerate",
           "--runThreadN", f'{self.n_core}',
           "--genomeDir", f'{outdir}',
           "--genomeFastaFiles", f'{fa_path}',
           "--sjdbGTFfile", f'{gtf_path}']

    cmdstr = " ".join(cmd)
    cp = subprocess.run(cmd)
    logging.debug(
        f"Ran cmd='{cmdstr}' returncode={cp.returncode} {type(cp.returncode)} ")
    if str(cp.returncode) == "0":
        pass

else:
    self.log.info(
        "... Genome Index has previously been generated with these GTF and FASTA files. Nothing to do.")


    
    
    #try:
    #    os.makedirs(outdir)
    #except FileExistsError:
    #    pass

    #spec2urls = {
    #    "human": {
    #        "fa": "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/GRCh38.primary_assembly.genome.fa.gz",
    #        "gtf": "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/gencode.v37.annotation.gtf.gz"
    #    },
    #    "mouse": {
    #        "fa": "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M26/GRCm39.primary_assembly.genome.fa.gz",
    #        "gtf": "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M26/gencode.vM26.annotation.gtf.gz"
    #    }
    #}


    #flag = 0
    # if the files are new, download them to the generic file names
    # and create an empty file to identify current gtf/fa version
    #fa_tail = fa_url.split("/")[-1]
    # this is a new fa file, let's overwrite the old genome.fa file
    #if not os.path.exists(f"{outdir}/{fa_tail}"):
    #    self.log.info(f"... Downloading FASTA file for {self.species}")
    #    os.system(f"wget -O {outdir}/genome.fa.gz {fa_url}")
    #    os.system(f"gunzip -f {outdir}/genome.fa.gz")
    #else:  # fa file given is the current version
    #    flag += 1

    #gtf_tail = gtf_url.split("/")[-1]
    #if not os.path.exists(f"{outdir}/{gtf_tail}"):
    #    self.log.info(f"... Downloading GTF file for {self.species}")
    #    os.system(f"wget -O {outdir}/annotation.gtf.gz {gtf_url}")
    #    os.system(f"gunzip -f {outdir}/annotation.gtf.gz")
    #else:  # gtf file given is the current version
    #    flag += 1

    # remove all gz objects in the directory - only keeps the most recent
    #for p in Path(outdir).glob("*.gz"):
    #    p.unlink()

    # touch a gz file to remember gtf/fa version
    #Path(f"{outdir}/{fa_tail}").touch()
    #Path(f"{outdir}/{gtf_tail}").touch()

    #return(f"{outdir}/annotation.gtf", f"{outdir}/genome.fa", outdir, flag)
    
    
    def _clean_up_tempdir(self, proj_id, solooutdir):
        newdirname = f'{self.cachedir}/{proj_id}/{base}'
        os.rename(solooutdir,newdirname)
        self.log.debug(f'moved {solooutdir} -> {newdirname}')
        
        starlog = base.replace('Solo.out','Log.final.out')
        newdirname = f'{self.cachedir}/{proj_id}/{starlog}'
        os.rename(f'{dirname}/{starlog}',newdirname)
        self.log.debug(f'moved {dirname}/{starlog} -> {newdirname}')     
    
    
    
    
    
    
    
    