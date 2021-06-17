


class Query(object):
    """
    
    
    """
    
    def _oldexecute(self):
        """
         Perform query, get ids, fetch for each id, parse XML response. 
         Put project and run info in project_metadata.tsv and project_runs.tsv
         Put completed project ids into query-donelist.txt
         
        """
        # url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&term=%22rna+seq%22[Strategy]+%22mus+musculus%22[Organism]+%22single+cell%22[Text Word]+%22brain%22[Text Word]&retstart=&retmax=50&retmode=json"
        self.log.info('querying SRA...')
        url = f"{self.sra_esearch}&term={self.search_term}&retmax={self.query_max}&retmode=json"
        self.log.debug(f"search url: {url}")
        r = requests.get(url)
        er = json.loads(r.content.decode('utf-8'))
        logging.debug(f"er: {er}")
        idlist = er['esearchresult']['idlist']
        logging.debug(f"got idlist: {idlist}")
        # filter ids by already done.
        donelist = readlist(self.uidfile)
        idlist = listdiff(idlist, donelist)
    
        if len(idlist) > 0:
            allrows = []
            allprojruns = []
            doneids = []
            for id in idlist:
                try:
                    url = f"{self.sra_efetch}&id={id}"
                    self.log.debug(f"fetch url={url}")
                    r = requests.post(url)
                    #self.log.debug(f'status code {r.status_code} type {type(r.status_code)}')
                    if r.status_code == 200:
                        rd = r.content.decode()
                        #logging.debug(f"data for id={id}: {rd}")
                        #rows = self._parse_experiment_pkg(rd)
                        (rows, proj_runs) = self.parse_experiment_package_set(rd)
                        allrows = itertools.chain(allrows, rows)
                        allprojruns = itertools.chain(allprojruns, proj_runs)
                        #allprojruns.append(proj_runs)
                        doneids.append(id)
                    else:
                        self.log.warn(f'bad HTTP response for NCBI uid {id}')
                    
                except Exception as ex:
                    self.log.error(f'problem with NCBI uid {id}')
                    logging.error(traceback.format_exc(None))
                
                finally:
                    self.log.debug(f"sleeping {self.query_sleep} secs between fetch calls...")
                    time.sleep(self.query_sleep)
    
    
            newdone = listmerge(donelist, doneids)
            self.log.info(f'updating uid done list...')
            writelist(self.uidfile, newdone)
            
            filepath = f'{self.metadir}/all_metadata.tsv'
            self.log.info(f'updating metadata df: {filepath}')
            adf = pd.DataFrame(allrows, columns = META_COLUMNS)
            #adf = self._impute_tech(adf) 
            self.log.debug(f'made all df: {adf}')           
            adf.drop_duplicates(inplace=True)
            merge_write_df(adf, filepath )
            #df["status"] = "tech_imputed"            
    
            filepath = f"{self.metadir}/project_runs.tsv"            
            self.log.info(f'updating proj-run df: {filepath}')
            self.log.debug(f'making dataframe from proj-runs: {allprojruns}')
            pdf = pd.DataFrame(allprojruns, columns = PROJ_RUN_COLUMNS)
            self.log.debug(f'made project-run df: {pdf}')   
            merge_write_df(pdf, filepath )
                        
            srplist = list(pdf.project.unique())
            #sralist = list(itertools.chain.from_iterable(sl))
            return srplist
            
            #self._split_df_by_project(df)   # saves dfs by project accession
        else:
            self.log.info('no new uids to process...')
            return []
    
    
    
    
    # added sample attributes
    def _parse_experiment_pkg(self, xmlstr):
        root = et.fromstring(xmlstr)
        self.log.debug(f"root={root}")
        rows = []
        for exp in root.iter("EXPERIMENT_PACKAGE"):
            for lcp in exp[0].iter("LIBRARY_CONSTRUCTION_PROTOCOL"):
                lcp = lcp.text
            SRXs = exp[0].get('accession')
            SRAs = exp[1].get('accession')
            SRPs = exp[3].get('accession')
            alias = exp[0].get('alias')
            # title = exp[3][1][0].text
            abstract = exp[3][1][2].text
            SRRs = []
            date = []
            taxon = []
            orgsm = []
            sample_attrib = {}
            for study in exp[3].iter("STUDY_TITLE"):
                title = study.text
            for study in exp[3].iter("STUDY_ABSTRACT"):
                abstract = study.text

            for run in exp.iter("RUN"):
                SRRs.append(run.attrib['accession'])
                date.append(run.attrib['published'])
                for mem in run.iter("Member"):
                    taxon.append(mem.attrib['tax_id'])
                    orgsm.append(mem.attrib['organism'])

            for sample in exp.iter("SAMPLE_ATTRIBUTE"):
                tag = sample[0].text
                value = sample[1].text
                sample_attrib[tag] = value

            row = [SRPs, SRXs, SRAs, SRRs, alias, date, taxon,
                   orgsm, lcp, title, abstract, sample_attrib]
            self.log.debug(f'got SRRs: {SRRs}')
            rows.append(row)
        return rows    


    
def parse_experiment_pkg(self, root):    
    for lcp in exp[0].iter("LIBRARY_CONSTRUCTION_PROTOCOL"):
        lcp = lcp.text
        SRXs = exp[0].get('accession')
        SRAs = exp[1].get('accession')
        SRPs = exp[3].get('accession')
        alias = exp[0].get('alias')
        # title = exp[3][1][0].text
        abstract = exp[3][1][2].text
        
        SRRs = []
        date = []
        taxon = []
        orgsm = []
        sample_attrib = {}
        for study in exp[3].iter("STUDY_TITLE"):
            title = study.text
        for study in exp[3].iter("STUDY_ABSTRACT"):
            abstract = study.text

        for run in exp.iter("RUN"):
            SRRs.append(run.attrib['accession'])
            date.append(run.attrib['published'])
            for mem in run.iter("Member"):
                taxon.append(mem.attrib['tax_id'])
                orgsm.append(mem.attrib['organism'])

        for sample in exp.iter("SAMPLE_ATTRIBUTE"):
            tag = sample[0].text
            value = sample[1].text
            sample_attrib[tag] = value

        row = [SRPs, SRXs, SRAs, SRRs, alias, date, taxon,
               orgsm, lcp, title, abstract, sample_attrib]
        self.log.debug(f'got SRRs: {SRRs}')
        rows.append(row)
    return rows
