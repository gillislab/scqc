


class Query(object):
    """
    
    
    """
    
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
