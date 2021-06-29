# given a set of sample and runs dataframe for a project, guess the batch based on sample attributes
import pandas as pd
import ast 

# run this prior to saving the dataframes to runs.tsv
# include another column header in runs.tsv for batch
def impute_batch(sdf, rdf):
    '''
    This should only be run at the project level dataframes, but just in case,
    Splits by project id first and assigns a set batches to each 
    project id 
    '''
    
    cols = rdf.columns.tolist()
    new_rdf = pd.DataFrame(columns= cols.append('batch'))
    for proj , df in sdf.groupby('proj_id') :

        samp2batch = pd.DataFrame({ 'samp_id' : df.samp_id ,
                        'batch': pd.factorize(df['attributes'])[0]})
        tm_rdf = rdf.loc[ rdf.proj_id == proj ,:] 
        # merge these batches with the runs 
        tm_rdf=tm_rdf.merge(samp2batch, how ="left", on="samp_id")

        new_rdf = pd.concat([new_rdf, tm_rdf])

    return new_rdf


def impute_tissue(sdf) :
    # tags = ['source_name','tissue','organ', 'cell_type']

    # look at source name - get unique 

    attr = pd.DataFrame({ 'attr': sdf.attributes.unique() })

    source_names = attr.applymap(lambda x: ast.literal_eval(x)['source_name'])
    # map these back to the sample id? 

    # map these to the runs? 
