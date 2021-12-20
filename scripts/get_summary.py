import glob
import os
import numpy as np
import scanpy as sc
import h5py
import pandas as pd
import matplotlib.pyplot as plt
from scipy import sparse
plt.style.use('seaborn-deep')

# dropped datasets: meninges, immune and vascular
# TODO fix NEMO header box field. 


DATADIR = '/home/ftp/data/MetaQC/'
# for_product_chart = pd.merge(ranked_result, pdf , left_index= True, right_on='proj_id', how ='left').reset_index(drop=True)
# for_product_chart = pd.merge(for_product_chart ,all_prs, right_index=True, left_on='proj_id', sort = 'left')

# allrows=[]
# for pid, df in global_dist.iloc[:,-3:].groupby('proj_id'):
#     ntech = df.tech.value_counts().to_dict()
#     ncell = df.shape[0]
#     date = df.pdat.mean()
#     row = [pid, ncell, ntech, date ]
#     allrows.append(row)
    
# tmpdf = pd.DataFrame(allrows,columns = ['proj_id','n_cells','n_techs','date'])

# for_product_chart = pd.merge(for_product_chart, tmpdf, on = 'proj_id',how='left')


CTS_LABS =[
    "total_counts",
    "n_genes_by_counts",
    "corr_to_mean",
    "gini",
    "ribo",
    "mt",
    "essential",
    "cell_cycle",
    "housekeeping",  
    "female",
    "male",  
    "highly_variable",
    "top_50_gene" ,
    "top_100_gene" ,
    "top_200_gene" ,
    "top_500_gene" 
]

# do we drop meninges? 

datasets2drop = [
    'SRP278583', # multi tissue T cells
    'SRP066963', # too few cells
    'SRP308387',# too few cells
    'SRP239491',# too few cells
    'SRP150863',# too few cells
    'SRP142629',# too few cells
    'SRP150630',# too few cells
    'SRP071876',# too few cells
    'SRP308826',# too few cells
    'SRP066314',# too few cells
    'SRP228572',# too few cells
    'SRP108034',# too few cells
    'SRP303200',# too few cells
    
    
    ]

def main(outdir = '/home/ftp/data/MetaQC',datasets2drop = []):
    global_dist, all_prs = get_global_distributions(datasets2drop = datasets2drop)
    global_bin_heights ,global_breaks ,proj_hist = get_global_hist_data(global_dist)
    proj_bin_heights = proj_hist.loc[0,:].set_index('proj_id')
    auc , ranked_auc = get_aurocs_with_global(global_bin_heights,proj_bin_heights)
    auc.to_csv(f'{outdir}/auc_wrt_global.tsv',sep="\t")
    ranked_auc.to_csv(f'{outdir}/ranked_auc_wrt_global.tsv',sep="\t")
    


    # save in h5file
    for stat_opt in proj_bin_heights.columns:
        
        df = pd.DataFrame()
        for pid in proj_bin_heights.index:
            df_col = pd.DataFrame(proj_bin_heights.loc[pid,stat_opt])
            df = pd.concat([df,df_col],axis=1)
        df.columns = proj_bin_heights.index

        with h5py.File(f'{outdir}/proj_bin_heights.hdf5','a') as f : 
            f.create_dataset(stat_opt, data =df)
            
    with h5py.File(f'{outdir}/proj_bin_heights.hdf5','a') as f : 
        f.create_dataset('projid_order', data =list(df.columns))
        f.create_dataset('global_bin_heights', data =global_bin_heights)
        f.create_dataset('breaks', data =global_breaks)
        f.create_dataset('stat_order', data =list(global_bin_heights.columns))


    # save the auroc between the project and global
    pass

def gini_coefficient_fast(X):
    """ 
        expects a CSR sparse matrix
        Sorting is O(n log n ) (here n is at most number of genes)
        loops over cells (m) instead of gene pairs. 
        Overall, at most O( m n log n)  but realistically, 
        density of 10% -> `effective n` is 0.1 * n
        
        only looks at nonzero elements
        Cells with no expression get a gini score of 0       
    """    
    # x = np.asarray(x)
    g = np.zeros(X.shape[0])    # ncells
    n = X.shape[1]          # ngenes
    for i in range(X.shape[0]): # loops for all cells
        # take the nonzero elements of the ith cell
        x = X[i,:]  
        x= x[:, x.indices].A.flatten()

        sorted_x = np.sort(x)   
        cumx = np.cumsum(sorted_x, dtype=float)

        if len(cumx) == 0 : # cell with zero expression - perfect equilibrium
            g[i] = 0
        else :
            g[i] =(n + 1 - 2 * np.sum(cumx) / cumx[-1]) / n
    return g


def get_global_distributions(datasets2drop=[] , saveit=False):
    global_dist = pd.DataFrame()
    h5paths = glob.glob(f'{DATADIR}*.h5ad')
    all_prs = pd.DataFrame()
    for h5path in h5paths:
        proj_id = os.path.basename(h5path).replace('.h5ad','')
        if proj_id not in datasets2drop:
            with h5py.File(h5path,'r') as f:
                print(h5path)
                ds_df = pd.DataFrame()
                
                for obsname in CTS_LABS:
                    ds_df[obsname]= f[f'obs/{obsname}'][()]
                    if obsname == 'total_counts' or obsname == 'n_genes_by_counts': 
                        ds_df[obsname] = np.log(ds_df[obsname])
                
                ds_df['pdat'] = f['uns']['header_box']['publish_date'][()]

                techs = f['uns']['header_box']['tech_count'].keys()
                tech_counts ={}
                for tech in techs:
                    tech_counts[tech] = f['uns']['header_box']['tech_count'][tech][()] 
                if len(tech_counts) == 1:
                    ds_df['tech'] = list(tech_counts.keys())[0]
                else : 
                    ds_df['tech'] = 'mixture'

                prec =  np.append( 1, f['uns']['MetaMarkers']['class_PR']['Precision'][()])
                prec = np.append(prec, 0 )

                recall =  np.append( 0,f['uns']['MetaMarkers']['class_PR']['Recall'][()])
                recall = np.append(recall, 1 )
                auprc_class = np.trapz(y=prec, x=recall)

                prec =  np.append( 1,f['uns']['MetaMarkers']['subclass_PR']['Precision'][()])
                prec = np.append(prec, 0 )
                recall =  np.append( 0,f['uns']['MetaMarkers']['subclass_PR']['Recall'][()])
                recall = np.append(recall, 1 )

                auprc_subclass = np.trapz(y=prec, x=recall)

                prec =  np.append( 1,f['uns']['MetaMarkers']['subclass_PR_hier']['Precision'][()])
                prec = np.append(prec, 0 )
                recall =  np.append( 0,f['uns']['MetaMarkers']['subclass_PR_hier']['Recall'][()])
                recall = np.append(recall, 1 )
                
                auprc_subclass_h = np.trapz(y=prec, x=recall)

                gini  = gini_coefficient_fast(sparse.csr_matrix( f['obs/total_counts'][()]))
                prs = pd.DataFrame({'class_pr':auprc_class, 'subclass_pr': auprc_subclass, 'subclass_pr_hier': auprc_subclass_h,'gini':gini },index = [proj_id])
                

        ds_df['proj_id'] = proj_id
        all_prs =all_prs.append(prs) 
        global_dist = global_dist.append(ds_df)
        print(global_dist.shape)


    global_dist['pdat'] = pd.to_datetime(global_dist.pdat.str.decode('utf-8'))
    global_dist = global_dist.reset_index(drop=True)
    if saveit:
        for obsname in CTS_LABS:
            hist_data = plt.hist(global_dist[obsname],50)

            with h5py.File(f'{DATADIR}global_histograms.hdf5','a') as f:
                f.create_dataset(f'{obsname}/frequency',data = hist_data[0], dtype=int)
                f.create_dataset(f'{obsname}/breaks',data = hist_data[1], dtype=float)

    return global_dist ,all_prs



def get_global_hist_data(global_dist,saveit=False):
    # density normalize by dataset
    hist_df = list()
    global_breaks  =global_dist.iloc[:,:-3].apply(lambda x : plt.hist(x,50,density=True ),axis =0 ) 

    for pid, df in  global_dist.groupby('proj_id'):
        print(pid)
        tmpdf = df.iloc[:,:-3].apply(lambda x : plt.hist(x,bins = global_breaks[x.name][1]  ,density=True),axis =0 )
        tmpdf['proj_id'] = pid
        
        hist_df.append(tmpdf)
    
    hist_df = pd.concat(hist_df,axis=0)
    densities = hist_df.loc[0,:]
    bin_heights = densities.iloc[:,:-1].apply(lambda x: x.mean(),axis=0)
    
    breaks_df = pd.DataFrame()
    for i in range(global_breaks.shape[1]):
        breaks_df[global_breaks.columns[i]] = global_breaks.iloc[1,i]

    if saveit :
        bin_heights.to_csv(f'{DATADIR}/global_histogram_heights.tsv',sep="\t")
        breaks_df.to_csv(f'{DATADIR}/global_histogram_breaks.tsv',sep="\t")
        tmp = hist_df.loc[0,:].set_index(hist_df.loc[0,:].proj_id)
        # tmp = tmp.apply(lambda x : str(x) , axis=0)
        tmp.to_csv(f'{DATADIR}/all_histogram_heights.tsv',sep="\t")

        # hist_df.loc[1,:].set_index(hist_df.loc[1,:].proj_id).to_csv(f'{DATADIR}/all_histogram_breaks.tsv',sep="\t")
    return  bin_heights ,breaks_df ,hist_df    # bin heights , breaks 


def get_aurocs_with_global(bin_heights,hist_df,saveit=False):
    # bin heights correspond to the global distributions
    # hist_df has distributions of each dataset
    
    global_cs = np.cumsum(bin_heights, axis =0 )
    result = pd.DataFrame()
    all_heights = hist_df.copy()
    for pid in all_heights.index:
        pid_heights = all_heights.loc[all_heights.index==pid,:]

        pid_df = pd.DataFrame()
        for obsname in pid_heights.columns: 
            pid_df[obsname] = pid_heights.loc[pid,obsname]
        pid_cs = np.cumsum(pid_df,axis=0)

        ranks = global_cs.append(pid_cs).apply(lambda x : x.rank(),axis =0 ).reset_index(drop=True)
        labels = pd.DataFrame({'label': ['global']*global_cs.shape[0] + ['project']*pid_cs.shape[0] })
        label_matrix =pd.get_dummies(labels.label)

        npos = label_matrix.sum(0)
        nneg = label_matrix.shape[0] - npos

        sum_pos = ranks.T @ label_matrix

        result[pid]=  1- ((sum_pos /npos - (npos+1)/2) /nneg).project
         
    ranked_result = result.rank(1) / result.shape[1]

    if saveit:
        ranked_result.T.to_csv(f'{DATADIR}/ranked_aurocs.tsv',sep="\t")
        result.T.to_csv(f'{DATADIR}/aurocs.tsv',sep="\t")

    return result.T , ranked_result.T
    



def plot_ncell_year(global_dist):

    global_dist['year'] = global_dist.pdat.apply( lambda x:x.year)
    cell_counts = global_dist['year'].value_counts().sort_index()
    plt.plot(cell_counts)
    plt.xlabel('Year')    
    plt.xlabel('Number of Cells')    
    plt.title('Brain Cells over time')
    plt.savefig(f'{DATADIR}ncell_year.png')
    plt.close()

    tmp = global_dist[['year','proj_id']].drop_duplicates()
    ds_count = tmp.year.value_counts().sort_index()
    plt.plot(ds_count)
    plt.xlabel('Year')    
    plt.xlabel('Number of Datasets')    
    plt.title('Brain Datasets over time')
    plt.savefig(f'{DATADIR}nds_year.png')
    plt.close()

def plot_ss_v_10x(global_dist):
    obsname = global_dist.columns[:17]
    global_dist['tech'] = global_dist['tech'].str.replace('v[0-9]','',regex=True)
    for obs in obsname:
        for tech , df in global_dist.groupby('tech'):
            if obs =='total_counts' or obs == 'n_genes_by_counts':
                dat = np.log(df[obs])
                xlab = f'log({obs})'
            else: 
                dat = df[obs]
                xlab='AUROC'
            if obs in obsname[:5]:
                xlab = obs
            plt.hist(dat , bins = 50, alpha = 0.3 ,label = tech,density = True)
            
        plt.legend()
        if obs == 'total_counts' or obs =='n_genes_by_counts':
            xlab = f'log({obs})'

        plt.xlabel(xlab)
        plt.ylabel('Density')
        # plt.title(obs)
        plt.savefig(f'/home/ftp/data/MetaQC/figures/{obs}_ss_10x.png')
        plt.close()
        

def get_density_dropout(h5paths,saveit=False) : 

    allrows = []
    for h5path in h5paths : 
        projid = os.path.basename(h5path.replace('.h5ad',''))
        print(projid)
        adata =sc.read_h5ad(h5path )
        pdropouts =(adata.var.n_cells_by_counts == 0).sum() / len(adata.var.n_cells_by_counts) * 100
        density = sum(adata.var.n_cells_by_counts) / (len(adata.var.n_cells_by_counts)* len(adata.obs.n_genes_by_counts))
        allrows.append([projid, pdropouts,density])
    
    df = pd.DataFrame(allrows, columns = ['proj_id','percent_dropout','density'])

    if saveit:
        df.to_csv(f'{DATADIR}/add2main_forPC.tsv')
    return df

def count_MM_assignments(h5paths,saveit=False):
    
    allrows = []
    for h5path in h5paths : 
        projid = os.path.basename(h5path.replace('.h5ad',''))
        print(projid)
        adata =sc.read_h5ad(h5path )
        pclass = adata.uns['MetaMarkers']['class_assign'].predicted.value_counts().idxmax()
        psubclass = adata.uns['MetaMarkers']['subclass_assign_hier'].predicted.value_counts().idxmax()

        allrows.append([projid, pclass,psubclass])
    
    df = pd.DataFrame(allrows, columns = ['proj_id','primary_class','primary_subclass'])
    
    if saveit:
        df.to_csv(f'{DATADIR}/add2main_forPC2.tsv')
    return df

def fix_prcs(h5paths ):

    for h5path in h5paths :
        try:
            print(h5path)
            adata = sc.read_h5ad(h5path)
            pr = MetaMarkers_PR(adata.uns['MetaMarkers']['class_enr'], class_pred = None)
            adata.uns['MetaMarkers']['class_PR'] = pr 

            pr = MetaMarkers_PR(adata.uns['MetaMarkers']['subclass_enr'], class_pred = None)
            adata.uns['MetaMarkers']['subclass_PR'] = pr 

            pr = MetaMarkers_PR(adata.uns['MetaMarkers']['subclass_enr_hier'], class_pred = None)
            adata.uns['MetaMarkers']['subclass_PR_hier'] = pr 

            adata.write_h5ad(h5path)
        except:
            pass


def get_dataset_gini(global_dist):
    ginis = pd.DataFrame()
    h5paths = glob.glob(f'{DATADIR}*.h5ad')
    for h5path in h5paths:
        proj_id = os.path.basename(h5path).replace('.h5ad','')
        with h5py.File(h5path,'r') as f:
            print(h5path)
            total_counts = sparse.csr_matrix(f['obs/total_counts'][()])
            log_total_counts = sparse.csr_matrix(np.log(f['obs/total_counts'][()]))
            
            tmp = pd.DataFrame({'gini':gini_coefficient_fast(total_counts),
                'log_gini':gini_coefficient_fast(log_total_counts)},index = [proj_id])
            
            ginis = ginis.append(tmp)

    return ginis
    

def plot_auprc_v_all(ginis,all_prs) :
    plt.plot(all_prs.class_pr, all_prs.subclass_pr_hier,'.')
    plt.plot(all_prs.class_pr, all_prs.subclass_pr,'.')
    plt.xlabel('Class AUPRC')
    plt.ylabel('Subclass AUPRC')
    plt.legend(['Hierarchical' ,'Classical'])

    
    for i in all_prs.columns:
        plt.plot(all_prs[i],all_prs.gini,'.')
        plt.xlabel(i)
        plt.ylabel('Gini')
        plt.savefig(f'{DATADIR}/figures/gini_{i}.png')
        plt.close()

    pass