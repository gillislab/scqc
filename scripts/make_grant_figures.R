library(rhdf5)


filter_PR <- function(pr_df){
    # used to make ggplot faster - give it fewer points
    # pr_df = as.data.frame(pr_df)

    # pr_df = data.frame('Precision' = prec, Recall = recall)
    n = nrow(pr_df)
    prec = pr_df$Precision
    recall = pr_df$Recall

    slopes = (prec[3:n] - prec[1:(n-2)]) / (recall[3:n] - recall[1:(n-2)])
    to_keep = slopes != 0 & slopes != -Inf & slopes != Inf

    
    pr_df = pr_df[ c(TRUE, to_keep, TRUE) ,]
    return(pr_df)
}


read_h5ad <- function(h5path, groups = c( 'obs','var','obsm','uns','obsp' ) ){
    if (FALSE){"
    Currently only reads 1D vectors from obs and var
    obs/var     - 1D vectors for the cell and gene respectively (N / M)

    TODO: include option to read the obsm/varm/obsp/varp/uns, store as lists.
    "}

    # 'var' , 'obs' etc
    output = vector(mode='list' ,length= length(groups) )
    names(output) = groups
    for (vn in groups){
        flag = 0
        tryCatch({
            outputlist = rhdf5::h5read(h5path, name = vn )
            flag = 1
        },error = function(e){}) 

        if (flag == 1 ){
            if (vn == 'obsm') {
                # ensure that the dimensions match
                # should be  - x nCells
                for (n in names(outputlist)) {
                    if (any(class(outputlist[[n]]) == 'list')) {
                        output[[vn]][[n]] = do.call(cbind.data.frame, outputlist[[n]] )
                    } else if (any(class(outputlist[[n]]) == 'matrix' )) { 
                        output[[vn]][[n]] = data.frame(t(outputlist[[n]]))
                    }
                }

            } else if (endsWith(vn, 'p')){
                #TODO
            } else if(vn =='uns') {
                for (n in names(outputlist)){   
                    # which ones should we convert to dataframes?
                    if (endsWith(n ,'_df') ){   # should only pick up samp, ext, run dfs
                        output[[vn]][[n]] = do.call(cbind.data.frame, outputlist[[n]] )
                    }else {
                        output[[vn]][[n]] = outputlist[[n]]
                    }
                }

            } else { # obs or var
                # transform categorical data
                cats = outputlist$'__categories'
                for (cat in names(cats)){
                    ind = outputlist[[cat]]
                    outputlist[[cat]] = outputlist$'__categories'[[cat]][ind+1]
                }

                l = lengths(outputlist)
                output[[vn]] = data.frame(outputlist[l == l['_index'] ])
                # rownames(output[[vn]]) =output[[vn]][['X_index']]
            }
        }
    }

    for (g in groups) {
        tryCatch({
            output[g] = taRifx::remove.factors(output[g])
        },error = function(e){})
    }
        


    rhdf5::h5closeAll()
    return(output)
}

get_auc  <-function(x,y){
    o = order(x)
    auc = sum(diff(x[o]) * zoo::rollmean(y[o], 2))
    return(auc)
}

plot_pr_curves <- function(pr_df ,filename){
    # df = data.frame(Recall = adata$uns[[id]]['Recall'] , Precision = adata$uns[[id]]['Precision'], Enrichment = adata$uns[[id]]['value'] )
    # pr_df is a list of dataframes.
    
    fig = ggplot2::ggplot(pr_df, ggplot2::aes(x = Recall, y= Precision,group = proj_id)) + 
            ggplot2::geom_line(color='grey20',size=0.5)+
            ggplot2::theme_classic() +
            ggplot2::guides(colour=ggplot2::guide_legend(title='Project'))+
            ggplot2::xlim(0,1)+ggplot2::ylim(0,1)

    ggplot2::ggsave(filename,width=3,height=3)
    
    return(fig)
}


plot_scatter <- function(val1, val2, val1name,val2name,filename = '~/test.png'){
    df = data.frame( x = val1, y = val2 )

    fig = ggplot2::ggplot(df , ggplot2::aes(x = x , y=y) ) + 
        ggplot2::geom_point() +
        ggplot2::theme_classic()+
        ggplot2::xlab(val1name) + 
        ggplot2::ylab(val2name)+
        ggplot2::xlim(0,1)+ggplot2::ylim(0,1)

        # ggplot2::geom_smooth(method = "lm",formula = y ~ x, se = TRUE)+
        # ggpubr::stat_cor(method="pearson",na.rm =TRUE, label.x.npc = 'center',label.y.npc = 'top')


    ggplot2::ggsave(filename,width=3,height=3)
    return(fig)
}

gini <- function(x ){
    n = length(x)          # ngenes
    # for i in range(x.shape[0]): # loops for all cells
        # take the nonzero elements of the ith cell

    sorted_x = sort(x)   
    cumx = cumsum(sorted_x)
    g =(n + 1 - 2 * sum(cumx) / cumx[length(x)]) / n

    return(g)
}
# data/hover/scqc/backup/20211026.1425.prebiccn/output
h5paths =dir('/data/hover/scqc/backup/20211026.1425.prebiccn/output',pattern = 'h5ad',full.names=TRUE)
i=1

for (h5path in h5paths){

    proj_id = sub('.h5ad','',basename(h5path))
    print(proj_id)
    adata = read_h5ad(h5path,c('uns','obs'))

    # a= adata$obs[c('X_index', 'cell_id','class_label','subclass_label' ,'gini','total_counts')]
    g = gini(adata$obs$total_counts)    # dataset specific
    
    class_pr = adata$uns[['MetaMarker_class_PR']]
    class_pr = data.frame(class_pr ,row.names='X_index'   )
    class_pr = filter_PR(class_pr)
    class_pr['proj_id'] = proj_id

    subclass_pr = adata$uns[['MetaMarker_subclass_PR']]
    subclass_pr = data.frame(subclass_pr ,row.names='X_index'   )
    subclass_pr = filter_PR(subclass_pr)
    subclass_pr['proj_id'] = proj_id

    subclass_pr_heir = adata$uns[['MetaMarker_subclass_PR_heir']]
    subclass_pr_heir = data.frame(subclass_pr_heir ,row.names='X_index'   )
    subclass_pr_heir = filter_PR(subclass_pr_heir)
    subclass_pr_heir['proj_id'] = proj_id


    # calculate the AUCs
    tmpauc_class = get_auc(class_pr$Recall, class_pr$Precision)
    tmpauc_subclass = get_auc(subclass_pr$Recall, subclass_pr$Precision)
    tmpauc_subclass_heir = get_auc(subclass_pr_heir$Recall, subclass_pr_heir$Precision)
    
    

    if (i ==1){
        allclass_pr = class_pr
        allsubclass_pr = subclass_pr
        allsubclass_pr_heir = subclass_pr_heir

        auc_class  = data.frame('proj_id' =proj_id,'gini'=g, 'AUC' = tmpauc_class ) 
        auc_subclass  = data.frame('proj_id' =proj_id,'gini'=g,'AUC' = tmpauc_subclass ) 
        auc_subclass_heir  = data.frame('proj_id' =proj_id,'gini'=g,'AUC' = tmpauc_subclass_heir ) 

        allgini = data.frame('proj_id' =proj_id, 'gini'= adata$obs$gini,'cell_id' =adata$obs$cell_id)
    } else{
        allclass_pr = rbind(allclass_pr,class_pr)
        allsubclass_pr = rbind(allsubclass_pr,subclass_pr)
        allsubclass_pr_heir = rbind(allsubclass_pr_heir,subclass_pr_heir)

        auc_class  = rbind(auc_class,c('proj_id' =proj_id,'gini'=g,'AUC' = tmpauc_class ) ) 
        auc_subclass  = rbind(auc_subclass,c('proj_id' =proj_id,'gini'=g,'AUC' = tmpauc_subclass ) ) 
        auc_subclass_heir  = rbind(auc_subclass_heir,c('proj_id' =proj_id,'gini'=g,'AUC' = tmpauc_subclass_heir ) ) 

        allgini = rbind(allgini, c('proj_id' =proj_id, 'gini'= adata$obs$gini,'cell_id' =adata$obs$cell_id))
    }

    i  = i + 1
}

aucs = data.frame(class_auc = as.numeric(auc_class$AUC) , subclass_auc = as.numeric(auc_subclass$AUC),subclass_auc_h = as.numeric(auc_subclass_heir$AUC))
auc2 = as.data.frame(data.table::melt (aucs))

ggplot(aucs,aes (x=class_auc,y = subclass_auc_h) )+ geom_point() + theme_classic() +xlab('Class AUPRC') + ylab('Subclass AUPRC')
ggsave('~/test.png',width =4,height=4)

ggplot(auc2,aes(x = value,fill= variable)) + geom_histogram() + theme_classic()



# plot the figures
fig1 = plot_pr_curves(allclass_pr,filename = '~/MM_class_PR.png')
fig2 = plot_pr_curves(allsubclass_pr,filename = '~/MM_subclass_PR.png')
fig3 = plot_pr_curves(allsubclass_pr_heir,filename = '~/MM_subclass_PR_heir.png')

fig4 = plot_scatter(as.numeric(auc_class$AUC), as.numeric(auc_class$gini),'AUPRC','Gini Coefficient',filename ='~/gini2AURPC_class.png')
fig5 = plot_scatter(as.numeric(auc_subclass$AUC), as.numeric(auc_subclass$gini),'AUPRC','Gini Coefficient',filename ='~/gini2AURPC_subclass.png')
fig6 = plot_scatter(as.numeric(auc_subclass_heir$AUC), as.numeric(auc_subclass_heir$gini),'AUPRC','Gini Coefficient',filename ='~/gini2AURPC_subclass_heir.png')
 





##########################################
# make  some toy data
cmap = RColorBrewer::brewer.pal(8,'Set2')
names(cmap)= paste0('c',1:8)

mat = matrix(runif(3*8,0,1 ),8,3)
rownames(mat) = paste0('c', 1:8)
colnames(mat) = paste0('Label ',1:3)
make_pr_curve(mat)
# plot the enrichment matrix
png('~/PR_enrichment.png'); 
gplots::heatmap.2(mat, Rowv = FALSE ,Colv=FALSE,
trace ='none',dendrogram='none',cexRow=1, cexCol=1 ,
RowSideColors=cmap,offsetRow=-32.5) ; 
dev.off()
# plot the reshaped enrichment matrix
mat2 = matrix(mat,8*3,2)
labmat = mat

for (r in rownames(mat)){
    for (c in colnames(mat)){
        labmat[r,c] = paste(c,cmap[r])
    }
}


labmat = matrix(labmat,8*3,1)
o = order(mat2[,1])
rownames(mat2) =labmat[o,1]

cols = sub('Label [0-9] ','',rownames(mat2))
labs = names(cmap[match(cols,cmap)])
png('~/PR_enrichment_flat.png')
gplots::heatmap.2(mat2[o,], 
    Rowv = FALSE ,Colv=FALSE,trace ='none',dendrogram='none',
    cexRow=1, cexCol=1 ,RowSideColors=cols,
    labRow= labs,offsetRow=-32.5) ; dev.off()


bw = c('#000000','#FFFFFF')

png('~/PR_enrichment_flat_2.png'); 
gplots::heatmap.2(mat2[o,], 
    Rowv = FALSE ,Colv=FALSE,trace ='none',dendrogram='none',
    cexRow=1, cexCol=1 ,RowSideColors=bw[as.numeric(duplicated(cols))+1],
    labRow= labs,offsetRow=32.5) ; dev.off()

    


make_pr_curve <-function(enr ){
    pr = data.frame(value = matrix(enr), cell =rownames(enr))
    
    ord = order(pr[,'value'],decreasing=TRUE)
    pr= pr[ord,]
    rownames(pr)=NULL
    pr['dup_cell'] = ! duplicated(pr$cell)
    pr['TP'] = cumsum(pr['dup_cell'])
    pr['P'] = as.numeric(rownames(pr)) 
    pr['Recall'] = pr$TP / dim(enr)[1]
    pr['Precision'] = pr$TP / pr$P
    ggplot2::ggplot(pr, ggplot2::aes(x = Recall,y=Precision) )+ggplot2::geom_line()+ ggplot2::theme_classic()+
    ggplot2::xlim(0, 1) +ggplot2::ylim(0, 1) 
    ggplot2::ggsave('~/test.png',width=3,height=3)
}

