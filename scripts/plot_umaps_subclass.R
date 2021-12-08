

SUBCLASS_COLORS = data.frame( 
    label = c(
        "Lamp5","Sncg","Vip","Sst Chodl","Sst","Pvalb" ,"L2-3 IT",
        "L4-5 IT","L5 IT","L6 IT","L6 IT Car3","L5 PT","L5-6 NP","L6 CT",
        "L6b",'Meis2','Oligo' , 'Astro','Endo', 'VLMC','SMC',  'Peri'  ,'Micro-PVM','unassigned','Other'),
    color =c(
        '#DA808C','#D633FF','#B864CC','#ECD704','#FF9900','#D93137','#C4EC04','#09CCC6','#50B2AD','#A19922', '#5100FF',
        '#0D5B78','#3E9E64','#2D8CB8','#53377D','#C60C0F','#2E3E39','#665C47', '#8D6C62','#697255', '#807059', '#665547','#94AF97','#D3D3D3','#000000'
    )
)


read_h5ad <- function(h5path, groups = c( 'obs','var' ) ){
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
                    # if (endsWith(n ,'_df') ){   # should only pick up samp, ext, run dfs
                    #     output[[vn]][[n]] = do.call(cbind.data.frame, outputlist[[n]] )
                    # }else {
                    output[[vn]][[n]] = outputlist[[n]]
                    # }
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


plot_umap_discrete <- function(adata, value = NULL,  colormap = NULL, ranges=NULL){
    if (FALSE) {"
    
    TODO specify by <normalization>_<variable>
    "}

    umap = adata$obsm$X_umap
    umap['v'] = value
    colnames(umap) = c('x','y','v')
    umap$v[umap$v =='NA'] = 'unassigned'
        
    # sort by number of colors - visually aesthetic
    a =names(sort(table(umap$v),decreasing=TRUE))
    umap=umap[order(match( umap$v,a)),]    
    # umap=umap[order(-as.numeric(factor(umap$v))),]


    # umap = data.frame( x = coords[1,] , y = coords[2,] , v = adata$obs[[obsname]])
    if (is.null(colormap)){
        n = length(unique(umap$v))
        colormap = data.frame(label = a, color = DISTINCT_COLORS[1:n] ,stringsAsFactors=FALSE)
    }
    
    umap = merge(umap, colormap ,by.x = 'v',by.y = 'label' ,sort=FALSE)
    # umap$v <- factor(umap$v , levels =a)
    if (is.null(ranges$x ) | is.null(ranges$y )   ){
        N = nrow(umap)
    } else{
        N = sum(umap$x > ranges$x[1] & umap$x < ranges$x[2] & umap$y > ranges$y[1] & umap$y < ranges$y[2] ) 
    }
    

    ptsize = (100/(N+1))^.5  + .1   
    ptshape =16
    

    fig = ggplot2::ggplot(umap, ggplot2::aes(x = x, y = y,col = v )) +
        ggplot2::geom_point(alpha = .9, shape = ptshape, size = ptsize) + #geom_text( data =annots, aes(x =tSNE_1, y =tSNE_2 ,label = Tissue)   )+
        ggplot2::theme_classic() + 
        ggplot2::scale_colour_manual(breaks=umap$v, values=umap$color)+
        # ggplot2::scale_color_gradient(low = "azure1", high = "dodgerblue4") + 
        ggplot2::xlab('UMAP 1') +
        ggplot2::ylab('UMAP 2') + 
        ggplot2::labs(colour = NULL) +
        ggplot2::theme(legend.position='left')+
        ggplot2::coord_cartesian(xlim = ranges$x, ylim = ranges$y) +
        ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size=3)))

    return(fig)
}


h5paths = dir(pattern = '.h5ad')
for (h5path in h5paths ){



    projid = gsub('.h5ad','',h5path)
    print(projid)
    adata =read_h5ad(h5path , c('obs','uns','obsm'))
    # pdropouts = sum(adata$var$n_cells_by_counts == 0) / length(adata$var$n_cells_by_counts)
    # density = sum(adata$obs$n_cells_by_counts) / (length(adata$var$n_cells_by_counts)* length(adata$obs$n_genes_by_counts))

    vals = adata$uns$MetaMarkers$subclass_assign_hier$predicted
    fig = plot_umap_discrete(adata, vals, SUBCLASS_COLORS,NULL)
    fig = fig   + ggplot2::theme_void()+ ggplot2::theme(legend.position = "none") 
    ggplot2::ggsave(filename = paste0('/data/public_labshare/resource/MetaQC/pcimage/',projid,'_umap.png' ) , dpi = 100, height = 5, width= 5  )
}


