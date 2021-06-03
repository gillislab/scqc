# unable to get MetaMarkers to download properly in a conda environment... Cairo dependency issue
library(rhdf5)
library(Matrix)
library(MetaMarkers) 
library(taRifx)
library(tibble)
library(dplyr)
library(SingleCellExperiment)
library(rlist)
library(ggplot2)
# library(reticulate)
# library(scater)
 


# cell ontology obo
# http://purl.obolibrary.org/obo/cl.obo


getMergedBICCN <- function(filepath='/home/johlee/data/biccn/full_biccn.rds'){
    # currently data is only for the mouse primary cortex. Whole brain can be added later.
    # colData should have: joint_class_label  joint_subclass_label joint_cluster_label
    return(readRDS(filepath))
}


# hierarchical annotations vignette
build_marker_sets_biccn <- function( outdir = "SupplementData/BICCN/"){

    # dir.create(outdir, recursive=  TRUE, showWarnings = FALSE)

    message("\nLoading merged BICCN data")
    full_biccn = getMergedBICCN() 
    cpm(full_biccn) = MetaMarkers::convert_to_cpm(assay(full_biccn))

    studies = unique(full_biccn$study_id )
    labelsets=  c("class", "subclass", "cluster" )

    # create meta marker sets for each of these labels 

            
    markers = vector(mode="list", length=length(labelsets))
    names(markers) = labelsets

    meta_markers = vector(mode="list", length=length(labelsets))
    names(meta_markers) = labelsets

    old_labelset=NULL
    for  (labelset in labelsets){
        markers[[labelset]] = list()
        # split by study 
        col_id = paste0("joint_", labelset, "_label") 
        for (study in studies) {
            message("Working on ", study ,", " ,labelset, " markers." )

            tmp = full_biccn[, full_biccn$study_id ==study  ]

            if (is.null(old_labelset) ){
                gl = rep("all",ncol(tmp))
            } else {
                old_col_id = paste0("joint_", old_labelset, "_label") 
                gl = colData(tmp)[,old_col_id ]
            }

            # print(unique(gl))
            

            markers[[labelset]] = rlist::list.append(markers[[labelset]], 
                MetaMarkers::compute_markers(
                    expression = cpm(tmp),   # messy but avoids RAM issue 
                    cell_type_labels =colData(tmp)[,col_id ] ,
                    group_labels = gl
                ) 
            )

        }
        old_labelset = labelset

        message("\nMaking Meta Markers for ", labelset, "\n")
        names(markers[[labelset]]) = studies
        # store meta markers in a list
        meta_markers[[labelset]] = MetaMarkers::make_meta_markers(markers[[labelset]], detailed_stats = TRUE)
        fname = paste0(outdir, "biccn_", labelset, "_marker_set.csv" )
        
        MetaMarkers::export_markers(meta_markers[[labelset]], fname , gzip = TRUE)
        gc()    # garbage collect
    }

    saveRDS(meta_markers ,file = paste0(outdir, "biccn_all_marker_sets.rds" ) )
    
    return(meta_markers)

}


plot_meta_markers <- function( metamarkers , sce ,labelset, outdir = "Figures/BICCN_markers") {
    # meta_markers should be a single df for a single label set   
    # sce is a SingleCellExperiment and should have a cpm field 
    # labelset is the name of the column to consider

    MetaMarkers::plot_pareto_summary( metamarkers, min_recurrence = 1 )
    ggsave( paste0(outdir, "/pareto_summary_",labelset , ".png")  )

    # plot expression levels for each marker on the pareto front
    cell_types = unique(metamarkers$cell_type )
    for (ct in  cell_types) {
        pareto_markers = get_pareto_markers(metamarkers,ct , min_recurrence=0)
        plot_marker_expression(cpm(sce), pareto_markers,  colData(sce)[, paste0("joint_",labelset,"_label") ]) # expression across all datasets in the sce object
        ct = gsub(" ", "_", ct)
        ct = gsub("/","-" ,ct)
        ggsave( paste0( outdir, "/marker_expression_",labelset, "_", ct,".png"))
    }

}


# deprecated
# metamarkers - First vignette "Computation of meta-analytic markers"
build_marker_sets <- function( label_set = "tissue", require_common_types=FALSE) { # tissue, subtissue, cell_ontology_class
    # load the data
    dropletData = get_tabula_data("SupplementData/allTabulaMurisDroplet.Rdata")
    FACSData = get_tabula_data("SupplementData/allTabulaMurisFACS.Rdata")

    if (label_set == "tissue") {
        FACSData
    }


    missing_ind = dropletData$coldata[,label_set] == "" | is.na(dropletData$coldata[,label_set] )
    dropletData$counts = dropletData$counts[ , ! missing_ind  ] 
    dropletData$coldata   = dropletData$coldata[ ! missing_ind , ]
    
    if (label_set =="tissue") {
        # replace the "Heart_and_Aorta" with "Heart"
        dropletData$coldata$tissue[  dropletData$coldata$tissue == "Heart_and_Aorta"  ] = "Heart"
    }

    missing_ind = FACSData$coldata[,label_set] == "" | is.na(FACSData$coldata[,label_set] )
    FACSData$counts = FACSData$counts[ , ! missing_ind  ] 
    FACSData$coldata   = FACSData$coldata[ ! missing_ind , ]
 

    # get the cell type annotations . 
    labelsDroplet = unique(dropletData$coldata[label_set])[,1]
    labelsFACS = unique(FACSData$coldata[label_set])[,1]

 


    # if true, we only take the intersection of cell types between the two datasets 
    if (require_common_types) {
        labels = intersect(labelsFACS,labelsDroplet) 
        
        indDroplet = dropletData$coldata[,label_set] %in% labels 
        indFACS = FACSData$coldata[,label_set] %in% labels 

        # filter the counts matrix
        dropletData$counts = dropletData$counts[,indDroplet  ]
        FACSData$counts = FACSData$counts[,indFACS  ]
        
        # filter the coldata
        dropletData$coldata = dropletData$coldata[indDroplet,  ]
        FACSData$coldata = FACSData$coldata[indFACS,  ]   

    } # otherwise do nothing
    
    # CPM normalize the counts 
    dropletData$cpm = MetaMarkers::convert_to_cpm(dropletData$counts)
    FACSData$cpm = MetaMarkers::convert_to_cpm(FACSData$counts)    

    # compute marker sets for each cell type
    droplet_markers = MetaMarkers::compute_markers(dropletData$cpm, dropletData$coldata[,label_set])
    FACS_markers = MetaMarkers::compute_markers(FACSData$cpm, FACSData$coldata[,label_set])

    # export metamarkers
    fname_droplet = paste0("SupplementData/TabulaMurisDroplet_" , label_set,"_markers.csv")
    fname_FACS = paste0("SupplementData/TabulaMurisFACS_" , label_set,"_markers.csv")
    fname_labelset = paste0("SupplementData/TabulaMurisAll_" , label_set,"_markers.csv")
    MetaMarkers::export_markers(droplet_markers, fname_droplet,gzip=TRUE)
    MetaMarkers::export_markers(FACS_markers, fname_FACS,gzip=TRUE)
    # meta_markers = lapply(markersetpaths, FUN=MetaMarkers::read_meta_markers )
    meta_markers = list(droplet = droplet_markers , FACS = FACS_markers )
    meta_markers = MetaMarkers::make_meta_markers(meta_markers, detailed_stats = TRUE)
    MetaMarkers::export_markers(meta_markers, fname_labelset , gzip = TRUE)

    return(meta_markers)
}


parse_STAR_output <- function(outpath ){

    # is there an mtx file in this directory? 
    if (any(grepl("mtx", dir(outpath) )) ) {
        tryCatch({
            genes = read.csv(paste0(outpath, "/features.tsv") ,sep = "\t", stringsAsFactors=FALSE,header=FALSE )
        },error=function(e){}, warning= function(w){} )

        tryCatch({        
            genes = read.csv(paste0(outpath, "/genes.tsv") ,sep = "\t", stringsAsFactors=FALSE,header=FALSE )
        },error=function(e){}, warning= function(w){} )

        mat = Matrix::readMM(paste0(outpath, "/matrix.mtx"))

        barcodes = read.csv(paste0(outpath, "/barcodes.tsv") ,sep = "\t", stringsAsFactors=FALSE,header=FALSE )

        rownames(mat) = genes[,2]       # the second column contains the gene symbols
        colnames(mat) = barcodes[,1]    # only one column containing the barcodes which identifies the cells

    } else if ( grepl( ".csv", outpath ) ) {    
        mat = as.matrix(read.delim(outpath, sep=",", stringsAsFactors=FALSE,row.names=1)  )
        mat = Matrix(mat, sparse=TRUE)

    # "SupplementData/FACS/Bladder-counts.csv"
    } else if ((any(grepl(".csv", dir(outpath) )) )) {  # otherwise, just read the csv file that STAR outputs
        # do nothing for now.
    }

    

    return(mat)

}


get_top_markers <- function(markersetpath="SupplementData/TabulaMurisAll_tissue_markers.csv.gz"  , max_rank  = 100 ){
    meta_markers = MetaMarkers::read_meta_markers(markersetpath)
    top_markers = dplyr::filter(meta_markers, rank <=max_rank)  
    return(top_markers)
}

# Assign cells hierarchically
assign_cell_type <- function(dataset, top_markers, group_assignment = NULL ,plot_it= FALSE ) {
    # given list of datasets as matrices, asign the cell type from the markerset - requires result from "build_marker_sets"
    # top markers should be a list of marker sets for hierarchical annotations
    # an example path - to be removed later
    # datasetpath = "SupplementData/droplet/Bladder-10X_P4_4" 

    # parse the data, score and assign cells using the top_marker

    # class scores
    ct_scores = MetaMarkers::score_cells(dataset, top_markers)
    ct_enrichment = MetaMarkers::compute_marker_enrichment(ct_scores)
    ct_pred = MetaMarkers::assign_cells(ct_scores , group_assignment = group_assignment ) 



    if (plot_it){
            
        sce = SingleCellExperiment::SingleCellExperiment(dataset)
        hvg = MetaNeighbor::variableGenes(sce , exp_labels = rep("None", ncol(sce) ) )
        umap = MetaMarkers::compute_umap(logcounts[hvg,])
        
        MetaMarkers::plot_assignments(ct_pred, umap_coordinates = umap[,2:3], enrichment_threshold = 1 )
        ggsave("Figures/umap.png")
    }


    return(ct_pred)
    
    # pr[[1]] = summarize_precision_recall(ct_enrichment, true_labels, seq(1,5,0.1)) %>%
    #     filter(get_cell_type(marker_set) == true_label)
}

# not needed
make_biccn_annotation_tree <- function(coldata ){
    # order 1 annotation tree
    coldata = coldata[,grep("*_label" ,colnames(coldata))]
    coldata = coldata[,grep("joint*" ,colnames(coldata))]
    
    classes= unique(coldata$joint_class_label)
    subclasses= unique(coldata$joint_subclass_label)
    clusters= unique(coldata$joint_cluster_label)

    is_a_df = data.frame(label=NULL, is_a=NULL)
    for (cls in classes){
        cd = coldata[coldata$joint_class_label == cls ,]
        cts = unique(cd$joint_subclass_label)
        df = data.frame( label = cts, is_a = cls)
        is_a_df = rbind(is_a_df, df)
    }

    for (scls in subclasses){
        cd = coldata[coldata$joint_subclass_label == scls ,]
        cts = unique(cd$joint_cluster_label)
        df = data.frame( label = cts, is_a = scls)
        is_a_df = rbind(is_a_df, df)
    }
    
    # i = match(is_a_df$label, unique(is_a_df$label) )
    # j = match(is_a_df$is_a, unique(is_a_df$is_a) )
    
    # mat = Matrix::sparseMatrix(i,j , dimnames = list( unique(is_a_df$label), unique(is_a_df$is_a)))

    return(is_a_df)

}


main <- function( acc="SRR11604218" ,STARout_direc= "STARout" ,biccnDirec= "SupplementData/BICCN",max_rank=100 ) {
    # need to figure out a good max_rank for each of the three label sets

    solo_out_path = paste0(STARout_direc ,"/", acc, "Solo.out/Gene/filtered")
    # solo_out_path ='/home/johlee/SCQC/SupplementData/FACS/Brain_Neurons-counts.csv'
    # solo_out_path = '/home/johlee/SCQC/SupplementData/FACS/Brain_Microglia-counts.csv'
    # solo_out_path = '/home/johlee/SCQC/SupplementData/FACS/Lung-counts.csv'
    dataset = parse_STAR_output(solo_out_path)   

    # load the markers - ideally, save the top markers instead of the whole thing
    marker_sets = readRDS(paste0(biccnDirec ,"/biccn_all_marker_sets.rds"))
    
    top_markers = lapply(marker_sets, function(x) {dplyr::filter(x, rank <=50) } )
    logcounts = log1p(MetaMarkers::convert_to_cpm(dataset))

    class_pred = assign_cell_type(logcounts,top_markers$class,  group_assignment=NULL, plot_it=FALSE )
    subclass_pred = assign_cell_type(logcounts,top_markers$subclass,  group_assignment=class_pred$predicted, plot_it=FALSE )
    cluster_pred = assign_cell_type(logcounts,top_markers$cluster,  group_assignment=subclass_pred$predicted, plot_it=FALSE )

    print(table(class_pred$predicted))
    print(table(subclass_pred$predicted) /sum(table(subclass_pred$predicted)))
}
