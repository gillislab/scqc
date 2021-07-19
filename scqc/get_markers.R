library(optparse)

# make sure we use the R in conda and not the system R
# configr::read.config(file="/home/johlee/git/scqc/etc/scqc.conf")

option_list = list(  
    optparse::make_option(
        c("-m","--mode"), type="character", default='annotate', 
        help="setup or annotate", metavar="mode"),

    optparse::make_option(
        c("-d", "--markerdir"), type="character", default='~/scqc/supplement_data/markersets/MoP', 
        help="Markerset Directory to read or save", metavar="markerdir"),

    optparse::make_option(
        c("-o", "--outprefix"), type="character", default="~/scqc/metamarker/", 
        help="Output prefix", metavar="outprefix"),

    optparse::make_option(
        c("-r", "--max_rank"), type="integer", default=100, 
        help="Maximum rank of marker genes. Ignored if setup.", metavar="max_rank"),

    # optparse::make_option(
    #     c("-p", "--rds_path"), type="character", default="~/scqc/supplement_data/markersets/biccn_MoP.rds", 
    #     help="rds path to SingleCellExperiment of merged data.", metavar="rds_path"),

    optparse::make_option(
        c("-s", "--subclass_markerset"), type="character", default=NULL, 
        help=" ", metavar="subclass_markerset"),

    optparse::make_option(
        c("-c", "--class_markerset"), type="character", default=NULL, 
        help=" ", metavar="class_markerset"),

    optparse::make_option(
        c("-p", "--h5path"), type="character", default=NULL, 
        help="STARsolo Gene output directory (required if annotate, ignored otherwise).", metavar="solo_out_dir")


)

opt_parser = optparse::OptionParser(option_list=option_list);
opt = optparse::parse_args(opt_parser);



# expects the user to merge data before hand 
# should be a single cell experiement object with 
get_merged_data <- function(filepath='~/scqc/supplement_data/markersets/biccn_MoP.rds'){
    # currently data is only for the mouse primary cortex. Whole brain can be added later.
    # colData should have: joint_class_label  joint_subclass_label joint_cluster_label
    return(readRDS(filepath))
}

# only done once
# hierarchical annotations vignette
build_marker_sets_biccn <- function( rds_path='~/scqc/supplement_data/markersets/biccn_MoP.rds' ,outdir = "~/scqc/supplement_data/markersets/MoP"){
    # dir.create(outdir, recursive=  TRUE, showWarnings = FALSE)


    if (! file.exists(paste0(outdir, "/all_marker_sets.rds" )) ){

        library(SingleCellExperiment)

        message("\nLoading merged BICCN data")
        full_biccn = get_merged_data(rds_path) 
        cpm(full_biccn) = MetaMarkers::convert_to_cpm(assay(full_biccn))

        studies = unique(full_biccn$study_id )
        labelsets=  c("class", "subclass" )

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


                markers[[labelset]] = rlist::list.append(markers[[labelset]], 
                    MetaMarkers::compute_markers(
                        expression = cpm(tmp),   
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
            fname = paste0(outdir, "/",labelset, "_marker_set.csv" )
            
            MetaMarkers::export_markers(meta_markers[[labelset]], fname , gzip = TRUE)
            gc()    # garbage collect
        }

        saveRDS(meta_markers ,file = paste0(outdir, "/all_marker_sets.rds" ) )

        return(meta_markers)

    } else {
        warning("Marker files have already been generated in this out directory. \n", 
                "... Extract data from all_marker_sets.rds")
    }


}




##### class annotate
parse_h5ad <- function( h5path){
    X = rhdf5::h5read(h5path, name ='X')
    # list data | indices | indptr
    var = rhdf5::h5read(h5path,'var/__categories/gene_symbol')    # genes
    ind = rhdf5::h5read(h5path,'var/gene_symbol')
    var = var[ind+1]
    obs = rhdf5::h5read(h5path, name ='obs/_index')    # cells 

    mat = Matrix::sparseMatrix(i = as.numeric(X$indices)+1, p = as.numeric(X$indptr) , x= as.numeric(X$data ) , dims = list(length(var),length(obs)))
    rownames(mat) =  var    
    colnames(mat) =  obs
    # outputs a gene x cell matrix
    # note that scanpy outputs a cell x gene matrix 
    return(mat)
}

#deprecated
# part 1 of annotate - deprecating - instead use the h5ad from gatherstats.py as input
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
        colnames(mat) = barcodes[,1]    # only one column containing the barcodes or run ids which identifies the cells

    } else if ( grepl( ".csv", outpath ) ) {    
        mat = as.matrix(read.delim(outpath, sep=",", stringsAsFactors=FALSE,row.names=1)  )
        mat = Matrix(mat, sparse=TRUE)

    # "SupplementData/FACS/Bladder-counts.csv"
    } else if ((any(grepl(".csv", dir(outpath) )) )) {  # otherwise, just read the csv file that STAR outputs
        # do nothing for now.
    }

    

    return(mat)

}

# filter the markers. 
get_top_markers <- function(markersetpath="~/scqc/supplement_data/markersets/MoP/class_marker_set.csv.gz"  , max_rank  = 100 ){
    meta_markers = MetaMarkers::read_meta_markers(markersetpath)
    top_markers = dplyr::filter(meta_markers, rank <=max_rank)  
    return(top_markers)
}

# Assign cells hierarchically
assign_cell_type <- function(dataset, top_markers, group_assignment = NULL ) {
    # given list of datasets as matrices, asign the cell type from the markerset - requires result from "build_marker_sets"
    # top markers should be a list of marker sets for hierarchical annotations
    # an example path - to be removed later
    # datasetpath = "SupplementData/droplet/Bladder-10X_P4_4" 

    # parse the data, score and assign cells using the top_marker

    ct_scores = MetaMarkers::score_cells(dataset, top_markers)
    ct_enrichment = MetaMarkers::compute_marker_enrichment(ct_scores)
    ct_pred = MetaMarkers::assign_cells(ct_scores , group_assignment = group_assignment ) 

    return(list(ct_pred , ct_scores, ct_enrichment))

}


annotate_execute <- function( h5path, class_ms ='class_marker_set.csv.gz',
    subclass_ms='subclass_marker_set.csv.gz', max_rank=100 ) {
        # need to figure out a good max_rank for each of the three label sets

    # marker_dir = "~/scqc/supplement_data/markersets/MoP"
    dataset = parse_h5ad(h5path)   
    # dataset should already be log-ed
    marker_sets = list(class=class_ms, subclass=subclass_ms)
    # marker_sets = list(class='class_marker_set.csv.gz' ,
    #                 subclass= 'subclass_marker_set.csv.gz')
    
    #require that the marker sets exist
    if (all(unlist(lapply(marker_sets, file.exists)))) {
        top_markers = lapply(marker_sets, get_top_markers, max_rank =max_rank)
        names(top_markers) = names(marker_sets)
        cpmcounts = MetaMarkers::convert_to_cpm(dataset) 

        class_stats = assign_cell_type(cpmcounts,top_markers$class,  group_assignment = NULL)
        class_pred = class_stats[[1]]
        class_scores = class_stats[[2]]
        class_enrichment = class_stats[[3]]

        subclass_stats = assign_cell_type(cpmcounts,top_markers$subclass,  group_assignment=class_pred$predicted)
        subclass_pred = subclass_stats[[1]]
        subclass_scores = subclass_stats[[2]]
        subclass_enrichment = subclass_stats[[3]]
        

        dfs2merge = list(class_scores_=class_scores,class_enrichment_=class_enrichment,
                subclass_scores_=subclass_scores,subclass_enrichment_=subclass_enrichment )
        i=1
        for (df in dfs2merge) {
            
            tmp = as.data.frame(t(df))
            colnames(tmp) = paste0(names(dfs2merge)[i] , colnames(tmp) )
            tmp['cell'] = rownames(tmp)
            if (i ==1){
                mergeddf = tmp
            } else{
                mergeddf = merge(mergeddf, tmp, by ='cell'  )
            }
            i = i+1 
        }

        colnames(class_pred) = paste0('class_',colnames(class_pred))
        class_pred$cell = rownames(class_pred)
        colnames(subclass_pred) = paste0('subclass_',colnames(subclass_pred))
        subclass_pred$cell = rownames(subclass_pred)
        
        preds = merge(class_pred, subclass_pred ,by='cell',sort=FALSE)
        mergeddf = merge(mergeddf, preds ,by='cell',sort=FALSE)

        # tried to save directly to h5ad file, but sc.read_h5ad doesn't recognize the type
        # save results to a temp file
        tmp_path = sub('.h5ad','.tsv',h5path)
        write.table(mergeddf, file=tmp_path,sep="\t" )
        # for (cn in colnames(preds)[2:ncol(preds)]){
        #     tryCatch({
        #         rhdf5::h5write(obj = preds[,cn] , file = h5path, name = paste0('obs/',cn) )
        #     }, error = function(e){warning(e)})
        # }
    }

}





#### driver
if (is.null(opt$mode)) {
    optparse::print_help(opt_parser)
} else if (opt$mode=="annotate"){
    # need to make sure I have a solo out directory
    if (is.null(opt$h5path)  ) {
        # no directory given, print help statement
    } else if (   file.exists(opt$h5path) ){
        annotate_execute(opt$h5path, opt$class_markerset, opt$subclass_markerset, opt$max_rank) 
    } else if ( ! file.exists(opt$h5path) ){
        # cant find directory 
        warning("no star output directory found. doing nothing")
    }

} else if (opt$mode =="setup" ){
    meta_markers= build_marker_sets_biccn(opt$rds_path , opt$markerdir )    
}




##### old functions


# not done in general
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
