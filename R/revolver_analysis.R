###########################################################
##                                                       ##
## Package: revolver                                     ##
## Version: 0.2.0                                        ##
## Date: 2013-10-2                                       ##
## Title: Repeated Evolution in Cancer (REVOLVER-updates)##
##                                                       ##
###########################################################
library(here)
library(revolver)
library(docopt)
library(dplyr)

'REVOLVER Analysis Script - Min Zhang, updated by Jan Rogel

Usage:
  revolver_analysis.R [options]

Options:
  --dataset=<FILE_PATH>     This is a required argument. Specify the path to the input file in .txt format - REQUIRED ARG
  --results=<FILE_PATH>     Speficy the outfolder to save the results in. - REQUIRED ARG
  --sspace=<VALUE>          Set the value for sspace.cutoff argument [default: 1000]
  --n_sampling=<VALUE>      [default: 500]
  --store_max=<VALUE>       [default: 200]
  --restarts=<VALUE>        [default: 10]
  --min_group_size=<VALUE>  [default: 3]
  --hc_method=<VALUE>       [default: ward]
  --split_method=<VALUE>    [default: cutreeHybrid]
  --tree_type               [defult:  clone_trees]
  -c --cohort               If set, looks for a constructed cohort in the ./ or results directory.
  -t --trees                If set, looks for computed trees in the ./ or results directory.
  -f --fit                  If set, looks for computed AND fitted trees in the ./ or results directory.
  -l --cluster              If set, looks for computed cohort, computed and fitted trees AND computed clusters. 
  -j --jackknife            If set, runs a quick non-intensive jackknife cross-validation.
  --version=0.1     Show version.

' -> doc

#---------------- Functions -----------
create_cohort <- function(file_path, results.dir = NULL, 
                          cohort_name = 'MEDUSA', 
                          load_cohort = F){
  if(load_cohort && missing(file_path)){
    return(readRDS(here(".meso_cohort.RDS"))) # Reads in from the root directory
  } else{
    
    dataframe <- data.table::fread(here(file_path), sep = "\t")
    
    dataframe$Misc    <- as.character(dataframe$Misc)
    dataframe$patientID <- as.character(dataframe$patientID)
    dataframe$variantID <- as.character(dataframe$variantID)
    dataframe$cluster <- as.character(dataframe$cluster)
    dataframe$is.driver <- as.logical(dataframe$is.driver)
    dataframe$is.clonal <- as.logical(dataframe$is.clonal)
    dataframe$CCF       <- as.character(dataframe$CCF)
    
    do.plots       <- TRUE
    cex.fit        <- 3
    
    ## construct revolver object
    meso.cohort <- revolver_cohort(
      dataset = dataframe,
      CCF_parser = revolver:::CCF_parser,
      annotation = cohort_name,
      ONLY.DRIVER = FALSE,
      MIN.CLUSTER.SIZE = 1
    )
    
    revolver_check_cohort(meso.cohort, stopOnError = FALSE)
    non.recurrent <- Stats_drivers(meso.cohort) %>%
      filter(N_tot == 1) %>%
      pull(variantID)
    
    meso.cohort <- remove_drivers(meso.cohort, non.recurrent)
    if(!missing(results.dir)){
      saveRDS(meso.cohort, here(results.dir, ".meso_cohort.RDS"))
    } else{
      saveRDS(meso.cohort, here(".meso_cohort.RDS"))
    }
    
    return(meso.cohort)
  }
  
}
create_directory <- function(results.path){
  dir.list <- c("Phylogenies", "Clustering", "Jackknife")
  
  out.dirs <- list()
  out.dirs$results <- results.path
  for(name in dir.list){
    dir.create(here(results.path, name), showWarnings = F, recursive = T)
    out.dirs[[name]] <- here(results.path, name)
    
  }
  
  return(out.dirs)
}

main <- function(meso.cohort = NULL, start.cohort = F, start.trees = F, start.fit = F, start.cluster = F, ...){ #... must match file_path, options.trees, options.fit, options.clustering, dirs
  arguments <- list(...)
  print(arguments)
  
  if(start.cohort){
    ## create meso cohort
    print(paste0("Creating a REVOLVER cohort from ", arguments$file_path))
    meso.cohort <- create_cohort(file_path = arguments$file_path, load_cohort = F)

    ## infer phylogeny tree of subclonal/clonal
    print("Computing mutation trees... ")
    meso.cohort <- compute_mutation_trees(meso.cohort,
                                          patients     = meso.cohort$patients,
                                          sspace.cutof = arguments$options.trees$sspace.cutof,
                                          n.sampling   = arguments$options.trees$n.sampling,
                                          store.max    = arguments$options.trees$store.max,
                                          overwrite    = F)
    saveRDS(meso.cohort, here(dirs$Phylogenies, ".meso_computed_mutation_trees.RDS"))
    
    ## fit repeat evolution tree
    meso.cohort <- revolver_fit(meso.cohort,
                                parallel = F,
                                n = 3,
                                initial.solution = arguments$options.fit$initial.solution)
    saveRDS(meso.cohort, here(dirs$Phylogenies,".meso_fitted_mutation_trees.RDS"))
    
    ## cluster the tree
    meso.cohort <- revolver_cluster(meso.cohort,
                                    split.method = arguments$options.clustering$split.method,
                                    hc.method = arguments$options.clustering$hc.method,
                                    min.group.size = arguments$options.clustering$min.group.size)
    saveRDS(meso.cohort, here(dirs$Clustering,".meso_clusters.RDS"))
    
    return(meso.cohort)
  } else 
    if(start.trees){
      ## infer phylogeny tree of subclonal/clonal
      meso.cohort <- compute_mutation_trees(meso.cohort,
                                            patients     = meso.cohort$patients,
                                            sspace.cutof = arguments$options.trees$sspace.cutof,
                                            n.sampling   = arguments$options.trees$n.sampling,
                                            store.max    = arguments$options.trees$store.max,
                                            overwrite    = F)
      saveRDS(meso.cohort, here(dirs$Phylogenies,".meso_computed_mutation_trees.RDS"))
      
      ## fit repeat evolution tree
      meso.fit <- revolver_fit(meso.cohort,
                               parallel = F,
                               n = 3,
                               initial.solution = arguments$options.fit$initial.solution)
      saveRDS(meso.cohort, here(dirs$Phylogenies,".meso_fitted_mutation_trees.RDS"))
      
      ## cluster the tree
      meso.cohort <- revolver_cluster(meso.cohort,
                                      split.method = arguments$options.clustering$split.method,
                                      hc.method = arguments$options.clustering$hc.method,
                                      min.group.size = arguments$options.clustering$min.group.size)
      saveRDS(meso.cohort, here(dirs$Clustering,".meso_clusters.RDS"))
      
      return(meso.cohort)
    } else 
      if(start.fit){
        ## fit repeat evolution tree
        meso.fit <- revolver_fit(meso.cohort,
                                 parallel = F,
                                 n = 3,
                                 initial.solution = arguments$options.fit$initial.solution)
        saveRDS(meso.cohort, here(dirs$Phylogenies,".meso_fitted_mutation_trees.RDS"))
        

        ## cluster the tree
        meso.cohort <- revolver_cluster(meso.cohort,
                                        split.method = arguments$options.clustering$split.method,
                                        hc.method = arguments$options.clustering$hc.method,
                                        min.group.size = arguments$options.clustering$min.group.size)
        saveRDS(meso.cohort, here(dirs$Clustering,".meso_clusters.RDS"))
        
        return(meso.cohort)
      } else 
        if(start.cluster){
          ## cluster the tree
          meso.cohort <- revolver_cluster(meso.cohort,
                                          split.method = arguments$options.clustering$split.method,
                                          hc.method = arguments$options.clustering$hc.method,
                                          min.group.size = arguments$options.clustering$min.group.size)
          saveRDS(meso.cohort, here(dirs$Clustering,".meso_clusters.RDS"))
          
          return(meso.cohort)
        }
}

#---------------- Arguments -----------
arguments <- docopt(doc) 
print(arguments)

results   <- arguments$results
file_path <- arguments$dataset

#---------------- Options -----------
options(crayon.enabled = T)
options(revolver.progressBar = T)
options(revolver.jamPDF = TRUE)
options.trees                <- list(sspace.cutof = as.integer(arguments$sspace), 
                                     n.sampling = as.integer(arguments$n_sampling), 
                                     store.max = as.integer(arguments$store_max), 
                                     overwrite = FALSE)
options.fit                  <- list(initial.solution = NA, 
                                     transitive.orderings = FALSE, 
                                     restarts = as.integer(arguments$restarts))
if(arguments$split_method == 'cutreeHybrid'){
  options.clustering    <- list(transitive.closure = FALSE, 
                                min.group.size = as.integer(arguments$min_group_size), 
                                hc.method = arguments$hc_method, 
                                cutoff.features_annotation = 1, 
                                split.method = arguments$split_method)
} else if(arguments$split_method == 'cutreeDynamic'){
  options.clustering   <- list(transitive.closure = FALSE, 
                               min.group.size = arguments$min_group_size, 
                               hc.method = arguments$hc_method, 
                               cutoff.features_annotation = 1, 
                               split.method = arguments$split_method)
}


# print("Setting the analysis to default options...")
# options.trees <- list(sspace.cutoff = 1000, n.sampling = 500, store.max = 200, overwrite = FALSE)
# options.fit <- list(initial.solution = NA, transitive.orderings = FALSE, restarts = 10)
# options.clustering <- list(transitive.closure = FALSE, min.group.size = 3, hc.method = 'ward', cutoff.features_annotation = 1, split.method = 'cutreeHybrid')
# options.clustering.withoutGL <- list(transitive.closure = FALSE, min.group.size = 3, hc.method = 'ward', cutoff.features_annotation = 1,split.method = 'cutreeDynamic')
# dataframe <- data.table::fread(here('data', 'medusa50.revolver.input.csv'), sep = "\t")


#------------ Directory Creation -----------
if(length(results) == 0){
  print("No results folder set, saving to default folder...")
  folder.output  <- "./"
} else if(length(results) == 1){
  print(paste0("Saving results to ", results))
  dirs <- create_directory(results)
}

#----------- Write Analysis Parameters -----------
writexl::write_xlsx(list("tree_parameters" = as.data.frame(options.trees), 
                         "model_fit_parameters" = as.data.frame(options.fit), 
                         "clustering_parameters" = as.data.frame(options.clustering)),
                    here(dirs$results, "training_parameters.xlsx"))

#------------ Script Main -----------
if(arguments$cohort){
  print("Loading from a constructed cohort...")
  meso.cohort <- create_cohort(load_cohort = T) 
  meso.cohort <- main(meso.cohort = meso.cohort, 
                         start.cohort = F, start.trees = T, start.fit = F, start.cluster = F,
                         options.fit = options.fit,
                         options.trees = options.trees, 
                         options.clustering = options.clustering,
                         file_path = file_path, dirs = dirs)

} else if(arguments$trees){
  # meso.cohort <- readRDS(here(dirs$results, '.meso_computed_clone_trees.RDS'))
  meso.cohort <- readRDS(here('.meso_computed_mutation_trees.RDS')) 
  meso.cohort <- main(meso.cohort = meso.cohort,
                         start.cohort = F, start.trees = F, start.fit = T, start.cluster = F,
                         options.fit = options.fit, 
                         options.clustering = options.clustering,
                         file_path = file_path, dirs = dirs)
  
} else if(arguments$fit){
  meso.cohort <- readRDS(here('.meso_fitted_mutation_trees.RDS'))
  meso.cohort <- main(meso.cohort = meso.cohort,
                         start.cohort = F, start.trees = F, start.fit = F, start.cluster = T,
                         options.clustering = options.clustering,
                         file_path = file_path, dirs = dirs)
  
} else if(arguments$cluster){
  meso.cohort <- readRDS(here(".meso_clusters.RDS"))
  
} else{
  print("Running whole analysis...")
  meso.cohort <- main(start.cohort = T, start.trees = F, start.fit = F, start.cluster = F, 
                         options.fit = options.fit,
                         options.trees = options.trees, 
                         options.clustering = options.clustering,
                         file_path = file_path, dirs = dirs)
}


pdf(file = here(dirs$result, "outplots.pdf"))
plot_DET_index(meso.cohort)
plot_clusters(meso.cohort)
plot_drivers_clonality(meso.cohort)
plot_drivers_graph(meso.cohort)
plot_drivers_occurrence(meso.cohort)
plot_penalty(meso.cohort)
plot_trajectories_per_cluster(meso.cohort, min_counts = 2)
plot_trajectories_per_cluster(meso.cohort, min_counts = 3)
plot_trajectories_per_cluster(meso.cohort, min_counts = 4)
plot_trajectories_per_cluster(meso.cohort, min_counts = 5)
dev.off()

if(arguments$jackknife){
  # --- cross-validate using jackknife approach
  meso.jackknife <- revolver_jackknife(meso.cohort)
  saveRDS(meso.jackknife, here(dirs$Jackknife, ".meso_jackknife.RDS"))
  
  # --- plot results for phylogenetic analysis, jacknife metrics
  pdf(file = here(dirs$Jackknife, "outplots.pdf"))
  plot_jackknife_cluster_stability(meso.jackknife)
  plot_jackknife_coclustering(meso.jackknife)
  plot_jackknife_trajectories_stability(meso.jackknife)
  dev.off()
}

