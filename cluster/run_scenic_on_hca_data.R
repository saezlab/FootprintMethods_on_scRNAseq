library(tidyverse)
library(SCENIC)

# utility functions
update_path = function(mat = NULL, technology = NULL) {
  mat %>% 
    data.frame() %>%
    rownames_to_column("id") %>%
    mutate(fileName = file.path("output/hca_data/scenic", technology, fileName)) %>%
    column_to_rownames("id") %>%
    as.matrix()
}

initialize_scenic_wrapper = function(technology, n_cores = 4,...) {
  scenic_options = initializeScenic(
    org="hgnc",
    datasetTitle = as.character(technology),
    dbDir="data/scenic",
    dbs = c("hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather",
            "hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather"),
    nCores=n_cores)
  
  if (!dir.exists(file.path("output/hca_data/scenic", technology))) {
    dir.create(file.path("output/hca_data/scenic", technology))
    dir.create(file.path("output/hca_data/scenic", technology, "int"))
    dir.create(file.path("output/hca_data/scenic", technology, "output"))
  }
  
  scenic_options@fileNames$int = update_path(scenic_options@fileNames$int, technology)
  scenic_options@fileNames$output = update_path(scenic_options@fileNames$output, technology)
  
  return(scenic_options)
}

run_scenic = function(expr, scenic_options, ...) {
  
  message(getDatasetInfo(scenic_options, "datasetTitle"))
  
  expr = expr %>% as.matrix()
  genesKept <- geneFiltering(expr, scenic_options)
  exprMat_filtered <- expr[genesKept, ]
  runCorrelation(exprMat_filtered, scenic_options)
  exprMat_filtered_log <- log2(exprMat_filtered+1)
  runGenie3(exprMat_filtered_log, scenic_options)

  ### Build and score the GRN
  exprMat_log <- log2(exprMat+1)
  runSCENIC_1_coexNetwork2modules(scenic_options)
  runSCENIC_2_createRegulons(scenic_options)
  
  return(NULL)
}

# load required data
# make sure that the working directory is set to root of project
# meta = readRDS("output/hca_data/meta/meta_df.rds")
expr_all = readRDS("output/hca_data/expression/raw_preprocessed.rds") %>%
  rename(expr = raw)

# set number of cores
n_cores = 4

# expr = expr_all %>% pluck(2,1)
# technology = expr_all %>% pluck(1,1)
# scenic_options = res %>% pluck(4,1)

y = expr_all %>%
  mutate(n_cores = n_cores) %>%
  transmute(technology, expr, scenic_options = pmap(., .f = initialize_scenic_wrapper)) %>%
  slice(1) %>%
  mutate(tmp = pmap(., .f = run_scenic))



