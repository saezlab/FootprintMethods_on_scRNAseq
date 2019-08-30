#' This functions runs DoRothEA (either version1 or version2) for a set of 
#' experiments. If required the quality of the experiments is checked and 
#' experiments which do not pass are excluded @seealso [quality_check()] 
#' 
#' @param df tidy data frame containing zscore normalized gene expression values
#'   for a given set of experiments and their corresponding meta information
#'   (source, organism, accession etc...).
#' @param qc logical indicating if a quality check should be performed 
#' @param missing_value chose between 0 and NA to replace missing values
#' @param regulon used regulon
#' @param v either "v1" or "v2" indication version1 or version2
#' @param ... other arguments passed to quality_check @seealso [quality_check()]
#' @return tidy data frame of DoRothEA results for a given set of experiments
#'   and their corresponding meta information 
run_dorothea = function(df, qc=F, qc_th = 0, missing_value=0, regulon,...) {

  if (qc == TRUE) {
    keep_ids = quality_check(df, qc_th = qc_th)
  } else {
    keep_ids = df %>%
      distinct(id) %>%
      pull()
  }
  
  if (length(keep_ids) == 0) {
    return(as_tibble(NULL))
  }
  
  meta_df = df %>%
    filter(id %in% keep_ids) %>%
    select(one_of("id", "accession", "tf", "platform", "info", "effect", 
                  "source","sign", "disease", "disease_name", "do_id", "run")) %>%
    distinct()
  
  meta_regulon = regulon %>%
    select(-c(target, mor, likelihood)) %>%
    distinct() %>%
    rename(d_tf = tf)
  
  expr = df %>%
    filter(id %in% keep_ids) %>%
    select(gene, id, expression) %>%
    spread(id, expression, fill=missing_value) %>% 
    na.omit() %>%
    data.frame(row.names=1, stringsAsFactors = F, check.names = F)
  
  viper_regulon = regulon %>%
    df2regulon()
    
  res = viper(eset = expr, regulon = viper_regulon, nes = T, method = 'none',
              minsize = 4, eset.filter = F) %>%
    data.frame(stringsAsFactors = F, check.names = F) %>%
    mutate(d_tf = rownames(.)) %>%
    gather(key=id, value="d_score", -d_tf) %>%
    as_tibble() %>%
    inner_join(., meta_df, by="id") %>%
    inner_join(., meta_regulon, by="d_tf")
    
  
  
}


#' This function checks the quality of a set of experiments. For each experiment
#' the zscore normalized gene expression value of the overexpressed/knockdown 
#' target gene should be above/below a specified threshold (qc_th, default=0)
#' 
#' @param df tidy data frame containing zscore normalized gene expression values
#'   for a given set of experiments and their corresponding meta information
#'   (source, organism, accession etc...).
#' @param qc_th float specifying the threshold. Experiments where the 
#'   overexpressed gene is below this threshold and knockdowned genes above this
#'   threshold will be discarded
#' @return vector of experiment ids which passed the quality check
quality_check = function(df, qc_th = 0) {
  keep = df %>%
    group_by(id) %>%
    filter(gene == tf) %>%
    arrange(effect, expression) %>%
    mutate(qc = case_when(effect == "knockdown" & expression <= qc_th ~ TRUE,
                          effect == "overexpression" & expression * -1 <= qc_th ~ TRUE,
                          TRUE ~ FALSE)) %>%
    filter(qc == T) %>%
    pull(id)
}


filter_common_tfs = function(df) {
  tmp = df %>%
    select(-activity)
  
  common_tfs = df %>%
    unnest(activity) %>%
    distinct(organism, name, confidence, d_tf) %>%
    group_by(organism, confidence, d_tf) %>%
    count() %>%
    filter(n == 2) %>%
    ungroup() %>%
    select(-n)
  
  df %>% unnest(activity) %>%
    inner_join(common_tfs, by = c("confidence", "organism", "d_tf")) %>%
    nest(-c(confidence, organism, name), .key="activity") %>%
    inner_join(tmp, by = c("confidence", "organism", "name"))
}
#' Title
#'
#' @param df
#'
#' @return
#' @export
#'
#' @examples
df2regulon = function(df) {
  regulon = df %>%
    split(.$tf) %>%
    map(function(dat) {
      tf = dat %>% distinct(tf) %>% pull()
      targets = setNames(dat$mor, dat$target)
      likelihood = dat$likelihood
      list(tfmode =targets, likelihood = likelihood)
    })
  return(regulon)
}


#' Title
#'
#' @param r
#'
#' @return
#' @export
#'
#' @examples
regulon2df = function(r) {
  res = r %>%
    map_df(.f = function(i) {
      tf_target = i$tfmode %>%
        enframe(name = "target", value="mor") %>%
        mutate(likelihood = i$likelihood)
    },
    .id = "tf")
  return(res)
}

#' @param df run_dorothea() output @seealso [run_dorothea()]
#' @param filter_tn logical flag indicating if unnescesary true negatives should
#'   be filtered out (unnescssary means that there are no true positives)
#' @return tidy data frame with meta information for each experiment and the
#'   response and the predictor value which are required for roc curve analysis
prepare_dorothea_for_roc = function(df, filter_tn = F) {
  res = df %>%
    mutate(response = case_when(d_tf == tf ~ 1, 
                                d_tf != tf ~ 0),
           predictor = d_score * sign)
  
  if (filter_tn == TRUE) {
    # Only TF which are perturbed and predicted are considered
    z = intersect(res$d_tf, res$tf)
    res = res %>%
      filter(d_tf %in% z, tf %in% z)
  } 
  res %>%
    select(-c(tf, d_score, accession, platform, info, effect, sign)) %>%
    rename(tf = d_tf)
}