# TPM-based, inference of propability vector, with replacement with variable libsize
sc_simulation = function(bulk_sample, sample, num_cells, libsize, ...) {
  s = bulk_sample %>%
    mutate(prob = tpm/sum(tpm))
  
  # NOTE: it can happen that the sampled count value is greater than the original
  map_df(1:num_cells, function(cell_num) {
    # negative values need to be excluded
    sample_size = abs(round(rnorm(1, libsize, libsize/2))) + 1 # it can happen that cells have very low read counts (~10) which yields to problems during normalization (negative size factors, actual 0)
    
    
    s %>%
      sample_n(sample_size, replace=T, weight=prob) %>%
      count(gene, name="count") %>%
      complete(gene, fill=list(count=0)) %>%
      mutate(cell = str_c(sample, 
                          "_cell_",
                          str_pad(cell_num, nchar(num_cells)+1, pad=0)))
    
  }) %>%
    spread(cell, count) %>%
    arrange(gene) %>%
    data.frame(row.names = 1, check.names = F, stringsAsFactors = F) %>%
    as.matrix() %>%
    Matrix(sparse=T)
}
construct_save_path = function(technology, type, class, num_cells, libsize, run) {
  # save path structure: number of cells, library size, simulation run
  file.path("output", "in_silico_benchmark","simulation", technology, 
            type,
            str_c(
              str_c(class,
                    type, 
                    num_cells, 
                    libsize,
                    run, sep = "_"),
              ".rds"))
}
merge_single_cells = function(sim_res) {
  sim_res %>%
    select(-sample) %>%
    nest(raw_counts) %>%
    mutate(raw_counts = data %>% map(function(data) {
      data %>% pull(raw_counts) %>% reduce(cbind)
    })) %>%
    select(-data)
}
normalize_sc_or_bulk = function(df, design_df) {
  df %>%
    left_join(design_df, by=c("id", "type", "num_cells")) %>%
    mutate(emat = future_pmap(., .f=function(raw_counts, type, design,...) {
      #message(id, " - ", type, " - ", libsize)
      if (type == "single_cells") {
        
        # filter out cells with library sizes < 100
        keep = colSums(raw_counts)[colSums(raw_counts) >= 100] %>% names()
        
        sce = SingleCellExperiment(list(counts=raw_counts[,keep, drop=F])) %>%
          computeSumFactors() %>% # computes size factor of 0 when library size is very low (e.g. 7)
          normalize() %>%
          exprs()
        
      } else if (type == "bulk_sample") {
        norm_count_matrix = raw_counts %>%
          DGEList() %>%
          calcNormFactors() %>%
          voom(design = design) %>%
          pluck("E")
      }
    }, .progress = T)) %>%
    select(-raw_counts, -design)
}
calc_contrast = function(df, design_df) {
  
  x = df %>% 
    left_join(design_df, by=c("id", "type", "num_cells")) %>%
    mutate(contrast = pmap(., .f = function(emat, design, ...) {
      # adjust design matrix as some cells might be lost due to filtering during normalization
      sub_design = design[colnames(emat),]
      stopifnot(rownames(sub_design) == colnames(emat))
      contrasts = makeContrasts(
        con = perturbed - control,
        levels = sub_design
      )
      
      lmFit(emat, sub_design) %>% 
        contrasts.fit(contrasts) %>%
        eBayes() %>%
        tidy() %>%
        select(gene, logfc = estimate)
    })) %>%
    select(-c(emat, design))
}

ss_calc_contrast = function(df, design_df, col_name) {
  
  x = df %>% 
    rename(emat = !!col_name) %>%
    left_join(design_df, by=c("id", "type", "num_cells")) %>%
    mutate(contrast = pmap(., .f = function(emat, design, ...) {
      # adjust design matrix as some cells might be lost due to filtering during normalization
      sub_design = design[colnames(emat),]
      stopifnot(rownames(sub_design) == colnames(emat))
      contrasts = makeContrasts(
        con = perturbed - control,
        levels = sub_design
      )
      
      lmFit(emat, sub_design) %>% 
        contrasts.fit(contrasts) %>%
        eBayes() %>%
        tidy() %>%
        select(gene, logfc = estimate)
    })) %>%
    select(-c(emat, design))
}

calc_viper_scores = function(df) {
  df %>%
    unnest(contrast) %>%
    nest(id, gene, logfc) %>%
    mutate(viper_result = data %>% map(run_viper, 
                                       regulon = human_dorothea, 
                                       value_name = "logfc", 
                                       id_name = "id")) %>%
    select(-data)
}
calc_progeny_scores = function(df) {
  models %>%
    select(footprints) %>%
    mutate(contrast = list(unnest(df, contrast))) %>%
    unnest(contrast) %>%
    nest(id, gene, logfc) %>%
    left_join(models, by="footprints") %>%
    mutate(progeny_result = pmap(., .f = function(data, M, ...) {
      run_progeny(data, M, value_name = "logfc", id_name = "id", permutation = 0)
    })) %>%
    select(-data, -M)
}

ss_calc_progeny_scores = function(emat,...) {
  x = models %>%
    select(footprints) %>%
    mutate(emat = list(emat)) %>%
    left_join(models, by="footprints") %>%
    mutate(ss_progeny_result = pmap(., .f = function(emat, M, ...) {
      M_matrix = M %>%
        spread(pathway, weight, fill=0) %>%
        data.frame(row.names = 1, check.names = F, stringsAsFactors = F)
      common_genes = intersect(rownames(emat), rownames(M_matrix))
      expr_matched = emat[common_genes,,drop=FALSE] %>%
        t()
      M_matched = M_matrix[common_genes,, drop=FALSE] %>%
        data.matrix()
      
      stopifnot(names(expr_matched) == rownames(M_matched))
      
      scores = expr_matched %*% M_matched
      progeny_scores = t(scale(scores))
      return(progeny_scores)
    })) %>%
    select(-emat, -M)
}


run_viper = function (E, regulon, gene_name = "gene", value_name = "expression", 
                      id_name = "sample", regulator_name = "tf", ...) {
  meta_data = E %>% select(-c(!!gene_name, !!value_name)) %>% 
    distinct()
  meta_regulon_data = regulon %>% select(-c(target, mor, likelihood)) %>% 
    distinct()
  emat = E %>% select(!!gene_name, !!id_name, !!value_name) %>% 
    spread(!!id_name, !!value_name, fill = 0) %>% drop_na() %>% 
    data.frame(row.names = 1, stringsAsFactors = F, check.names = F)
  viper_regulon = regulon %>% df2regulon(regulator_name = regulator_name)
  activity_scores = viper(eset = emat, regulon = viper_regulon, 
                          nes = T, method = "none", minsize = 4, eset.filter = F, 
                          adaptive.size = F) %>% data.frame(stringsAsFactors = F, 
                                                            check.names = F) %>% rownames_to_column(var = regulator_name) %>% 
    gather(key = !!id_name, value = "activity", -!!regulator_name) %>% 
    as_tibble() %>% inner_join(., meta_data, by = id_name) %>% 
    inner_join(., meta_regulon_data, by = regulator_name)
  return(activity_scores)
}
run_progeny = function (E, M, gene_name = "gene", value_name = "expression",
                        id_name = "sample", permutation = 10000, ...) {
  plan(multiprocess)
  E = E %>% mutate_if(is.factor, as.character)
  if (permutation > 0) {
    null_model = future_map_dfr(1:permutation, .progress = T,
                                function(p) {
                                  E %>% group_by(!!!syms(id_name)) %>% sample_frac() %>%
                                    ungroup() %>% mutate(`:=`(!!gene_name, E[[gene_name]])) %>%
                                    run_progeny(M, gene_name = gene_name, value_name = value_name,
                                                id_name = id_name, permutation = 0)
                                }) %>% group_by(!!!syms(id_name), pathway) %>% summarise(m = mean(activity),
                                                                                         s = sd(activity)) %>% ungroup()
  }
  meta_data = E %>% select(-c(!!gene_name, !!value_name)) %>%
    distinct()
  emat = E %>% select(!!gene_name, !!id_name, !!value_name) %>%
    spread(!!id_name, !!value_name, fill = 0) %>% drop_na() %>%
    data.frame(row.names = 1, stringsAsFactors = F, check.names = F)
  model = M %>% spread(pathway, weight, fill = 0) %>% data.frame(row.names = 1,
                                                                 check.names = F, stringsAsFactors = F)
  common_genes = intersect(rownames(emat), rownames(model))
  emat_matched = emat[common_genes, , drop = FALSE] %>% t()
  model_matched = model[common_genes, , drop = FALSE] %>%
    data.matrix()
  stopifnot(names(emat_matched) == rownames(model_matched))
  progeny_scores = emat_matched %*% model_matched %>% data.frame(stringsAsFactors = F,
                                                                 check.names = F) %>% rownames_to_column(id_name) %>%
    gather(key = pathway, value = activity, -!!id_name) %>%
    as_tibble() %>% inner_join(meta_data, by = id_name) %>%
    group_by(pathway) %>%
    mutate(activity = scale(activity)) %>%
    ungroup()
  if (permutation > 0) {
    progeny_z_scores = progeny_scores %>% inner_join(null_model,
                                                     by = c(id_name, "pathway")) %>% mutate(activity = (activity -
                                                                                                          m)/s) %>% select(!!id_name, pathway, activity)
    return(progeny_z_scores)
  }
  else {
    return(progeny_scores)
  }
}

sc_progeny = function(expr, M, ...) {
  common_genes = intersect(rownames(expr), rownames(M))
  expr_matched = expr[common_genes,,drop=FALSE] %>%
    t()
  M_matched = M[common_genes,, drop=FALSE] %>%
    data.matrix()
  
  stopifnot(names(expr_matched) == rownames(M_matched))
  
  scores = expr_matched %*% M_matched
  res = t(scale(scores))
  return(res)
}

df2regulon <- function(df, regulator_name="tf") {
  regulon = df %>% split(.[regulator_name]) %>% map(function(dat) {
    targets = setNames(dat$mor, dat$target)
    likelihood = dat$likelihood
    list(tfmode = targets, likelihood = likelihood)
  })
  return(regulon)
}
