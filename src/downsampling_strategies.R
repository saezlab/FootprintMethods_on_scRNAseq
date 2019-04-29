generate_single_cells = function(bulk_sample, capture_prob, num_cells, ...) {
  sample = bulk_sample %>%
    filter(tpm != 0) %>%
    mutate(index = 1:n(), prob = tpm/sum(tpm), cs = cumsum(tpm/sum(tpm)))
  
  kkk = map_df(1:num_cells, function(cell_num) {
    contruct_single_cell_vec(sample, capture_prob) %>%
      mutate(cell = as.character(cell_num))
  })
}

generate_single_cells2 = function(bulk_sample, capture_prob, num_cells, ...) {
  sample = bulk_sample %>%
    transmute(gene, prob = tpm/sum(tpm))
  
  sample_size = capture_prob * 1e+6
  
  map_df(1:num_cells, function(cell_num) {
    sample %>%
      sample_n(sample_size, replace=T, weight=prob) %>%
      count(gene, name = "tpm") %>%
      complete(gene, fill=list(tpm = 0)) %>%
      mutate(cell = as.character(cell_num))
  })
}

generate_single_cells3 = function(bulk_sample, sample, num_cells, ...) {
  s = bulk_sample %>%
    mutate(prob = count/sum(count))
  
  libsize = sum(s$count)
  
  # it can happen that the sampled count value is greater than the original
  map_df(1:num_cells, function(cell_num) {
    sample_size = 1e10
    while (sample_size > libsize) {
      # negative values need to be excluded
      sample_size = abs(round(rnorm(1, 2e4, 1e4)))
    }
    s %>%
      sample_n(sample_size, replace=T, weight=prob) %>%
      count(gene, name="count") %>%
      complete(gene, fill=list(count=0)) %>%
      mutate(cell = str_c(sample, 
                          "_cell_",
                          str_pad(cell_num, nchar(num_cells)+1, pad=0)))
      
  }) 
}

generate_single_cells4 = function(bulk_sample, sample, num_cells, ...) {
  s = bulk_sample %>%
    mutate(prob = count/sum(count))
  
  libsize = sum(s$count)
  
  # it can happen that the sampled count value is greater than the original
  map_df(1:num_cells, function(cell_num) {
    sample_size = 1e10
    while (sample_size > libsize) {
      # negative values need to be excluded
      sample_size = abs(round(rnorm(1, 2e4, 1e4)))
    }
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

contruct_single_cell_vec = function(sample, capture_prob, ...) {
  num_of_non_zero_genes = round(capture_prob * nlevels(sample$gene))
  ix = unique(findInterval(runif(num_of_non_zero_genes, 0, 1), sample$cs) + 1)
  
  while (length(ix) != num_of_non_zero_genes) {
    
    missing = num_of_non_zero_genes - length(ix)
    new_ix = unique(findInterval(runif(missing, 0, 1), sample$cs) + 1)
    
    ix = unique(c(ix, new_ix))
  }
  
  single_cell = sample %>%
    transmute(gene, tpm = case_when(index %in% ix ~ tpm,
                                    TRUE ~ 0)) %>%
    complete(gene, fill=list(tpm = 0)) %>%
    mutate(scaling_factor = capture_prob * 1000000 / sum(tpm)) %>%
    mutate(tpm = case_when(scaling_factor > 0 ~ tpm * scaling_factor,
                           scaling_factor < 0 ~ tpm / scaling_factor,
                           scaling_factor == 0 ~ tpm)) %>%
    select(gene, tpm) %>%
    mutate(gene = as.character(gene))
  return(single_cell)
}

run_ss_progeny = function(df, M, type, ...) {
  if (type == "bulk_sample") {
    expr = df %>% 
      mutate(id = "bulk") %>%
      mutate(tpm = log2(tpm+1)) %>%
      spread(id, tpm, fill=0) %>%
      data.frame(row.names=1, check.names = F, stringsAsFactors = F)
  } else if (type == "single_cells") {
    expr = df %>%
      select(gene, tpm, cell) %>%
      mutate(tpm = log2(tpm+1)) %>%
      spread(cell, tpm, fill=0) %>%
      data.frame(row.names=1, check.names = F, stringsAsFactors = F)
  }
  
  common_genes = intersect(rownames(expr), rownames(M))
  expr_matched = expr[common_genes,,drop=FALSE] %>%
    t()
  M_matched = M[common_genes,, drop=FALSE] %>%
    data.matrix()
  
  stopifnot(names(expr_matched) == rownames(M_matched))
  
  scores = expr_matched %*% M_matched %>%
    data.frame(stringsAsFactors = F, check.names = F) %>%
    rownames_to_column("id") %>%
    gather(key = pathway, value="activity", -id) %>%
    as_tibble() 
  
  return(scores)
}