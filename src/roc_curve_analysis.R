#' This function is to scale feature(pathway/TF) scores before applying them to ROC curve
#' analyis. Before scaling the data frame is grouped by the organism. It can be
#' scaled only feature or sample wise or both. Please note that if you want to 
#' scale both feature and sample wise, feature wise scaling is done first.
#' 
#' @param df tidy data frame with response, predictor values and meta 
#'   information for each experiment. The to be scaled values should be in the
#'   "predictor" column
#' @param feature name of the feature (mostly either "pathway" or "tf")
#' @param scale_feat logical flag indicating if the scaling should be feature wise
#' @param scale_sp logical flag indication in the scaling should be sample wise  
#' @return same data frame as input but with scaled predictor values
scale_this = function(df, feature="pathway", scale_feat=T, scale_sp=T) {
  if (scale_feat == TRUE & scale_sp == FALSE) {
    message("Scaling only feature wise...")
    dat = df %>%
      group_by_(feature) %>%
      mutate(predictor = scale(predictor))
  } else if(scale_feat == FALSE & scale_sp == TRUE) {
    message("Scaling only sample wise...")
    dat = df %>%
      group_by(id) %>%
      mutate(predictor = scale(predictor))
  } else if (scale_feat == TRUE & scale_sp == TRUE) {
    message("Scaling feature and sample wise...")
    dat = df %>%
      group_by_(feature) %>%
      mutate(predictor = scale(predictor)) %>%
      ungroup() %>%
      group_by(id) %>%
      mutate(predictor = scale(predictor)) 
  } else if (scale_feat == FALSE & scale_sp == FALSE) {
      message("The function scale_this() does nothing...")
      dat = df
  }
  res = dat %>% 
    ungroup() %>%
    mutate(scale_feat = as.logical(scale_feat),
           scale_sp = as.logical(scale_sp))
}


#' Performs ROC curve analysis and calculates confidence intervalls if required
#' 
#' @param df tidy data frame containing "response" and "predictor" as column
#'   names. "response" contains the true values (0 or 1) and "predict" contains
#'   the predicted values.
#' @param calc_ci flag indicating in confidence intervalls should be calculated
#'   for auc values
#' @return tidy data frame with false positive rate, true positive rate, area
#'   under the curve and the corresponding meta information. If confidence
#'   intervalls have been computed than also the confidence intervall, the lower
#'   bound, upper bound and the significant level is contained
calc_roc_auc = function(df, calc_ci=FALSE) {
  if (sum(df$response) != 0) {
    r = df %>%
      roc(response = .$response, predictor = .$predictor, direction = "<")
    
    dat = tibble(tpr = r$sensitivities,
                 fpr = 1-r$specificities,
                 auc = r$auc,
                 n = sum(r$original.response == 1),
                 roc_obj = list(r)) %>%
      arrange(fpr, tpr)
    
    if (calc_ci == TRUE) {
      res = calc_ci_auc(r)
      dat = dat %>%
        mutate(sig_level = res$sig_level,
               ci = res$ci)
    } else {
      return(dat)
    }
  } else {
    return(as_tibble(NULL))
  }
}

#' This function calculates confidence intervalls for area under the curve (auc)
#' values for different intervalls (90, 95, 99, 99.9 and 99.99 %)
#' 
#' @param r output (object) from pROC::roc() function
#' @return tibble with confidence intervall, upper bound, lower bound  and
#'   significance level (symbolic)
calc_ci_auc = function(r) {
  sig_levels = tibble(pval = c(0.1, 0.05, 0.01, 0.001),
                      symbol = c("+", "*", "**", "***"))
  
  sig_test = map_df(c(0.9, 0.95, 0.99, 0.999), function(conf_level) {
    x = ci(r, of="auc", conf.level=conf_level)
    list(ci = conf_level, lb = x[1], ub = x[3])
  })
  
  if (r$auc >= 0.5) {
    sig_symbol = sig_test %>%
      filter(lb >= 0.5) %>%
      slice(which.min(lb)) %>%
      mutate(
        sig_level = with(sig_levels,
                         as.character(symbol)[match(round(1-ci,8), pval)])
        ) %>%
      select(sig_level) %>%
      mutate(ci = list(sig_test))
  } else if (r$auc < 0.5) {
    sig_symbol = sig_test %>%
      filter(ub < 0.5) %>%
      slice(which.max(ub)) %>%
      mutate(
        sig_level = with(sig_levels, 
                         as.character(symbol)[match(round(1-ci,8), pval)])
        ) %>%
      select(sig_level) %>%
      mutate(ci = list(sig_test))
  }
  if (nrow(sig_symbol) == 0) {
    sig_symbol = tibble(sig_level = "", ci = list(sig_test))
  } else {
    return(sig_symbol)
  }
}

#' This function tests if auc values from two different ROC curves differ
#' significantly
#' 
#' @param df
#' @return 
compare_roc_auc = function(df, group = "organism") {
  # test if a ROC curve analysis has been performed for the specific group
  if (nrow(df) > 0) {
    roc_list = df %>%
      select_(group, "roc_obj") %>%
      group_by_(group) %>%
      slice(1) %>%
      pull(roc_obj)
    
    if (length(roc_list) == 2) {
      res = roc.test(roc_list[[1]], roc_list[[2]]) %>%
        tidy() %>%
        as_tibble()
      
    } else {
      res = tibble(estimate1 = NA,
                   estimate2 = NA,
                   statistic = NA,
                   p.value = NA,
                   parameter = NA,
                   method = NA,
                   alternative = NA)
    }
  } else {
    res = tibble(estimate1 = NA,
                 estimate2 = NA,                                     
                 statistic = NA,
                 p.value = NA,
                 parameter = NA,
                 method = NA,
                 alternative = NA)
  }
}


#' This function corrects p values for multiple hypothesis testing
#' 
#' @param df
#' @return
correct_p_values = function(df) {
  g = group_vars(df)
  
  nested_roc = df %>% 
    select_(., .dots=c(g, "roc"))
  
  if (length(g) == 2) {
    res = df %>%
      ungroup() %>%
      unnest(stats) %>%
      select(-roc, -organism) %>%
      distinct()  %>%
      mutate(adj_p_val = p.adjust(p.value, method="fdr")) %>%
      mutate(symbol = case_when(
        adj_p_val <= 0.1 & adj_p_val > 0.05 ~ "+",
        adj_p_val <= 0.05 & adj_p_val > 0.01 ~ "*",
        adj_p_val <= 0.01 & adj_p_val > 0.001 ~ "**",
        adj_p_val <= 0.001 & adj_p_val > 0.0001 ~ "***", 
        adj_p_val <= 0.0001 & adj_p_val > 0.00001 ~ "****",
        adj_p_val < 0.00001 ~ "*****",
        TRUE ~ "n.s.")) %>%
      nest(estimate1, estimate2, statistic, p.value, parameter, method, 
           alternative, adj_p_val, symbol, .key="stats") %>%
      full_join(., nested_roc) %>%
      group_by_(.dots=g)
    
  } else if (length(g) == 1) {
    res = df %>%
      ungroup() %>%
      unnest(stats) %>%
      select(-roc) %>%
      distinct()  %>%
      mutate(adj_p_val = p.adjust(unique(p.value), method="fdr")) %>%
      mutate(symbol = case_when(
        adj_p_val <= 0.1 & adj_p_val > 0.05 ~ "+",
        adj_p_val <= 0.05 & adj_p_val > 0.01 ~ "*",
        adj_p_val <= 0.01 & adj_p_val > 0.001 ~ "**",
        adj_p_val <= 0.001 & adj_p_val > 0.0001 ~ "***", 
        adj_p_val <= 0.0001 & adj_p_val > 0.00001 ~ "****", 
        TRUE ~ "n.s.")) %>%
      nest(estimate1, estimate2, statistic, p.value, parameter, method, 
           alternative, adj_p_val, symbol, .key="stats") %>%
      full_join(., nested_roc) %>%
      group_by_(.dots=g)
    
  } else {
    stop(paste(
      "Only pathway/tf and organism are allowed to be grouping variables.",
      "Your groups are currently: ", g))
  }
}

#' This function compares all ROC curves with each other. P values are adjusted
#' using false discovery rate
#' 
#' @param df tidy data frame with a setup column containing different setups 
#'   (various regulons/organisms/methods)
#' @return 1x1 nested data frame containing tidy output from roc.test()
compare_all_auc = function(df) {
  stats = combn(df$setup, 2, FUN=function(g) {
    roc1 = filter(df, setup == g[1]) %>%
      unnest(roc) %>%
      pull(roc_obj) %>%
      unique()
    
    roc2 = filter(df, setup == g[2]) %>%
      unnest(roc) %>%
      pull(roc_obj) %>%
      unique()
    
    test_res = roc.test(roc1[[1]], roc2[[1]], paired=F) %>%
      tidy() %>% 
      as_tibble() %>%
      mutate(group1 = g[1], group2 = g[2])
    
  }, simplify = F) %>%
    bind_rows() %>%
    mutate(adj_p_val = p.adjust(p.value, method="fdr")) %>%
    mutate(symbol = case_when(
      adj_p_val <= 0.1 & adj_p_val > 0.05 ~ "+",
      adj_p_val <= 0.05 & adj_p_val > 0.01 ~ "*",
      adj_p_val <= 0.01 & adj_p_val > 0.001 ~ "**",
      adj_p_val <= 0.001 & adj_p_val > 0.0001 ~ "***", 
      adj_p_val <= 0.0001 & adj_p_val > 0.00001 ~ "****",
      adj_p_val < 0.00001 ~ "*****",
      TRUE ~ "n.s.")) 
  return(stats)
}

#' Plots ROC curves for different groups
#' 
#' @param df output from calc_roc_auc() function  @seealso [calc_roc_auc()]
plot_roc_curve = function(df) {
  g = group_vars(df)
  
  anno_df = df %>%
    select(-c(tpr, fpr, roc_obj, ci)) %>%
    distinct() %>%
    mutate(x = 0.75, y = case_when(organism == "human" ~ 0.25,
                                   organism == "mouse" ~0.05),
           lab = if (exists("sig_level")) {
             paste0(round(auc,3), " (", n, ") ", sig_level)
           } else {
             paste0(round(auc,3), " (", n, ")")}
           )

  p = ggplot(df , aes(x=fpr, y=tpr, color=organism)) + 
    geom_line() +
    geom_abline(linetype=2) +
    coord_fixed() +
    geom_text(data = anno_df, aes(x,y,label=lab, color=organism),
              inherit.aes = F, size=2) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  if (length(g) == 1) {
    sig = df %>%
      ungroup() %>%
      unnest(stats) %>%
      select(symbol) %>%
      distinct()
    p = p +
      annotate("text", x=0.1, y=1, label=pull(sig))
  } else if (length(g) == 2) {
    sig = df %>%
      ungroup() %>%
      unnest(stats) %>%
      select(-c(roc, estimate1, estimate2, statistic, p.value, parameter,
                method, alternative, organism)) %>%
      unique() %>%
      mutate(lab = case_when(
        is.na(adj_p_val) ~ paste(pathway, "(NA)"),
        !is.na(adj_p_val) ~ paste0(pathway," (", symbol, ")")))
    
    grid_labels = setNames(sig$lab, sig$pathway)
    
    p = p + facet_wrap(as.formula(paste("~", g[2])),
                       labeller=labeller(.default = grid_labels))
  } else if (length(g) == 3) {
    p = p + facet_grid(as.formula(paste(g[2], "~", g[3])))
  } else {
    p
  }
}




calc_roc_curve_backup = function(df, downsampling=F, times = 1000, sig_perm=F, nperm = 1000) {
  # ensure that either downsampling or sig_perm are TRUE but not both an the same time
  if (downsampling == TRUE & sig_perm == TRUE) {
    stop("The args downsampling and sig_perm cannot be TRUE at the same time!")
  }
  
  tn = df %>% filter(response == 0)
  tp = df %>% filter(response == 1)
  
  feature_coverage = df %>% 
    select(one_of("pathway", "tf")) %>% 
    distinct() %>% 
    nrow()
  
  if (downsampling == T) {
    # number of true positives
    num_tp = nrow(tp)
    
    res = map_df(seq(from=1, to=times, by=1), .f= function(i) {
      df_sub = sample_n(tn, num_tp, replace=TRUE) %>% 
        bind_rows(tp)
      
      
      r_sub = roc.curve(scores.class0 = df_sub$predictor, 
                      weights.class0 = df_sub$response,
                      curve=T) 
      res_tmp = r_sub$curve %>%
        as_tibble() %>%
        setNames(., c("fpr", "tpr", "th")) %>%
        mutate(auc = r_sub$auc,
               type = r_sub$type,
               n = sum(df_sub$response),
               tp = nrow(tp),
               tn = nrow(tn),
               coverage = feature_coverage) %>%
        mutate_("run" = i)
      })
  } else if (sig_perm == T) {
    res_perm = map_df(seq(from=1, to=nperm, by=1), .f=function(j) {
      # shuffle predicted values
      shuffled_predictor = sample(df$predictor)
      
      r_perm = roc.curve(scores.class0 = shuffled_predictor, 
                         weights.class0 = df$response,
                         curve=T) 
      
      res_perm = r_perm$curve %>%
        as_tibble() %>%
        setNames(., c("fpr", "tpr", "th")) %>%
        mutate(auc = r_perm$auc,
               type = r_perm$type,
               n = sum(df$response),
               tp = nrow(tp),
               tn = nrow(tn),
               coverage = feature_coverage) %>%
        mutate_("run" = j) %>%
        mutate(stat = "permutation")
    })
    
    r_orig = roc.curve(scores.class0 = df$predictor,
                       weights.class0 = df$response,
                       curve=T)
    res_orig = r_orig$curve %>%
      as_tibble() %>%
      setNames(., c("fpr", "tpr", "th")) %>%
      mutate(auc = r_orig$auc,
             type = r_orig$type,
             n = sum(df$response),
             tp = nrow(tp),
             tn = nrow(tn),
             coverage = feature_coverage) %>%
      mutate(run = 0) %>%
      mutate(stat = "original") %>%
      arrange(fpr, tpr)
    
    res = bind_rows(res_orig, res_perm)
  } else {
    r = roc.curve(scores.class0 = df$predictor,
                  weights.class0 = df$response,
                  curve=T)
    res = r$curve %>%
      as_tibble() %>%
      setNames(., c("fpr", "tpr", "th")) %>%
      mutate(auc = r$auc,
             type = r$type,
             n = sum(df$response),
             tp = nrow(tp),
             tn = nrow(tn),
             coverage = feature_coverage) %>%
      arrange(fpr, tpr)
  }
  return(res)
}

#' This function calculates precision recall curves
#' 
#' @param df tidy data frame containing predicted values in column predictor and
#'   true values in column response
#' @return tidy data frame containing recall, precision, auc and number of TP
calc_pr_curve = function(df) {
  if (sum(df$response) == 0) {
    return(as_tibble(NULL))
  } 
  
  if (nrow(distinct(df, response)) != 2) {
    return(as_tibble(NULL))
  }
  
  feature_coverage = df %>% 
    select(one_of("pathway", "tf")) %>% 
    distinct() %>% 
    nrow()
  
  tn = df %>% filter(response == 0)
  tp = df %>% filter(response == 1)
  
  r = pr.curve(scores.class0 = df$predictor,
               weights.class0 = df$response,
               curve=T)
  res = r$curve %>%
    as_tibble() %>%
    setNames(., c("recall", "precision", "th")) %>%
    mutate(auc = r$auc.davis.goadrich,
           type = r$type,
           n = sum(df$response),
           tp = nrow(tp),
           tn = nrow(tn),
           coverage = feature_coverage,
           run = unique(df$run)) %>%
    arrange(recall, desc(precision))
}


calc_roc_curve = function(df, downsampling=F, times = 1000) {
  if (sum(df$response) == 0) {
    return(as_tibble(NULL))
  } 
  
  if (nrow(distinct(df, response)) != 2) {
    return(as_tibble(NULL))
  }
  
  tn = df %>% filter(response == 0)
  tp = df %>% filter(response == 1)
  
  feature_coverage = df %>% 
    select(one_of("pathway", "tf")) %>% 
    distinct() %>% 
    nrow()
  
  if (downsampling == T) {
    # number of true positives
    num_tp = nrow(tp)
    
    res = map_df(seq(from=1, to=times, by=1), .f= function(i) {
      df_sub = sample_n(tn, num_tp, replace=TRUE) %>% 
        bind_rows(tp)
      
      r_sub = roc(response = df_sub$response, predictor = df_sub$predictor, direction = "<")
      
      res_sub = tibble(tpr = r_sub$sensitivities,
                   fpr = 1-r_sub$specificities,
                   th = r_sub$thresholds,
                   auc = r_sub$auc,
                   n = sum(df$response),
                   tp = nrow(tp),
                   tn = nrow(tn),
                   coverage = feature_coverage) %>%
        mutate_("run" = i)

    })
  } else {
    r = roc(response = df$response, predictor = df$predictor, direction = "<")
    
    res = tibble(tpr = r$sensitivities,
                 fpr = 1-r$specificities,
                 th = r$thresholds,
                 auc = r$auc,
                 n = sum(df$response),
                 tp = nrow(tp),
                 tn = nrow(tn),
                 coverage = feature_coverage,
                 run = unique(df$run)) %>%
      arrange(fpr, tpr)
    
    ci = calc_ci_auc(r)
    res = res %>%
      mutate(sig_level = ci$sig_level,
             ci = ci$ci)
  }
  return(res)
}

get_roc_object = function(df) {
  if (sum(df$response) == 0) {
    return(as_tibble(NULL))
  } 
  
  if (nrow(distinct(df, response)) != 2) {
    return(as_tibble(NULL))
  }
  
  roc(response = df$response, predictor = df$predictor, direction = "<")
}

calc_mcc = function(df) {
  p = df %>%
    filter(response == 1) %>%
    nrow()
  
  n = df %>%
    filter(response == 0) %>%
    nrow()
  
  r = roc(response = df$response, predictor = df$predictor, direction = "<")
  res = tibble(sensitivity = r$sensitivities,
               specificity = r$specificities,
               p = p,
               n = n) %>%
    arrange(sensitivity, desc(specificity)) %>%
    mutate(tp = sensitivity * p,
           fn = p - tp,
           tn = specificity * n,
           fp = n - tn) %>%
    mutate(denom = as.double(tp + fp) * (tp + fn) * (tn + fp) * (tn + fn),
           denom = case_when(denom == 0 ~ 1,
                            TRUE ~ denom),
           mcc = ((tp * tn) - (fp * fn))/sqrt(denom))
  return(res)
}

calc_bac = function(df) {
  p = df %>%
    filter(response == 1) %>%
    nrow()
  
  n = df %>%
    filter(response == 0) %>%
    nrow()
  
  r = roc(response = df$response, predictor = df$predictor, direction = "<")
  res = tibble(sensitivity = r$sensitivities,
               specificity = r$specificities,
               p = p,
               n = n)  %>%
    arrange(sensitivity, desc(specificity)) %>%
    mutate(tp = sensitivity * p,
           tn = specificity * n) %>%
    mutate(bac = 0.5*(tp/p + tn/n))
    
  return(res)
  
}
