#------ NAUDOJAMOS BIBLIOTEKOS --------
library(dplyr)
library(jsonlite)
library(ggplot2)
library(tidyr)
library(readr)
library(tidyverse)
library(ggstatsplot)
library(psych) 
library(car)
library(lmtest)
library(GGally)
library(grid)
library(pROC)
library(vegan)
library(boot)
library(QuantPsyc)
library(ROCR)
library(gridExtra)
library(MLmetrics)
library(uwot)
library(corrplot)
library(reshape2)
library(randomForest)
library(tibble)
library(patchwork)
library(caret)
library(glmnet)
library(FSelector)
library(scales)
library(cluster)
library(factoextra)
library(rcompanion)
library(ggcorrplot)
library(mclust)
library(e1071)
library(pheatmap)
library(readxl)
library(stringr)
library(GGally)
library(plotly)
library(SaturnCoefficient)
library(rstatix)
library(ggtext)
library(Hmisc)
library(aricode)
library(clue)
library(jsonlite)
library(purrr)
library(xgboost)
library(kernlab)
library(klaR)
library(dbscan)
#--------------DUOMENŲ NUSKAITYMAS, POŽYMIAI IR MEDIKŲ METRIKOS ---------
# Nuskaitomi teksturos pozymiai tik globali struktūra be fazių----------------
json_data <- fromJSON("/Users/arinaperzu/Desktop/MRT darbas/MRT2/src/features/patient_features_phases_1_to_25.json")
df_raw <- as.data.frame(t(sapply(json_data, unlist)))

df_raw <- df_raw |>
  rownames_to_column("Patient_ID") |>
  mutate(Class = ifelse(grepl("^AS", Patient_ID), 1, 0))

# Struktūra, kad galima būtų dirbti su visomis zonomis ir palyginti
extract_structure <- function(df, structure_name) {
  
  structure_pattern <- paste0("_", structure_name, "$")
  feature_cols <- colnames(df)[
    stringr::str_detect(colnames(df), structure_pattern)
  ]
  selected_df <- df |>
    dplyr::select(
      Patient_ID,
      dplyr::all_of(feature_cols),
      Class
    )
  
  clean_names <- stringr::str_remove(feature_cols, structure_pattern)
  colnames(selected_df)[
    match(feature_cols, colnames(selected_df))
  ] <- clean_names
  
  selected_df
}
structures <- c("global", "bottom", "mid", "top")

#------------ POŽYMIŲ PAVADINIMŲ TRUMPINIAI --------------------------------

feature_short_names <- c(
  Fourier_Energy = "F_En",
  Fourier_Entropy = "F_E",
  Fourier_Mean = "F_M",
  Fourier_Variance = "F_V",
  
  Fractal_mean = "FR_M",
  Fractal_variance = "FR_V",
  
  GLCM_ASM = "GLCM_ASM",
  GLCM_contrast = "GLCM_Con",
  GLCM_correlation = "GLCM_Cor",
  GLCM_dissimilarity = "GLCM_D",
  GLCM_energy = "GLCM_En",
  GLCM_homogeneity = "GLCM_Hom",
  
  GLRLM_GrayLevelNonUniformity = "GLRLM_GLNU",
  GLRLM_GrayLevelNonUniformityNormalized = "GLRLM_GLNUN",
  GLRLM_GrayLevelVariance = "GLRLM_GLV",
  GLRLM_HighGrayLevelRunEmphasis = "GLRLM_HGLE",
  GLRLM_LongRunEmphasis = "GLRLM_LRE",
  GLRLM_LongRunHighGrayLevelEmphasis = "GLRLM_LRHGE",
  GLRLM_LongRunLowGrayLevelEmphasis = "GLRLM_LRLGE",
  GLRLM_LowGrayLevelRunEmphasis = "GLRLM_LGLE",
  GLRLM_RunEntropy = "GLRLM_RE",
  GLRLM_RunLengthNonUniformity = "GLRLM_RLNU",
  GLRLM_RunLengthNonUniformityNormalized = "GLRLM_RLNUN",
  GLRLM_RunPercentage = "GLRLM_RP",
  GLRLM_RunVariance = "GLRLM_RV",
  GLRLM_ShortRunEmphasis = "GLRLM_SRE",
  GLRLM_ShortRunHighGrayLevelEmphasis = "GLRLM_SRHGLE",
  GLRLM_ShortRunLowGrayLevelEmphasis = "GLRLM_SRLGLE",
  
  HOG_mean = "HOG_M",
  HOG_variance = "HOG_V",
  
  LBP_contrast = "LBP_Con",
  LBP_uniformity = "LBP_U",
  
  Wavelet_level_1_HH_mean = "WL1_HH_M",
  Wavelet_level_1_HH_sd = "WL1_HH_Sd",
  Wavelet_level_1_HL_mean = "WL1_HL_M",
  Wavelet_level_1_HL_sd = "WL1_HL_Sd",
  Wavelet_level_1_LH_mean = "WL1_LH_M",
  Wavelet_level_1_LH_sd = "WL1_LH_Sd",
  
  Wavelet_level_2_HH_mean = "WL2_HH_M",
  Wavelet_level_2_HH_sd = "WL2_HH_Sd",
  Wavelet_level_2_HL_mean = "WL2_HL_M",
  Wavelet_level_2_HL_sd = "WL2_HL_Sd",
  Wavelet_level_2_LH_mean = "WL2_LH_M",
  Wavelet_level_2_LH_sd = "WL2_LH_Sd",
  Wavelet_level_2_LL_mean = "WL2_LL_M",
  Wavelet_level_2_LL_sd = "WL2_LL_Sd"
)

rename_features_short <- function(data, short_names) {
  
  old_names <- names(short_names)
  new_names <- unname(short_names)
  names(data) <- ifelse(
    names(data) %in% old_names,
    new_names[match(names(data), old_names)],
    names(data)
  )
  data
}

datasets <- structures |>
  purrr::set_names() |>
  purrr::map(\(x) extract_structure(df_raw, x))
names(datasets)

datasets <- datasets |>
  purrr::map(\(data) rename_features_short(data, feature_short_names))

#-------------- ATSKIRI NENORMALIZUOTI DUOMENŲ RINKINIAI ---------
df_global <- datasets$global
df_bottom <- datasets$bottom
df_mid <- datasets$mid
df_top <- datasets$top

# Grafikams ir lentelėms, kad būtų patogu lyginti tarpusavyje zonas
df_all_structures <- datasets |>
  purrr::imap_dfr(\(data, structure_name) {
    data |>
      dplyr::mutate(Structure = structure_name)
  })

# patikrinimas
df_all_structures |>
  dplyr::count(Structure, Class)

# ----------- DUOMENU MIN_MAX NORMALIZAVIMAS -------------
minmax_scale <- function(x) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

datasets_scaled <- datasets |>
  purrr::map(\(data) {
    data |>
      dplyr::mutate(
        dplyr::across(
          .cols = -c(Patient_ID, Class),
          .fns = minmax_scale
        )
      )
  })

# atskiros lenteles
df_global_scaled <- datasets_scaled$global
df_bottom_scaled <- datasets_scaled$bottom
df_mid_scaled    <- datasets_scaled$mid
df_top_scaled    <- datasets_scaled$top

# bendra lentele normalizuotiems duomenims
df_all_structures_scaled <- datasets_scaled |>
  purrr::imap_dfr(\(data, structure_name) {
    data |>
      dplyr::mutate(Structure = structure_name)
  })

# Klinikiniu metriku nuskaitymas is Exel lenteles ---------------------
# Imame tik ID ir kelias metrikas pagal pavadinimus
pacientu_lentele_su_med_metrikomis <- read_excel("/Users/arinaperzu/Desktop/MRT darbas/MRT Rkodas/med.xlsx")

pacientai <- pacientu_lentele_su_med_metrikomis |>
  dplyr::select(ID, k_sex, h_CVF, 'm_Nat T1', m_ECV, e_GLS)|>
  dplyr::rename(
    T1 = `m_Nat T1`,
    CVF = h_CVF,
    sex = k_sex,
    ECV = m_ECV,
    GLS = e_GLS
  )

# Duomenu apjungimas -----------------
# Paimam numerius, kad galima butu sujungti, turime tik sergantiems pacientams med metrikas, todel ziurime tik AS
add_clinical_metrics <- function(data, pacientai) {
  data |>
    dplyr::mutate(
      ID = dplyr::if_else(
        stringr::str_starts(Patient_ID, "AS"),
        readr::parse_number(Patient_ID),
        NA_real_
      )
    ) |>
    dplyr::left_join(pacientai, by = "ID")
}

# medicininiu metriku normalizavimas 
pacientai_scaled <- pacientai |>
  dplyr::mutate(
    dplyr::across(
      .cols = c(CVF, T1, ECV, GLS),
      .fns = minmax_scale
    )
  )

# Medicininių metrikų apjungimas su tekstūros požymiais
datasets_with_clinical <- datasets |>
  purrr::map(\(data) add_clinical_metrics(data, pacientai))

datasets_scaled_with_clinical <- datasets_scaled |>
  purrr::map(\(data) add_clinical_metrics(data, pacientai))

datasets_scaled_with_clinical_scaled <- datasets_scaled |>
  purrr::map(\(data) add_clinical_metrics(data, pacientai_scaled))

# bendros lentelės
df_all_structures_with_clinical <- datasets_with_clinical |>
  purrr::imap_dfr(\(data, structure_name) {
    data |>
      dplyr::mutate(Structure = structure_name)
  })

df_all_structures_scaled_with_clinical <- datasets_scaled_with_clinical |>
  purrr::imap_dfr(\(data, structure_name) {
    data |>
      dplyr::mutate(Structure = structure_name)
  })

df_all_structures_scaled_with_clinical_scaled <- datasets_scaled_with_clinical_scaled |>
  purrr::imap_dfr(\(data, structure_name) {
    data |>
      dplyr::mutate(Structure = structure_name)
  })


# ---------- REIKSMINGU POZYMIU ATRANKA ----------------------------------
#------------ Mann Whitney testas -------------------------
run_mann_whitney <- function(data, class_col = "Class", id_col = "Patient_ID") {
  
  feature_cols <- setdiff(names(data), c(id_col, class_col))
  feature_cols <- feature_cols[
    purrr::map_lgl(data[feature_cols], is.numeric)
  ]
  results <- purrr::map_dfr(feature_cols, function(feature_name) {
    
    group0 <- data |>
      dplyr::filter(.data[[class_col]] == 0) |>
      dplyr::pull(dplyr::all_of(feature_name))
    
    group1 <- data |>
      dplyr::filter(.data[[class_col]] == 1) |>
      dplyr::pull(dplyr::all_of(feature_name))
    
    group0 <- group0[!is.na(group0)]
    group1 <- group1[!is.na(group1)]
    
    if (
      length(group0) < 2 ||
      length(group1) < 2 ||
      length(unique(c(group0, group1))) < 2
    ) {
      p_value <- NA_real_
    } else {
      test_result <- wilcox.test(
        group0,
        group1,
        exact = FALSE
      )
      p_value <- test_result$p.value
    }
    
    tibble::tibble(
      Feature = feature_name,
      N_Control = length(group0),
      N_AS = length(group1),
      Median_Control = median(group0, na.rm = TRUE),
      Median_AS = median(group1, na.rm = TRUE),
      Mean_Control = mean(group0, na.rm = TRUE),
      Mean_AS = mean(group1, na.rm = TRUE),
      Median_Difference = median(group1, na.rm = TRUE) - median(group0, na.rm = TRUE),
      Mean_Difference = mean(group1, na.rm = TRUE) - mean(group0, na.rm = TRUE),
      p_value = p_value
    )
  }) |>
    dplyr::mutate(
      p_adj = p.adjust(p_value, method = "BH"),
      Significant = dplyr::case_when(
        is.na(p_adj) ~ NA,
        p_adj < 0.05 ~ TRUE,
        TRUE ~ FALSE
      )
    ) |>
    dplyr::arrange(p_adj, p_value)
  results
}

select_mw_features <- function(mw_results, top_n = 10, alpha = 0.05) {
  mw_results |>
    dplyr::filter(!is.na(p_adj)) |>
    dplyr::filter(p_adj < alpha) |>
    dplyr::arrange(p_adj, p_value) |>
    dplyr::slice_head(n = top_n) |>
    dplyr::pull(Feature)
}

#------------ Spearmeno koreliacija -------------------------
spearman_filter_features <- function(data, features, cutoff = 0.85) {
  
  if (length(features) <= 1) {
    return(features)
  }
  feature_data <- data |>
    dplyr::select(dplyr::all_of(features))
  cor_matrix <- cor(
    feature_data,
    method = "spearman",
    use = "pairwise.complete.obs"
  )
  cor_matrix[is.na(cor_matrix)] <- 0
  remove_idx <- caret::findCorrelation(
    cor_matrix,
    cutoff = cutoff,
    names = FALSE
  )
  if (length(remove_idx) == 0) {
    return(features)
  }
  features[-remove_idx]
}

#------------ Lasso regresija -------------------------
run_lasso_selection <- function(data, class_col = "Class", id_col = "Patient_ID", top_n = 10) {
  model_data <- data |>
    dplyr::select(-dplyr::all_of(id_col))
  
  x <- model.matrix(
    as.formula(paste(class_col, "~ .")),
    data = model_data
  )[, -1]
  
  y <- model_data[[class_col]]
  set.seed(123)
  
  lasso_model <- cv.glmnet(
    x = x,
    y = y,
    family = "binomial",
    alpha = 1,
    standardize = FALSE
  )
  coef_df <- coef(lasso_model, s = "lambda.1se") |>
    as.matrix() |>
    as.data.frame() |>
    tibble::rownames_to_column(var = "Feature")
  
  colnames(coef_df)[2] <- "Coefficient"
  selected <- coef_df |>
    dplyr::filter(Feature != "(Intercept)") |>
    dplyr::filter(Coefficient != 0) |>
    dplyr::mutate(Abs_Coefficient = abs(Coefficient)) |>
    dplyr::arrange(dplyr::desc(Abs_Coefficient)) |>
    dplyr::slice_head(n = top_n)
  
  # Jeigu lambda.1se nieko neatrinko, naudojame lambda.min
  if (nrow(selected) == 0) {
    coef_df <- coef(lasso_model, s = "lambda.min") |>
      as.matrix() |>
      as.data.frame() |>
      tibble::rownames_to_column(var = "Feature")
    
    colnames(coef_df)[2] <- "Coefficient"
    selected <- coef_df |>
      dplyr::filter(Feature != "(Intercept)") |>
      dplyr::filter(Coefficient != 0) |>
      dplyr::mutate(Abs_Coefficient = abs(Coefficient)) |>
      dplyr::arrange(dplyr::desc(Abs_Coefficient)) |>
      dplyr::slice_head(n = top_n)
  }
  selected
}

#------------ Random forest -------------------------
run_rf_selection <- function(data, class_col = "Class", id_col = "Patient_ID", top_n = 10) {
  model_data <- data |>
    dplyr::select(-dplyr::all_of(id_col)) |>
    dplyr::mutate(
      Class = factor(.data[[class_col]], levels = c(0, 1), labels = c("Control", "AS"))
    )
  
  set.seed(123)
  rf_model <- randomForest(
    as.formula(paste(class_col, "~ .")),
    data = model_data,
    importance = TRUE,
    ntree = 500
  )
  importance_scores <- as.data.frame(importance(rf_model)) |>
    tibble::rownames_to_column(var = "Feature")
  
  selected <- importance_scores |>
    dplyr::arrange(dplyr::desc(MeanDecreaseGini)) |>
    dplyr::slice_head(n = top_n)
  selected
}

#------------ XGBOOST POŽYMIŲ ATRANKA -------------------------
run_xgb_selection <- function(data, class_col = "Class", id_col = "Patient_ID", top_n = 10) {
  model_data <- data |>
    dplyr::select(-dplyr::all_of(id_col))
  
  feature_cols <- setdiff(names(model_data), class_col)
  x <- model_data |>
    dplyr::select(dplyr::all_of(feature_cols)) |>
    as.matrix()
  y <- model_data[[class_col]]
  y <- as.numeric(as.character(y))
  dtrain <- xgboost::xgb.DMatrix(data = x, label = y)
  
  params <- list(
    objective = "binary:logistic",
    eval_metric = "logloss",
    max_depth = 2,
    eta = 0.05,
    subsample = 0.8,
    colsample_bytree = 0.8,
    min_child_weight = 1
  )
  
  set.seed(123)
  xgb_cv <- xgboost::xgb.cv(
    params = params,
    data = dtrain,
    nrounds = 300,
    nfold = 5,
    stratified = TRUE,
    early_stopping_rounds = 20,
    verbose = 0
  )
  
  # iteraciju skaicius
  if (!is.null(xgb_cv$best_iteration) && length(xgb_cv$best_iteration) > 0) {
    best_nrounds <- xgb_cv$best_iteration
  } else {
    eval_log <- xgb_cv$evaluation_log
    
    if ("test_logloss_mean" %in% names(eval_log)) {
      best_nrounds <- which.min(eval_log$test_logloss_mean)
    } else {
      best_nrounds <- 50
    }
  }
  # Papildoma apsauga
  if (is.null(best_nrounds) || length(best_nrounds) == 0 || is.na(best_nrounds) || best_nrounds < 1) {
    best_nrounds <- 50
  }
  
  set.seed(123)
  
  xgb_model <- xgboost::xgb.train(
    params = params,
    data = dtrain,
    nrounds = best_nrounds,
    verbose = 0
  )
  importance_scores <- xgboost::xgb.importance(
    feature_names = feature_cols,
    model = xgb_model
  )
  if (nrow(importance_scores) == 0) {
    return(
      tibble::tibble(
        Feature = character(),
        Gain = numeric(),
        Cover = numeric(),
        Frequency = numeric()
      )
    )
  }
  selected <- importance_scores |>
    dplyr::slice_head(n = top_n)
  selected
}

#------------ PILNA POŽYMIŲ ATRANKA VIENAI STRUKTŪRAI -------------------------
run_feature_selection_for_structure <- function(
    data,
    structure_name,
    class_col = "Class",
    id_col = "Patient_ID",
    top_n = 10,
    alpha = 0.05,
    cor_cutoff = 0.85
) {
  
  message("Atliekama požymių atranka struktūrai: ", structure_name)
  
  # Mann-Whitney
  mw_results <- run_mann_whitney(
    data = data,
    class_col = class_col,
    id_col = id_col
  )
  mw_features <- select_mw_features(
    mw_results = mw_results,
    top_n = top_n,
    alpha = alpha
  )
  
  # Spirmeno koreliacija
  mw_features_filtered <- spearman_filter_features(
    data = data,
    features = mw_features,
    cutoff = cor_cutoff
  )
  
  # LASSO
  lasso_results <- run_lasso_selection(
    data = data,
    class_col = class_col,
    id_col = id_col,
    top_n = top_n
  )
  lasso_features <- lasso_results$Feature
  
  lasso_features_filtered <- spearman_filter_features(
    data = data,
    features = lasso_features,
    cutoff = cor_cutoff
  )
  
  # Random Forest
  rf_results <- run_rf_selection(
    data = data,
    class_col = class_col,
    id_col = id_col,
    top_n = top_n
  )
  rf_features <- rf_results$Feature
  
  rf_features_filtered <- spearman_filter_features(
    data = data,
    features = rf_features,
    cutoff = cor_cutoff
  )
  
  # XGBoost
  xgb_results <- run_xgb_selection(
    data = data,
    class_col = class_col,
    id_col = id_col,
    top_n = top_n
  )
  xgb_features <- xgb_results$Feature
  xgb_features_filtered <- spearman_filter_features(
    data = data,
    features = xgb_features,
    cutoff = cor_cutoff
  )
  
  # Visi atrinkti požymiai pagal metodus
  selected_by_method <- tibble::tibble(
    Structure = structure_name,
    Method = c(
      rep("Mann_Whitney", length(mw_features_filtered)),
      rep("LASSO", length(lasso_features_filtered)),
      rep("Random_Forest", length(rf_features_filtered)),
      rep("XGBoost", length(xgb_features_filtered))
    ),
    Feature = c(
      mw_features_filtered,
      lasso_features_filtered,
      rf_features_filtered,
      xgb_features_filtered
    )
  )
  
  # Stabilumas pagal metodus
  feature_votes <- selected_by_method |>
    dplyr::count(Structure, Feature, name = "N_methods") |>
    dplyr::arrange(dplyr::desc(N_methods), Feature)
  
  # Požymiai, kuriuos atrinko bent 2 metodai
  stable_features <- feature_votes |>
    dplyr::filter(N_methods >= 2) |>
    dplyr::pull(Feature)
  
  if (length(stable_features) == 0) {
    stable_features <- feature_votes |>
      dplyr::slice_head(n = min(top_n, n())) |>
      dplyr::pull(Feature)
  }
  
  # Galutinis duomenų rinkinys
  final_data <- data |>
    dplyr::select(
      dplyr::all_of(c(id_col, class_col)),
      dplyr::all_of(stable_features)
    )
  list(
    structure = structure_name,
    mann_whitney_results = mw_results,
    mann_whitney_features = mw_features_filtered,
    lasso_results = lasso_results,
    lasso_features = lasso_features_filtered,
    rf_results = rf_results,
    rf_features = rf_features_filtered,
    xgb_results = xgb_results,
    xgb_features = xgb_features_filtered,
    selected_by_method = selected_by_method,
    feature_votes = feature_votes,
    stable_features = stable_features,
    final_data = final_data
  )
}

#------------ POŽYMIŲ ATRANKA VISOMS STRUKTŪROMS -------------------------
feature_selection_results <- purrr::imap(
  datasets_scaled,
  ~ run_feature_selection_for_structure(
    data = .x,
    structure_name = .y,
    top_n = 10,
    alpha = 0.05,
    cor_cutoff = 0.8
  )
)

#------------ ATRINKTI DUOMENYS ATSKIROMS STRUKTŪROMS -------------------------
df_global_selected <- feature_selection_results$global$final_data
df_bottom_selected <- feature_selection_results$bottom$final_data
df_mid_selected <- feature_selection_results$mid$final_data
df_top_selected <- feature_selection_results$top$final_data


# -------------- GLOBAL STRUKTŪROS KLASTERIZAVIMAS --------------
# -----------K-means, HDBSCAN, Hierarchinis metoda -------------
set.seed(123)

df_cluster_global <- df_global_selected |>
  dplyr::distinct(Patient_ID, .keep_all = TRUE) |>
  tidyr::drop_na()

# Tik atrinkti tekstūros požymiai
x_global <- df_cluster_global |>
  dplyr::select(-Patient_ID, -Class)

x_global <- x_global |>
  dplyr::mutate(
    dplyr::across(
      dplyr::everything(),
      as.numeric
    )
  )
x_mat <- as.matrix(x_global)

true_class <- df_cluster_global$Class
patient_ids <- df_cluster_global$Patient_ID

# ---------- Pagalbinės funkcijos ---------------------------------------
calculate_mean_silhouette <- function(x, clusters) {
  clusters <- as.integer(as.factor(clusters))
  if (length(unique(clusters)) < 2) {
    return(NA_real_)
  }
  if (length(unique(clusters)) >= length(clusters)) {
    return(NA_real_)
  }
  dist_x <- dist(x)
  sil <- cluster::silhouette(clusters, dist_x)
  mean(sil[, "sil_width"], na.rm = TRUE)
}

calculate_elbow_wss <- function(x, k_range = 2:10) {
  purrr::map_dfr(k_range, function(k) {
    set.seed(123)
    km <- kmeans(
      x = x,
      centers = k,
      nstart = 50,
      iter.max = 100
    )
    tibble::tibble(
      k = k,
      WSS = km$tot.withinss
    )
  })
}


calculate_kmeans_silhouette <- function(x, k_range = 2:10) {
  purrr::map_dfr(k_range, function(k) {
    set.seed(123)
    km <- kmeans(
      x = x,
      centers = k,
      nstart = 50,
      iter.max = 100
    )
    tibble::tibble(
      k = k,
      Mean_Silhouette = calculate_mean_silhouette(x, km$cluster)
    )
  })
}


calculate_hclust_silhouette <- function(x, k_range = 2:10, method = "ward.D2") {
  dist_x <- dist(x)
  hc <- hclust(dist_x, method = method)
  purrr::map_dfr(k_range, function(k) {
    clusters <- cutree(hc, k = k)
    tibble::tibble(
      k = k,
      Mean_Silhouette = calculate_mean_silhouette(x, clusters)
    )
  })
}

# ---------- k parinkimas K-means ir hierarchiniam metodui --------------
k_range <- 2:10

# K-means elbow
kmeans_elbow <- calculate_elbow_wss(
  x = x_mat,
  k_range = k_range
)

# K-means silhouette
kmeans_silhouette <- calculate_kmeans_silhouette(
  x = x_mat,
  k_range = k_range
)
best_k_kmeans <- kmeans_silhouette |>
  dplyr::filter(!is.na(Mean_Silhouette)) |>
  dplyr::arrange(dplyr::desc(Mean_Silhouette)) |>
  dplyr::slice(1) |>
  dplyr::pull(k)

# Hierarchinis silhouette
hclust_silhouette <- calculate_hclust_silhouette(
  x = x_mat,
  k_range = k_range,
  method = "ward.D2"
)
best_k_hclust <- hclust_silhouette |>
  dplyr::filter(!is.na(Mean_Silhouette)) |>
  dplyr::arrange(dplyr::desc(Mean_Silhouette)) |>
  dplyr::slice(1) |>
  dplyr::pull(k)

best_k_kmeans
best_k_hclust

# ---------- k parinkimo grafikai ---------------------------------------
p_elbow <- ggplot(kmeans_elbow, aes(x = k, y = WSS)) +
  geom_line() +
  geom_point(size = 2) +
  scale_x_continuous(breaks = k_range) +
  theme_minimal() +
  labs(
    title = "K-means elbow metodas: global struktūra",
    x = "Klasterių skaičius k",
    y = "Total within-cluster sum of squares"
  )

p_kmeans_sil <- ggplot(kmeans_silhouette, aes(x = k, y = Mean_Silhouette)) +
  geom_line() +
  geom_point(size = 2) +
  scale_x_continuous(breaks = k_range) +
  theme_minimal() +
  labs(
    title = "K-means vidutinė silueto reikšmė",
    x = "Klasterių skaičius k",
    y = "Vidutinis silhouette"
  )

p_hclust_sil <- ggplot(hclust_silhouette, aes(x = k, y = Mean_Silhouette)) +
  geom_line() +
  geom_point(size = 2) +
  scale_x_continuous(breaks = k_range) +
  theme_minimal() +
  labs(
    title = "Hierarchinio klasterizavimo vidutinė silueto reikšmė",
    x = "Klasterių skaičius k",
    y = "Vidutinis silhouette"
  )

p_elbow
p_kmeans_sil
p_hclust_sil

# ----------  Klasterizavimo metodai -------------------------------------
# K-means
set.seed(123)
kmeans_model <- kmeans(
  x = x_mat,
  centers = 5,
  nstart = 50,
  iter.max = 100
)
clusters_kmeans <- kmeans_model$cluster

# ----------- Hierarchinis klasterizavimas---------------
dist_global <- dist(x_mat)

hclust_model <- hclust(
  d = dist_global,
  method = "ward.D2"
)

clusters_hclust <- cutree(
  hclust_model,
  k = 5
)


#------------- HDBSCAN ----------------
set.seed(123)
hdbscan_model <- dbscan::hdbscan(
  x = x_mat,
  minPts = 5
)
clusters_hdbscan <- hdbscan_model$cluster
table(clusters_hdbscan)

# ---------- Vidutinė silueto reikšmė -----------------------------------
silhouette_results <- tibble::tibble(
  Method = c("K-means", "HDBSCAN", "Hierarchical"),
  Mean_Silhouette = c(
    calculate_mean_silhouette(x_mat, clusters_kmeans),
    calculate_mean_silhouette(x_mat, clusters_hdbscan),
    calculate_mean_silhouette(x_mat, clusters_hclust)
  ),
  N_clusters = c(
    length(unique(clusters_kmeans)),
    length(unique(clusters_hdbscan)),
    length(unique(clusters_hclust))
  )
)
silhouette_results

# ---------- ARI ir NMI tarp tikrosios klasės ir metodo klasterių ----------
cluster_assignments <- tibble::tibble(
  Patient_ID = patient_ids,
  Class = true_class,
  Kmeans = clusters_kmeans,
  HDBSCAN = clusters_hdbscan,
  Hierarchical = clusters_hclust
)
cluster_assignments

compare_class_with_clusters <- function(true_labels, cluster_labels, method_name) {
  tibble::tibble(
    Method = method_name,
    ARI = mclust::adjustedRandIndex(true_labels, cluster_labels),
    NMI = aricode::NMI(true_labels, cluster_labels)
  )
}

class_cluster_agreement_results <- dplyr::bind_rows(
  compare_class_with_clusters(
    true_labels = cluster_assignments$Class,
    cluster_labels = cluster_assignments$Kmeans,
    method_name = "K-means"
  ),
  compare_class_with_clusters(
    true_labels = cluster_assignments$Class,
    cluster_labels = cluster_assignments$HDBSCAN,
    method_name = "HDBSCAN"
  ),
  compare_class_with_clusters(
    true_labels = cluster_assignments$Class,
    cluster_labels = cluster_assignments$Hierarchical,
    method_name = "Hierarchical"
  )
)

class_cluster_agreement_results

# ---------- UMAP projekcija --------------------------------------------
set.seed(123)
umap_global <- uwot::umap(
  X = x_mat,
  n_neighbors = 15,
  min_dist = 0.1,
  n_components = 2,
  metric = "euclidean",
  scale = FALSE
)

umap_df <- tibble::tibble(
  Patient_ID = patient_ids,
  Class = factor(true_class, levels = c(0, 1), labels = c("Kontrolinė", "Patologinė")),
  UMAP1 = umap_global[, 1],
  UMAP2 = umap_global[, 2],
  Kmeans = factor(clusters_kmeans),
  HDBSCAN = factor(clusters_hdbscan),
  Hierarchical = factor(clusters_hclust)
)

# ---------- UMAP grafikai ---------------------------------------------
custom_colors <- c(
  "1" = "orange",
  "2" = "royalblue1",
  "3" = "royalblue4",
  "4" = "violetred4",
  "5" = "violet",
  "6" = "grey60"   
)

p_umap_kmeans <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Kmeans, shape = Class)) +
  geom_point(size = 2.5, alpha = 0.85) +
  theme_minimal() +
  scale_color_manual(values = custom_colors) +
  scale_shape_manual(
    values = c(
      "Kontrolinė" = 17,  
      "Patologinė" = 16       
    )
  ) +
  theme_minimal() +
  theme(
    legend.position = "none"
  ) +
  labs(
    title = "UMAP: K-means klasteriai",
    color = "Klasteris",
    shape = "Grupė"
  )

p_umap_hdbscan <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = HDBSCAN)) +
  geom_point(size = 2.5, alpha = 0.85) +
  scale_color_manual(values = custom_colors) +
  theme_minimal() +
  labs(
    title = "UMAP: HDBSCAN klasteriai",
    color = "Klasteris"
  )

p_umap_hclust <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Hierarchical, shape = Class)) +
  geom_point(size = 2.5, alpha = 0.85) +
  scale_color_manual(values = custom_colors) +
  scale_shape_manual(
    values = c(
      "Kontrolinė" = 17,  
      "Patologinė" = 16       
    )
  ) +
  theme_minimal() +
  labs(
    title = "UMAP: hierarchinio metodo klasteriai",
    color = "Klasteris",
    shape = "Grupė"
  )

p_umap_kmeans
p_umap_hdbscan
p_umap_hclust

p_umap_kmeans_clean <- p_umap_kmeans +
  labs(
    title = expression(bold(italic(k)) * bold(" vidurkių")),
    subtitle = NULL
  ) +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5))

p_umap_hclust_clean <- p_umap_hclust +
  labs(
    title = "Hierarchinis klasterizavimas",
    subtitle = NULL
  ) +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5))


p_umap_kmeans_hclust <- p_umap_kmeans_clean + p_umap_hclust_clean +
  plot_annotation(
    title = "UMAP projekcija: K-means ir hierarchinio klasterizavimo palyginimas",
    subtitle = "UMAP parametrai: n = 15, min_d = 0.1, metric = euclidean",
    theme = theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5)
    )
  )

p_umap_kmeans_hclust

# ----------  Projekcijos kokybės metrikos -------------------------------
get_rank_matrix <- function(dist_matrix) {
  n <- nrow(dist_matrix)
  rank_matrix <- matrix(NA_real_, nrow = n, ncol = n)
  for (i in seq_len(n)) {
    rank_matrix[i, ] <- rank(dist_matrix[i, ], ties.method = "first")
  }
  rank_matrix
}

calculate_trustworthiness <- function(x_high, x_low, k = 10) {
  n <- nrow(x_high)
  if (k >= n / 2) {
    stop("k turi būti mažesnis nei n / 2")
  }
  
  dist_high <- as.matrix(dist(x_high))
  dist_low <- as.matrix(dist(x_low))
  rank_high <- get_rank_matrix(dist_high)
  trust_sum <- 0
  for (i in seq_len(n)) {
    neighbors_high <- order(dist_high[i, ])[2:(k + 1)]
    neighbors_low <- order(dist_low[i, ])[2:(k + 1)]
    unexpected_neighbors <- setdiff(neighbors_low, neighbors_high)
    if (length(unexpected_neighbors) > 0) {
      trust_sum <- trust_sum + sum(rank_high[i, unexpected_neighbors] - k)
    }
  }
  trustworthiness <- 1 - (
    2 / (n * k * (2 * n - 3 * k - 1))
  ) * trust_sum
  trustworthiness
}


calculate_continuity <- function(x_high, x_low, k = 10) {
  n <- nrow(x_high)
  if (k >= n / 2) {
    stop("k turi būti mažesnis nei n / 2")
  }
  dist_high <- as.matrix(dist(x_high))
  dist_low <- as.matrix(dist(x_low))
  rank_low <- get_rank_matrix(dist_low)
  continuity_sum <- 0
  
  for (i in seq_len(n)) {
    neighbors_high <- order(dist_high[i, ])[2:(k + 1)]
    neighbors_low <- order(dist_low[i, ])[2:(k + 1)]
    missing_neighbors <- setdiff(neighbors_high, neighbors_low)
    
    if (length(missing_neighbors) > 0) {
      continuity_sum <- continuity_sum + sum(rank_low[i, missing_neighbors] - k)
    }
  }
  continuity <- 1 - (
    2 / (n * k * (2 * n - 3 * k - 1))
  ) * continuity_sum
  continuity
}

# ---------- UMAP trustworthiness ir continuity -------------------------
umap_projection_metrics <- tibble::tibble(
  Projection = "UMAP",
  k_neighbors = 15,
  Trustworthiness = calculate_trustworthiness(
    x_high = x_mat,
    x_low = umap_global,
    k = 15
  ),
  Continuity = calculate_continuity(
    x_high = x_mat,
    x_low = umap_global,
    k = 15
  )
)

umap_projection_metrics

# ---------- Galutinė metodų palyginimo lentelė -------------------------
silhouette_results 
class_cluster_agreement_results
umap_projection_metrics

calculate_hdbscan_silhouette_no_noise <- function(x, clusters) {
  keep <- clusters != 0
  if (sum(keep) < 3) {
    return(NA_real_)
  }
  calculate_mean_silhouette(
    x = x[keep, , drop = FALSE],
    clusters = clusters[keep]
  )
}

hdbscan_grid <- purrr::map_dfr(3:6, function(minPts_value) {
  model <- dbscan::hdbscan(
    x = x_mat,
    minPts = minPts_value
  )
  tibble::tibble(
    minPts = minPts_value,
    N_clusters = length(setdiff(unique(model$cluster), 0)),
    N_noise = sum(model$cluster == 0),
    Mean_Silhouette_all = calculate_mean_silhouette(x_mat, model$cluster),
    Mean_Silhouette_no_noise = calculate_hdbscan_silhouette_no_noise(
      x_mat,
      model$cluster
    )
  )
})
hdbscan_grid

# ---------- Pacientai pagal klasterius -----------------
cluster_class_df <- tibble::tibble(
  Patient_ID = patient_ids,
  Class = factor(
    true_class,
    levels = c(0, 1),
    labels = c("Kontrolinė", "Patologinė")
  ),
  Kmeans = factor(clusters_kmeans),
  Hierarchical = factor(clusters_hclust)
)

cluster_class_df


class_colors <- c(
  "Kontrolinė" = "royalblue2",
  "Patologinė" = "violetred3"
)

# ---------- K-means: pacientų skaičius klasteriuose ----------
kmeans_cluster_counts <- cluster_class_df |>
  dplyr::count(Kmeans, Class, name = "N")

# ---------- Hierarchinis: pacientų skaičius klasteriuose ----------
hclust_cluster_counts <- cluster_class_df |>
  dplyr::count(Hierarchical, Class, name = "N")


kmeans_cluster_totals <- kmeans_cluster_counts |>
  dplyr::group_by(Kmeans) |>
  dplyr::summarise(
    Total = sum(N),
    .groups = "drop"
  )

hclust_cluster_totals <- hclust_cluster_counts |>
  dplyr::group_by(Hierarchical) |>
  dplyr::summarise(
    Total = sum(N),
    .groups = "drop"
  )

max_cluster_n <- max(
  c(
    kmeans_cluster_totals$Total,
    hclust_cluster_totals$Total
  ),
  na.rm = TRUE
)

# ---------- K-means grafikas ----------
p_kmeans_cluster_counts_stacked <- ggplot(
  kmeans_cluster_counts,
  aes(x = Kmeans, y = N, fill = Class)
) +
  geom_col(
    width = 0.7,
    color = "white"
  ) +
  geom_text(
    aes(label = N),
    position = position_stack(vjust = 0.5),
    size = 3.5,
    color = "white",
    fontface = "bold"
  ) +
  scale_fill_manual(values = class_colors) +
  scale_y_continuous(
    limits = c(0, max_cluster_n * 1.08),
    breaks = scales::pretty_breaks(n = 6),
    expand = expansion(mult = c(0, 0))
  )+
  theme_minimal() +
  labs(
    title = expression(bold(italic(k)) * bold(" vidurkių")),
    x = "Klasteris",
    y = "Pacientų skaičius",
    fill = "Grupė"
  ) +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title.x = element_text(size = 12, color = "black"),
    axis.title.y = element_text(size = 12, color = "black")
  )

p_kmeans_cluster_counts_stacked

# ---------- Hierarchinio klasterizavimo grafikas ----------
p_hclust_cluster_counts_stacked <- ggplot(
  hclust_cluster_counts,
  aes(x = Hierarchical, y = N, fill = Class)
) +
  geom_col(
    width = 0.7,
    color = "white"
  ) +
  geom_text(
    aes(label = N),
    position = position_stack(vjust = 0.5),
    size = 3.5,
    color = "white",
    fontface = "bold"
  ) +
  scale_fill_manual(values = class_colors) +
  scale_y_continuous(
    limits = c(0, max_cluster_n * 1.08),
    breaks = scales::pretty_breaks(n = 6),
    expand = expansion(mult = c(0, 0))
  )+
  theme_minimal() +
  labs(
    title = "Hierarchinis klasterizavimas",
    x = "Klasteris",
    y = NULL,
    fill = "Grupė"
  ) +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.text.y = element_blank(),
    axis.title.x = element_text(size = 12, color = "black"),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank()
  )

p_hclust_cluster_counts_stacked

# ---------- Bendras palyginimo grafikas ----------
p_cluster_counts_comparison <- p_kmeans_cluster_counts_stacked +
  p_hclust_cluster_counts_stacked +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom"
  )

p_cluster_counts_comparison <- p_cluster_counts_comparison +
  plot_annotation(
    title = "Kontrolinės ir patologinės grupių pasiskirstymas tekstūros klasteriuose",
    theme = theme(
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5)
    )
  )

p_cluster_counts_comparison

# ---------- Pacientu sarasas --------------------
cluster_class_df <- tibble::tibble(
  Patient_ID = patient_ids,
  Class = factor(
    true_class,
    levels = c(0, 1),
    labels = c("Sveiki", "AS")
  ),
  Kmeans = factor(clusters_kmeans),
  Hierarchical = factor(clusters_hclust)
)

kmeans_patient_list <- cluster_class_df |>
  dplyr::arrange(Kmeans, Class, Patient_ID) |>
  dplyr::group_by(Kmeans) |>
  dplyr::summarise(
    N_patients = dplyr::n(),
    Patient_IDs = paste(Patient_ID, collapse = ", "),
    .groups = "drop"
  )

kmeans_patient_list

hclust_patient_list <- cluster_class_df |>
  dplyr::arrange(Hierarchical, Class, Patient_ID) |>
  dplyr::group_by(Hierarchical) |>
  dplyr::summarise(
    N_patients = dplyr::n(),
    Patient_IDs = paste(Patient_ID, collapse = ", "),
    .groups = "drop"
  )

hclust_patient_list

kmeans_patient_list_by_class <- cluster_class_df |>
  dplyr::arrange(Kmeans, Class, Patient_ID) |>
  dplyr::group_by(Kmeans, Class) |>
  dplyr::summarise(
    N_patients = dplyr::n(),
    Patient_IDs = paste(Patient_ID, collapse = ", "),
    .groups = "drop"
  )

kmeans_patient_list_by_class

hclust_patient_list_by_class <- cluster_class_df |>
  dplyr::arrange(Hierarchical, Class, Patient_ID) |>
  dplyr::group_by(Hierarchical, Class) |>
  dplyr::summarise(
    N_patients = dplyr::n(),
    Patient_IDs = paste(Patient_ID, collapse = ", "),
    .groups = "drop"
  )

hclust_patient_list_by_class

#  ------------- Klasteriu analize pagal kruskall wallis ----------------
texture_cluster_df <- df_cluster_global |>
  dplyr::mutate(
    Kmeans = factor(clusters_kmeans),
    Hierarchical = factor(clusters_hclust),
    Class = factor(
      Class,
      levels = c(0, 1),
      labels = c("Sveiki", "AS")
    )
  )

texture_features <- texture_cluster_df |>
  dplyr::select(-Patient_ID, -Class, -Kmeans, -Hierarchical) |>
  names()

texture_features

# ------------ KRUSKAL-WALLIS FUNKCIJA --------------
run_kruskal_by_cluster <- function(data, features, cluster_col) {
  purrr::map_dfr(features, function(feature_name) {
    test_data <- data |>
      dplyr::select(
        Feature_value = dplyr::all_of(feature_name),
        Cluster = dplyr::all_of(cluster_col)
      ) |>
      tidyr::drop_na()
    
    n_clusters <- dplyr::n_distinct(test_data$Cluster)
    if (nrow(test_data) < 3 || n_clusters < 2) {
      return(
        tibble::tibble(
          Feature = feature_name,
          Cluster_method = cluster_col,
          statistic = NA_real_,
          df = NA_real_,
          p_value = NA_real_
        )
      )
    }
    test_result <- kruskal.test(
      Feature_value ~ Cluster,
      data = test_data
    )
    tibble::tibble(
      Feature = feature_name,
      Cluster_method = cluster_col,
      statistic = unname(test_result$statistic),
      df = unname(test_result$parameter),
      p_value = test_result$p.value
    )
  }) |>
    dplyr::mutate(
      p_adj = p.adjust(p_value, method = "BH"),
      Significant = dplyr::case_when(
        is.na(p_adj) ~ NA_character_,
        p_adj < 0.001 ~ "***",
        p_adj < 0.01 ~ "**",
        p_adj < 0.05 ~ "*",
        TRUE ~ "ns"
      )
    ) |>
    dplyr::arrange(p_adj, p_value)
}

kw_texture_kmeans <- run_kruskal_by_cluster(
  data = texture_cluster_df,
  features = texture_features,
  cluster_col = "Kmeans"
)
kw_texture_kmeans

kw_texture_hclust <- run_kruskal_by_cluster(
  data = texture_cluster_df,
  features = texture_features,
  cluster_col = "Hierarchical"
)
kw_texture_hclust

kw_texture_all <- dplyr::bind_rows(
  kw_texture_kmeans,
  kw_texture_hclust
)
kw_texture_all

significant_texture_features <- kw_texture_all |>
  dplyr::filter(!is.na(p_adj)) |>
  dplyr::filter(p_adj < 0.05)
significant_texture_features

# ---------- PORINIAI PALYGINIMAI ------------
run_dunn_by_cluster <- function(data, features, cluster_col) {
  purrr::map_dfr(features, function(feature_name) {
    test_data <- data |>
      dplyr::select(
        Feature_value = dplyr::all_of(feature_name),
        Cluster = dplyr::all_of(cluster_col)
      ) |>
      tidyr::drop_na()
    
    if (
      nrow(test_data) < 3 ||
      dplyr::n_distinct(test_data$Cluster) < 2
    ) {
      return(NULL)
    }
    dunn_result <- rstatix::dunn_test(
      data = test_data,
      Feature_value ~ Cluster,
      p.adjust.method = "BH"
    )
    
    dunn_result |>
      dplyr::mutate(
        Feature = feature_name,
        Cluster_method = cluster_col
      ) |>
      dplyr::select(
        Feature,
        Cluster_method,
        group1,
        group2,
        n1,
        n2,
        statistic,
        p,
        p.adj,
        p.adj.signif
      )
  })
}

sig_texture_kmeans <- kw_texture_kmeans |>
  dplyr::filter(!is.na(p_adj), p_adj < 0.05) |>
  dplyr::pull(Feature)

dunn_texture_kmeans <- run_dunn_by_cluster(
  data = texture_cluster_df,
  features = sig_texture_kmeans,
  cluster_col = "Kmeans"
)
dunn_texture_kmeans

sig_texture_hclust <- kw_texture_hclust |>
  dplyr::filter(!is.na(p_adj), p_adj < 0.05) |>
  dplyr::pull(Feature)

dunn_texture_hclust <- run_dunn_by_cluster(
  data = texture_cluster_df,
  features = sig_texture_hclust,
  cluster_col = "Hierarchical"
)
dunn_texture_hclust

dunn_texture_all <- dplyr::bind_rows(
  dunn_texture_kmeans,
  dunn_texture_hclust
)
dunn_texture_all

# --------- MEDIKŲ METRIKOS + KLASTERIAI ----------
clinical_cluster_df <- cluster_class_df |>
  dplyr::mutate(
    ID = dplyr::if_else(
      stringr::str_starts(Patient_ID, "AS"),
      readr::parse_number(Patient_ID),
      NA_real_
    )
  ) |>
  dplyr::left_join(
    pacientai,
    by = "ID"
  )
clinical_cluster_df

clinical_features <- c("CVF", "T1", "ECV", "GLS")

kw_clinical_kmeans <- run_kruskal_by_cluster(
  data = clinical_cluster_df,
  features = clinical_features,
  cluster_col = "Kmeans"
)
kw_clinical_kmeans

kw_clinical_hclust <- run_kruskal_by_cluster(
  data = clinical_cluster_df,
  features = clinical_features,
  cluster_col = "Hierarchical"
)
kw_clinical_hclust

kw_clinical_all <- dplyr::bind_rows(
  kw_clinical_kmeans,
  kw_clinical_hclust
)
kw_clinical_all

sig_clinical_kmeans <- kw_clinical_kmeans |>
  dplyr::filter(!is.na(p_adj), p_adj < 0.05) |>
  dplyr::pull(Feature)

dunn_clinical_kmeans <- run_dunn_by_cluster(
  data = clinical_cluster_df,
  features = sig_clinical_kmeans,
  cluster_col = "Kmeans"
)

dunn_clinical_kmeans
sig_clinical_hclust <- kw_clinical_hclust |>
  dplyr::filter(!is.na(p_adj), p_adj < 0.05) |>
  dplyr::pull(Feature)

dunn_clinical_hclust <- run_dunn_by_cluster(
  data = clinical_cluster_df,
  features = sig_clinical_hclust,
  cluster_col = "Hierarchical"
)
dunn_clinical_hclust

dunn_clinical_all <- dplyr::bind_rows(
  dunn_clinical_kmeans,
  dunn_clinical_hclust
)
dunn_clinical_all

# ----------- APRAŠOMOJI STATISTIKA --------------
summarise_features_by_cluster <- function(data, features, cluster_col) {
  data |>
    dplyr::select(
      Cluster = dplyr::all_of(cluster_col),
      dplyr::all_of(features)
    ) |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(features),
      names_to = "Feature",
      values_to = "Value"
    ) |>
    tidyr::drop_na(Value) |>
    dplyr::group_by(Cluster, Feature) |>
    dplyr::summarise(
      N = dplyr::n(),
      Median = median(Value, na.rm = TRUE),
      Q1 = quantile(Value, 0.25, na.rm = TRUE),
      Q3 = quantile(Value, 0.75, na.rm = TRUE),
      IQR = IQR(Value, na.rm = TRUE),
      Mean = mean(Value, na.rm = TRUE),
      SD = sd(Value, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      Median_IQR = paste0(
        round(Median, 3),
        " [",
        round(Q1, 3),
        "; ",
        round(Q3, 3),
        "]"
      ),
      Cluster_method = cluster_col
    ) |>
    dplyr::select(
      Cluster_method,
      Cluster,
      Feature,
      N,
      Median,
      Q1,
      Q3,
      IQR,
      Mean,
      SD,
      Median_IQR
    )
}

desc_texture_kmeans <- summarise_features_by_cluster(
  data = texture_cluster_df,
  features = texture_features,
  cluster_col = "Kmeans"
)

desc_texture_hclust <- summarise_features_by_cluster(
  data = texture_cluster_df,
  features = texture_features,
  cluster_col = "Hierarchical"
)

desc_texture_all <- dplyr::bind_rows(
  desc_texture_kmeans,
  desc_texture_hclust
)
desc_texture_all

desc_clinical_kmeans <- summarise_features_by_cluster(
  data = clinical_cluster_df,
  features = clinical_features,
  cluster_col = "Kmeans"
)

desc_clinical_hclust <- summarise_features_by_cluster(
  data = clinical_cluster_df,
  features = clinical_features,
  cluster_col = "Hierarchical"
)

desc_clinical_all <- dplyr::bind_rows(
  desc_clinical_kmeans,
  desc_clinical_hclust
)
desc_clinical_all

# ----------- Staciakampes diagramos ------------
texture_long_kmeans <- texture_cluster_df |>
  dplyr::select(Kmeans, dplyr::all_of(texture_features)) |>
  tidyr::pivot_longer(
    cols = dplyr::all_of(texture_features),
    names_to = "Feature",
    values_to = "Value"
  ) |>
  tidyr::drop_na(Value)

p_texture_boxplot_kmeans <- ggplot(
  texture_long_kmeans,
  aes(x = Kmeans, y = Value, fill = Kmeans)
) +
  geom_boxplot(
    alpha = 0.75,
    outlier.alpha = 0.6,
    outlier.size = 1.5,
    color = "black"
  ) +
  scale_fill_manual(values = custom_colors) +
  facet_wrap(
    ~ Feature,
    scales = "free_y"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 11, face = "bold", color = "black"),
    axis.text.x = element_text(size = 11, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title.x = element_text(size = 12, color = "black"),
    axis.title.y = element_text(size = 12, color = "black"),
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5)
  ) +
  labs(
    title = "Tekstūros požymių pasiskirstymas pagal K-means klasterius",
    subtitle = "Global struktūra, atrinkti tekstūros požymiai",
    x = "K-means klasteris",
    y = "Požymio reikšmė"
  )
p_texture_boxplot_kmeans

texture_long_hclust <- texture_cluster_df |>
  dplyr::select(Hierarchical, dplyr::all_of(texture_features)) |>
  tidyr::pivot_longer(
    cols = dplyr::all_of(texture_features),
    names_to = "Feature",
    values_to = "Value"
  ) |>
  tidyr::drop_na(Value)

p_texture_boxplot_hclust <- ggplot(
  texture_long_hclust,
  aes(x = Hierarchical, y = Value, fill = Hierarchical)
) +
  geom_boxplot(
    alpha = 0.75,
    outlier.alpha = 0.6,
    outlier.size = 1.5,
    color = "black"
  ) +
  scale_fill_manual(values = custom_colors) +
  facet_wrap(
    ~ Feature,
    scales = "free_y"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 11, face = "bold", color = "black"),
    axis.text.x = element_text(size = 11, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title.x = element_text(size = 12, color = "black"),
    axis.title.y = element_text(size = 12, color = "black"),
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5)
  ) +
  labs(
    title = "Tekstūros požymių pasiskirstymas pagal hierarchinius klasterius",
    subtitle = "Global struktūra, atrinkti tekstūros požymiai",
    x = "Hierarchinis klasteris",
    y = "Požymio reikšmė"
  )
p_texture_boxplot_hclust


# AS PACIENTŲ KLASTERIZAVIMAS PAGAL MEDICININES METRIKAS ------------
clinical_cluster_as_df <- cluster_class_df |>
  dplyr::filter(Class == "AS") |>
  dplyr::mutate(
    ID = readr::parse_number(Patient_ID)
  ) |>
  dplyr::left_join(
    pacientai,
    by = "ID"
  ) |>
  tidyr::drop_na(CVF, T1, ECV, GLS)


kw_clinical_as_kmeans <- run_kruskal_by_cluster(
  data = clinical_cluster_as_df,
  features = clinical_features,
  cluster_col = "Kmeans"
)

kw_clinical_as_hclust <- run_kruskal_by_cluster(
  data = clinical_cluster_as_df,
  features = clinical_features,
  cluster_col = "Hierarchical"
)

clinical_cluster_as_df |>
  dplyr::count(Kmeans)

clinical_cluster_as_df |>
  dplyr::count(Hierarchical)

run_kruskal_by_cluster <- function(data, features, cluster_col) {
  purrr::map_dfr(features, function(feature_name) {
    test_data <- data |>
      dplyr::transmute(
        Feature_value = .data[[feature_name]],
        Cluster = .data[[cluster_col]]
      ) |>
      tidyr::drop_na()
    n_clusters <- dplyr::n_distinct(test_data$Cluster)
    
    if (nrow(test_data) < 3 || n_clusters < 2) {
      return(
        tibble::tibble(
          Feature = feature_name,
          Cluster_method = cluster_col,
          statistic = NA_real_,
          df = NA_real_,
          p_value = NA_real_
        )
      )
    }
    test_result <- kruskal.test(
      Feature_value ~ Cluster,
      data = test_data
    )
    tibble::tibble(
      Feature = feature_name,
      Cluster_method = cluster_col,
      statistic = unname(test_result$statistic),
      df = unname(test_result$parameter),
      p_value = test_result$p.value
    )
  }) |>
    dplyr::mutate(
      p_adj = p.adjust(p_value, method = "BH"),
      Significant = dplyr::case_when(
        is.na(p_adj) ~ NA_character_,
        p_adj < 0.001 ~ "***",
        p_adj < 0.01 ~ "**",
        p_adj < 0.05 ~ "*",
        TRUE ~ "ns"
      )
    ) |>
    dplyr::arrange(p_adj, p_value)
}

kw_clinical_as_kmeans <- run_kruskal_by_cluster(
  data = clinical_cluster_as_df,
  features = clinical_features,
  cluster_col = "Kmeans"
)
kw_clinical_as_kmeans

kw_clinical_as_hclust <- run_kruskal_by_cluster(
  data = clinical_cluster_as_df,
  features = clinical_features,
  cluster_col = "Hierarchical"
)
kw_clinical_as_hclust

kw_clinical_as_all <- dplyr::bind_rows(
  kw_clinical_as_kmeans,
  kw_clinical_as_hclust
)
kw_clinical_as_all

# -------- DUNN ---------
run_dunn_by_cluster <- function(data, features, cluster_col) {
  purrr::map_dfr(features, function(feature_name) {
    test_data <- data |>
      dplyr::transmute(
        Feature_value = .data[[feature_name]],
        Cluster = .data[[cluster_col]]
      ) |>
      tidyr::drop_na()
    
    if (
      nrow(test_data) < 3 ||
      dplyr::n_distinct(test_data$Cluster) < 2
    ) {
      return(NULL)
    }
    
    dunn_result <- rstatix::dunn_test(
      data = test_data,
      Feature_value ~ Cluster,
      p.adjust.method = "BH"
    )
    
    dunn_result |>
      dplyr::mutate(
        Feature = feature_name,
        Cluster_method = cluster_col
      ) |>
      dplyr::select(
        Feature,
        Cluster_method,
        group1,
        group2,
        n1,
        n2,
        statistic,
        p,
        p.adj,
        p.adj.signif
      )
  })
}

sig_clinical_as_kmeans <- kw_clinical_as_kmeans |>
  dplyr::filter(!is.na(p_adj), p_adj < 0.05) |>
  dplyr::pull(Feature)

dunn_clinical_as_kmeans <- run_dunn_by_cluster(
  data = clinical_cluster_as_df,
  features = sig_clinical_as_kmeans,
  cluster_col = "Kmeans"
)

dunn_clinical_as_kmeans

sig_clinical_as_hclust <- kw_clinical_as_hclust |>
  dplyr::filter(!is.na(p_adj), p_adj < 0.05) |>
  dplyr::pull(Feature)

dunn_clinical_as_hclust <- run_dunn_by_cluster(
  data = clinical_cluster_as_df,
  features = sig_clinical_as_hclust,
  cluster_col = "Hierarchical"
)
dunn_clinical_as_hclust

dunn_clinical_as_all <- dplyr::bind_rows(
  dunn_clinical_as_kmeans,
  dunn_clinical_as_hclust
)
dunn_clinical_as_all

dunn_clinical_as_significant <- dunn_clinical_as_all |>
  dplyr::filter(!is.na(p.adj), p.adj < 0.05)
dunn_clinical_as_significant

# -------- aprasomoji statistika --------
summarise_features_by_cluster <- function(data, features, cluster_col) {
  data |>
    dplyr::select(
      Cluster = dplyr::all_of(cluster_col),
      dplyr::all_of(features)
    ) |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(features),
      names_to = "Feature",
      values_to = "Value"
    ) |>
    tidyr::drop_na(Value) |>
    dplyr::group_by(Cluster, Feature) |>
    dplyr::summarise(
      N = dplyr::n(),
      Median = median(Value, na.rm = TRUE),
      Q1 = quantile(Value, 0.25, na.rm = TRUE),
      Q3 = quantile(Value, 0.75, na.rm = TRUE),
      IQR = IQR(Value, na.rm = TRUE),
      Mean = mean(Value, na.rm = TRUE),
      SD = sd(Value, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      Median_IQR = paste0(
        round(Median, 3),
        " [",
        round(Q1, 3),
        "; ",
        round(Q3, 3),
        "]"
      ),
      Cluster_method = cluster_col
    ) |>
    dplyr::select(
      Cluster_method,
      Cluster,
      Feature,
      N,
      Median,
      Q1,
      Q3,
      IQR,
      Mean,
      SD,
      Median_IQR
    )
}

desc_clinical_as_kmeans <- summarise_features_by_cluster(
  data = clinical_cluster_as_df,
  features = clinical_features,
  cluster_col = "Kmeans"
)
desc_clinical_as_kmeans

desc_clinical_as_hclust <- summarise_features_by_cluster(
  data = clinical_cluster_as_df,
  features = clinical_features,
  cluster_col = "Hierarchical"
)
desc_clinical_as_hclust

# ---------- grafikai --------
clinical_long_kmeans <- clinical_cluster_as_df |>
  dplyr::select(Kmeans, dplyr::all_of(clinical_features)) |>
  tidyr::pivot_longer(
    cols = dplyr::all_of(clinical_features),
    names_to = "Feature",
    values_to = "Value"
  ) |>
  tidyr::drop_na(Value)

clinical_long_kmeans <- clinical_cluster_as_df |>
  dplyr::filter(Kmeans != 2) |>
  dplyr::mutate(Kmeans = droplevels(Kmeans)) |>
  dplyr::select(Kmeans, dplyr::all_of(clinical_features)) |>
  tidyr::pivot_longer(
    cols = dplyr::all_of(clinical_features),
    names_to = "Feature",
    values_to = "Value"
  ) |>
  tidyr::drop_na(Value)

p_clinical_as_boxplot_kmeans <- ggplot(
  clinical_long_kmeans,
  aes(x = Kmeans, y = Value, fill = Kmeans)
) +
  geom_boxplot(
    alpha = 0.75,
    outlier.alpha = 0.6,
    outlier.size = 1.5,
    color = "black"
  ) +
  scale_fill_manual(
    values = custom_colors,
    drop = TRUE
  ) +
  facet_wrap(
    ~ Feature,
    scales = "free_y"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 11, face = "bold", color = "black"),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title.x = element_text(size = 12, color = "black"),
    axis.title.y = element_text(size = 12, color = "black"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5)
  ) +
  labs(
    title = "AS pacientų medicininės metrikos pagal K-means tekstūros klasterius",
    subtitle = "Klasteriai gauti iš pilnos imties pagal tekstūros požymius",
    x = "K-means klasteris",
    y = "Reikšmė"
  )
p_clinical_as_boxplot_kmeans

clinical_long_hclust <- clinical_cluster_as_df |>
  dplyr::select(Hierarchical, dplyr::all_of(clinical_features)) |>
  tidyr::pivot_longer(
    cols = dplyr::all_of(clinical_features),
    names_to = "Feature",
    values_to = "Value"
  ) |>
  tidyr::drop_na(Value)

p_clinical_as_boxplot_hclust <- ggplot(
  clinical_long_hclust,
  aes(x = Hierarchical, y = Value, fill = Hierarchical)
) +
  geom_boxplot(
    alpha = 0.75,
    outlier.alpha = 0.6,
    outlier.size = 1.5,
    color = "black"
  ) +
  scale_fill_manual(values = custom_colors) +
  facet_wrap(
    ~ Feature,
    scales = "free_y"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 11, face = "bold", color = "black"),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title.x = element_text(size = 12, color = "black"),
    axis.title.y = element_text(size = 12, color = "black"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5)
  ) +
  labs(
    title = "AS pacientų medicininės metrikos pagal hierarchinius tekstūros klasterius",
    subtitle = "Klasteriai gauti iš pilnos imties pagal tekstūros požymius",
    x = "Hierarchinis klasteris",
    y = "Medicininės metrikos reikšmė"
  )
p_clinical_as_boxplot_hclust


# NUOKRYPIAI NUO NORMOS PAGAL TEKSTŪROS KLASTERIUS ------------
add_clinical_deviations <- function(data) {
  data |>
    dplyr::mutate(
      # CVF: virš 4 blogiau
      dev_CVF = dplyr::if_else(CVF > 4, CVF - 4, 0),
      
      # GLS: jei didesnis nei -20, blogiau
      dev_GLS = dplyr::if_else(GLS > -20, GLS + 20, 0),
      
      # T1 pagal lytį
      T1_upper = dplyr::case_when(
        sex == 1 ~ 1045,  # vyrai: 978 + 67
        sex == 2 ~ 969,   # moterys: 931 + 38
        TRUE ~ 1015       # bendras variantas: 956 + 59
      ),
      dev_T1 = dplyr::if_else(T1 > T1_upper, T1 - T1_upper, 0),
      
      # ECV pagal lytį
      ECV_upper = dplyr::case_when(
        sex == 1 ~ 29,    # vyrai: 26 + 3
        sex == 2 ~ 23,    # moterys: 22 + 1
        TRUE ~ 27         # bendra riba: 24 + 3
      ),
      dev_ECV = dplyr::if_else(ECV > ECV_upper, ECV - ECV_upper, 0)
    )
}

clinical_cluster_as_dev_df <- clinical_cluster_as_df |>
  add_clinical_deviations()
clinical_cluster_as_dev_df

clinical_cluster_as_dev_df |>
  dplyr::select(
    Patient_ID, sex, Kmeans, Hierarchical,
    CVF, GLS, T1, ECV,
    dev_CVF, dev_GLS, dev_T1, dev_ECV
  )
deviation_features <- c("dev_CVF", "dev_GLS", "dev_T1", "dev_ECV")

summarise_deviations_by_cluster <- function(data, cluster_col) {
  data |>
    dplyr::select(
      Cluster = dplyr::all_of(cluster_col),
      dplyr::all_of(deviation_features)
    ) |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(deviation_features),
      names_to = "Metric",
      values_to = "Deviation"
    ) |>
    tidyr::drop_na(Deviation) |>
    dplyr::group_by(Metric, Cluster) |>
    dplyr::summarise(
      n = dplyr::n(),
      mean = mean(Deviation, na.rm = TRUE),
      sd = sd(Deviation, na.rm = TRUE),
      median = median(Deviation, na.rm = TRUE),
      Q1 = quantile(Deviation, 0.25, na.rm = TRUE),
      Q3 = quantile(Deviation, 0.75, na.rm = TRUE),
      min = min(Deviation, na.rm = TRUE),
      max = max(Deviation, na.rm = TRUE),
      mad = mad(Deviation, constant = 1, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      Cluster_method = cluster_col,
      Median_IQR = paste0(
        round(median, 3),
        " [",
        round(Q1, 3),
        "; ",
        round(Q3, 3),
        "]"
      )
    ) |>
    dplyr::select(
      Cluster_method,
      Metric,
      Cluster,
      n,
      mean,
      sd,
      median,
      Q1,
      Q3,
      min,
      max,
      mad,
      Median_IQR
    )
}

# ------ k-means ------
cluster_summary_dev_kmeans <- summarise_deviations_by_cluster(
  data = clinical_cluster_as_dev_df,
  cluster_col = "Kmeans"
)
cluster_summary_dev_kmeans

# ------- hierarchinis -------
cluster_summary_dev_hclust <- summarise_deviations_by_cluster(
  data = clinical_cluster_as_dev_df,
  cluster_col = "Hierarchical"
)
cluster_summary_dev_hclust

cluster_summary_dev_all <- dplyr::bind_rows(
  cluster_summary_dev_kmeans,
  cluster_summary_dev_hclust
)
cluster_summary_dev_all


clinical_cluster_as_dev_df <- clinical_cluster_as_dev_df |>
  dplyr::filter(Kmeans != 2) |>
  dplyr::mutate(Kmeans = droplevels(Kmeans))

cluster_summary_dev_kmeans <- summarise_deviations_by_cluster(
  data = clinical_cluster_as_dev_df,
  cluster_col = "Kmeans"
)

# ------ Kruskall wallis -------
kw_deviation_kmeans <- run_kruskal_by_cluster(
  data = clinical_cluster_as_dev_df,
  features = deviation_features,
  cluster_col = "Kmeans"
)
kw_deviation_kmeans

kw_deviation_hclust <- run_kruskal_by_cluster(
  data = clinical_cluster_as_dev_df,
  features = deviation_features,
  cluster_col = "Hierarchical"
)
kw_deviation_hclust

sig_deviation_kmeans <- kw_deviation_kmeans |>
  dplyr::filter(!is.na(p_adj), p_adj < 0.05) |>
  dplyr::pull(Feature)

dunn_deviation_kmeans <- run_dunn_by_cluster(
  data = clinical_cluster_as_dev_df,
  features = sig_deviation_kmeans,
  cluster_col = "Kmeans"
)
dunn_deviation_kmeans

sig_deviation_hclust <- kw_deviation_hclust |>
  dplyr::filter(!is.na(p_adj), p_adj < 0.05) |>
  dplyr::pull(Feature)

dunn_deviation_hclust <- run_dunn_by_cluster(
  data = clinical_cluster_as_dev_df,
  features = sig_deviation_hclust,
  cluster_col = "Hierarchical"
)
dunn_deviation_hclust

# -------- MAD --------
plot_deviation_median_mad <- function(summary_data, cluster_method_label) {
  summary_data |>
    dplyr::filter(Metric %in% c("dev_CVF", "dev_GLS")) |>
    dplyr::mutate(
      Metric_label = dplyr::recode(
        Metric,
        dev_CVF = "CVF",
        dev_GLS = "GLS"
      )
    ) |>
    ggplot2::ggplot(
      ggplot2::aes(
        x = factor(Cluster),
        y = median
      )
    ) +
    ggplot2::geom_hline(
      yintercept = 0,
      linetype = "dashed",
      color = "black"
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = median - mad,
        ymax = median + mad
      ),
      width = 0.15,
      color = "grey40"
    ) +
    ggplot2::geom_point(
      size = 2.5,
      color = "violetred3"
    ) +
    ggplot2::facet_wrap(
      ~ Metric_label,
      scales = "free_y"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::labs(
      title = paste0(
        "AS pacientų medianinis nuokrypis nuo normos pagal ",
        cluster_method_label,
        " tekstūros klasterius"
      ),
      subtitle = "Taškai rodo medianą, paklaidos juostos – MAD",
      x = "Tekstūros klasteris",
      y = "Nuokrypis nuo normos"
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        face = "bold",
        size = 12
      ),
      plot.subtitle = ggplot2::element_text(
        hjust = 0.5,
        size = 10
      ),
      strip.text = ggplot2::element_text(
        size = 10,
        face = "bold"
      ),
      axis.text = ggplot2::element_text(color = "black"),
      axis.title = ggplot2::element_text(color = "black")
    )
}

p_deviation_kmeans <- plot_deviation_median_mad(
  summary_data = cluster_summary_dev_kmeans,
  cluster_method_label = "K-means"
)
p_deviation_kmeans

p_deviation_hclust <- plot_deviation_median_mad(
  summary_data = cluster_summary_dev_hclust,
  cluster_method_label = "hierarchinius"
)
p_deviation_hclust


# SPEARMANO KORELIACIJA TARP TEKSTŪROS POŽYMIŲ IR KLINIKINIŲ METRIKŲ ---------
spearman_df <- df_global_selected |>
  dplyr::filter(Class == 1) |>
  dplyr::mutate(
    ID = readr::parse_number(Patient_ID)
  ) |>
  dplyr::left_join(
    pacientai,
    by = "ID"
  ) |>
  tidyr::drop_na(CVF, T1, ECV, GLS)

texture_features_corr <- spearman_df |>
  dplyr::select(-Patient_ID, -Class, -ID, -sex, -CVF, -T1, -ECV, -GLS) |>
  names()
clinical_features_corr <- c("CVF", "T1", "ECV", "GLS")

run_spearman_texture_clinical <- function(data, texture_features, clinical_features) {
  expand.grid(
    Texture_feature = texture_features,
    Clinical_metric = clinical_features,
    stringsAsFactors = FALSE
  ) |>
    tibble::as_tibble() |>
    dplyr::rowwise() |>
    dplyr::mutate(
      N = sum(
        complete.cases(
          data[, c(Texture_feature, Clinical_metric)]
        )
      ),
      rho = cor(
        data[[Texture_feature]],
        data[[Clinical_metric]],
        method = "spearman",
        use = "complete.obs"
      ),
      p_value = cor.test(
        data[[Texture_feature]],
        data[[Clinical_metric]],
        method = "spearman",
        exact = FALSE
      )$p.value
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      p_adj = p.adjust(p_value, method = "BH"),
      Significant = dplyr::case_when(
        p_adj < 0.001 ~ "***",
        p_adj < 0.01  ~ "**",
        p_adj < 0.05  ~ "*",
        TRUE ~ "ns"
      )
    ) |>
    dplyr::arrange(p_adj, p_value)
}

spearman_texture_clinical_all <- run_spearman_texture_clinical(
  data = spearman_df,
  texture_features = texture_features_corr,
  clinical_features = clinical_features_corr
)
spearman_texture_clinical_all

spearman_texture_clinical_significant <- spearman_texture_clinical_all |>
  dplyr::filter(p_adj < 0.05) |>
  dplyr::mutate(
    rho = round(rho, 3),
    p_value = round(p_value, 4),
    p_adj = round(p_adj, 4)
  )

spearman_texture_clinical_significant