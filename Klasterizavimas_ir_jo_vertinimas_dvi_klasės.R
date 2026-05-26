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
library(forcats)

#--------------DUOMENU NUSKAITYMAS, POŽYMIAI IR MEDIKŲ METRIKOS ---------
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

#------------ POZYMIŲ PAVADINIMŲ TRUMPINIAI --------------------------------
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

# Grafikams ir lentelems, kad butu patogu lyginti tarpusavyje zonas
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

#------------ XGBOOST -------------------------
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
  
  # apsauga
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
  
  # Galutinis duomenu rinkinys
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

#------------ POZYMIU ATRANKA VISOMS STRUKTUROMS -------------------------
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

#------------ ATRINKTI DUOMENYS ATSKIROMS STRUKTUROMS -------------------------
df_global_selected <- feature_selection_results$global$final_data
df_bottom_selected <- feature_selection_results$bottom$final_data
df_mid_selected <- feature_selection_results$mid$final_data
df_top_selected <- feature_selection_results$top$final_data

#------------ BENDRA ATRINKTU POZYMIU LENTELE -------------------------
all_selected_by_method <- purrr::map_dfr(
  feature_selection_results,
  "selected_by_method"
)

all_feature_votes <- purrr::map_dfr(
  feature_selection_results,
  "feature_votes"
)

all_stable_features <- all_feature_votes |>
  dplyr::filter(N_methods >= 2) |>
  dplyr::arrange(Structure, dplyr::desc(N_methods), Feature)

all_selected_by_method
all_feature_votes
all_stable_features

structure_labels <- c(
  global = "Visa struktūra",
  bottom = "Bazinė zona",
  mid    = "Vidurinė zona",
  top    = "Apikalinė zona"
)

add_structure_labels <- function(data) {
  data |>
    dplyr::mutate(
      Structure_LT = dplyr::recode(
        Structure,
        !!!structure_labels
      ),
      Structure_LT = factor(
        Structure_LT,
        levels = c(
          "Visa struktūra",
          "Bazinė zona",
          "Vidurinė zona",
          "Apikalinė zona"
        )
      )
    )
}

all_feature_votes <- add_structure_labels(all_feature_votes)
all_selected_by_method <- add_structure_labels(all_selected_by_method)
all_stable_features <- add_structure_labels(all_stable_features)

#------------ POŽYMIŲ STABILUMAS TARP STRUKTŪRŲ -------------------------
feature_stability_table <- all_stable_features |>
  dplyr::distinct(Structure, Feature) |>
  dplyr::mutate(Selected = 1) |>
  tidyr::pivot_wider(
    names_from = Structure,
    values_from = Selected,
    values_fill = 0
  ) |>
  dplyr::mutate(
    N_structures = global + bottom + mid + top
  ) |>
  dplyr::arrange(dplyr::desc(N_structures), Feature)

feature_stability_table

#------------ POŽYMIAI, ATRINKTI BENT 2 METODAIS -------------------------
all_feature_votes_filtered <- all_feature_votes |>
  dplyr::filter(N_methods >= 2)

ggplot(
  all_feature_votes_filtered,
  aes(
    x = Structure_LT,
    y = reorder(Feature, N_methods),
    fill = N_methods
  )
) +
  geom_tile(color = "white") +
  labs(
    title = "Atrinkti stabilūs požymiai",
    x = "Struktūra",
    y = "Požymis"
  ) +
  theme_minimal() +
  scale_fill_gradient(
    low = "orange1",
    high = "violetred4",
    name = "Metodų\nskaičius",
    breaks = 1:4
  )+
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.title = element_text(size = 12),
    axis.text.x = element_text(angle = 30, hjust = 1, color = "black"),
    axis.text.y = element_text(size = 8, color = "black"),
    panel.grid = element_blank()
  )
# Pozymiu populiarumas tarp strukturu
feature_stability_across_structures <- feature_stability_table |>
  dplyr::filter(N_structures >= 2)

ggplot(feature_stability_across_structures, aes(x = reorder(Feature, N_structures), y = N_structures)) +
  geom_col(fill = "violetred4", width = 0.6) +
  coord_flip() +
  labs(
    title = "Požymių stabilumas tarp struktūrų",
    x = "Požymis",
    y = "Keliose struktūrose atrinktas"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 9, color = "black")
  )

# Grizimas prie pradiniu duomenu ----------
selected_datasets_raw <- purrr::imap(
  datasets,
  function(data, structure_name) {
    
    selected_features <- feature_selection_results[[structure_name]]$stable_features
    data |>
      dplyr::select(
        Patient_ID,
        Class,
        dplyr::all_of(selected_features)
      )
  }
)

df_global_selected_raw <- selected_datasets_raw$global
df_bottom_selected_raw <- selected_datasets_raw$bottom
df_mid_selected_raw    <- selected_datasets_raw$mid
df_top_selected_raw    <- selected_datasets_raw$top

# ATRINKTI DUOMENYS PCA / UMAP ANALIZEI ----------------
selected_datasets_scaled <- purrr::imap(
  datasets_scaled,
  function(data, structure_name) {
    selected_features <- feature_selection_results[[structure_name]]$stable_features
    data |>
      dplyr::select(
        Patient_ID,
        Class,
        dplyr::all_of(selected_features)
      )
  }
)

selected_feature_counts <- purrr::imap_dfr(
  selected_datasets_scaled,
  function(data, structure_name) {
    tibble::tibble(
      Structure = structure_name,
      Structure_LT = structure_labels[structure_name],
      N_features = ncol(data) - 2,
      N_patients = nrow(data)
    )
  }
)
selected_feature_counts


# BENDRA DUOMENŲ PARUOŠIMO FUNKCIJA --------------
prepare_embedding_data <- function(data) {
  
  data_clean <- data |>
    dplyr::mutate(
      Class = factor(
        Class,
        levels = c(0, 1),
        labels = c("Kontrolė", "Patologija")
      )
    ) |>
    tidyr::drop_na()
  
  feature_cols <- setdiff(names(data_clean), c("Patient_ID", "Class"))
  
  x <- data_clean |>
    dplyr::select(dplyr::all_of(feature_cols)) |>
    as.data.frame()
  
  non_zero_var <- sapply(x, function(col) stats::sd(col, na.rm = TRUE) > 0)
  x <- x[, non_zero_var, drop = FALSE]
  
  x_scaled <- as.matrix(x)
  
  meta <- data_clean |>
    dplyr::select(Patient_ID, Class)
  
  list(
    x = x_scaled,
    meta = meta,
    feature_cols = colnames(x)
  )
}


# TRUSTWORTHINESS IR CONTINUITY -------------
calculate_knn_ranks <- function(distance_matrix) {
  n <- nrow(distance_matrix)
  rank_matrix <- matrix(NA_integer_, nrow = n, ncol = n)
  for (i in seq_len(n)) {
    order_i <- order(distance_matrix[i, ])
    order_i <- order_i[order_i != i]
    rank_matrix[i, order_i] <- seq_along(order_i)
  }
  rank_matrix
}


calculate_trustworthiness_continuity <- function(x_high, x_low, k = 5) {
  x_high <- as.matrix(x_high)
  x_low  <- as.matrix(x_low)
  
  if (nrow(x_high) != nrow(x_low)) {
    stop("x_high and x_low must have the same number of rows.")
  }
  
  if (k < 1) {
    stop("k must be at least 1.")
  }
  n <- nrow(x_high)
  
  if (n <= 3) {
    return(
      tibble::tibble(
        k = k,
        Trustworthiness = NA_real_,
        Continuity = NA_real_
      )
    )
  }
  max_k <- floor((n - 1) / 2)
  if (k > max_k) {
    warning("k was reduced to floor((n - 1) / 2).")
    k <- max_k
  }
  
  dist_high <- as.matrix(dist(x_high))
  dist_low  <- as.matrix(dist(x_low))
  rank_high <- calculate_knn_ranks(dist_high)
  rank_low  <- calculate_knn_ranks(dist_low)
  trust_sum <- 0
  cont_sum <- 0
  
  for (i in seq_len(n)) {
    high_neighbors <- which(rank_high[i, ] <= k)
    low_neighbors  <- which(rank_low[i, ] <= k)
    unexpected_neighbors <- setdiff(low_neighbors, high_neighbors)
    missing_neighbors    <- setdiff(high_neighbors, low_neighbors)
    trust_sum <- trust_sum +
      sum(rank_high[i, unexpected_neighbors] - k, na.rm = TRUE)
    cont_sum <- cont_sum +
      sum(rank_low[i, missing_neighbors] - k, na.rm = TRUE)
  }
  
  denominator <- n * k * (2 * n - 3 * k - 1)
  tibble::tibble(
    k = k,
    Trustworthiness = 1 - (2 / denominator) * trust_sum,
    Continuity = 1 - (2 / denominator) * cont_sum
  )
}

calculate_stress <- function(x_high, x_low) {
  x_high <- as.matrix(x_high)
  x_low  <- as.matrix(x_low)
  
  if (nrow(x_high) != nrow(x_low)) {
    stop("x_high and x_low must have the same number of rows.")
  }
  if (nrow(x_high) <= 2) {
    return(
      tibble::tibble(
        Stress = NA_real_
      )
    )
  }
  
  dist_high <- as.matrix(dist(x_high))
  dist_low  <- as.matrix(dist(x_low))
  upper_idx <- upper.tri(dist_high)
  d_high <- dist_high[upper_idx]
  d_low  <- dist_low[upper_idx]
  numerator <- sum((d_high - d_low)^2, na.rm = TRUE)
  denominator <- sum(d_high^2, na.rm = TRUE)
  stress_value <- sqrt(numerator / denominator)
  tibble::tibble(
    Stress = stress_value
  )
}

add_structure_info <- function(data, structure_name, method_name = NULL) {
  data <- data |>
    dplyr::mutate(
      Structure = structure_name,
      Structure_LT = structure_labels[structure_name]
    )
  
  if (!is.null(method_name)) {
    data <- data |>
      dplyr::mutate(Method = method_name)
  }
  data
}

# PCA -----------------
run_pca_embedding <- function(data, structure_name) {
  prepared <- prepare_embedding_data(data)
  pca_model <- prcomp(
    prepared$x,
    center = FALSE,
    scale. = FALSE
  )
  
  explained_var <- summary(pca_model)$importance[2, 1:2] * 100
  embedding <- as.data.frame(pca_model$x[, 1:2]) |>
    dplyr::rename(
      Dim1 = PC1,
      Dim2 = PC2
    ) |>
    dplyr::bind_cols(prepared$meta) |>
    dplyr::mutate(
      Method = "PCA",
      Structure = structure_name,
      Structure_LT = structure_labels[structure_name],
      Dim1_label = paste0("PC1 (", round(explained_var[1], 1), "%)"),
      Dim2_label = paste0("PC2 (", round(explained_var[2], 1), "%)")
    )
  
  list(
    embedding = embedding,
    model = pca_model,
    original_x = prepared$x,
    low_x = pca_model$x[, 1:2],
    explained_var = explained_var
  )
}

# UMAP ---------------
run_umap_embedding <- function(
    data,
    structure_name,
    n_neighbors = 15,
    min_dist = 0.01,
    metric = "euclidean",
    seed = 123,
    method_name = "UMAP"
) {
  
  prepared <- prepare_embedding_data(data)
  n <- nrow(prepared$x)
  
  n_neighbors <- min(n_neighbors, n - 1)
  
  set.seed(seed)
  umap_matrix <- uwot::umap(
    prepared$x,
    n_neighbors = n_neighbors,
    min_dist = min_dist,
    metric = metric,
    n_components = 2,
    scale = FALSE,
    verbose = FALSE
  )
  
  embedding <- as.data.frame(umap_matrix) |>
    dplyr::rename(
      Dim1 = V1,
      Dim2 = V2
    ) |>
    dplyr::bind_cols(prepared$meta) |>
    dplyr::mutate(
      Method = method_name,
      Structure = structure_name,
      Structure_LT = structure_labels[structure_name],
      Dim1_label = "UMAP1",
      Dim2_label = "UMAP2",
      n_neighbors = n_neighbors,
      min_dist = min_dist,
      metric = metric
    )
  
  list(
    embedding = embedding,
    original_x = prepared$x,
    low_x = umap_matrix,
    n_neighbors = n_neighbors,
    min_dist = min_dist,
    metric = metric
  )
}

# UMAP PARAMETRŲ PARINKIMAS -------------
tune_umap_parameters <- function(
    data,
    structure_name,
    n_neighbors_grid = c(5, 10, 15),
    min_dist_grid = c(0.01, 0.05, 0.1),
    metric_grid = c("euclidean"),
    k_quality = 5,
    seed = 123
) {
  
  prepared <- prepare_embedding_data(data)
  x <- prepared$x
  n <- nrow(x)
  param_grid <- tidyr::expand_grid(
    n_neighbors = n_neighbors_grid,
    min_dist = min_dist_grid,
    metric = metric_grid
  ) |>
    dplyr::filter(n_neighbors < n)
  
  purrr::pmap_dfr(
    param_grid,
    function(n_neighbors, min_dist, metric) {
      set.seed(seed)
      umap_matrix <- tryCatch(
        uwot::umap(
          x,
          n_neighbors = n_neighbors,
          min_dist = min_dist,
          metric = metric,
          n_components = 2,
          scale = FALSE,
          verbose = FALSE
        ),
        error = function(e) NULL
      )
      
      if (is.null(umap_matrix)) {
        return(
          tibble::tibble(
            Structure = structure_name,
            Structure_LT = structure_labels[structure_name],
            n_neighbors = n_neighbors,
            min_dist = min_dist,
            metric = metric,
            k = k_quality,
            Trustworthiness = NA_real_,
            Continuity = NA_real_,
            Quality_mean = NA_real_
          )
        )
      }
      
      quality <- calculate_trustworthiness_continuity(
        x_high = x,
        x_low = umap_matrix,
        k = k_quality
      )
      
      quality |>
        dplyr::mutate(
          Structure = structure_name,
          Structure_LT = structure_labels[structure_name],
          n_neighbors = n_neighbors,
          min_dist = min_dist,
          metric = metric,
          Quality_mean = mean(c(Trustworthiness, Continuity), na.rm = TRUE)
        )
    }
  ) |>
    dplyr::arrange(
      dplyr::desc(Quality_mean),
      dplyr::desc(Trustworthiness),
      dplyr::desc(Continuity)
    )
}

# ----------- PCA VISOMS STRUKTŪROMS ----------------
pca_results <- purrr::imap(
  selected_datasets_scaled,
  run_pca_embedding
)

pca_embeddings <- purrr::map_dfr(
  pca_results,
  "embedding"
)

# --------- DEFAULT UMAP VISOMS STRUKTŪROMS -------------
umap_results_default <- purrr::imap(
  selected_datasets_scaled,
  function(data, structure_name) {
    run_umap_embedding(
      data = data,
      structure_name = structure_name,
      n_neighbors = 15,
      min_dist = 0.01,
      metric = "euclidean",
      seed = 123,
      method_name = "UMAP_default"
    )
  }
)

umap_embeddings_default <- purrr::map_dfr(
  umap_results_default,
  "embedding"
)

# ------------- UMAP PARAMETRŲ PARINKIMAS VISOMS STRUKTŪROMS --------------
umap_tuning_results <- purrr::imap_dfr(
  selected_datasets_scaled,
  function(data, structure_name) {
    tune_umap_parameters(
      data = data,
      structure_name = structure_name,
      n_neighbors_grid = c(5,10,15),
      min_dist_grid = c(0.01,0.05,0.1),
      metric_grid = c("euclidean"),
      k_quality = 5,
      seed = 123
    )
  }
)

best_umap_params_by_structure <- umap_tuning_results |>
  dplyr::group_by(Structure, Structure_LT) |>
  dplyr::slice_max(
    order_by = Quality_mean,
    n = 1,
    with_ties = FALSE
  ) |>
  dplyr::ungroup()

best_umap_params_by_structure

umap_params_subtitle <- best_umap_params_by_structure |>
  dplyr::mutate(
    param_text = paste0(
      "n_neighbors = ", n_neighbors,
      ", min_dist = ", min_dist,
      ", metric = ", metric
    )
  ) |>
  dplyr::pull(param_text) |>
  paste(collapse = "; ")

umap_params_subtitle

# ------------ TUNED UMAP VISOMS STRUKTŪROMS ---------------
umap_results_tuned <- purrr::imap(
  selected_datasets_scaled,
  function(data, structure_name) {
    
    params <- best_umap_params_by_structure |>
      dplyr::filter(Structure == structure_name)
    
    run_umap_embedding(
      data = data,
      structure_name = structure_name,
      n_neighbors = params$n_neighbors[1],
      min_dist = params$min_dist[1],
      metric = params$metric[1],
      seed = 123,
      method_name = "UMAP_tuned"
    )
  }
)

umap_embeddings_tuned <- purrr::map_dfr(
  umap_results_tuned,
  "embedding"
)

umap_embeddings_tuned_labeled <- umap_embeddings_tuned |>
  dplyr::mutate(
    Facet_label = paste0(
      Structure_LT,
      "\n",
      "n = ", n_neighbors,
      ", min_d = ", min_dist,
      ", metric = ", metric
    )
  )

# ---------- PROJEKCIJŲ KOKYBĖS METRIKOS ----------------
calculate_embedding_quality <- function(results_list, method_name, k = 5) {
  purrr::imap_dfr(
    results_list,
    
    function(result, structure_name) {
      trust_cont <- calculate_trustworthiness_continuity(
        x_high = result$original_x,
        x_low = result$low_x,
        k = k
      )
      stress <- calculate_stress(
        x_high = result$original_x,
        x_low = result$low_x
      )
      dplyr::bind_cols(trust_cont, stress) |>
        dplyr::mutate(
          Method = method_name,
          Structure = structure_name,
          Structure_LT = structure_labels[structure_name]
        )
    }
  )
}

embedding_result_sets <- list(
  PCA = pca_results,
  UMAP_default = umap_results_default,
  UMAP_tuned = umap_results_tuned
)

embedding_quality_metrics <- purrr::imap_dfr(
  embedding_result_sets,
  function(results_list, method_name) {
    calculate_embedding_quality(
      results_list = results_list,
      method_name = method_name,
      k = 5
    )
  }
) |>
  dplyr::select(
    Structure,
    Structure_LT,
    Method,
    k,
    Trustworthiness,
    Continuity,
    Stress
  ) |>
  dplyr::arrange(
    Structure_LT,
    Method
  )

embedding_quality_metrics

# -------------- PCA Sklaida --------------
pca_variance_by_structure <- purrr::imap_dfr(
  pca_results,
  function(result, structure_name) {
    tibble::tibble(
      Structure = structure_name,
      Structure_LT = structure_labels[structure_name],
      PC1_percent = result$explained_var[1],
      PC2_percent = result$explained_var[2],
      PC1_PC2_total_percent = sum(result$explained_var[1:2])
    )
  }
)

pca_variance_by_structure

pca_variance_long <- pca_variance_by_structure |>
  tidyr::pivot_longer(
    cols = c(PC1_percent, PC2_percent),
    names_to = "Component",
    values_to = "Explained_percent"
  ) |>
  dplyr::mutate(
    Component = dplyr::recode(
      Component,
      PC1_percent = "PC1",
      PC2_percent = "PC2"
    )
  )

pca_variance_plot <- ggplot(
  pca_variance_long,
  aes(x = Structure_LT, y = Explained_percent, fill = Component)
) +
  geom_col(position = "stack", width = 0.5) +
  geom_text(
    aes(label = paste0(round(Explained_percent, 1), "%")),
    position = position_stack(vjust = 0.5),
    color = "white",
    size = 3
  ) +
  scale_fill_manual(
    values = c(
      "PC1" = "violetred1",
      "PC2" = "violetred4"
    )
  ) +
  labs(
    title = "PCA pirmųjų dviejų komponenčių paaiškinamoji variacija",
    x = "Struktūra",
    y = "Paaiškinamoji variacija (%)",
    fill = "Komponentė"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text.x = element_text(angle=30, hjust = 1, color = "black", size=10),
    axis.text.y = element_text(color = "black", size=10)
  )
pca_variance_plot


# ----------- BENDRA PROJEKCIJŲ GRAFIKO FUNKCIJA ---------------
plot_embedding_facets <- function(
    embedding_data,
    plot_title,
    x_lab,
    y_lab,
    add_ellipse = FALSE
) {
  
  facet_var <- if ("Facet_label" %in% names(embedding_data)) {
    "Facet_label"
  } else {
    "Structure_LT"
  }
  
  p <- ggplot(
    embedding_data,
    aes(x = Dim1, y = Dim2, color = Class)
  ) +
    geom_point(size = 2, alpha = 0.85) +
    facet_wrap(
      stats::as.formula(paste("~", facet_var)),
      scales = "free"
    ) +
    scale_color_manual(
      values = c(
        "Kontrolė" = "royalblue2",
        "Patologija" = "violetred3"
      ),
      name = "Grupė"
    ) +
    labs(
      title = plot_title,
      x = x_lab,
      y = y_lab
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      strip.text = element_text(face = "plain", size = 10),
      axis.text = element_text(color = "black"),
      legend.position = "right"
    )
  
  if (add_ellipse) {
    p <- p +
      stat_ellipse(
        aes(group = Class),
        type = "norm",
        linewidth = 0.7,
        alpha = 0.7,
        show.legend = FALSE
      )
  }
  p
}

# ---------------- PCA GRAFIKAI ----------------
pca_plot <- plot_embedding_facets(
  embedding_data = pca_embeddings,
  plot_title = "PCA projekcija pagal atrinktus tekstūros požymius",
  x_lab = "PC1",
  y_lab = "PC2",
  add_ellipse = FALSE
)
pca_plot

extract_pca_loadings <- function(pca_results, top_n = 5) {
  purrr::imap_dfr(
    pca_results,
    function(result, structure_name) {
      loadings <- as.data.frame(result$model$rotation[, 1:2]) |>
        tibble::rownames_to_column("Feature") |>
        dplyr::rename(
          PC1_loading = PC1,
          PC2_loading = PC2
        ) |>
        dplyr::mutate(
          Structure = structure_name,
          Structure_LT = structure_labels[structure_name],
          PC1_abs = abs(PC1_loading),
          PC2_abs = abs(PC2_loading)
        )
      pc1_top <- loadings |>
        dplyr::arrange(dplyr::desc(PC1_abs)) |>
        dplyr::slice_head(n = top_n) |>
        dplyr::mutate(Component = "PC1", Abs_loading = PC1_abs)
      
      pc2_top <- loadings |>
        dplyr::arrange(dplyr::desc(PC2_abs)) |>
        dplyr::slice_head(n = top_n) |>
        dplyr::mutate(Component = "PC2", Abs_loading = PC2_abs)
      dplyr::bind_rows(pc1_top, pc2_top)
    }
  ) |>
    dplyr::select(
      Structure,
      Structure_LT,
      Component,
      Feature,
      PC1_loading,
      PC2_loading,
      Abs_loading
    ) |>
    dplyr::arrange(
      Structure_LT,
      Component,
      dplyr::desc(Abs_loading)
    )
}
pca_top_loadings <- extract_pca_loadings(
  pca_results = pca_results,
  top_n = 5
)

pca_top_loadings

pca_top_loadings_plot <- pca_top_loadings |>
  dplyr::mutate(
    Feature = forcats::fct_reorder(Feature, Abs_loading)
  )

ggplot(
  pca_top_loadings_plot,
  aes(
    x = Feature,
    y = Abs_loading,
    fill = Structure_LT
  )
) +
  geom_col(width = 0.65) +
  coord_flip() +
  facet_wrap(~ Component, scales = "free_y") +
  scale_fill_manual(
    values = c(
      "Visa struktūra" = "gray40",
      "Bazinė zona" = "royalblue2",
      "Vidurinė zona" = "orange2",
      "Apikalinė zona" = "violetred3"
    )
  ) +
  labs(
    title = "Svarbiausi tekstūros požymių svoriai PCA komponentėse",
    x = "Požymis",
    y = "Absoliuti svorio reikšmė",
    fill = "Struktūra"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text = element_text(color = "black"),
    strip.text = element_text(face = "bold", size = 11),
    legend.position = "bottom"
  )

# --------------- UMAP GRAFIKAI -----------
umap_plot_default <- plot_embedding_facets(
  embedding_data = umap_embeddings_default,
  plot_title = "UMAP projekcija su pradiniais parametrais",
  x_lab = "UMAP1",
  y_lab = "UMAP2",
  add_ellipse = FALSE
)

umap_plot_tuned <- plot_embedding_facets(
  embedding_data = umap_embeddings_tuned_labeled,
  plot_title = "UMAP projekcija su parinktais parametrais",
  x_lab = "UMAP1",
  y_lab = "UMAP2",
  add_ellipse = FALSE
)

umap_plot_default
umap_plot_tuned

# ---------- GRUPIŲ ATSKYRIMO ĮVERTINIMAS PROJEKCIJOJE: SILUETAS -----------
calculate_class_separation <- function(embedding_data) {
  embedding_data |>
    dplyr::group_by(Method, Structure, Structure_LT) |>
    dplyr::group_modify(function(df, key) {
      
      if (length(unique(df$Class)) < 2 || nrow(df) < 4) {
        return(
          tibble::tibble(
            Silhouette_mean = NA_real_
          )
        )
      }
      
      coords <- df |>
        dplyr::select(Dim1, Dim2) |>
        as.matrix()
      class_numeric <- as.numeric(df$Class)
      
      sil <- cluster::silhouette(
        class_numeric,
        dist(coords)
      )
      tibble::tibble(
        Silhouette_mean = mean(sil[, "sil_width"], na.rm = TRUE)
      )
    }) |>
    dplyr::ungroup() |>
    dplyr::arrange(dplyr::desc(Silhouette_mean))
}

all_embeddings_for_separation <- dplyr::bind_rows(
  pca_embeddings,
  umap_embeddings_default,
  umap_embeddings_tuned
)

class_separation_results <- calculate_class_separation(
  all_embeddings_for_separation
)
class_separation_results

# ------------ GALUTINĖ PROJEKCIJŲ PALYGINIMO LENTELĖ --------------
projection_summary <- embedding_quality_metrics |>
  dplyr::mutate(
    Projection_quality_mean = mean(c(Trustworthiness, Continuity), na.rm = TRUE),
    Stress_quality = 1 - Stress,
    Projection_quality_with_stress = mean(
      c(Trustworthiness, Continuity, Stress_quality),
      na.rm = TRUE
    ),
    .by = c(Structure, Structure_LT, Method)
  ) |>
  dplyr::left_join(
    class_separation_results,
    by = c("Structure", "Structure_LT", "Method")
  ) |>
  dplyr::arrange(
    Structure_LT,
    dplyr::desc(Projection_quality_with_stress),
    dplyr::desc(Silhouette_mean)
  )

projection_summary



