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

# medicininiu metriku apjungimas su teksturos pozymiais
datasets_with_clinical <- datasets |>
  purrr::map(\(data) add_clinical_metrics(data, pacientai))

datasets_scaled_with_clinical <- datasets_scaled |>
  purrr::map(\(data) add_clinical_metrics(data, pacientai))

datasets_scaled_with_clinical_scaled <- datasets_scaled |>
  purrr::map(\(data) add_clinical_metrics(data, pacientai_scaled))

# bendros lenteles
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

#------------ XGBOOST  -------------------------
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

#------------ PILNA POZYMIU ATRANKA VIENAI STRUKTURAI -------------------------
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
  
  # Pozymiai, kuriuos atrinko bent 2 metodai
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

#------------ POZYMIU STABILUMAS TARP STRUKTURU -------------------------
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

#------------ POZYMIAI, ATRINKTI BENT 2 METODAIS -------------------------
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

# ----------- KLASIFIKAVIMAS SU GAUTAIS REIKSMINGAIS POZYMIAIS ----------------
# Duomenu rinkiniai --------
selected_datasets <- list(
  global = df_global_selected_raw,
  bottom = df_bottom_selected_raw,
  mid    = df_mid_selected_raw,
  top    = df_top_selected_raw
)

structure_labels <- c(
  global = "Visa struktūra",
  bottom = "Bazinė zona",
  mid    = "Vidurinė zona",
  top    = "Apikalinė zona"
)

set.seed(123)
# -----------Duomenu paruošimas klasifikavimui ------------
prepare_classification_data <- function(data) {
  
  data |>
    dplyr::select(-Patient_ID) |>
    dplyr::mutate(
      Class = factor(
        Class,
        levels = c(0, 1),
        labels = c("Control", "AS")
      )
    ) |>
    tidyr::drop_na()
}

selected_datasets_prepared <- selected_datasets |>
  purrr::map(prepare_classification_data)

# Klasiu pasiskirstymo patikrinimas
class_distribution_selected <- selected_datasets_prepared |>
  purrr::imap_dfr(\(data, structure_name) {
    data |>
      dplyr::count(Class) |>
      dplyr::mutate(
        Structure = structure_name,
        Structure_LT = structure_labels[structure_name]
      )
  })
class_distribution_selected

# ------------------ KLASIFIKAVIMO METRIKU VERTINIMAS -------------
macro_f1_score <- function(obs, pred, lev = NULL) {
  
  if (is.null(lev)) {
    lev <- levels(obs)
  }
  obs <- factor(obs, levels = lev)
  pred <- factor(pred, levels = lev)
  
  f1_values <- purrr::map_dbl(lev, function(class_name) {
    
    tp <- sum(pred == class_name & obs == class_name, na.rm = TRUE)
    fp <- sum(pred == class_name & obs != class_name, na.rm = TRUE)
    fn <- sum(pred != class_name & obs == class_name, na.rm = TRUE)
    
    precision <- ifelse((tp + fp) == 0, 0, tp / (tp + fp))
    recall    <- ifelse((tp + fn) == 0, 0, tp / (tp + fn))
    ifelse(
      (precision + recall) == 0,
      0,
      2 * precision * recall / (precision + recall)
    )
  })
  
  mean(f1_values, na.rm = TRUE)
}


classification_summary <- function(data, lev = NULL, model = NULL) {
  negative_class <- "Control"
  positive_class <- "AS"
  
  data$obs  <- factor(data$obs, levels = c(negative_class, positive_class))
  data$pred <- factor(data$pred, levels = c(negative_class, positive_class))
  cm_table <- table(
    Prediction = data$pred,
    Reference = data$obs
  )
  
  TP <- cm_table[positive_class, positive_class]
  TN <- cm_table[negative_class, negative_class]
  FP <- cm_table[positive_class, negative_class]
  FN <- cm_table[negative_class, positive_class]
  
  sensitivity <- ifelse((TP + FN) == 0, NA_real_, TP / (TP + FN))
  specificity <- ifelse((TN + FP) == 0, NA_real_, TN / (TN + FP))
  balanced_accuracy <- mean(
    c(sensitivity, specificity),
    na.rm = TRUE
  )
  
  macro_f1 <- macro_f1_score(
    obs = data$obs,
    pred = data$pred,
    lev = c(negative_class, positive_class)
  )
  auc_value <- NA_real_
  
  if (positive_class %in% colnames(data)) {
    auc_value <- tryCatch(
      {
        roc_obj <- pROC::roc(
          response = data$obs,
          predictor = data[[positive_class]],
          levels = c(negative_class, positive_class),
          direction = "<",
          quiet = TRUE
        )
        as.numeric(pROC::auc(roc_obj))
      },
      error = function(e) NA_real_
    )
  }
  
  out <- c(
    Macro_F1 = macro_f1,
    Balanced_Accuracy = balanced_accuracy,
    Sensitivity = sensitivity,
    Specificity = specificity,
    AUC = auc_value
  )
  out[is.nan(out)] <- NA_real_
  out
}

# --------------- Kryzminis validavimas -----------------
set.seed(123)

cv_control <- caret::trainControl(
  method = "repeatedcv",
  number = 5,
  repeats = 10,
  classProbs = TRUE,
  summaryFunction = classification_summary,
  savePredictions = "final",
  sampling = "up",
  verboseIter = FALSE
)

#--------------- Modeliu treniravimas ------------------------
train_models_for_structure <- function(data, structure_name) {
  
  message("Treniruojami modeliai struktūrai: ", structure_name)
  
  models <- list()
  
  # ---------------- SVM Linear ----------------
  set.seed(123)
  models$SVM_Linear <- caret::train(
    Class ~ .,
    data = data,
    method = "svmLinear",
    trControl = cv_control,
    metric = "Balanced_Accuracy",
    preProcess = c("zv", "nzv", "range"),
    tuneLength = 5
  )
  
  # ---------------- SVM RBF ----------------
  set.seed(123)
  models$SVM_RBF <- caret::train(
    Class ~ .,
    data = data,
    method = "svmRadial",
    trControl = cv_control,
    metric = "Balanced_Accuracy",
    preProcess = c("zv", "nzv", "range"),
    tuneLength = 5
  )
  
  # ---------------- Random Forest ----------------
  set.seed(123)
  models$Random_Forest <- caret::train(
    Class ~ .,
    data = data,
    method = "rf",
    trControl = cv_control,
    metric = "Balanced_Accuracy",
    tuneLength = 5,
    ntree = 500
  )
  
  # ---------------- Naive Bayes ----------------
  nb_grid <- expand.grid(
    fL = c(0, 1),
    usekernel = c(TRUE),
    adjust = c(1, 1.5, 2)
  )
  
  set.seed(123)
  models$Naive_Bayes <- caret::train(
    Class ~ .,
    data = data,
    method = "nb",
    trControl = cv_control,
    metric = "Balanced_Accuracy",
    preProcess = c("zv", "nzv", "range"),
    tuneGrid = nb_grid
  )
  
  # ---------------- XGBoost ----------------
  xgb_grid <- expand.grid(
    nrounds = c(30, 50),
    max_depth = c(1, 2),
    eta = c(0.05, 0.10),
    gamma = c(0, 0.5, 1),
    colsample_bytree = c(0.7, 0.8, 0.9),
    min_child_weight = c(3, 5),
    subsample = c(0.8)
  )
  
  set.seed(123)
  models$XGBoost <- caret::train(
    Class ~ .,
    data = data,
    method = "xgbTree",
    trControl = cv_control,
    metric = "Balanced_Accuracy",
    tuneGrid = xgb_grid,
    verbose = FALSE,
    verbosity = 0
  )
  
  # ---------------- k-NN ----------------
  set.seed(123)
  models$KNN <- caret::train(
    Class ~ .,
    data = data,
    method = "knn",
    trControl = cv_control,
    metric = "Balanced_Accuracy",
    preProcess = c("zv", "nzv", "range"),
    tuneLength = 10
  )
  models
}

# Modeliu treniravimas visoms strukturoms ------------
classification_models <- purrr::imap(
  selected_datasets_prepared,
  ~ train_models_for_structure(
    data = .x,
    structure_name = .y
  )
)

# Metriku istraukimas is geriausio hiperparametru rinkinio -----------------
extract_best_metrics <- function(model, structure_name, model_name) {
  
  best_tune <- model$bestTune
  results <- model$results
  for (param_name in names(best_tune)) {
    results <- results |>
      dplyr::filter(.data[[param_name]] == best_tune[[param_name]])
  }
  
  results |>
    dplyr::slice(1) |>
    dplyr::mutate(
      Structure = structure_name,
      Structure_LT = structure_labels[structure_name],
      Model = model_name
    ) |>
    dplyr::select(
      Structure,
      Structure_LT,
      Model,
      Macro_F1,
      Balanced_Accuracy,
      Sensitivity,
      Specificity,
      AUC,
      dplyr::everything()
    )
}

classification_metrics <- purrr::imap_dfr(
  classification_models,
  function(models_for_structure, structure_name) {
    purrr::imap_dfr(
      models_for_structure,
      function(model, model_name) {
        extract_best_metrics(
          model = model,
          structure_name = structure_name,
          model_name = model_name
        )
      }
    )
  }
)

classification_metrics <- classification_metrics |>
  dplyr::arrange(
    Structure,
    dplyr::desc(Balanced_Accuracy),
    dplyr::desc(Macro_F1),
    dplyr::desc(AUC)
  ) |>
  dplyr::select(
    Structure,
    Structure_LT,
    Model,
    Macro_F1,
    Sensitivity,
    Specificity,
    AUC,
    Balanced_Accuracy
  )
classification_metrics

# Geriausias modelis kiekvienai strukturai -------------------
klasifik <- classification_metrics |>
  dplyr::select(
    Structure,
    Structure_LT,
    Model,
    Macro_F1,
    Sensitivity,
    Specificity,
    AUC,
    Balanced_Accuracy
  ) |>
  dplyr::arrange(
    Structure_LT,
    dplyr::desc(Balanced_Accuracy),
    dplyr::desc(Macro_F1),
    dplyr::desc(AUC)
  )

classification_metrics_clean <- classification_metrics |>
  dplyr::group_by(Structure, Structure_LT) |>
  dplyr::select(
    Structure_LT,
    Model,
    Macro_F1,
    Sensitivity,
    Specificity,
    AUC,
    Balanced_Accuracy
  ) |>
  dplyr::arrange(
    Structure_LT,
    dplyr::desc(Balanced_Accuracy),
    dplyr::desc(Macro_F1),
    dplyr::desc(AUC)
  )|>
  dplyr::slice(1) |>
  dplyr::ungroup()


classification_metrics_clean


# Paziurime parametrus -----------------
classification_models$top$Naive_Bayes$bestTune
classification_models$bottom$Naive_Bayes$bestTune
classification_models$mid$Naive_Bayes$bestTune
classification_models$global$KNN$bestTune


# ------------- Klasifikavimas vienam pacientui --------------
get_patient_level_confusion_matrix <- function(
    model,
    data_with_ids,
    positive_class = "AS",
    threshold = 0.5
) {
  
  preds <- model$pred
  best_tune <- model$bestTune
  for (param_name in names(best_tune)) {
    if (param_name %in% names(preds)) {
      preds <- preds |>
        dplyr::filter(.data[[param_name]] == best_tune[[param_name]])
    }
  }
  
  id_lookup <- data_with_ids |>
    dplyr::mutate(rowIndex = dplyr::row_number()) |>
    dplyr::select(rowIndex, Patient_ID)
  
  # Viena prognoze vienam pacientui
  patient_predictions <- preds |>
    dplyr::left_join(id_lookup, by = "rowIndex") |>
    dplyr::group_by(rowIndex, Patient_ID, obs) |>
    dplyr::summarise(
      Mean_AS_Probability = mean(.data[[positive_class]], na.rm = TRUE),
      N_predictions = dplyr::n(),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      pred = ifelse(
        Mean_AS_Probability >= threshold,
        positive_class,
        "Control"
      ),
      pred = factor(pred, levels = c("Control", "AS")),
      obs = factor(obs, levels = c("Control", "AS"))
    )
  
  cm <- caret::confusionMatrix(
    data = patient_predictions$pred,
    reference = patient_predictions$obs,
    positive = positive_class
  )
  
  list(
    patient_predictions = patient_predictions,
    confusion_matrix = cm
  )
}

prepare_classification_data_with_id <- function(data) {
  data |>
    dplyr::mutate(
      Class = factor(
        Class,
        levels = c(0, 1),
        labels = c("Control", "AS")
      )
    ) |>
    tidyr::drop_na()
}

selected_datasets_with_id <- selected_datasets |>
  purrr::map(prepare_classification_data_with_id)


# ------- Pritaikymas klasifikavimo matricos ---------
# --------- Top -------------
top_nb_patient_cm <- get_patient_level_confusion_matrix(
  model = classification_models$top$Naive_Bayes,
  data_with_ids = selected_datasets_with_id$top
)

top_nb_patient_cm$confusion_matrix$table
top_nb_patient_cm$confusion_matrix

top_pred <- top_nb_patient_cm$patient_predictions |>
  dplyr::arrange(Patient_ID)

# -------------- Bottom -------------
bottom_nb_patient_cm <- get_patient_level_confusion_matrix(
  model = classification_models$bottom$Naive_Bayes,
  data_with_ids = selected_datasets_with_id$bottom
)

bottom_nb_patient_cm$confusion_matrix$table
bottom_nb_patient_cm$confusion_matrix

bottom_pred <- bottom_nb_patient_cm$patient_predictions |>
  dplyr::arrange(Patient_ID)

# ---------- Mid ------------
mid_nb_patient_cm <- get_patient_level_confusion_matrix(
  model = classification_models$mid$Naive_Bayes,
  data_with_ids = selected_datasets_with_id$mid
)

mid_nb_patient_cm$confusion_matrix$table
mid_nb_patient_cm$confusion_matrix

mid_pred <- mid_nb_patient_cm$patient_predictions |>
  dplyr::arrange(Patient_ID)


# ---------- Global --------------
global_knn_patient_cm <- get_patient_level_confusion_matrix(
  model = classification_models$global$KNN,
  data_with_ids = selected_datasets_with_id$global
)

global_knn_patient_cm$confusion_matrix$table
global_knn_patient_cm$confusion_matrix

global_pred <- global_knn_patient_cm$patient_predictions |>
  dplyr::arrange(Patient_ID)


# ------------- Bendra lentele pacientams --------------
extract_patient_level_metrics <- function(cm_object, structure_name, structure_label, model_name) {
  
  cm <- cm_object$confusion_matrix
  tibble::tibble(
    Structure = structure_name,
    Structure_LT = structure_label,
    Model = model_name,
    Accuracy = as.numeric(cm$overall["Accuracy"]),
    Balanced_Accuracy = as.numeric(cm$byClass["Balanced Accuracy"]),
    Sensitivity = as.numeric(cm$byClass["Sensitivity"]),
    Specificity = as.numeric(cm$byClass["Specificity"]),
    Pos_Pred_Value = as.numeric(cm$byClass["Pos Pred Value"]),
    Neg_Pred_Value = as.numeric(cm$byClass["Neg Pred Value"])
  )
}

patient_level_metrics <- dplyr::bind_rows(
  extract_patient_level_metrics(
    top_nb_patient_cm,
    "top",
    "Apikalinė zona",
    "Naive_Bayes"
  ),
  extract_patient_level_metrics(
    bottom_nb_patient_cm,
    "bottom",
    "Bazinė zona",
    "Naive_Bayes"
  ),
  extract_patient_level_metrics(
    mid_nb_patient_cm,
    "mid",
    "Vidurinė zona",
    "Naive_Bayes"
  ),
  extract_patient_level_metrics(
    global_knn_patient_cm,
    "global",
    "Visa struktūra",
    "KNN"
  )
)

patient_level_metrics


# ---------- Apmokytas modelis, duomenu vertinimas ----------
# Sis etapas skirtas ne validavimui, o galutiniu modeliu apmokymui

train_final_model_on_all_data <- function(data, old_model, model_name) {
  best_tune <- old_model$bestTune
  
  ctrl_none <- caret::trainControl(
    method = "none",
    classProbs = TRUE
  )
  if (model_name == "SVM_Linear") {
    
    final_model <- caret::train(
      Class ~ .,
      data = data,
      method = "svmLinear",
      trControl = ctrl_none,
      preProcess = c("zv", "nzv", "range"),
      tuneGrid = best_tune
    )
  } else if (model_name == "SVM_RBF") {
    final_model <- caret::train(
      Class ~ .,
      data = data,
      method = "svmRadial",
      trControl = ctrl_none,
      preProcess = c("zv", "nzv", "range"),
      tuneGrid = best_tune
    )
    
  } else if (model_name == "Random_Forest") {
    final_model <- caret::train(
      Class ~ .,
      data = data,
      method = "rf",
      trControl = ctrl_none,
      tuneGrid = best_tune,
      ntree = 500
    )
    
  } else if (model_name == "Naive_Bayes") {
    final_model <- caret::train(
      Class ~ .,
      data = data,
      method = "nb",
      trControl = ctrl_none,
      preProcess = c("zv", "nzv", "range"),
      tuneGrid = best_tune
    )
    
  } else if (model_name == "XGBoost") {
    final_model <- caret::train(
      Class ~ .,
      data = data,
      method = "xgbTree",
      trControl = ctrl_none,
      tuneGrid = best_tune,
      verbose = FALSE,
      verbosity = 0
    )
    
  } else if (model_name == "KNN") {
    final_model <- caret::train(
      Class ~ .,
      data = data,
      method = "knn",
      trControl = ctrl_none,
      preProcess = c("zv", "nzv", "range"),
      tuneGrid = best_tune
    )
  } else {
    stop("Neatpažintas modelis: ", model_name)
  }
  final_model
}

best_model_per_structure <- classification_metrics_clean

final_models <- purrr::pmap(
  list(
    structure_name = best_model_per_structure$Structure,
    model_name = best_model_per_structure$Model
  ),
  function(structure_name, model_name) {
    data <- selected_datasets_prepared[[structure_name]]
    old_model <- classification_models[[structure_name]][[model_name]]
    
    train_final_model_on_all_data(
      data = data,
      old_model = old_model,
      model_name = model_name
    )
  }
)

names(final_models) <- best_model_per_structure$Structure

final_models$global
final_models$bottom
final_models$mid
final_models$top

#new_patient_global
#predict(final_models$global, newdata = new_patient_global)
#predict(final_models$global, newdata = new_patient_global, type = "prob")
# predict(final_models$top, newdata = new_patient_top)
# predict(final_models$mid, newdata = new_patient_mid)
# predict(final_models$bottom, newdata = new_patient_bottom)
#classification_metrics
#patient_level_metrics

# --------- Pozymiai, kurie atrinkti daugiau nei 2 strukturose -----------
features_repeated_3_structures <- feature_stability_table |>
  dplyr::filter(N_structures >= 2) |>
  dplyr::pull(Feature)

features_repeated_3_structures
length(features_repeated_3_structures)

combined_structures_data <- purrr::imap_dfr(
  datasets,
  function(data, structure_name) {
    
    data |>
      dplyr::select(
        Patient_ID,
        Class,
        dplyr::all_of(features_repeated_3_structures)
      ) |>
      dplyr::mutate(
        Structure = structure_name
      )
  }
)

combined_structures_data |>
  dplyr::count(Structure, Class)

 
# Paruosimas klasifikavimui visos strukturos ----------------
combined_structures_prepared <- combined_structures_data |>
  dplyr::mutate(
    Class = factor(
      Class,
      levels = c(0, 1),
      labels = c("Control", "AS")
    ),
    Structure = factor(
      Structure,
      levels = c("global", "bottom", "mid", "top")
    )
  ) |>
  tidyr::drop_na()

combined_structures_prepared |>
  dplyr::count(Structure, Class)

combined_structures_prepared |>
  dplyr::distinct(Patient_ID, Class) |>
  dplyr::count(Class)

# Grouped 5-fold CV pagal Patient_ID -----------
create_group_folds <- function(data, group_col, class_col, k = 5, seed = 123) {
  
  set.seed(seed)
  patient_level <- data |>
    dplyr::distinct(
      Patient_ID = .data[[group_col]],
      Class = .data[[class_col]]
    )
  
  folds_patients <- caret::createFolds(
    patient_level$Class,
    k = k,
    returnTrain = FALSE
  )
  index <- list()
  indexOut <- list()
  
  for (f in seq_along(folds_patients)) {
    validation_patients <- patient_level$Patient_ID[folds_patients[[f]]]
    validation_rows <- which(data[[group_col]] %in% validation_patients)
    training_rows <- setdiff(seq_len(nrow(data)), validation_rows)
    fold_name <- paste0("Fold", f)
    index[[fold_name]] <- training_rows
    indexOut[[fold_name]] <- validation_rows
  }
  list(
    index = index,
    indexOut = indexOut
  )
}

group_folds_5 <- create_group_folds(
  data = combined_structures_prepared,
  group_col = "Patient_ID",
  class_col = "Class",
  k = 5,
  seed = 123
)

# Patikrinimas, ar nera Patient_ID tarp train ir validation -----------
check_group_folds <- function(data, folds, group_col = "Patient_ID") {
  
  purrr::imap_dfr(
    folds$index,
    function(train_rows, fold_name) {
      val_rows <- folds$indexOut[[fold_name]]
      train_ids <- unique(data[[group_col]][train_rows])
      val_ids <- unique(data[[group_col]][val_rows])
      tibble::tibble(
        Fold = fold_name,
        N_train_patients = length(train_ids),
        N_validation_patients = length(val_ids),
        N_overlap_patients = length(intersect(train_ids, val_ids))
      )
    }
  )
}

check_group_folds(
  data = combined_structures_prepared,
  folds = group_folds_5,
  group_col = "Patient_ID"
)

# Modeliui Patient_ID nenaudojame
combined_model_data <- combined_structures_prepared |>
  dplyr::select(-Patient_ID)

# Train control jungtiniam modeliui ---------------
cv_control_combined <- caret::trainControl(
  method = "cv",
  number = 5,
  index = group_folds_5$index,
  indexOut = group_folds_5$indexOut,
  classProbs = TRUE,
  summaryFunction = classification_summary,
  savePredictions = "final",
  sampling = "up",
  verboseIter = FALSE
)


# Modeliu treniravimas jungtiniam visu strukturu rinkiniui
train_models_combined_structures <- function(data) {
  models <- list()
  
  # ---------------- SVM Linear ----------------
  set.seed(123)
  models$SVM_Linear <- caret::train(
    Class ~ .,
    data = data,
    method = "svmLinear",
    trControl = cv_control_combined,
    metric = "Balanced_Accuracy",
    preProcess = c("zv", "nzv", "range"),
    tuneLength = 5
  )
  
  # ---------------- SVM RBF ----------------
  set.seed(123)
  models$SVM_RBF <- caret::train(
    Class ~ .,
    data = data,
    method = "svmRadial",
    trControl = cv_control_combined,
    metric = "Balanced_Accuracy",
    preProcess = c("zv", "nzv", "range"),
    tuneLength = 5
  )
  
  # ---------------- Random Forest ----------------
  set.seed(123)
  models$Random_Forest <- caret::train(
    Class ~ .,
    data = data,
    method = "rf",
    trControl = cv_control_combined,
    metric = "Balanced_Accuracy",
    tuneLength = 5,
    ntree = 500
  )
  
  # ---------------- Naive Bayes ----------------
  nb_grid <- expand.grid(
    fL = c(0, 1),
    usekernel = c(TRUE),
    adjust = c(1, 1.5, 2)
  )
  
  set.seed(123)
  models$Naive_Bayes <- caret::train(
    Class ~ .,
    data = data,
    method = "nb",
    trControl = cv_control_combined,
    metric = "Balanced_Accuracy",
    preProcess = c("zv", "nzv", "range"),
    tuneGrid = nb_grid
  )
  
  # ---------------- XGBoost ----------------
  xgb_grid <- expand.grid(
    nrounds = c(30, 50),
    max_depth = c(1, 2),
    eta = c(0.05, 0.10),
    gamma = c(0, 0.5, 1),
    colsample_bytree = c(0.7, 0.8, 0.9),
    min_child_weight = c(3, 5),
    subsample = c(0.8)
  )
  
  set.seed(123)
  models$XGBoost <- caret::train(
    Class ~ .,
    data = data,
    method = "xgbTree",
    trControl = cv_control_combined,
    metric = "Balanced_Accuracy",
    tuneGrid = xgb_grid,
    verbose = FALSE,
    verbosity = 0
  )
  
  # ---------------- k-NN ----------------
  set.seed(123)
  models$KNN <- caret::train(
    Class ~ .,
    data = data,
    method = "knn",
    trControl = cv_control_combined,
    metric = "Balanced_Accuracy",
    preProcess = c("zv", "nzv", "range"),
    tuneLength = 10
  )
  models
}

combined_structure_models <- train_models_combined_structures(
  combined_model_data
)

# Metriku istraukimas jungtiniams modeliams -----------------
extract_combined_metrics <- function(model, model_name) {
  
  best_tune <- model$bestTune
  results <- model$results
  for (param_name in names(best_tune)) {
    results <- results |>
      dplyr::filter(.data[[param_name]] == best_tune[[param_name]])
  }
  results |>
    dplyr::slice(1) |>
    dplyr::mutate(
      Dataset = "Visos struktūros kartu",
      Model = model_name
    ) |>
    dplyr::select(
      Dataset,
      Model,
      Macro_F1,
      Sensitivity,
      Specificity,
      AUC,
      Balanced_Accuracy,
      dplyr::everything()
    )
}

combined_structure_metrics <- purrr::imap_dfr(
  combined_structure_models,
  extract_combined_metrics
) |>
  dplyr::arrange(
    dplyr::desc(Balanced_Accuracy),
    dplyr::desc(Macro_F1),
    dplyr::desc(AUC)
  ) |>
  dplyr::select(
    Dataset,
    Model,
    Macro_F1,
    Sensitivity,
    Specificity,
    AUC,
    Balanced_Accuracy
  )
combined_structure_metrics

best_combined_structure_model <- combined_structure_metrics |>
  dplyr::slice(1)
best_combined_structure_model

# Paciento lygmens confusion matrix jungtiniam modeliui ------------
get_combined_patient_level_confusion_matrix <- function(
    model,
    data_with_ids,
    positive_class = "AS",
    threshold = 0.5
) {
  preds <- model$pred
  best_tune <- model$bestTune
  
  for (param_name in names(best_tune)) {
    if (param_name %in% names(preds)) {
      preds <- preds |>
        dplyr::filter(.data[[param_name]] == best_tune[[param_name]])
    }
  }
  id_lookup <- data_with_ids |>
    dplyr::mutate(rowIndex = dplyr::row_number()) |>
    dplyr::select(rowIndex, Patient_ID, Structure)
  
  patient_predictions <- preds |>
    dplyr::left_join(id_lookup, by = "rowIndex") |>
    dplyr::group_by(Patient_ID, obs) |>
    dplyr::summarise(
      Mean_AS_Probability = mean(.data[[positive_class]], na.rm = TRUE),
      N_predictions = dplyr::n(),
      N_structures = dplyr::n_distinct(Structure),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      pred = ifelse(
        Mean_AS_Probability >= threshold,
        positive_class,
        "Control"
      ),
      pred = factor(pred, levels = c("Control", "AS")),
      obs = factor(obs, levels = c("Control", "AS"))
    )
  cm <- caret::confusionMatrix(
    data = patient_predictions$pred,
    reference = patient_predictions$obs,
    positive = positive_class
  )
  list(
    patient_predictions = patient_predictions,
    confusion_matrix = cm
  )
}

best_combined_model_name <- best_combined_structure_model$Model[1]

combined_patient_cm <- get_combined_patient_level_confusion_matrix(
  model = combined_structure_models[[best_combined_model_name]],
  data_with_ids = combined_structures_prepared,
  positive_class = "AS",
  threshold = 0.5
)

combined_patient_cm$confusion_matrix$table
combined_patient_cm$confusion_matrix


combined_patient_predictions <- combined_patient_cm$patient_predictions |> 
  dplyr::arrange(Patient_ID)
combined_patient_predictions

extract_patient_level_metrics <- function(cm_object, dataset_name, model_name) {
  cm <- cm_object$confusion_matrix
  
  tibble::tibble(
    Dataset = dataset_name,
    Model = model_name,
    Accuracy = as.numeric(cm$overall["Accuracy"]),
    Balanced_Accuracy = as.numeric(cm$byClass["Balanced Accuracy"]),
    Sensitivity = as.numeric(cm$byClass["Sensitivity"]),
    Specificity = as.numeric(cm$byClass["Specificity"]),
    Pos_Pred_Value = as.numeric(cm$byClass["Pos Pred Value"]),
    Neg_Pred_Value = as.numeric(cm$byClass["Neg Pred Value"])
  )
}

combined_patient_level_metrics <- extract_patient_level_metrics(
  cm_object = combined_patient_cm,
  dataset_name = "Visos struktūros kartu",
  model_name = best_combined_model_name
)
combined_patient_level_metrics

combined_for_comparison <- combined_structure_metrics |>
  dplyr::mutate(
    Structure_LT = "Visos struktūros kartu"
  ) |>
  dplyr::select(
    Structure_LT,
    Model,
    Macro_F1,
    Sensitivity,
    Specificity,
    AUC,
    Balanced_Accuracy
  )

comparison_structures_and_combined <- dplyr::bind_rows(
  classification_metrics |>
    dplyr::select(
      Structure_LT,
      Model,
      Macro_F1,
      Sensitivity,
      Specificity,
      AUC,
      Balanced_Accuracy
    ),
  combined_for_comparison
)

comparison_structures_and_combined |>
  dplyr::arrange(
    Structure_LT,
    dplyr::desc(Balanced_Accuracy),
    dplyr::desc(Macro_F1),
    dplyr::desc(AUC)
  )



