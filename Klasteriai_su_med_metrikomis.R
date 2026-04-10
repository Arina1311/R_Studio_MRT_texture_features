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
library(purrr)
library(ggtext)
library(Hmisc)
library(aricode)
library(clue)

#--------------DUOMENŲ NUSKAITYMAS, POŽYMIAI IR MEDIKŲ METRIKOS ---------
# Nuskaitomi teksturos pozymiai ----------------
json_data <- fromJSON("/Users/arinaperzu/Desktop/MRT darbas/MRT2/src/features/patient_features_phases_1_to_25.json")
df <- as.data.frame(t(sapply(json_data, unlist)))
head(df)

# Imame tik pagal agreguotus MRT vaizdus
df_su_teksturos_pozymiais <- df[, grepl("global$", colnames(df))]|>
  mutate(Patient_ID = rownames(df)) |>
  dplyr::select(Patient_ID, everything()) |>
  mutate(Class = ifelse(grepl("^AS", Patient_ID), 1, 0))

lentele_su_pacientu_ID <- df_su_teksturos_pozymiais

# Lentele be pacientu ID
df_su_teksturos_pozymiais <- df_su_teksturos_pozymiais[,-1]
table(df_su_teksturos_pozymiais$Class)
# 13 sveiku ir 74 su pazeistu raumeniu

# Klinikiniu metriku nuskaitymas is Exel lenteles ---------------------
# Imame tik ID ir kelias metrikas pagal pavadinimus
pacientu_lentele_su_med_metrikomis <- read_excel("/Users/arinaperzu/Desktop/MRT darbas/MRT Rkodas/pataisyta lentelė.xlsx")

pacientai <- pacientu_lentele_su_med_metrikomis |>
  dplyr::select(ID, h_CVF, 'm_Nat T1', m_ECV, e_GLS)

# Duomenu apjungimas -----------------
# Paimam numerius, kad galima butu sujungti, turime tik sergantiems pacientams med metrikas, todel ziurime tik AS
duomenys_klasteriams <- lentele_su_pacientu_ID |>
  mutate( ID = if_else(str_starts(Patient_ID, "AS"),
          parse_number(Patient_ID), NA_real_)) |>
  left_join(pacientai, by = "ID")

#------------ PRISKIRIAME TRUMPINIUS POŽYMIŲ PAVADINIMAMS --------------------------------
df_su_teksturos_pozymiais <- df_su_teksturos_pozymiais %>%
  rename(
    F_En  = Fourier_Energy_global,
    F_E   = Fourier_Entropy_global,
    F_M   = Fourier_Mean_global,
    F_V   = Fourier_Variance_global,
    FR_M  = Fractal_mean_global,
    FR_V  = Fractal_variance_global,
    GLCM_ASM   = GLCM_ASM_global,
    GLCM_Con   = GLCM_contrast_global,
    GLCM_Cor  = GLCM_correlation_global,
    GLCM_D   = GLCM_dissimilarity_global,
    GLCM_En   = GLCM_energy_global,
    GLCM_Hom   = GLCM_homogeneity_global,
    GLRLM_GLNU   = GLRLM_GrayLevelNonUniformity_global,
    GLRLM_GLNUN  = GLRLM_GrayLevelNonUniformityNormalized_global,
    GLRLM_GLV   = GLRLM_GrayLevelVariance_global,
    GLRLM_HGLE  = GLRLM_HighGrayLevelRunEmphasis_global,
    GLRLM_LRE   = GLRLM_LongRunEmphasis_global,
    GLRLM_LRHGE = GLRLM_LongRunHighGrayLevelEmphasis_global,
    GLRLM_LRLGE = GLRLM_LongRunLowGrayLevelEmphasis_global,
    GLRLM_LGLE  = GLRLM_LowGrayLevelRunEmphasis_global,
    GLRLM_RE    = GLRLM_RunEntropy_global,
    GLRLM_RLNU  = GLRLM_RunLengthNonUniformity_global,
    GLRLM_RLNUN = GLRLM_RunLengthNonUniformityNormalized_global,
    GLRLM_RP    = GLRLM_RunPercentage_global,
    GLRLM_RV    = GLRLM_RunVariance_global,
    GLRLM_SRE   = GLRLM_ShortRunEmphasis_global,
    GLRLM_SRHGLE= GLRLM_ShortRunHighGrayLevelEmphasis_global,
    GLRLM_SRLGLE= GLRLM_ShortRunLowGrayLevelEmphasis_global,
    HOG_M = HOG_mean_global,
    HOG_V = HOG_variance_global,
    LBP_Con = LBP_contrast_global,
    LBP_U = LBP_uniformity_global,
    WL1_HH_M = Wavelet_level_1_HH_mean_global,
    WL1_HH_Sd = Wavelet_level_1_HH_sd_global,
    WL1_HL_M = Wavelet_level_1_HL_mean_global,
    WL1_HL_Sd = Wavelet_level_1_HL_sd_global,
    WL1_LH_M = Wavelet_level_1_LH_mean_global,
    WL1_LH_Sd = Wavelet_level_1_LH_sd_global,
    WL2_HH_M = Wavelet_level_2_HH_mean_global,
    WL2_HH_Sd = Wavelet_level_2_HH_sd_global,
    WL2_HL_M = Wavelet_level_2_HL_mean_global,
    WL2_HL_Sd = Wavelet_level_2_HL_sd_global,
    WL2_LH_M = Wavelet_level_2_LH_mean_global,
    WL2_LH_Sd = Wavelet_level_2_LH_sd_global,
    WL2_LL_M = Wavelet_level_2_LL_mean_global,
    WL2_LL_Sd = Wavelet_level_2_LL_sd_global
  )

# Standartizuoti teksturos pozymiai
teksturos_pozymiai <- df_su_teksturos_pozymiais %>%
  dplyr::select(-Class)

standartizuoti_teksturos_pozymiai <- as.data.frame(scale(teksturos_pozymiai))

standartizuoti_teksturos_pozymiai <- dplyr::bind_cols(
  Class = df_su_teksturos_pozymiais$Class,
  standartizuoti_teksturos_pozymiai
)
# ---------- INFORMATYVIAUSIU POZYMIU GAVIMAS ----------------------------------
# Vidurkiu analize ----------------
analyze_features <- function(df, class_col = "Class", alpha = 0.01) {
  feature_cols <- setdiff(names(df), class_col)
  
  # Pagalbinė funkcija normalumo testui, patikriname ar reiksmiu uztenka, pasaliname NA
  get_shapiro_p <- function(x) {
    x <- x[!is.na(x)]
    if (length(x) < 10 || length(unique(x)) < 2) {
      return(NA_real_)
    }
    shapiro.test(x)$p.value
  }
  
  # Pagalbine funkcija vieno pozymio analizei, iskirstome pagal grupes ir tikriname normaluma
  analyze_single_feature <- function(feature_name) {
    group0 <- df[df[[class_col]] == 0, feature_name, drop = TRUE]
    group1 <- df[df[[class_col]] == 1, feature_name, drop = TRUE]
    
    norm0 <- get_shapiro_p(group0)
    norm1 <- get_shapiro_p(group1)
    
    is_normal <- !is.na(norm0) &&
      !is.na(norm1) &&
      norm0 > 0.05 &&
      norm1 > 0.05
    
    mean0 <- mean(group0, na.rm = TRUE)
    mean1 <- mean(group1, na.rm = TRUE)
    
    group0_clean <- group0[!is.na(group0)]
    group1_clean <- group1[!is.na(group1)]
    combined <- c(group0_clean, group1_clean)
    
    # Tikriname grupiu atsiskyrima, jei normaliai pasiskirste ir jei ne
    if (length(combined) == 0 || length(unique(combined)) == 1) {
      p_val <- 1
    } else if (is_normal) {
      p_val <- t.test(group0_clean, group1_clean)$p.value
    } else {
      p_val <- wilcox.test(group0_clean, group1_clean)$p.value
    }
    
    data.frame(
      Feature = feature_name,
      Normality = ifelse(is_normal, "Yes", "No"),
      Mean_Group0 = mean0,
      Mean_Group1 = mean1,
      p_value = p_val,
      Difference = round(abs(mean0 - mean1), 3),
      stringsAsFactors = FALSE
    )
  }
  
  # Visu pozymiu analize
  results <- bind_rows(lapply(feature_cols, analyze_single_feature)) |>
    arrange(p_value)
  # Tik reiksmingi pozymiai
  sig_features <- results|>
    filter(p_value < alpha)
  # Duomenys tolimesniam modeliavimui
  filtered_data <- df |>
    dplyr::select(all_of(c(class_col, sig_features$Feature)))
  
  return(list(
    results_table = results,
    significant_features = sig_features,
    filtered_data = filtered_data
  ))
}

# Teksturos pozymiu rezultatai pagal t.test ir Mann Whitnio testus
rezultatai <- analyze_features(df_su_teksturos_pozymiais)
rezultatai_su_visa_statistika <- rezultatai$results_table
rezultatai_tik_reiksmingi<- rezultatai$significant_features
duomenys_su_informatyviausiais_pozymiais <- rezultatai$filtered_data

# Imame didziausius skirtumus tarp grupiu 
top_10_skirtumas <- rezultatai_tik_reiksmingi |>
  arrange(desc(Difference)) |>
  head(10)

# Surenkame visus reiksmingus pozymius i bendra lentele 
pozymiai_reiksmingi_bendras <- top_10_skirtumas |>
  dplyr::select(Mean_difference = Feature)


# Atsitiktiniu misku metodas ----------------
set.seed(123)
df_su_teksturos_pozymiais$Class <- as.factor(df_su_teksturos_pozymiais$Class)
rf_model <- randomForest(Class ~ ., data = df_su_teksturos_pozymiais, importance = TRUE, ntree = 500)
importance_scores <- as.data.frame(importance(rf_model))
importance_scores$Feature <- rownames(importance_scores)

colnames(importance_scores)[colnames(importance_scores) == "%IncMSE"] <- "MeanDecreaseAccuracy"
colnames(importance_scores)[colnames(importance_scores) == "IncNodePurity"] <- "MeanDecreaseGini"

# Top 10 pagal tiksluma ir Gini koeficienta
top10_MeanDecreaseAccuracy <- importance_scores |>
  arrange(desc(MeanDecreaseAccuracy)) |>
  slice(1:10)

top10_MeanDecreaseGini <- importance_scores |>
  arrange(desc(MeanDecreaseGini)) |>
  slice(1:10)

pozymiai_reiksmingi_bendras$Random_Forest <- top10_MeanDecreaseGini$Feature

# Regresija su L1 reguliarizacija ----------
set.seed(123)
x <- model.matrix(Class ~ . - 1, data = standartizuoti_teksturos_pozymiai)
y <- standartizuoti_teksturos_pozymiai$Class

lasso_model <- cv.glmnet(x = x, y = y, family = "binomial", alpha = 1, standardize = FALSE)

coef_df <- coef(lasso_model, s = "lambda.min") |>
  as.matrix() |>
  as.data.frame() |>
  tibble::rownames_to_column(var = "feature")

colnames(coef_df)[2] <- "coefficient"

# Atrenkame tik nenulinius koeficientus
coef_df_filtered <- coef_df |>
  dplyr::filter(feature != "(Intercept)", coefficient != 0)

# 10 pirmųjų požymių
pozymiai_lasso <- coef_df_filtered |>
  mutate(abs_coefficient = abs(coefficient)) |>
  arrange(desc(abs_coefficient)) |>
  slice_head(n = 10)

pozymiai_reiksmingi_bendras$Regression <- NA
pozymiai_reiksmingi_bendras$Regression[1:nrow(pozymiai_lasso)] <- pozymiai_lasso$feature

# Pozymiu reiksmingumas pagal INFORMATION GAIN -------
info_gain <- information.gain(Class ~ ., df_su_teksturos_pozymiais) |>
  as.data.frame() |>
  tibble::rownames_to_column(var = "feature") |>
  dplyr::select(feature, attr_importance) |>
  dplyr::arrange(desc(attr_importance))

pozymiai_univar <- info_gain |>
  dplyr::slice_head(n = 10)

pozymiai_reiksmingi_bendras$Info_Gain <- pozymiai_univar$feature

#---------------- POZYMIU PO PRITAIKYTU METODU ATRINKIMAS ------------
reiksmingiausi_pozymiai_bendras <- pozymiai_reiksmingi_bendras |>
  pivot_longer(
    cols = everything(),
    names_to = "metodas",
    values_to = "pozymis"
  ) |>
  filter(!is.na(pozymis)) |>
  count(pozymis, sort = TRUE) |>
  filter(n >= 2)

reiksmingiausi_pozymiai_bendras

teksturos_pozymiai_reiksmingiausi_pilni_duomenys <- standartizuoti_teksturos_pozymiai |>
  dplyr::select(all_of(reiksmingiausi_pozymiai_bendras$pozymis))

# Apskaiciuojame koreliacijas tarp pozymiu ir pasaliname labai stipriai koreliuojancius
koreliacijos <- cor(teksturos_pozymiai_reiksmingiausi_pilni_duomenys, method = "spearman")
stipri_koreliacija <- findCorrelation(koreliacijos, cutoff = 0.8)

teksturos_pozymiai_reiksmingiausi_pilni_duomenys <- teksturos_pozymiai_reiksmingiausi_pilni_duomenys[, -stipri_koreliacija]
colnames(teksturos_pozymiai_reiksmingiausi_pilni_duomenys)

# Apskaiciuojam nauja koreliaciju matrica po filtravimo
kor_matrica <- cor(teksturos_pozymiai_reiksmingiausi_pilni_duomenys, method = "spearman")
ggcorrplot(kor_matrica,
           method = "square",
           type = "lower",     
           lab = TRUE,
           lab_size = 3.5,
           title = "Koreliacijų matrica tarp išrinktų požymių",
           show.legend = TRUE,
           legend.title = "Spearmano\nkoreliacija",
           colors = c("#6D9EC1", "white","#E64B35"),
           tl.col = "black",
           tl.cex = 9   )+
  theme(
    plot.title = element_text(hjust = 0, face = "bold"),
    legend.title = element_text( size = 10)
  )


Atrinkti_pozymiai <- standartizuoti_teksturos_pozymiai |>
  dplyr::select(Class, all_of(colnames(teksturos_pozymiai_reiksmingiausi_pilni_duomenys)))

Atrinkti_pozymiai <- standartizuoti_teksturos_pozymiai |>
  dplyr::select(Class, F_En, GLCM_Cor, GLRLM_LRHGE, GLRLM_RE, GLRLM_RLNUN, LBP_Con)
# --------------- ATRINKTU POZYMIU ANALIZE IR STATISTIKOS ------------
# Bendras vaizdas pozymiu ir koreliaciju
ggpairs(
  data = Atrinkti_pozymiai |> dplyr::select(-Class),
  upper = list(continuous = wrap("cor", method = "spearman", size = 4)),
  lower = list(continuous = wrap("points", alpha = 0.7, size = 1)),
  diag = list(continuous = wrap("densityDiag")),
  title = "Korelograma pagal atrinktus požymius"
) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16)
  )

# ----------- POZYMIU ATSISKYRIMAS PAGAL UMAP PROJEKCIJA ---------------
set.seed(42)

X_umap_sel <- dplyr::select(Atrinkti_pozymiai, -Class)
umap_sel_result <- umap(X_umap_sel, n_neighbors = 15, min_dist = 0.05, metric = "euclidean")

umap_sel_df <- as.data.frame(umap_sel_result)
colnames(umap_sel_df) <- c("UMAP1", "UMAP2")
umap_sel_df$Class <- Atrinkti_pozymiai$Class
umap_sel_df$Patient_ID <- lentele_su_pacientu_ID$Patient_ID

# Interaktyvus grafikas
p <- ggplot(umap_sel_df, aes(
  x = UMAP1,
  y = UMAP2,
  color = as.factor(Class),
  text = paste0("Pacientas: ", Patient_ID)
)) +
  theme_minimal(base_size = 12) +
  geom_point(alpha = 0.85, size = 2.8) +
  scale_color_manual(values = c("0" = "#4DBBD5", "1" = "#E64B35")) +
  labs(
    title = "UMAP dimensijos mažinimas - atrinkti požymiai<br><sup>Parametrai: n = 15, min_dist = 0.05, metric = Euclidean</sup>",
    color = "Klasė",
    x = "UMAP1",
    y = "UMAP2"
  )

ggplotly(p, tooltip = "text") %>%
  layout(
    title = list(
      text = "UMAP dimensijos mažinimas - atrinkti požymiai<br><sup>Parametrai: n = 15, min_dist = 0.05, metric = Euclidean</sup>",
      font = list(size = 16)  
    ),
    xaxis = list(title = list(text = "UMAP1", font = list(size = 12))),
    yaxis = list(title = list(text = "UMAP2", font = list(size = 12)))
  )

# Bandymas klasterizuoti duomenis naudojant K-means --------------
X <- dplyr::select(Atrinkti_pozymiai, -Class)

# Silueto metoda - parodo 4 klasterius
fviz_nbclust(X, kmeans, method = "silhouette") +
  labs(title = "Silueto metodas klasterių skaičiui parinkti",
       x = "Klasterių skaičius",
       y = "Vidutinė silueto reikšmė") +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 11, color = "black")
  )

# Alkunes metodas - parodo 3 ar 7 klasterius
fviz_nbclust(X, kmeans, method = "wss") +
  labs(title = "Elbow metodas klasterių skaičiui parinkti",
       x = "Klasterių skaičius",
       y = "Vidinė klasterių dispersija") +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 11, color = "black")
  )

set.seed(42)
k <- 5
kmeans_model <- kmeans(X, centers = k)
d <- dist(X)
sil <- silhouette(kmeans_model$cluster, d)
summary(sil)

umap_result <- umap(X, n_neighbors = 30, min_dist = 0.1, metric = "manhattan", n_epochs = 500)
umap_result <- as.data.frame(umap_result)
colnames(umap_result) <- c("UMAP1", "UMAP2")

# Projekcijos kokybes vertinimas
projection_quality <- calculatesSaturnContinuityTrustworthiness(
  original_matrix = as.matrix(X),
  umap_output_layout = as.matrix(umap_result[,1:2]),
  VERBOSE = FALSE
)
print(projection_quality)

umap_result$Cluster <- as.factor(kmeans_model$cluster)
umap_result$Patient_ID <- lentele_su_pacientu_ID$Patient_ID
umap_df <- umap_result

# Klateriu vizualizavimas
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Cluster)) +
  geom_point(size = 2.5, alpha = 0.8) +
  labs(
    title = "Clustering using k-means",
    subtitle = "Clustering: k-means (k = 5)\nVisualization: UMAP (n_neighbors = 10, min_dist = 0.05, Euclidean distance)",
    color = "Cluster",
    x = "UMAP1",
    y = "UMAP2"
  ) +
  scale_color_manual(
    values = c(
      "royalblue",  # Cluster 1 – green
      "violetred4",  # Cluster 2 – orange
      "skyblue1",  # Cluster 3 – purple
      "royalblue4",  # Cluster 4 – pink
      "orchid4"   # Cluster 5 – olive/green
    )
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(hjust = 0, face = "bold", size = 12),
    plot.subtitle = element_text(hjust = 0, size = 9, color = "gray20"),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 10, color = "black")
  )

# Vizualizacija
library(plotly)
p <- ggplot(umap_df, aes(
  x = UMAP1,
  y = UMAP2,
  color = as.factor(Cluster),
  text = paste0("Pacientas: ", Patient_ID,
                "<br>Klasteris: ", Cluster)
)) +
  geom_point(size = 2.5, alpha = 0.8) +
  labs(
    title = "Klasterizavimas naudojant k-means metodą<br><sup>Parametrai k-means: k = 5; Parametrai UMAP: n = 5, min_dist = 0.01, metric = Euclidean</sup>",
    color = "Klasteris",
    x = "UMAP1",
    y = "UMAP2"
  ) +
  theme_minimal(base_size = 12)

ggplotly(p, tooltip = "text") %>%
  layout(
    title = list(
      text = "Klasterizavimas naudojant k-means metodą<br><sup>Parametrai k-means: k = 5; Parametrai UMAP: n = 5, min_dist = 0.01, metric = Euclidean</sup>",
      font = list(size = 16)
    ),
    xaxis = list(title = list(text = "UMAP1", font = list(size = 12))),
    yaxis = list(title = list(text = "UMAP2", font = list(size = 12))),
    legend = list(title = list(text = "Klasteris", font = list(size = 12)))
  )

Atrinkti_pozymiai$Cluster <- umap_result$Cluster 
Atrinkti_pozymiai$group1 <- ifelse(Atrinkti_pozymiai$Class == 0, 1, 2)

# Paziurime su PCA vizualiai -------------
set.seed(42)
pca_res <- prcomp(X, center = TRUE, scale. = TRUE)

vizualizavimui_df <- data.frame(
  PC1 = pca_res$x[, 1],
  PC2 = pca_res$x[, 2],
  Cluster = as.factor(kmeans_model$cluster), 
  Class = as.factor(Atrinkti_pozymiai$Class),
  Patient_ID = lentele_su_pacientu_ID$Patient_ID
)

# Vizualizacija
ggplot(vizualizavimui_df, aes(x = PC1, y = PC2, color = Cluster, shape = Class)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(aes(group = Cluster), linetype = 2, alpha = 0.5) + 
  scale_shape_manual(values = c(16, 17), 
                     labels = c("Sveikas (0)", "Pažeistas (1)")) +
  scale_color_brewer(palette = "Set1") + 
  theme_minimal() +
  labs(
    title = "MRT tekstūros klasterizavimas (K-means)",
    subtitle = "Projekcija į PCA erdvę",
    x = paste0("PC1 (", round(summary(pca_res)$importance[2,1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(summary(pca_res)$importance[2,2] * 100, 1), "%)"),
    color = "Klasteris",
    shape = "Tikroji klasė"
  )

# PACIENTU DUOMENYS PAGAL KLASTERIUS ----------------
Atrinkti_pozymiai_su_paciento_id <- Atrinkti_pozymiai
Atrinkti_pozymiai_su_paciento_id$pacientas <- lentele_su_pacientu_ID$Patient_ID

# Iskirstome pacientus pagal klasterius
Pacientai_pagal_klasterius <- Atrinkti_pozymiai_su_paciento_id %>%
  group_by(Cluster) %>%
  summarise(
    pacientai = paste(sort(pacientas), collapse = ", "),
    n = dplyr::n(),
    .groups = "drop"
  ) %>%
  arrange(Cluster)

# Paziurime kiek pacientu sukrenta i klasteri, koks procentas sveiku ir serganciu
Klasterio_analize <- Atrinkti_pozymiai %>%
  group_by(Cluster, Class) %>%
  summarise(count = n(), .groups = "drop") %>%
  arrange(Cluster, Class) %>%
  group_by(Cluster) %>%
  mutate(
    suma = sum(count),
    procentas = round(100 * count / suma, 1)
  ) %>%
  ungroup() %>%
  dplyr::select(-suma) %>%
  arrange(Cluster, Class)

Klasterio_analize$grupe <- ifelse(Klasterio_analize$Class == 0, "Sveikas", "Pažeistas")

# Lenteles vizualizavimas
ggplot(Klasterio_analize, aes(x = as.factor(Cluster), y = count, fill = grupe)) +
  geom_bar(
    stat = "identity",
    position = position_dodge2(width = 0.6, preserve = "single"),
    width = 0.7,
    alpha = 0.8
  ) +
  geom_text(
    aes(label = paste0(procentas, "%")),
    position = position_dodge2(width = 0.7),
    vjust = -0.2,
    size = 4,
    fontface = "bold",
    color = "black"
  ) +
  scale_fill_manual(
    values = c(
      "Pažeistas" = "darkred",
      "Sveikas"   = "skyblue1"
    ),
    name = "Raumuo"
  ) +
  labs(
    title = "Pacientų pasiskirstymas tarp klasteriu",
    x = "Klasteris",
    y = "Pacientų skaičius"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0, face = "bold", size = 15),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10, color = "black")
  )

# Statistikos kiekvienam klasteriui -----------
Klasteriu_statistika <- Atrinkti_pozymiai |> 
  dplyr::select(-c(Class, group1)) |>
  group_by(Cluster) |>
  summarise(across(where(is.numeric), list(
    mean = ~mean(.),
    median = ~median(.),
    sd = ~sd(.),
    min = ~min(.),
    max = ~max(.)
  ), .names = "{.col}_{.fn}")) |>
  pivot_longer(
    cols = -Cluster, 
    names_to = "Metric_Stat", 
    values_to = "Value"
  ) |>
  separate(Metric_Stat, into = c("Metrika", "Stat"), sep = "_(?=[^_]+$)") |>
  pivot_wider(
    names_from = Stat, 
    values_from = Value
  ) |>
  arrange(Metrika, Cluster)
print(Klasteriu_statistika)

#------------- POZYMIU PASISKIRSTYMAS TARP KLASTERIU ----------------
Atrinkti_pozymiai_long <- Atrinkti_pozymiai |>
  dplyr::select(-c(group1)) |>
  pivot_longer(cols = -c(Cluster, Class), names_to = "pozymis", values_to = "reiksme")

custom_colors <- c("1" = "royalblue",  
                   "2" = "violetred4",  
                   "3" = "skyblue1",  
                   "4" = "royalblue4",
                   "5" ="orchid4" )  

ggplot(Atrinkti_pozymiai_long, aes(x = factor(Cluster), y = reiksme, fill = factor(Cluster))) +
  geom_boxplot(alpha = 0.8, outlier.color = "black", outlier.size = 1.2) +
  facet_wrap(~ pozymis, scales = "free", ncol = 3) +
  labs(
    x = "Klasteris",
    y = "Požymio reikšmė",
    title = "Tekstūros požymių pasiskirstymas pagal klasterius",
    fill = "Klasteris"
  ) +
  scale_fill_manual(values = custom_colors) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0, face = "bold", size=14),
    strip.text = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    legend.position = "none"
  )


# Kruskal Wallis testas teksturos pozymiu skirtumui tarp klasteriu ----------------
Rezultatai_Kruskal_Wallis <- sapply(names(Atrinkti_pozymiai)[sapply(Atrinkti_pozymiai, is.numeric) & !names(Atrinkti_pozymiai) %in% c("Cluster", "Class", "group1")], function(var) {
  test <- kruskal.test(reformulate("factor(Cluster)", response = var), data = Atrinkti_pozymiai)
  return(test$p.value)
})

# Lentele su rezultatais
df_long <- Atrinkti_pozymiai |>
  pivot_longer(
    cols = where(is.numeric) & !c(Cluster, Class), 
    names_to = "Metrika", 
    values_to = "Reiksme"
  ) |>
  mutate(Cluster = as.factor(Cluster))

kw_res <- data.frame(
  Metrika = names(Rezultatai_Kruskal_Wallis),
  p_value = as.numeric(Rezultatai_Kruskal_Wallis)
) |>
  mutate(p_BH = p.adjust(p_value, method = "BH")) |> 
  filter(p_BH < 0.05) |>
  arrange(p_BH)
sig_metrics <- kw_res$Metrika

# Dunn testas (porinis palyginimas)
posthoc_dunn <- df_long |>
  filter(Metrika %in% sig_metrics) |>
  group_by(Metrika) |>
  dunn_test(Reiksme ~ Cluster, p.adjust.method = "BH")

# Wilcoxon porinis palyginimas + Efekto dydis (r)
eff_pairs_texturos_pozymiai <- df_long |>
  filter(Metrika %in% sig_metrics) |>
  group_by(Metrika) |>
  pairwise_wilcox_test(Reiksme ~ Cluster, p.adjust.method = "BH") |>
  left_join(
    df_long |>
      filter(Metrika %in% sig_metrics) |>
      group_by(Metrika) |>
      wilcox_effsize(Reiksme ~ Cluster),
    by = c("Metrika", "group1", "group2")
  ) |>
  select(Metrika, group1, group2, p.adj, effsize, magnitude) |>
  arrange(Metrika, p.adj) |>
  filter(p.adj < 0.05) 

print(eff_pairs_texturos_pozymiai)


# ---------------- PACIENTU PRIJUNGIMAS -------------------
nauja_pacientai_klasteriai <- duomenys_klasteriams %>% 
  select(Patient_ID, Class, h_CVF, 'm_Nat T1', m_ECV, e_GLS) %>% 
  mutate(Cluster = Atrinkti_pozymiai$Cluster)

metric_cols <- nauja_pacientai_klasteriai %>%
  select(where(is.numeric)) %>%
  names() %>%
  setdiff(c("Class","Cluster"))

# NA reiksmes pasalinime kiekvienai metrikai atskirai
nauja_pacientai_klasteriai_long <- nauja_pacientai_klasteriai %>%
  select(Cluster, all_of(metric_cols)) %>%
  pivot_longer(cols = -Cluster,
               names_to = "Metrika",
               values_to = "Reiksme") %>%
  filter(!is.na(Reiksme)) %>%          
  mutate(Cluster = factor(Cluster))

# Medicininiu metriku pasiskirstymas tarp klasteriu ----------------
summary_pagal_med_metrikas <- nauja_pacientai_klasteriai_long %>%
  group_by(Metrika, Cluster) %>%
  summarise(
    n = dplyr::n(),                
    mean = mean(Reiksme),
    median = median(Reiksme),
    sd = sd(Reiksme),
    q1 = quantile(Reiksme, 0.25),
    q3 = quantile(Reiksme, 0.75),
    min = min(Reiksme),
    max = max(Reiksme),
    .groups = "drop"
  ) %>%
  arrange(Metrika, Cluster)
print(summary_pagal_med_metrikas)

# Vizualizavimas pasiskirstymo -------------
nauja_pacientai_klasteriai_long$Cluster <- factor(
  nauja_pacientai_klasteriai_long$Cluster,
  levels = c("1", "2", "3", "4", "5")
)

# 5 klasteris turi tik 1 serganti pacienta, todel jam priskiriamas taskas
ggplot(nauja_pacientai_klasteriai_long,
       aes(x = Cluster, y = Reiksme, fill = Cluster)) +
  geom_boxplot(
    data = subset(nauja_pacientai_klasteriai_long, Cluster != "5"),
    alpha = 0.75, outlier.alpha = 0.7, width = 0.7
  ) +
  geom_point(
    data = subset(nauja_pacientai_klasteriai_long, Cluster == "5"),
    size = 2.7,
    shape = 21,
    color = "black"
  ) +
  facet_wrap(~ Metrika, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = custom_colors) +
  labs(
    title = "Metrikų pasiskirstymas pagal klasterius",
    x = "Klasteris", y = "Reikšmė"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0, face = "bold", size = 15),
    strip.text = element_text(size = 10, face = "bold"),
    legend.position = "none"
  )

# Analizuojame ar medicinines metrikos skiriasi tarp klasteriu -----------
df_long <- nauja_pacientai_klasteriai_long %>%
  filter(!is.na(Reiksme)) %>%
  mutate(
    Cluster = factor(Cluster),
    Reiksme = as.numeric(Reiksme)   
  )

kw_res <- df_long |>
  group_by(Metrika) |>
  kruskal_test(Reiksme ~ Cluster) |>   
  ungroup() |>
  mutate(p_BH = p.adjust(p, method = "BH")) |>
  arrange(p_BH)
kw_res

# POST-HOC Dunn 
sig_metrics <- kw_res |> filter(p_BH < 0.05) |> pull(Metrika)
posthoc_dunn <- df_long |>
  filter(Metrika %in% sig_metrics) |>
  group_by(Metrika) |>
  dunn_test(Reiksme ~ Cluster, p.adjust.method = "BH") |>
  arrange(Metrika, p.adj) |>
  ungroup()
posthoc_dunn

# Poriniai Wilcoxon + efekto dydis (r)
eff_pairs <- df_long %>%
  filter(Metrika %in% sig_metrics) %>%
  group_by(Metrika) %>%
  pairwise_wilcox_test(Reiksme ~ Cluster, p.adjust.method = "BH") %>%
  mutate(p.signif = rstatix::p_format(p.adj)) %>%
  left_join(
    df_long %>% group_by(Metrika) %>% wilcox_effsize(Reiksme ~ Cluster, paired = FALSE),
    by = c("Metrika","group1","group2")
  ) %>%
  arrange(Metrika, p.adj)
eff_pairs

# ------------------- PACIENTAI SU TEKSTUROS POZYMIAIS ! --------------------------
nauja_pacientai_features <- nauja_pacientai_klasteriai %>%
  mutate(F_En = Atrinkti_pozymiai$F_En,
         GLCM_Cor = Atrinkti_pozymiai$FR_M,
         GLRLM_LRHGE = Atrinkti_pozymiai$GLRLM_LRHGE,
         GLRLM_RE = Atrinkti_pozymiai$GLRLM_RE,
         GLRLM_RLNU = Atrinkti_pozymiai$GLRLM_RLNU,
         LBP_Con = Atrinkti_pozymiai$LBP_Con)

nauja_pacientai_features$Cluster <- as.factor(nauja_pacientai_features$Cluster)

# Paimame tik medicinines metrikas
med_cols <- names(nauja_pacientai_features)[3:6]
feature_cols <- names(nauja_pacientai_features)[8:12]

# Nagrinejame koreliacijas kiekviename klasteryje -------------------
cross_cor_cluster <- function(d, med_cols, feature_cols, cluster_id) {
  X <- as.matrix(d[, med_cols,  drop = FALSE])
  Y <- as.matrix(d[, feature_cols, drop = FALSE])
  R <- suppressWarnings(cor(X, Y, method = "spearman", 
                            use = "pairwise.complete.obs"))
  as_tibble(R, rownames = "med_metric") |>
  pivot_longer(-med_metric, names_to = "feature", values_to = "rho") |>
  mutate(Cluster = cluster_id)
}

corr_long <- nauja_pacientai_features |>
  filter(!is.na(Cluster)) |>
  group_split(Cluster) |>
  map_dfr(~ cross_cor_cluster(.x, med_cols, feature_cols, unique(.x$Cluster)))

# Ziurime giliaus i koreliacijas klasteriu viduje -----------
corr_with_p <- nauja_pacientai_features |>
  filter(!is.na(Cluster)) |>
  group_by(Cluster) |>
  group_modify(~{
    dd <- .x
    expand_grid(med_metric = med_cols, feature = feature_cols) |>
      mutate(
        test = pmap(list(med_metric, feature), function(m, f) {
          x <- dd[[m]]; y <- dd[[f]]
          ok <- is.finite(x) & is.finite(y)
          if (sum(ok) >= 3) suppressWarnings(cor.test(x[ok], y[ok], method = "spearman")) else NULL
        }),
        rho = map_dbl(test, ~ if (is.null(.x)) NA_real_ else unname(.x$estimate)),
        p   = map_dbl(test, ~ if (is.null(.x)) NA_real_ else .x$p.value)
      ) |>
      mutate(p_BH = p.adjust(p, method = "BH")) 
  }) |>
  ungroup()

# Reiksmingu koreliaciju lentele ----------
top_corr <- corr_with_p |>
  filter(!is.na(rho)) |>
  arrange(Cluster, desc(abs(rho))) |>
  filter(!is.na(p_BH), p_BH < 0.1) |>
  group_by(Cluster) |>
  slice_head(n = 30) |>
  ungroup()
print(top_corr)

# Silumos zemelapis vizualizavimui
# Filtruojam klasteri 5 nes ten tik vienas pacientas
corr_long_f <- corr_long |> filter(as.character(Cluster) != "5")
corr_with_p_f <- corr_with_p |> filter(as.character(Cluster) != "5")

top_feats_by_cluster <- corr_long_f |>
  group_by(Cluster, feature) |>
  summarise(max_abs_rho = max(abs(rho), na.rm = TRUE), .groups = "drop") |>
  group_by(Cluster) |>
  slice_max(max_abs_rho, n = 20, with_ties = FALSE) |>
  ungroup()

plot_df <- corr_long_f |>
  semi_join(top_feats_by_cluster, by = c("Cluster","feature")) |>
  left_join(
    corr_with_p_f |> 
      dplyr::select(Cluster, med_metric, feature, p_BH),
      by = c("Cluster","med_metric","feature")
  ) |>
    mutate(
    label  = if_else(is.na(rho), "", sprintf("%.2f", rho)),
    txt_col = case_when(
      is.na(rho) ~ NA_character_,
      abs(rho) >= 0.5 ~ "white",
      TRUE ~ "black"
    )
  )

# Koreliaciju vizualizavimas
ggplot(plot_df, aes(x = feature, y = med_metric, fill = rho)) +
  geom_tile() +
  geom_point(
    data = subset(plot_df, !is.na(p_BH) & p_BH < 0.05),
    shape = 19, size = 1.6, stroke = 0.2, colour = "black"
  ) +
  geom_text(aes(label = label, color = txt_col), size = 3, na.rm = TRUE, show.legend = FALSE) +
  scale_color_identity() + 
  facet_wrap(~ Cluster, scales = "free_x", nrow = 2, ncol = 2) +
  scale_fill_gradient2(
    limits = c(-1, 1), na.value = "grey90",
    name = "ρ"   
  ) +
  labs(
    title = "Spearman'o koreliacijos klasteriuose (be klasterio 5)",
    x = "Išskirtas požymis", y = "Medikų metrika"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title   = element_text(face = "bold"),
    plot.subtitle = ggtext::element_markdown(),
    axis.text.x  = element_text(angle = 45, hjust = 1),
    panel.spacing = unit(1, "lines")
  )

# Bendros koreliacijos ------------
num_df <- nauja_pacientai_features %>% select(where(is.numeric),-Class)
res <- rcorr(as.matrix(num_df), type = "spearman")

ggcorrplot(res$r, 
           p.mat = res$P,          
           method = "square",       
           type = "lower",          
           insig = "blank",         
           lab = TRUE,              
           lab_size = 3,
           colors = c("#6D9EC1", "white", "#E46726"), 
           title = "Spearmano koreliacijos",
           ggtheme = theme_minimal())

#----------------KLASTERIZAVIMAS PAGAL MEDICININES METRIKAS-----------
nauja_pacientai_klasteriai2 <- nauja_pacientai_klasteriai |>
  select(-1, -2, -last_col()) |>
  drop_na()

X_pacientai <- scale(nauja_pacientai_klasteriai2) 
fviz_nbclust(X_pacientai, kmeans, method = "wss")
fviz_nbclust(X_pacientai, kmeans, method = "silhouette")

set.seed(123)
k_best <- 4
km <- kmeans(X_pacientai, centers = k_best, nstart = 50)
distancija <- dist(X_pacientai)
sil <- silhouette(km$cluster, distancija)
summary(sil)

nauja_pacientai_klasteriai_final <- nauja_pacientai_klasteriai2 |>
  mutate(Cluster = km$cluster)

# ----------------KLASTERIZAVIMO PALYGINIMAS-------------
bendri_pacientai <- cbind(nauja_pacientai_klasteriai, Atrinkti_pozymiai_su_paciento_id)
colnames(bendri_pacientai)[7] <- "Cluster_1"   

texture_vars <- c("F_En", "GLCM_Cor", "GLRLM_LRHGE",
                  "GLRLM_RE", "GLRLM_RLNUN", "LBP_Con")

medic_vars <- c("m_Nat T1", "m_ECV", "e_GLS")

nauja_lentele_pacientai <- bendri_pacientai |>
  select(Patient_ID, "F_En", "GLCM_Cor", "GLRLM_LRHGE",
         "GLRLM_RE", "GLRLM_RLNUN", "LBP_Con","m_Nat T1", "m_ECV", "e_GLS",  Cluster)

pacientai_be_NA <- nauja_lentele_pacientai |>
  drop_na()

X_texture <- pacientai_be_NA |>
  select(texture_vars) 

set.seed(123)
fviz_nbclust(X_texture, kmeans, method = "wss")
fviz_nbclust(X_texture, kmeans, method = "silhouette")

# pritaikome klasterizavima
km_texture <- kmeans(X_texture, centers = 4)

distancija <- dist(X_texture)
sil <- silhouette(km_texture$cluster, distancija)
summary(sil)
pacientai_be_NA$Cluster_texture <- as.factor(km_texture$cluster)

# paimam klinikines metrikas ir normalizuojame
X_medic <- pacientai_be_NA %>%
  select('m_Nat T1', m_ECV, e_GLS) %>% scale()

set.seed(123)
fviz_nbclust(X_medic, kmeans, method = "wss")
fviz_nbclust(X_medic, kmeans, method = "silhouette")

km_medic <- kmeans(X_medic, centers = 4)

pacientai_be_NA$Cluster_medic <- as.factor(km_medic$cluster)

# Surenkame visa informacija u klasteriais i viena lentele
cluster_table <- table(pacientai_be_NA$Cluster_texture,
                       pacientai_be_NA$Cluster_medic)

cluster_table
prop.table(cluster_table, margin = 1) * 100

# ------------- PALYGINIMAS SU STATISTINEMIS METRIKOMIS ---------
data.frame(
  ARI_med_text = adjustedRandIndex(pacientai_be_NA$Cluster_texture,
                          pacientai_be_NA$Cluster_medic),
  NMI_med_text = NMI(pacientai_be_NA$Cluster_texture,
            pacientai_be_NA$Cluster_medic),
  AMI_med_text = AMI(pacientai_be_NA$Cluster_texture,
            pacientai_be_NA$Cluster_medic),
  ARI_base_text = adjustedRandIndex(pacientai_be_NA$Cluster_texture,
                          pacientai_be_NA$Cluster),
  NMI_base_text = NMI(pacientai_be_NA$Cluster_texture,
            pacientai_be_NA$Cluster),
  AMI_base_text = AMI(pacientai_be_NA$Cluster_texture,
            pacientai_be_NA$Cluster),
  ARI_base_med = adjustedRandIndex(pacientai_be_NA$Cluster,
                          pacientai_be_NA$Cluster_medic),
  NMI_base_med = NMI(pacientai_be_NA$Cluster,
            pacientai_be_NA$Cluster_medic),
  AMI_base_med = AMI(pacientai_be_NA$Cluster,
            pacientai_be_NA$Cluster_medic)
)

# Vizualizavimas su UMAP
set.seed(41)
custom_colors <- c("royalblue", 
                   "violetred4", 
                   "skyblue1", 
                   "royalblue4")

umap_config <- umap::umap.defaults
umap_config$n_neighbors <- 10
umap_config$min_dist <- 0.05

# UMAP skaiciavimas: Teksturos erdve ---------
umap_tex <- umap::umap(X_texture, config = umap_config)
umap_tex_coords <- as.data.frame(umap_tex$layout)
colnames(umap_tex_coords) <- c("UMAP1_tex", "UMAP2_tex")

# UMAP skaiciavimas: Medicinine erdve ---------
umap_med <- umap::umap(X_medic, config = umap_config)
umap_med_coords <- as.data.frame(umap_med$layout)
colnames(umap_med_coords) <- c("UMAP1_med", "UMAP2_med")

plot_df <- cbind(
  pacientai_be_NA, 
  umap_tex_coords, 
  umap_med_coords
) %>%
  mutate(
    Cluster_texture = as.factor(Cluster_texture),
    Cluster_medic = as.factor(Cluster_medic)
  )

# Vizualizavimas
plot_umap <- function(df, x, y, color_var, title, subtitle = NULL, legend_title = "Cluster") {
  ggplot(df, aes_string(x = x, y = y, color = color_var)) +
    geom_point(size = 2.5, alpha = 0.8) +
    scale_color_manual(values = custom_colors) +
    theme_minimal(base_size = 11) +
    labs(
      title = title,
      subtitle = subtitle,
      color = legend_title,
      x = "UMAP 1",
      y = "UMAP 2"
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      legend.title = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )
}

p1 <- plot_umap(plot_df, "UMAP1_tex", "UMAP2_tex", "Cluster_texture", 
                "UMAP projekcija tekstūros požymių erdvėje")
p1
p2 <- plot_umap(plot_df, "UMAP1_med", "UMAP2_med", "Cluster_medic", 
                "UMAP projekcija klinikinių metrikų erdvėje")
p2


# KORELIACIJU NAGRINEJIMAS MAZOJE IMTYJE ----------
cor_mat <- matrix(NA,
                  nrow = length(texture_vars),
                  ncol = length(medic_vars),
                  dimnames = list(texture_vars, medic_vars))

for (i in texture_vars) {
  for (j in medic_vars) {
    cor_mat[i, j] <- cor(pacientai_be_NA[[i]],
                         pacientai_be_NA[[j]],
                         method = "spearman",
                         use = "complete.obs")
  }
}

round(cor_mat, 2)

# Skirtumu skaiciavimas -------------
df_long_medic <- pacientai_be_NA |>
  select(Cluster_texture, all_of(medic_vars)) |>
  pivot_longer(cols = -Cluster_texture, names_to = "variable", values_to = "value") |>
  drop_na()

# Kruskal-Wallis testas + Efekto dydis 
res_kw <- df_long_medic |>
  group_by(variable) |>
  kruskal_test(value ~ Cluster_texture) |>
  add_significance() |>
  left_join(
    df_long_medic |> group_by(variable) |> kruskal_effsize(value ~ Cluster_texture),
    by = "variable"
  ) |>
  mutate(p_adj_fdr = p.adjust(p, method = "BH")) |>
  arrange(p_adj_fdr)

# POST-HOC Dunn testas
res_posthoc <- df_long_medic |>
  filter(variable %in% (res_kw |> filter(p_adj_fdr < 0.05) |> pull(variable))) |>
  group_by(variable) |>
  dunn_test(value ~ Cluster_texture, p.adjust.method = "BH")
res_posthoc

# Medicininiu metriku pasiskirstymas ---------
pacientai_be_NA |>
  pivot_longer(cols = all_of(medic_vars), names_to = "metric", values_to = "value") |>
  ggplot(aes(x = Cluster_texture, y = value, fill = Cluster_texture)) +
  geom_boxplot(outlier.alpha = 0.4, alpha =0.8) +
  facet_wrap(~ metric, scales = "free_y") +
  theme_minimal(base_size = 10) +
  scale_fill_manual(values = custom_colors) +
  labs(
    title = "Medicininiu metriku pasiskirstymas pagal klasterius",
    x = "Klasteris",
    y = "Reikšmė"
  ) +
  theme(legend.position = "none") +
  theme(
    plot.title = element_text(hjust = 0, face = "bold", size = 12),
    strip.text = element_text(size = 10, face = "bold")
  )

# Koreliaciju vizualizavimas
pheatmap(cor_mat, color = colorRampPalette(c("#4575B4", "white", "#D73027"))(100),
         display_numbers = TRUE, 
         main = "Spermano bendros koreliacijos tarp metrikų")

# Tekstūros požymių skirtumas -----------
df_long_tex <- pacientai_be_NA |>
  select(Cluster_texture, all_of(texture_vars)) |>
  pivot_longer(cols = -Cluster_texture, names_to = "variable", values_to = "value") 

# Kruskal-Wallis 
res_tex_kw <- df_long_tex |>
  group_by(variable) |>
  kruskal_test(value ~ Cluster_texture) |>
  left_join(
    df_long_tex |> group_by(variable) |> kruskal_effsize(value ~ Cluster_texture),
    by = "variable"
  ) |>
  mutate(p_adj_fdr = p.adjust(p, method = "BH")) |>
  arrange(desc(effsize)) 

res_tex_posthoc <- df_long_tex |>
  filter(variable %in% res_tex_kw$variable) |> 
  group_by(variable) |>
  dunn_test(value ~ Cluster_texture, p.adjust.method = "BH")

sig_features <- res_tex_kw |> 
  filter(p_adj_fdr < 0.05) |> 
  pull(variable)

res_tex_posthoc_sig <- df_long_tex |>
  filter(variable %in% sig_features) |>
  group_by(variable) |>
  dunn_test(value ~ Cluster_texture, p.adjust.method = "BH")

pacientai_be_NA |>
  pivot_longer(cols = all_of(sig_features), names_to = "feature", values_to = "value") |>
  ggplot(aes(x = Cluster_texture, y = value, fill = Cluster_texture)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  facet_wrap(~ feature, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = custom_colors) + 
  theme_minimal(base_size = 12) +
  labs(
    title = "Tekstūros požymių pasiskirstymas tarp klasterių",
    x = "Klasteris",
    y = "Reikšmė"
  ) +
  theme(legend.position = "none", strip.text = element_text(face = "bold"))


# REZULTATŲ SUMMARY LENTELĖS ------------------
# Medicinines metrikos
sum_tbl <- pacientai_be_NA %>%
  mutate(Cluster = as.factor(Cluster_texture)) %>%   
  pivot_longer(cols = all_of(medic_vars),
               names_to = "Metrika",
               values_to = "value") %>%
  group_by(Metrika, Cluster) %>%
  summarise(
    n = sum(!is.na(value)),
    mean = mean(value, na.rm = TRUE),
    median = median(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    q1 = quantile(value, 0.25, na.rm = TRUE),
    q3 = quantile(value, 0.75, na.rm = TRUE),
    min = min(value, na.rm = TRUE),
    max = max(value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Metrika, Cluster)
print(sum_tbl)

# Tekstūros požymiai ------------
sum_tbl_texture <- pacientai_be_NA %>%
  mutate(Cluster = as.factor(Cluster_texture)) %>%   # grupuojam pagal texture klasterius
  pivot_longer(cols = all_of(texture_vars),
               names_to = "Metrika",
               values_to = "value") %>%
  group_by(Metrika, Cluster) %>%
  summarise(
    n = sum(!is.na(value)),
    mean = mean(value, na.rm = TRUE),
    median = median(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    q1 = quantile(value, 0.25, na.rm = TRUE),
    q3 = quantile(value, 0.75, na.rm = TRUE),
    min = min(value, na.rm = TRUE),
    max = max(value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Metrika, Cluster)
print(sum_tbl_texture)


# Pacientai
text_lentele <- pacientai_be_NA %>%
  group_by(Cluster_texture) %>%
  summarise(
    Pacientai = paste(Patient_ID, collapse = ", "),
    N = n(),
    .groups = "drop"
  )

med_lentele <- pacientai_be_NA %>%
  group_by(Cluster_medic) %>%
  summarise(
    Pacientai = paste(Patient_ID, collapse = ", "),
    N = n(),
    .groups = "drop"
  )


