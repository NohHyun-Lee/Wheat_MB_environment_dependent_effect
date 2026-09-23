### library
library(readxl)
library(brms)
library(ggplot2)
library(dplyr)
### data
df_for_more_analysis <- read_xlsx("D:/Microbiome/000.codes/FHB_microbiome_code/FHB_analysis/2025/MB_code_2025/Output/3. obs_pred_graph/2025/2025_v12/df_for_more_analysis_good_20_bad_100.xlsx")
correlation_df <-read_xlsx(paste0("D:/Microbiome/000.codes/FHB_causal_infer_Korea/output/", 
                                  "v1.1", "/correlation_plot/correlation_plot_", 
                                  "v1.1",  "(wth_cor_res_df2).xlsx")) 
network_res <- read_xlsx("D:/Microbiome/000.codes/FHB_causal_infer_Korea/output/v1.1/network/network_edges_genus_high_abudance_v1.1.xlsx")
MB_max_filtered_data <- read_xlsx(path = "D:/Microbiome/000.data/FHB_microbiome_data/FHB_intensity_and_wth_data/For_phyloseq/2025/Filtered_genus_ITSfull_v4.1.xlsx")

### data preproccess
network_res$from <- stringr::str_replace(network_res$from, "g__", "")
network_res$to <- stringr::str_replace(network_res$to, "g__", "")


MB_max_filtered_data <- MB_max_filtered_data[-which(MB_max_filtered_data$Taxa == "unclassified"),]
MB_max_filtered_microbes <- MB_max_filtered_data$Taxa
length(MB_max_filtered_microbes)
names(df_for_more_analysis)


### data preprocessing
# change 0 to 0.001 (because of beta distridution)
min(df_for_CI_selected$incidence[which(df_for_CI_selected$incidence!= 0)])
df_for_CI_selected$incidence <- ifelse(df_for_CI_selected$incidence == 0, 0.00001, df_for_CI_selected$incidence)
df_for_CI_selected$incidence <- ifelse(df_for_CI_selected$incidence == 1, 0.9999, df_for_CI_selected$incidence)
if(max(df_for_CI_selected$incidence) > 1.1){
  df_for_CI_selected$incidence <- df_for_CI_selected$incidence / 100
  print(df_for_CI_selected$incidence)
}else{
  print(df_for_CI_selected$incidence)
}
df_for_CI_selected$temp_flower_sampling_10days
df_scaled <- df_for_CI_selected
df_scaled$temp_flower_sampling_10days

df_scaled[vars_for_scale] <- scale(df_scaled[vars_for_scale])
df_scaled$temp_flower_sampling_10days

# for loop -----------------------------------------------------------------------------------------------------
library(brms)
microbes <- c("Alternaria", "Epicoccum", "Hannaella", "Periconia", "Cladosporium", "Papiliotrema")
network_microbes <- c("Alternaria", "Epicoccum", "Hannaella", "Periconia", "Cladosporium", "Papiliotrema")#MB_max_filtered_microbes
#######====================================================================================================================================================
imsi_save_dir <- "P:/000.codes/MB_causal_inference/output" #"D:/Microbiome/000.codes/FHB_microbiome_code/FHB_analysis/2025/MB_code_2025/Output/3. obs_pred_graph/2025/2025_v12/causal_inference/ITE_ATE_CATE"
ITE_save_version <- "v5.3"
ITE_save_version2 <- "GZ_interaction_v5.3"
high_weather_quantile <- 0.75
low_weather_quantile <- 0.25
microbes_q <- seq(0.1,0.9,0.1)
temp_seq <- c(13,14,15,16)
rhum_seq <- c(60,65,70,75,80)
prcp_seq <- c(80,120,160,200)

incidence_cor <- correlation_df %>% dplyr::filter(inc_or_M == "incidence") %>% dplyr::filter(p_val <= 0.05)
cor_wth_with_inc <- c("temp_flower_sampling_10days", "rhum_flower_sampling_10days", "prcp_after_sampling")  #incidence_cor$wth #wth variables which has correlation with incidence

for(ITE_save_version in c(ITE_save_version)){
  print(ITE_save_version)
  print(ITE_save_version2)
  results <- list()
  ATE_results <- data.frame(matrix(nrow = nrow(df_scaled), ncol = 0))
  CATE_high_results <- data.frame(matrix(nrow = nrow(df_scaled), ncol = 0))
  CATE_low_results <- data.frame(matrix(nrow = nrow(df_scaled), ncol = 0))
  
  for (mi in c(1:length(microbes))) {
    # microbe <- microbes[1]
    #mi=1
    microbe <- microbes[mi]
    print(microbe)
    # # 모델 정의
    ### make model ---------------------------------------------------------------------------------------------------------------------------
    ## select X variable -------------------------------------------------------
    base_terms <- paste0(
      "s(", microbe, ", k=3)"
    )
    cor_terms <- "s(temp_flower_sampling_10days, k=3) + s(rhum_flower_sampling_10days, k=3)"
    
    
    ## network ------------------------------------------------------------------
    network_filter_M <- network_res %>% dplyr::filter(from == microbe | to == microbe)
    imsi_microbes <- network_microbes[-which(network_microbes == microbe)]
    from_num <- which(network_filter_M$from %in% imsi_microbes)
    to_num <- which(network_filter_M$to %in% imsi_microbes)
    
    # interact Ms with imsi microbe
    imsi_selected_microbes <- vector()
    if(length(from_num) != 0){
      imsi_selected_microbes <- c(network_filter_M$from[from_num],imsi_selected_microbes)
    }
    if(length(to_num) != 0){
      imsi_selected_microbes <- c(network_filter_M$to[to_num],imsi_selected_microbes)
    }
    imsi_unique_selected_microbes <- unique(imsi_selected_microbes)
    
    # make formula
    if(length(imsi_unique_selected_microbes) != 0){
      interaction_term <- ""
      for(un in imsi_unique_selected_microbes){
        #un <- imsi_unique_selected_microbes[1]
        one_interaction <- paste0("t2(", microbe, ",", un, ", k=2)")
        interaction_term <- paste0(one_interaction, "+", interaction_term)
      }
    }else{
      interaction_term <- ""
    }
    interaction_term <- sub("\\+\\s*$", "", interaction_term)
    ## final formula -----------------------------------------------------------
    X_formula <- paste0(base_terms, "+", cor_terms, "+", interaction_term)
    X_formula <- sub("\\+\\s*$", "", X_formula)
    
    formula_str <- as.formula(paste0("incidence ~ ", X_formula))
    
    model <- brm(
      formula = formula_str,
      family = Beta(),
      data = df_scaled,
      prior = c(
        prior(normal(0,1), class="b"),
        prior(exponential(1), class="sds")
      ),
      chains = 4,
      iter = 6000,
      warmup = 3000,
      cores = 4,
      control = list(adapt_delta=0.995, max_treedepth=15),
      silent = 2
    )
    
    cat(paste0(microbe, " model finished"))
    
    ##### ATE #####
    
    new_high <- df_scaled
    new_low  <- df_scaled
    
    new_high[[microbe]] <- quantile(df_scaled[[microbe]],0.9)
    new_low[[microbe]]  <- quantile(df_scaled[[microbe]],0.1)
    
    pred_high <- posterior_epred(model,newdata=new_high)
    pred_low  <- posterior_epred(model,newdata=new_low)
    
    ATE_draws <- pred_high - pred_low
    ATE_draw_mean <- rowMeans(ATE_draws)
    ITE_draw_mean <- colMeans(ATE_draws)
    
    # save ITE
    pred_high_df <- as.data.frame(pred_high)
    pred_low_df <- as.data.frame(pred_low)
    pred_high_minus_low <- pred_high_df-pred_low_df
    
    names(pred_high_minus_low) <- names(pred_low_df) <- names(pred_high_df) <- paste0(df_for_more_analysis$wth_ID, "_", df_for_more_analysis$year)
    
    if(!dir.exists(paste0(imsi_save_dir, "/", ITE_save_version))){
      dir.create(paste0(imsi_save_dir, "/", ITE_save_version))
    }
    
    writexl::write_xlsx(pred_low_df, path = paste0(imsi_save_dir, "/", ITE_save_version, "/ITE_results_pred_low_",microbe, "_",ITE_save_version2,".xlsx"))
    writexl::write_xlsx(pred_high_df, path = paste0(imsi_save_dir, "/", ITE_save_version, "/ITE_results_pred_high_",microbe, "_",ITE_save_version2,".xlsx"))
    writexl::write_xlsx(pred_high_minus_low, path = paste0(imsi_save_dir, "/", ITE_save_version, "/ITE_results_pred_high_minus_low_",microbe, "_",ITE_save_version2,".xlsx"))
    
    print(paste0(microbe, " ITE finished"))
    
      
    
    ### CATE -----------------------------------------------------------------------
    # 데이터 복제
    new_cate_high_wth_T1 <- df_scaled
    new_cate_high_wth_T0 <- df_scaled
    
    new_cate_low_wth_T1  <- df_scaled
    new_cate_low_wth_T0  <- df_scaled
    
    # 해당 미생물의 90% quantile vs 10% quantile
    new_cate_high_wth_T1[[microbe]] <- quantile(df_scaled[[microbe]], 0.9, na.rm = TRUE)
    new_cate_high_wth_T0[[microbe]]  <- quantile(df_scaled[[microbe]], 0.1, na.rm = TRUE)
    
    new_cate_low_wth_T1[[microbe]] <- quantile(df_scaled[[microbe]], 0.9, na.rm = TRUE)
    new_cate_low_wth_T0[[microbe]]  <- quantile(df_scaled[[microbe]], 0.1, na.rm = TRUE)
    
    new_cate_high_wth_T1$temp_flower_sampling_10days <-  new_cate_high_wth_T0$temp_flower_sampling_10days <-
      quantile(df_scaled$temp_flower_sampling_10days, high_weather_quantile)
    
    new_cate_high_wth_T1$rhum_flower_sampling_10days <- new_cate_high_wth_T0$rhum_flower_sampling_10days <-
      quantile(df_scaled$rhum_flower_sampling_10days, high_weather_quantile)
    
    new_cate_high_wth_T1$prcp_after_sampling <- new_cate_high_wth_T0$prcp_after_sampling <- 
      quantile(df_scaled$prcp_after_sampling, high_weather_quantile)
    
    new_cate_low_wth_T1$temp_flower_sampling_10days <-  new_cate_low_wth_T0$temp_flower_sampling_10days <-
      quantile(df_scaled$temp_flower_sampling_10days, low_weather_quantile)
    
    new_cate_low_wth_T1$rhum_flower_sampling_10days <- new_cate_low_wth_T0$rhum_flower_sampling_10days <-
      quantile(df_scaled$rhum_flower_sampling_10days, low_weather_quantile)
    
    new_cate_low_wth_T1$prcp_after_sampling <- new_cate_low_wth_T0$prcp_after_sampling <- 
      quantile(df_scaled$prcp_after_sampling, low_weather_quantile)
    
    ## 예측
    # high pred
    pred_high_wth_T1 <- posterior_epred(model, newdata = new_cate_high_wth_T1)
    pred_high_wth_T0 <- posterior_epred(model, newdata = new_cate_high_wth_T0)
    CATE_high_draws <- pred_high_wth_T1 - pred_high_wth_T0
    
    #low pred
    pred_low_wth_T1 <- posterior_epred(model, newdata = new_cate_low_wth_T1)
    pred_low_wth_T0 <- posterior_epred(model, newdata = new_cate_low_wth_T0)
    CATE_low_draws <- pred_low_wth_T1 - pred_low_wth_T0
    
    ## save CITE
    #CATE_high
    pred_high_wth_T1_df <- as.data.frame(pred_high_wth_T1)
    pred_high_wth_T0_df <- as.data.frame(pred_high_wth_T0)
    pred_high_wth_T1_minus_T0 <- pred_high_wth_T1_df - pred_high_wth_T0_df
    
    names(pred_high_wth_T1_minus_T0) <- names(pred_high_wth_T1_df) <- names(pred_high_wth_T0_df) <- paste0(df_for_more_analysis$wth_ID, "_", df_for_more_analysis$year)
    
    writexl::write_xlsx(pred_high_wth_T0_df, path = paste0(imsi_save_dir, "/", ITE_save_version, "/CITE_results_high_wth_T0_",microbe, "_",ITE_save_version2,".xlsx"))
    writexl::write_xlsx(pred_high_wth_T1_df, path = paste0(imsi_save_dir, "/", ITE_save_version, "/CITE_results_high_wth_T1_",microbe, "_",ITE_save_version2,".xlsx"))
    writexl::write_xlsx(pred_high_wth_T1_minus_T0, path = paste0(imsi_save_dir, "/", ITE_save_version, "/CITE_results_high_wth_T1_minus_T0_",microbe, "_",ITE_save_version2,".xlsx"))
    
    #CATE_low
    pred_low_wth_T1_df <- as.data.frame(pred_low_wth_T1)
    pred_low_wth_T0_df <- as.data.frame(pred_low_wth_T0)
    pred_low_wth_T1_minus_T0 <- pred_low_wth_T1_df - pred_low_wth_T0_df
    
    names(pred_low_wth_T1_minus_T0) <- names(pred_low_wth_T1_df) <- names(pred_low_wth_T0_df) <- paste0(df_for_more_analysis$wth_ID, "_", df_for_more_analysis$year)
    
    writexl::write_xlsx(pred_low_wth_T0_df, path = paste0(imsi_save_dir, "/", ITE_save_version, "/CITE_results_low_wth_T0_",microbe, "_",ITE_save_version2,".xlsx"))
    writexl::write_xlsx(pred_low_wth_T1_df, path = paste0(imsi_save_dir, "/", ITE_save_version, "/CITE_results_low_wth_T1_",microbe, "_",ITE_save_version2,".xlsx"))
    writexl::write_xlsx(pred_low_wth_T1_minus_T0, path = paste0(imsi_save_dir, "/", ITE_save_version, "/CITE_results_low_wth_T1_minus_T0_",microbe, "_",ITE_save_version2,".xlsx"))
    
    
    # CATE
    CATE_high <- colMeans(pred_high_wth_T1_minus_T0)
    CATE_high_results <- cbind(CATE_high_results, CATE_high)
    names(CATE_high_results)[mi] <- paste0(microbe, "_CATE_high")
    
    CATE_low <- colMeans(pred_low_wth_T1_minus_T0)
    CATE_low_results <- cbind(CATE_low_results, CATE_low)
    names(CATE_low_results)[mi] <- paste0(microbe, "_CATE_low")
    
    ##### observation 평균 → effect per draw #####
    
    ATE_draw_mean <- rowMeans(ATE_draws)
    CATE_high_draw_mean <- rowMeans(CATE_high_draws)
    CATE_low_draw_mean <- rowMeans(CATE_low_draws)
    
    ITE_draw_mean <- colMeans(ATE_draws)
    CITE_high_draw_mean <- colMeans(CATE_high_draws)
    CITE_low_draw_mean <- colMeans(CATE_low_draws)
    ##### Microbe × Weather surface #####
    microbe_seq <- as.vector(c(quantile(df_scaled[[microbe]], 0.1, na.rm = TRUE), quantile(df_scaled[[microbe]], 0.9, na.rm = TRUE)))
    
    CATE_temp_draws_mean <- data.frame(matrix(nrow = 12000, ncol = 0))
    CITE_temp_draws_mean <- data.frame(matrix(nrow = nrow(df_scaled), ncol = 0))
    for(w in c(1:length(temp_seq))){
      for(m_l in microbe_seq){
        for(m_h in microbe_seq){
          #--------------------------
          if(m_h > m_l){
            new_imsi_cate_temp_T1 <- df_scaled
            new_imsi_cate_temp_T0 <- df_scaled
            
            new_imsi_cate_temp_T1[[microbe]] <- m_h
            new_imsi_cate_temp_T0[[microbe]]  <- m_l
            
            new_imsi_cate_temp_T1$temp_flower_sampling_10days <-  
              new_imsi_cate_temp_T0$temp_flower_sampling_10days <- temp_seq[w]
            
            pred_imsi_temp_wth_T1 <- posterior_epred(model, newdata = new_imsi_cate_temp_T1)
            pred_imsi_temp_wth_T0 <- posterior_epred(model, newdata = new_imsi_cate_temp_T0)
            CATE_imsi_temp_draws <- pred_imsi_temp_wth_T1 - pred_imsi_temp_wth_T0
            
            CATE_imsi_temp_draws_mean <- rowMeans(CATE_imsi_temp_draws)  
            CATE_temp_draws_mean <- cbind(CATE_temp_draws_mean, CATE_imsi_temp_draws_mean)
            names(CATE_temp_draws_mean)[ncol(CATE_temp_draws_mean)] <- paste0("temp_",temp_seq[w]) #paste0("mh_", which(microbe_seq == m_h), "ml_",which(microbe_seq == m_l))
            
            CITE_imsi_temp_draws_mean <- colMeans(CATE_imsi_temp_draws)
            CITE_temp_draws_mean <- cbind(CITE_temp_draws_mean, CITE_imsi_temp_draws_mean)
            names(CITE_temp_draws_mean)[ncol(CITE_temp_draws_mean)] <- paste0("temp_",temp_seq[w])
            
            cat("temp = ",temp_seq[w] ,"m_l = ", m_l, "m_h = ", m_h, "//")
          }else{
            cat("temp = ",temp_seq[w] ,"m_l = ", m_l, "m_h = ", m_h, "//")
          }
          #--------------------------
        }
      }
    }
    
    CATE_rhum_draws_mean<- data.frame(matrix(nrow = 12000, ncol = 0))
    CITE_rhum_draws_mean <- data.frame(matrix(nrow = nrow(df_scaled), ncol = 0))
    for(w in c(1:length(rhum_seq))){
      for(m_l in microbe_seq){
        for(m_h in microbe_seq){
          #--------------------------
          if(m_h > m_l){
            new_imsi_cate_rhum_T1 <- df_scaled
            new_imsi_cate_rhum_T0 <- df_scaled
            
            new_imsi_cate_rhum_T1[[microbe]] <- m_h
            new_imsi_cate_rhum_T0[[microbe]]  <- m_l
            
            new_imsi_cate_rhum_T1$rhum_flower_sampling_10days <-  
              new_imsi_cate_rhum_T0$rhum_flower_sampling_10days <- rhum_seq[w]
            
            pred_imsi_rhum_wth_T1 <- posterior_epred(model, newdata = new_imsi_cate_rhum_T1)
            pred_imsi_rhum_wth_T0 <- posterior_epred(model, newdata = new_imsi_cate_rhum_T0)
            CATE_imsi_rhum_draws <- pred_imsi_rhum_wth_T1 - pred_imsi_rhum_wth_T0
            
            CATE_imsi_rhum_draws_mean <- rowMeans(CATE_imsi_rhum_draws)  
            
            CATE_rhum_draws_mean <- cbind(CATE_rhum_draws_mean, CATE_imsi_rhum_draws_mean)
            names(CATE_rhum_draws_mean)[ncol(CATE_rhum_draws_mean)] <- paste0("rhum_",rhum_seq[w])# paste0("mh_", which(microbe_seq == m_h), "ml_",which(microbe_seq == m_l))
            
            CITE_imsi_rhum_draws_mean <- colMeans(CATE_imsi_rhum_draws)
            CITE_rhum_draws_mean <- cbind(CITE_rhum_draws_mean, CITE_imsi_rhum_draws_mean)
            names(CITE_rhum_draws_mean)[ncol(CITE_rhum_draws_mean)] <- paste0("rhum_",rhum_seq[w])
            
            cat("rhum = ",rhum_seq[w] ,"m_l = ", m_l, "m_h = ", m_h, "//")
          }else{
            cat("rhum = ",rhum_seq[w] ,"m_l = ", m_l, "m_h = ", m_h, "//")
          }
          #--------------------------
        }
      }
    }
    
    CATE_prcp_draws_mean <- data.frame(matrix(nrow = 12000, ncol = 0))
    CITE_prcp_draws_mean <- data.frame(matrix(nrow = nrow(df_scaled), ncol = 0))
    for(w in c(1:length(prcp_seq))){
      for(m_l in microbe_seq){
        for(m_h in microbe_seq){
          #--------------------------
          if(m_h > m_l){
            new_imsi_cate_prcp_T1 <- df_scaled
            new_imsi_cate_prcp_T0 <- df_scaled
            
            new_imsi_cate_prcp_T1[[microbe]] <- m_h
            new_imsi_cate_prcp_T0[[microbe]]  <- m_l
            
            new_imsi_cate_prcp_T1$prcp_after_sampling <-  
              new_imsi_cate_prcp_T0$prcp_after_sampling <- prcp_seq[w]
            
            pred_imsi_prcp_wth_T1 <- posterior_epred(model, newdata = new_imsi_cate_prcp_T1)
            pred_imsi_prcp_wth_T0 <- posterior_epred(model, newdata = new_imsi_cate_prcp_T0)
            CATE_imsi_prcp_draws <- pred_imsi_prcp_wth_T1 - pred_imsi_prcp_wth_T0
            
            CATE_imsi_prcp_draws_mean <- rowMeans(CATE_imsi_prcp_draws)  
            
            CATE_prcp_draws_mean <- cbind(CATE_prcp_draws_mean, CATE_imsi_prcp_draws_mean)
            names(CATE_prcp_draws_mean)[ncol(CATE_prcp_draws_mean)] <- paste0("prcp_",prcp_seq[w])# paste0("mh_", which(microbe_seq == m_h), "ml_",which(microbe_seq == m_l))
            
            CITE_imsi_prcp_draws_mean <- colMeans(CATE_imsi_prcp_draws)
            CITE_prcp_draws_mean <- cbind(CITE_prcp_draws_mean, CITE_imsi_prcp_draws_mean)
            names(CITE_prcp_draws_mean)[ncol(CITE_prcp_draws_mean)] <- paste0("prcp_",prcp_seq[w])
            
            
            cat("prcp = ",prcp_seq[w] ,"m_l = ", m_l, "m_h = ", m_h, "//")
          }else{
            cat("prcp = ",prcp_seq[w] ,"m_l = ", m_l, "m_h = ", m_h, "//")
          }
          #--------------------------
        }
      }
    }
    
    
    
    # 결과 저장
    results[[microbe]] <- list(
      
      model=model,
      
      ATE_draw=ATE_draw_mean,
      CATE_high_draw=CATE_high_draw_mean,
      CATE_low_draw=CATE_low_draw_mean, 
      
      ITE_draw_mean = ITE_draw_mean, 
      CITE_high_draw_mean = CITE_high_draw_mean,
      CITE_low_draw_mean = CITE_low_draw_mean,
      
      CATE_temp_draws_mean= CATE_temp_draws_mean,
      CATE_rhum_draws_mean= CATE_rhum_draws_mean,
      CATE_prcp_draws_mean= CATE_prcp_draws_mean, 
      
      CITE_temp_draws_mean= CITE_temp_draws_mean,
      CITE_rhum_draws_mean= CITE_rhum_draws_mean,
      CITE_prcp_draws_mean= CITE_prcp_draws_mean
    )
    
    cat("// Done:", microbe, "\n")
    
    ATE_CATE_draw_mean <- data.frame(ATE_draw_mean = ATE_draw_mean, 
                                     CATE_high_draw_mean = CATE_high_draw_mean, 
                                     CATE_low_draw_mean = CATE_low_draw_mean)
    ITE_CITE_draw_mean <- data.frame(ITE_draw_mean = ITE_draw_mean, 
                                     CITE_high_draw_mean = CITE_high_draw_mean, 
                                     CITE_low_draw_mean = CITE_low_draw_mean)
    
    writexl::write_xlsx(ATE_CATE_draw_mean, path = paste0(imsi_save_dir, "/", ITE_save_version, "/ATE_CATE_draw_mean_",microbe, "_",ITE_save_version2,".xlsx"))
    writexl::write_xlsx(ITE_CITE_draw_mean, path = paste0(imsi_save_dir, "/", ITE_save_version, "/ITE_CITE_draw_mean_",microbe, "_",ITE_save_version2,".xlsx"))
    writexl::write_xlsx(CATE_temp_draws_mean, path = paste0(imsi_save_dir, "/", ITE_save_version, "/CATE_temp_draws_mean_",microbe, "_",ITE_save_version2,".xlsx"))
    writexl::write_xlsx(CATE_rhum_draws_mean, path = paste0(imsi_save_dir, "/", ITE_save_version, "/CATE_rhum_draws_mean_",microbe, "_",ITE_save_version2,".xlsx"))
    writexl::write_xlsx(CATE_prcp_draws_mean, path = paste0(imsi_save_dir, "/", ITE_save_version, "/CATE_prcp_draws_mean_",microbe, "_",ITE_save_version2,".xlsx"))
  }
  writexl::write_xlsx(ATE_results, path = paste0(imsi_save_dir, "/", ITE_save_version, "/ATE_results_",ITE_save_version2,".xlsx"))
  writexl::write_xlsx(CATE_high_results, path = paste0(imsi_save_dir, "/", ITE_save_version, "/CATE_high_results_",ITE_save_version2,".xlsx"))
  writexl::write_xlsx(CATE_low_results, path = paste0(imsi_save_dir, "/", ITE_save_version, "/CATE_low_results_",ITE_save_version2,".xlsx"))
}

ATE_results
CATE_high_results
CATE_low_results

colMeans(ATE_results)
colMeans(CATE_high_results)
colMeans(CATE_low_results)

# graph
ITE_plot_df <- data.frame()
for(microbe in microbes){
  ITE <- results[[microbe]]$ITE_draw_mean
  CITE_high <- results[[microbe]]$CITE_high_draw_mean
  CITE_low  <- results[[microbe]]$CITE_low_draw_mean
  
  ite_tmp <- data.frame(
    microbe = microbe,
    effect_type = c("ITE","CITE_high","CITE_low"),
    mean = c(mean(ITE),
             mean(CITE_high),
             mean(CITE_low)),
    lower = c(quantile(ITE,0.025),
              quantile(CITE_high,0.025),
              quantile(CITE_low,0.025)),
    upper = c(quantile(ITE,0.975),
              quantile(CITE_high,0.975),
              quantile(CITE_low,0.975))
  )
  
  ITE_plot_df <- rbind(ITE_plot_df,ite_tmp)
}

ITE_g <- ggplot(ITE_plot_df,
                aes(x = mean,
                    y = microbe,
                    color = effect_type)) +
  
  geom_vline(xintercept = 0,
             linetype = "dashed",
             color = "grey40") +
  
  geom_errorbarh(aes(xmin = lower,
                     xmax = upper),
                 height = 0.2,
                 position = position_dodge(width = 0.5),
                 size = 1.2) +
  
  geom_point(size = 4,
             shape = "x",
             position = position_dodge(width = 0.5)) +
  
  scale_color_manual(values = c(
    ITE = "black",
    CITE_high = "#d73027",
    CITE_low = "#4575b4"
  )) +
  
  labs(
    x = "ITE",
    y = "",
    color = ""
  ) +
  theme_bw(base_size = 15) + 
  theme(
    plot.title = element_text(size = 10, face = 'bold'),
    axis.title.x = element_text(size = 15, hjust = 0.5, face = 'bold'),
    axis.title.y = element_text(size = 15, hjust = 0.5, face = 'bold'),
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.4,
                               size = 15, face = 'bold', color = 'black'),
    axis.text.y = element_text(size = 15, face = 'bold', color = 'black')
  )
ITE_g

if(!dir.exists(paste0(imsi_save_dir, "/", ITE_save_version,"/graph/"))){
  dir.create(paste0(imsi_save_dir, "/", ITE_save_version,"/graph/"))
}

ggsave(
  plot = ITE_g,
  file = paste0(
    imsi_save_dir, "/", ITE_save_version,"/graph/",
    "ITE_plot_",
    ITE_save_version2,
    ".png"
  ),
  width = 15,
  height = 10,
  units = c("cm")
)

colMeans(results[["Alternaria"]]$CATE_temp_draws_mean)
colMeans(results[[microbe]]$CATE_temp_draws_mean)

## CATE ---------------------------------------
# temp graph
CITE_temp_plot_df <- data.frame()
for(microbe in microbes){
  ITE <- results[[microbe]]$ITE_draw_mean
  
  imsi_M_CITE_df <- results[[microbe]]$CITE_temp_draws_mean
  CITE_temp_13 <- imsi_M_CITE_df$temp_13 #results[[microbe]]$CITE_high_draw_mean
  CITE_temp_14  <- imsi_M_CITE_df$temp_14 #results[[microbe]]$CITE_low_draw_mean
  CITE_temp_15  <- imsi_M_CITE_df$temp_15
  CITE_temp_16  <- imsi_M_CITE_df$temp_16
  
  ite_tmp <- data.frame(
    microbe = microbe,
    effect_type = c("ITE","temp_13","temp_14", "temp_15","temp_16"),
    mean = c(mean(ITE),
             mean(CITE_temp_13),
             mean(CITE_temp_14), 
             mean(CITE_temp_15),
             mean(CITE_temp_16)
    ),
    lower = c(quantile(ITE,0.025),
              quantile(CITE_temp_13,0.025),
              quantile(CITE_temp_14,0.025),
              quantile(CITE_temp_15,0.025),
              quantile(CITE_temp_16,0.025)
    ),
    upper = c(quantile(ITE,0.975),
              quantile(CITE_temp_13,0.975),
              quantile(CITE_temp_14,0.975),
              quantile(CITE_temp_15,0.975),
              quantile(CITE_temp_16,0.975))
  )
  
  CITE_temp_plot_df <- rbind(CITE_temp_plot_df,ite_tmp)
}

if(max(abs(CITE_temp_plot_df$upper)) < 1){
  CITE_temp_plot_df[,c(3:5)] <- CITE_temp_plot_df[,c(3:5)]*100
}

CITE_temp_plot_df2 <- CITE_temp_plot_df %>% filter(effect_type != "ITE")

CITE_temp_g <- ggplot(CITE_temp_plot_df2,
                      aes(x = mean,
                          y = microbe,
                          color = effect_type)) +
  
  geom_vline(xintercept = 0,
             linetype = "dashed",
             color = "grey40") +
  
  geom_errorbarh(aes(xmin = lower,
                     xmax = upper),
                 height = 0.2,
                 position = position_dodge(width = 0.5),
                 size = 1.2) +
  
  geom_point(size = 4,
             shape = "x",
             position = position_dodge(width = 0.5)) +
  
  scale_color_manual(values = c(
    ITE = "black",
    temp_13 = "#FF8000",
    temp_14 = "#F25656", 
    temp_15 = "#EE1515", 
    temp_16 = "#AE0505"
  )) +
  
  labs(
    x = "ITE (%)",
    y = "",
    color = ""
  ) +
  theme_bw(base_size = 15) + 
  theme(
    plot.title = element_text(size = 10, face = 'bold'),
    axis.title.x = element_text(size = 15, hjust = 0.5, face = 'bold'),
    axis.title.y = element_text(size = 15, hjust = 0.5, face = 'bold'),
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.4,
                               size = 15, face = 'bold', color = 'black'),
    axis.text.y = element_text(size = 15, face = 'bold', color = 'black'),
    legend.position = "bottom",
    legend.direction = "horizontal"
  ) +
  guides(color = guide_legend(nrow = 2))
CITE_temp_g
ggsave(
  plot = CITE_temp_g,
  file = paste0(
    imsi_save_dir, "/", ITE_save_version,"/graph/",
    "CITE_temp_plot_",
    ITE_save_version2,
    "(%).png"
  ),
  width = 15,
  height = 15,
  units = c("cm")
)


# rhum graph -------------------------
CITE_rhum_plot_df <- data.frame()
for(microbe in microbes){
  ITE <- results[[microbe]]$ITE_draw_mean
  
  imsi_M_CITE_df <- results[[microbe]]$CITE_rhum_draws_mean
  CITE_rhum_60  <- imsi_M_CITE_df$rhum_60 #results[[microbe]]$CITE_high_draw_mean
  CITE_rhum_65  <- imsi_M_CITE_df$rhum_65 #results[[microbe]]$CITE_low_draw_mean
  CITE_rhum_70  <- imsi_M_CITE_df$rhum_70
  CITE_rhum_75  <- imsi_M_CITE_df$rhum_75
  CITE_rhum_80  <- imsi_M_CITE_df$rhum_80
  
  ite_tmp <- data.frame(
    microbe = microbe,
    effect_type = c("ITE","rhum_60","rhum_65", "rhum_70","rhum_75", "rhum_80"),
    mean = c(mean(ITE),
             mean(CITE_rhum_60),
             mean(CITE_rhum_65), 
             mean(CITE_rhum_70),
             mean(CITE_rhum_75),
             mean(CITE_rhum_80)
    ),
    lower = c(quantile(ITE,0.025),
              quantile(CITE_rhum_60,0.025),
              quantile(CITE_rhum_65,0.025),
              quantile(CITE_rhum_70,0.025),
              quantile(CITE_rhum_75,0.025),
              quantile(CITE_rhum_80,0.025)
    ),
    upper = c(quantile(ITE,0.975),
              quantile(CITE_rhum_60,0.975),
              quantile(CITE_rhum_65,0.975),
              quantile(CITE_rhum_70,0.975),
              quantile(CITE_rhum_75,0.975),
              quantile(CITE_rhum_80,0.975))
  )
  
  CITE_rhum_plot_df <- rbind(CITE_rhum_plot_df,ite_tmp)
}

if(max(abs(CITE_rhum_plot_df$upper)) < 1){
  CITE_rhum_plot_df[,c(3:5)] <- CITE_rhum_plot_df[,c(3:5)]*100
}
CITE_rhum_plot_df2 <- CITE_rhum_plot_df %>% filter(effect_type != "ITE")
CITE_rhum_g <- ggplot(CITE_rhum_plot_df2,
                      aes(x = mean,
                          y = microbe,
                          color = effect_type)) +
  
  geom_vline(xintercept = 0,
             linetype = "dashed",
             color = "grey40") +
  
  geom_errorbarh(aes(xmin = lower,
                     xmax = upper),
                 height = 0.2,
                 position = position_dodge(width = 0.5),
                 size = 1.2) +
  
  geom_point(size = 4,
             shape = "x",
             position = position_dodge(width = 0.5)) +
  
  scale_color_manual(values = c(
    ITE = "black",
    rhum_60 = "#B6C46F",
    rhum_65 = "#8CC46F", 
    rhum_70 = "#5BAF72", 
    rhum_75 = "#26833F",
    rhum_80 = "#336600"
  )) +
  
  labs(
    x = "ITE (%)",
    y = "",
    color = ""
  ) +
  theme_bw(base_size = 15) + 
  theme(
    plot.title = element_text(size = 10, face = 'bold'),
    axis.title.x = element_text(size = 15, hjust = 0.5, face = 'bold'),
    axis.title.y = element_text(size = 15, hjust = 0.5, face = 'bold'),
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.4,
                               size = 15, face = 'bold', color = 'black'),
    axis.text.y = element_text(size = 15, face = 'bold', color = 'black'),
    legend.position = "bottom",
    legend.direction = "horizontal"
  ) +
  guides(color = guide_legend(nrow = 2))
CITE_rhum_g
ggsave(
  plot = CITE_rhum_g,
  file = paste0(
    imsi_save_dir, "/", ITE_save_version,"/graph/",
    "CITE_rhum_plot_",
    ITE_save_version2,
    "(%).png"
  ),
  width = 15,
  height = 15,
  units = c("cm")
)

## CATE ---------------------------------------
# prcp graph
CITE_prcp_plot_df <- data.frame()
for(microbe in microbes){
  ITE <- results[[microbe]]$ITE_draw_mean
  
  imsi_M_CITE_df <- results[[microbe]]$CITE_prcp_draws_mean
  CITE_prcp_80 <- imsi_M_CITE_df$prcp_80 #results[[microbe]]$CITE_high_draw_mean
  CITE_prcp_120  <- imsi_M_CITE_df$prcp_120 #results[[microbe]]$CITE_low_draw_mean
  CITE_prcp_160  <- imsi_M_CITE_df$prcp_160
  CITE_prcp_200  <- imsi_M_CITE_df$prcp_200
  
  ite_tmp <- data.frame(
    microbe = microbe,
    effect_type = c("ITE","prcp_80","prcp_120", "prcp_160","prcp_200"),
    mean = c(mean(ITE),
             mean(CITE_prcp_80),
             mean(CITE_prcp_120), 
             mean(CITE_prcp_160),
             mean(CITE_prcp_200)
    ),
    lower = c(quantile(ITE,0.025),
              quantile(CITE_prcp_80,0.025),
              quantile(CITE_prcp_120,0.025),
              quantile(CITE_prcp_160,0.025),
              quantile(CITE_prcp_200,0.025)
    ),
    upper = c(quantile(ITE,0.975),
              quantile(CITE_prcp_80,0.975),
              quantile(CITE_prcp_120,0.975),
              quantile(CITE_prcp_160,0.975),
              quantile(CITE_prcp_200,0.975))
  )
  
  CITE_prcp_plot_df <- rbind(CITE_prcp_plot_df,ite_tmp)
}

CITE_prcp_plot_df$effect_type <- factor(CITE_prcp_plot_df$effect_type, levels = c("ITE", "prcp_80", "prcp_120", "prcp_160", "prcp_200"))

if(max(abs(CITE_prcp_plot_df$upper)) < 1){
  CITE_prcp_plot_df[,c(3:5)] <- CITE_prcp_plot_df[,c(3:5)]*100
}
CITE_prcp_plot_df2 <- CITE_prcp_plot_df %>% filter(effect_type != "ITE")
CITE_prcp_g <- ggplot(CITE_prcp_plot_df2,
                      aes(x = mean,
                          y = microbe,
                          color = effect_type)) +
  
  geom_vline(xintercept = 0,
             linetype = "dashed",
             color = "grey40") +
  
  geom_errorbarh(aes(xmin = lower,
                     xmax = upper),
                 height = 0.2,
                 position = position_dodge(width = 0.5),
                 size = 1.2) +
  
  geom_point(size = 4,
             shape = "x",
             position = position_dodge(width = 0.5)) +
  
  scale_color_manual(values = c(
    ITE = "black",
    prcp_80 = "#8CCCF0",
    prcp_120 = "#619EC1", 
    prcp_160 = "#1C6189", 
    prcp_200 = "#0B3E5B"
  )) +
  
  labs(
    x = "ITE (%)",
    y = "",
    color = ""
  ) +
  theme_bw(base_size = 15) + 
  theme(
    plot.title = element_text(size = 10, face = 'bold'),
    axis.title.x = element_text(size = 15, hjust = 0.5, face = 'bold'),
    axis.title.y = element_text(size = 15, hjust = 0.5, face = 'bold'),
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.4,
                               size = 15, face = 'bold', color = 'black'),
    axis.text.y = element_text(size = 15, face = 'bold', color = 'black'),
    legend.position = "bottom",
    legend.direction = "horizontal"
  ) +
  guides(color = guide_legend(nrow = 2))
CITE_prcp_g
ggsave(
  plot = CITE_prcp_g,
  file = paste0(
    imsi_save_dir, "/", ITE_save_version,"/graph/",
    "CITE_prcp_plot_",
    ITE_save_version2,
    "(%).png"
  ),
  width = 15,
  height = 15,
  units = c("cm")
)

t.test(results[["Alternaria"]]$CITE_rhum_draws_mean$rhum_65, results[["Alternaria"]]$CITE_rhum_draws_mean$rhum_60)
t.test(results[["Alternaria"]]$CITE_rhum_draws_mean$rhum_75, results[["Alternaria"]]$CITE_rhum_draws_mean$rhum_60)
t.test(results[["Epicoccum"]]$CITE_rhum_draws_mean$rhum_75, results[["Epicoccum"]]$CITE_rhum_draws_mean$rhum_60) 
t.test(results[["Epicoccum"]]$CITE_prcp_draws_mean$prcp_80, results[["Epicoccum"]]$CITE_prcp_draws_mean$prcp_120) 


### save --------------------------------------------------
results
CITE_temp_plot_df
CITE_rhum_plot_df
CITE_prcp_plot_df
writexl::write_xlsx(CITE_temp_plot_df, path = file.path(imsi_save_dir, ITE_save_version, "CITE_temp_plot_df(%).xlsx") )
writexl::write_xlsx(CITE_rhum_plot_df, path = file.path(imsi_save_dir, ITE_save_version, "CITE_rhum_plot_df(%).xlsx") )
writexl::write_xlsx(CITE_prcp_plot_df, path = file.path(imsi_save_dir, ITE_save_version, "CITE_prcp_plot_df(%).xlsx") )
saveRDS(results, file = file.path(imsi_save_dir, ITE_save_version, "results.rds"))

### analysis
df_analysis <- df_scaled
results <- readRDS(file.path(imsi_save_dir, ITE_save_version, "results.rds"))
results[["Alternaria"]]$ITE_draw_mean
results[["Epicoccum"]]$ITE_draw_mean

results[["Alternaria"]]$CITE_temp_draws_mean$temp_13
plot(results[["Alternaria"]]$ITE_draw_mean, results[["Alternaria"]]$CITE_temp_draws_mean$temp_16)
plot(results[["Epicoccum"]]$ITE_draw_mean, results[["Epicoccum"]]$CITE_temp_draws_mean$temp_16)
plot(results[["Epicoccum"]]$ITE_draw_mean, results[["Epicoccum"]]$CITE_rhum_draws_mean$rhum_80)

plot(results[["Epicoccum"]]$CITE_rhum_draws_mean$rhum_70, results[["Epicoccum"]]$CITE_rhum_draws_mean$rhum_80)


df_analysis$Alt_ITE <- results[["Alternaria"]]$ITE_draw_mean
df_analysis$Epi_ITE <- results[["Epicoccum"]]$ITE_draw_mean
df_analysis$Han_ITE <- results[["Hannaella"]]$ITE_draw_mean
df_analysis$Per_ITE <- results[["Periconia"]]$ITE_draw_mean
df_analysis$Cla_ITE <- results[["Cladosporium"]]$ITE_draw_mean
df_analysis$Pap_ITE <- results[["Papiliotrema"]]$ITE_draw_mean

library(ggplot2)
library(cowplot)
scatter_function_col_incidence <- function(Data, X_var, Y_var, Point_color, Title, X_name, Y_name){
  
  ggplot(Data, aes(x = {{ X_var }}, y = {{ Y_var }})) +
    geom_point(
      aes(color = {{ Point_color }}),
      size = 3
    ) +
    labs(
      title = Title, 
      x = paste0("\n", X_name), 
      y = paste0(Y_name, "\n"),
      color = "Incidence"
    ) +
    scale_color_gradient2(
      low = "#6DD284",
      mid = "#DF3030",
      high = "#000000", 
      midpoint = 0.75,
      limits = c(0, 1)
    ) +
    theme_classic() + 
    theme(
      plot.title = element_text(size = 10, face = "bold"),
      axis.title.x = element_text(size = 15, hjust = 0.5, face = "bold"),
      axis.title.y = element_text(size = 15, hjust = 0.5, face = "bold"),
      axis.text.x = element_text(
        angle = 0, hjust = 0.5, vjust = 0.4,
        size = 15, face = "bold", color = "black"
      ),
      axis.text.y = element_text(
        size = 15, face = "bold", color = "black"
      )
    )
}

scatter_function_col_temp <- function(Data, X_var, Y_var, Point_color, Title, X_name, Y_name){
  
  ggplot(Data, aes(x = {{ X_var }}, y = {{ Y_var }})) +
    geom_point(
      aes(color = {{ Point_color }}),
      size = 3
    ) +
    labs(
      title = Title, 
      x = paste0("\n", X_name), 
      y = paste0(Y_name, "\n"),
      color = "Temperature"
    ) +
    scale_color_gradient2(
      low = "#FFFFCC",
      mid = "#FF0000",
      high = "black", 
      midpoint = 14
    ) +
    theme_classic() + 
    theme(
      plot.title = element_text(size = 10, face = "bold"),
      axis.title.x = element_text(size = 15, hjust = 0.5, face = "bold"),
      axis.title.y = element_text(size = 15, hjust = 0.5, face = "bold"),
      axis.text.x = element_text(
        angle = 0, hjust = 0.5, vjust = 0.4,
        size = 15, face = "bold", color = "black"
      ),
      axis.text.y = element_text(
        size = 15, face = "bold", color = "black"
      )
    )
}

scatter_function_col_rhum <- function(Data, X_var, Y_var, Point_color, Title, X_name, Y_name){
  
  ggplot(Data, aes(x = {{ X_var }}, y = {{ Y_var }})) +
    geom_point(
      aes(color = {{ Point_color }}),
      size = 3
    ) +
    labs(
      title = Title, 
      x = paste0("\n", X_name), 
      y = paste0(Y_name, "\n"),
      color = "Relative\nhumidity"
    ) +
    scale_color_gradient(
      low = "#BEE6E9",
      high = "#02134A"
    ) +
    theme_classic() + 
    theme(
      plot.title = element_text(size = 10, face = "bold"),
      axis.title.x = element_text(size = 15, hjust = 0.5, face = "bold"),
      axis.title.y = element_text(size = 15, hjust = 0.5, face = "bold"),
      axis.text.x = element_text(
        angle = 0, hjust = 0.5, vjust = 0.4,
        size = 15, face = "bold", color = "black"
      ),
      axis.text.y = element_text(
        size = 15, face = "bold", color = "black"
      )
    )
}

scatter_function_col_M <- function(
    Data, X_var, Y_var, Point_color, 
    Title, X_name, Y_name, Point_color_name){
  
  midpoint_value <- mean(
    Data[[rlang::as_name(rlang::ensym(Point_color))]], 
    na.rm = TRUE
  )
  
  ggplot(Data, aes(x = {{ X_var }}, y = {{ Y_var }})) +
    geom_point(
      aes(color = {{ Point_color }}),
      size = 3
    ) +
    labs(
      title = Title, 
      x = paste0("\n", X_name), 
      y = paste0(Y_name, "\n"),
      color = Point_color_name
    ) +
    scale_color_gradient2(
      low = "darkblue",
      mid = "white",
      high = "#CC0000",
      midpoint = 0#midpoint_value
    ) +
    theme_classic() + 
    theme(
      plot.title = element_text(size = 10, face = "bold"),
      axis.title.x = element_text(size = 15, hjust = 0.5, face = "bold"),
      axis.title.y = element_text(size = 15, hjust = 0.5, face = "bold"),
      axis.text.x = element_text(
        angle = 0, hjust = 0.5, vjust = 0.4,
        size = 15, face = "bold", color = "black"
      ),
      axis.text.y = element_text(
        size = 15, face = "bold", color = "black"
      )
    )
}

scatter_function_inc_vs_M <- function(Data, X_var, Y_var, Title, X_name, Y_name){
  
  # Convert variable names to strings
  x_var <- rlang::as_string(rlang::ensym(X_var))
  y_var <- rlang::as_string(rlang::ensym(Y_var))
  
  # Pearson correlation
  cor_result <- cor.test(
    Data[[x_var]],
    Data[[y_var]],
    method = "pearson"
  )
  
  r_value <- cor_result$estimate
  p_value <- cor_result$p.value
  
  # Title
  Title_cor <- paste0(
    Title,
    " (Pearson r = ", round(r_value, 2),
    ", p = ", format.pval(p_value, digits = 2, eps = 0.001),
    ")"
  )
  
  ggplot(Data, aes(x = {{ X_var }}, y = {{ Y_var }})) +
    geom_point(
      color = "#202020",
      size = 3
    ) +
    geom_smooth(
      method = "lm",
      formula = y ~ x,
      se = TRUE,
      color = "#0A0ACE",
      linewidth = 1
    ) +
    labs(
      title = Title_cor,
      x = paste0("\n", X_name),
      y = paste0(Y_name, "\n")
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(size = 10, face = "bold"),
      axis.title.x = element_text(
        size = 15, hjust = 0.5, face = "bold"
      ),
      axis.title.y = element_text(
        size = 15, hjust = 0.5, face = "bold"
      ),
      axis.text.x = element_text(
        angle = 0,
        hjust = 0.5,
        vjust = 0.4,
        size = 15,
        face = "bold",
        color = "black"
      ),
      axis.text.y = element_text(
        size = 15,
        face = "bold",
        color = "black"
      )
    )
}

scatter_function_col_red <- function(Data, X_var, Y_var, Point_color, Title, X_name, Y_name){
  
  ggplot(Data, aes(x = {{ X_var }}, y = {{ Y_var }})) +
    geom_point(
      # aes(color = "black"),
      color = "black",
      size = 3
    ) +
    geom_smooth(
      method = "gam",
      formula = y ~ s(x, k = 4),
      se = TRUE,
      color = "#F94545",
      linewidth = 1.3
    ) +
    labs(
      title = Title, 
      x = "Temp flower (℃)" ,#paste0(X_name), 
      y = paste0(Y_name)
    ) +
    theme_classic() + 
    theme(
      plot.title = element_text(size = 10, face = "bold"),
      axis.title.x = element_text(size = 15, hjust = 0.5, face = "bold"),
      axis.title.y = element_text(size = 15, hjust = 0.5, face = "bold"),
      axis.text.x = element_text(
        angle = 0, hjust = 0.5, vjust = 0.4,
        size = 15, face = "bold", color = "black"
      ),
      axis.text.y = element_text(
        size = 15, face = "bold", color = "black"
      ),
      plot.margin = ggplot2::margin(15, 15, 15, 15)
    )
}
scatter_function_col_blue <- function(Data, X_var, Y_var, Point_color, Title, X_name, Y_name){
  
  ggplot(Data, aes(x = {{ X_var }}, y = {{ Y_var }})) +
    geom_point(
      # aes(color = "black"),
      color = "black",
      size = 3
    ) +
    geom_smooth(
      method = "gam",
      formula = y ~ s(x, k = 4),
      se = TRUE,
      color = "#218721",
      linewidth = 1.3
    ) +
    labs(
      title = Title, 
      x = "Rhum flower (%)" ,# paste0(X_name), 
      y = paste0(Y_name)
    ) +
    theme_classic() + 
    theme(
      plot.title = element_text(size = 10, face = "bold"),
      axis.title.x = element_text(size = 15, hjust = 0.5, face = "bold"),
      axis.title.y = element_text(size = 15, hjust = 0.5, face = "bold"),
      axis.text.x = element_text(
        angle = 0, hjust = 0.5, vjust = 0.4,
        size = 15, face = "bold", color = "black"
      ),
      axis.text.y = element_text(
        size = 15, face = "bold", color = "black"
      ),
      plot.margin = ggplot2::margin(15, 15, 15, 15) # 상, 우, 하, 좌 여백 추가 (단위: pt)
    )
}

Temp_Alt_inc_g <- scatter_function_col_red(Data = df_analysis, X_var = temp_flower_sampling_10days, Y_var = Alt_ITE, Point_color = incidence,
                               Title = "Temp : Alternaria", X_name = "Temp_flower", Y_name = "Alternaria ITE")
Temp_Epi_inc_g <- scatter_function_col_red(Data = df_analysis, X_var = temp_flower_sampling_10days, Y_var = Epi_ITE, Point_color = incidence,
                               Title = "Temp : Epicoccum", X_name = "Temp_flower", Y_name = "Epicoccum ITE")
Temp_Han_inc_g <- scatter_function_col_red(Data = df_analysis, X_var = temp_flower_sampling_10days, Y_var = Han_ITE, Point_color = incidence,
                               Title = "Temp : Hannaella", X_name = "Temp_flower", Y_name = "Hannaella ITE")
Temp_Per_inc_g <- scatter_function_col_red(Data = df_analysis, X_var = temp_flower_sampling_10days, Y_var = Per_ITE, Point_color = incidence,
                               Title = "Temp : Periconia", X_name = "Temp_flower", Y_name = "Periconia ITE")
Temp_Cla_inc_g <- scatter_function_col_red(Data = df_analysis, X_var = temp_flower_sampling_10days, Y_var = Cla_ITE, Point_color = incidence,
                                   Title = "Temp : Cladosporium", X_name = "Temp_flower", Y_name = "Cladosporium ITE")
Temp_Pap_inc_g <- scatter_function_col_red(Data = df_analysis, X_var = temp_flower_sampling_10days, Y_var = Pap_ITE, Point_color = incidence,
                                   Title = "Temp : Papiliotrema", X_name = "Temp_flower", Y_name = "Papiliotrema ITE")

plot(df_analysis$temp_flower_sampling_10days,df_analysis$Cla_ITE)
plot(df_analysis$temp_flower_sampling_10days,df_analysis$Han_ITE)
max(df_analysis$Per_ITE)

Temp_col_inc_g <- plot_grid(
  Temp_Alt_inc_g,
  Temp_Cla_inc_g,
  Temp_Epi_inc_g,
  Temp_Han_inc_g,
  Temp_Pap_inc_g,
  Temp_Per_inc_g,
  ncol = 3,
  nrow = 2,
  align = "hv"
)

ggsave(plot = Temp_Alt_inc_g, filename = paste0(imsi_save_dir,"/", ITE_save_version, "/graph/", "Temp_Alt_inc_g" ,".png"), height = 4, width = 5)
ggsave(plot = Temp_Epi_inc_g, filename = paste0(imsi_save_dir,"/", ITE_save_version, "/graph/", "Temp_Epi_inc_g" ,".png"), height = 4, width = 5)
ggsave(plot = Temp_Han_inc_g, filename = paste0(imsi_save_dir,"/", ITE_save_version, "/graph/", "Temp_Han_inc_g" ,".png"), height = 4, width = 5)
ggsave(plot = Temp_Per_inc_g, filename = paste0(imsi_save_dir,"/", ITE_save_version, "/graph/", "Temp_Per_inc_g" ,".png"), height = 4, width = 5)
ggsave(plot = Temp_Cla_inc_g, filename = paste0(imsi_save_dir,"/", ITE_save_version, "/graph/", "Temp_Cla_inc_g" ,".png"), height = 4, width = 5)
ggsave(plot = Temp_Pap_inc_g, filename = paste0(imsi_save_dir,"/", ITE_save_version, "/graph/", "Temp_Pap_inc_g" ,".png"), height = 4, width = 5)
ggsave(plot = Temp_col_inc_g, filename = paste0(imsi_save_dir,"/", ITE_save_version, "/graph/", "Temp_all_inc_g" ,".png"), height = 6, width = 11)


rhum_Alt_inc_g <- scatter_function_col_blue(Data = df_analysis, X_var = rhum_flower_sampling_10days, Y_var = Alt_ITE, Point_color = incidence,
                               Title = "Rhum : Alternaria", X_name = "rhum_flower", Y_name = "Alternaria ITE")
rhum_Epi_inc_g <- scatter_function_col_blue(Data = df_analysis, X_var = rhum_flower_sampling_10days, Y_var = Epi_ITE, Point_color = incidence,
                               Title = "Rhum : Epicoccum", X_name = "rhum_flower", Y_name = "Epicoccum ITE")
rhum_Han_inc_g <- scatter_function_col_blue(Data = df_analysis, X_var = rhum_flower_sampling_10days, Y_var = Han_ITE, Point_color = incidence,
                               Title = "Rhum : Hannaella", X_name = "rhum_flower", Y_name = "Hannaella ITE")
rhum_Per_inc_g <- scatter_function_col_blue(Data = df_analysis, X_var = rhum_flower_sampling_10days, Y_var = Per_ITE, Point_color = incidence,
                               Title = "Rhum : Periconia", X_name = "rhum_flower", Y_name = "Periconia ITE")
rhum_Cla_inc_g <- scatter_function_col_blue(Data = df_analysis, X_var = rhum_flower_sampling_10days, Y_var = Cla_ITE, Point_color = incidence,
                                   Title = "Rhum : Cladosporium", X_name = "rhum_flower", Y_name = "Cladosporium ITE")
rhum_Pap_inc_g <- scatter_function_col_blue(Data = df_analysis, X_var = rhum_flower_sampling_10days, Y_var = Pap_ITE, Point_color = incidence,
                                   Title = "Rhum : Papiliotrema", X_name = "rhum_flower", Y_name = "Papiliotrema ITE")
rhum_col_inc_g <- plot_grid(
  rhum_Alt_inc_g,
  rhum_Cla_inc_g,
  rhum_Epi_inc_g,
  rhum_Han_inc_g,
  rhum_Pap_inc_g,
  rhum_Per_inc_g,
  ncol = 3,
  nrow = 2,
  align = "hv"
)
ggsave(plot = rhum_Alt_inc_g, filename = paste0(imsi_save_dir,"/", ITE_save_version, "/graph/", "rhum_Alt_inc_g" ,".png"), height = 4, width = 5)
ggsave(plot = rhum_Epi_inc_g, filename = paste0(imsi_save_dir,"/", ITE_save_version, "/graph/", "rhum_Epi_inc_g" ,".png"), height = 4, width = 5)
ggsave(plot = rhum_Han_inc_g, filename = paste0(imsi_save_dir,"/", ITE_save_version, "/graph/", "rhum_Han_inc_g" ,".png"), height = 4, width = 5)
ggsave(plot = rhum_Per_inc_g, filename = paste0(imsi_save_dir,"/", ITE_save_version, "/graph/", "rhum_Per_inc_g" ,".png"), height = 4, width = 5)
ggsave(plot = rhum_Cla_inc_g, filename = paste0(imsi_save_dir,"/", ITE_save_version, "/graph/", "rhum_Cla_inc_g" ,".png"), height = 4, width = 5)
ggsave(plot = rhum_Pap_inc_g, filename = paste0(imsi_save_dir,"/", ITE_save_version, "/graph/", "rhum_Pap_inc_g" ,".png"), height = 4, width = 5)
ggsave(plot = rhum_col_inc_g, filename = paste0(imsi_save_dir,"/", ITE_save_version, "/graph/", "rhum_all_inc_g" ,".png"), height = 6, width = 11)


Temp_Alt_col_rhum_g <- scatter_function_col_rhum(Data = df_analysis, X_var = temp_flower_sampling_10days, Y_var = Alt_ITE, Point_color = rhum_flower_sampling_10days,
                                   Title = "Temp : Alternaria", X_name = "Temp_flower", Y_name = "Alternaria ITE")
Temp_Epi_col_rhum_g <- scatter_function_col_rhum(Data = df_analysis, X_var = temp_flower_sampling_10days, Y_var = Epi_ITE, Point_color = rhum_flower_sampling_10days,
                                   Title = "Temp : Epicoccum", X_name = "Temp_flower", Y_name = "Epicoccum ITE")
Temp_Han_col_rhum_g <- scatter_function_col_rhum(Data = df_analysis, X_var = temp_flower_sampling_10days, Y_var = Han_ITE, Point_color = rhum_flower_sampling_10days,
                                   Title = "Temp : Hannaella", X_name = "Temp_flower", Y_name = "Hannaella ITE")
Temp_Per_col_rhum_g <- scatter_function_col_rhum(Data = df_analysis, X_var = temp_flower_sampling_10days, Y_var = Per_ITE, Point_color = rhum_flower_sampling_10days,
                                   Title = "Temp : Periconia", X_name = "Temp_flower", Y_name = "Periconia ITE")
Temp_Cla_col_rhum_g <- scatter_function_col_rhum(Data = df_analysis, X_var = temp_flower_sampling_10days, Y_var = Cla_ITE, Point_color = rhum_flower_sampling_10days,
                                   Title = "Temp : Cladosporium", X_name = "Temp_flower", Y_name = "Cladosporium ITE")
Temp_Pap_col_rhum_g <- scatter_function_col_rhum(Data = df_analysis, X_var = temp_flower_sampling_10days, Y_var = Pap_ITE, Point_color = rhum_flower_sampling_10days,
                                   Title = "Temp : Papiliotrema", X_name = "Temp_flower", Y_name = "Papiliotrema ITE")
Temp_col_rhum_g <- plot_grid(
  Temp_Alt_col_rhum_g,
  Temp_Epi_col_rhum_g,
  Temp_Han_col_rhum_g,
  Temp_Per_col_rhum_g,
  Temp_Cla_col_rhum_g,
  Temp_Pap_col_rhum_g,
  ncol = 3,
  nrow = 2,
  align = "hv"
)
ggsave(plot = Temp_Alt_col_rhum_g, filename = paste0(imsi_save_dir,"/", ITE_save_version, "/graph/", "Temp_Alt_col_rhum_g" ,".png"), height = 4, width = 5)
ggsave(plot = Temp_Epi_col_rhum_g, filename = paste0(imsi_save_dir,"/", ITE_save_version, "/graph/", "Temp_Epi_col_rhum_g" ,".png"), height = 4, width = 5)
ggsave(plot = Temp_Han_col_rhum_g, filename = paste0(imsi_save_dir,"/", ITE_save_version, "/graph/", "Temp_Han_col_rhum_g" ,".png"), height = 4, width = 5)
ggsave(plot = Temp_Per_col_rhum_g, filename = paste0(imsi_save_dir,"/", ITE_save_version, "/graph/", "Temp_Per_col_rhum_g" ,".png"), height = 4, width = 5)
ggsave(plot = Temp_Cla_col_rhum_g, filename = paste0(imsi_save_dir,"/", ITE_save_version, "/graph/", "Temp_Cla_col_rhum_g" ,".png"), height = 4, width = 5)
ggsave(plot = Temp_Pap_col_rhum_g, filename = paste0(imsi_save_dir,"/", ITE_save_version, "/graph/", "Temp_Pap_col_rhum_g" ,".png"), height = 4, width = 5)
ggsave(plot = Temp_col_rhum_g, filename = paste0(imsi_save_dir,"/", ITE_save_version, "/graph/", "Temp_all_rhum_g" ,".png"), height = 8, width = 16)



rhum_Alt_col_temp_g <- scatter_function_col_temp(Data = df_analysis, X_var = rhum_flower_sampling_10days, Y_var = Alt_ITE, Point_color = temp_flower_sampling_10days,
                                   Title = "Rhum : Alternaria", X_name = "rhum_flower", Y_name = "Alternaria ITE")
rhum_Epi_col_temp_g <- scatter_function_col_temp(Data = df_analysis, X_var = rhum_flower_sampling_10days, Y_var = Epi_ITE, Point_color = temp_flower_sampling_10days,
                                   Title = "Rhum : Epicoccum", X_name = "rhum_flower", Y_name = "Epicoccum ITE")
rhum_Han_col_temp_g <- scatter_function_col_temp(Data = df_analysis, X_var = rhum_flower_sampling_10days, Y_var = Han_ITE, Point_color = temp_flower_sampling_10days,
                                   Title = "Rhum : Hannaella", X_name = "rhum_flower", Y_name = "Hannaella ITE")
rhum_Per_col_temp_g <- scatter_function_col_temp(Data = df_analysis, X_var = rhum_flower_sampling_10days, Y_var = Per_ITE, Point_color = temp_flower_sampling_10days,
                                   Title = "Rhum : Periconia", X_name = "rhum_flower", Y_name = "Periconia ITE")
rhum_Cla_col_temp_g <- scatter_function_col_temp(Data = df_analysis, X_var = rhum_flower_sampling_10days, Y_var = Cla_ITE, Point_color = temp_flower_sampling_10days,
                                   Title = "Rhum : Cladosporium", X_name = "rhum_flower", Y_name = "Cladosporium ITE")
rhum_Pap_col_temp_g <- scatter_function_col_temp(Data = df_analysis, X_var = rhum_flower_sampling_10days, Y_var = Per_ITE, Point_color = temp_flower_sampling_10days,
                                   Title = "Rhum : Papiliotrema", X_name = "rhum_flower", Y_name = "Papiliotrema ITE")
rhum_col_temp_g <- plot_grid(
  rhum_Alt_col_temp_g,
  rhum_Epi_col_temp_g,
  rhum_Han_col_temp_g,
  rhum_Per_col_temp_g,
  rhum_Cla_col_temp_g,
  rhum_Pap_col_temp_g,
  ncol = 3,
  nrow = 2,
  align = "hv"
)
ggsave(plot = rhum_Alt_col_temp_g, filename = paste0(imsi_save_dir,"/", ITE_save_version, "/graph/", "rhum_Alt_col_temp_g" ,".png"), height = 4, width = 5)
ggsave(plot = rhum_Epi_col_temp_g, filename = paste0(imsi_save_dir,"/", ITE_save_version, "/graph/", "rhum_Epi_col_temp_g" ,".png"), height = 4, width = 5)
ggsave(plot = rhum_Han_col_temp_g, filename = paste0(imsi_save_dir,"/", ITE_save_version, "/graph/", "rhum_Han_col_temp_g" ,".png"), height = 4, width = 5)
ggsave(plot = rhum_Per_col_temp_g, filename = paste0(imsi_save_dir,"/", ITE_save_version, "/graph/", "rhum_Per_col_temp_g" ,".png"), height = 4, width = 5)
ggsave(plot = rhum_Cla_col_temp_g, filename = paste0(imsi_save_dir,"/", ITE_save_version, "/graph/", "rhum_Cla_col_temp_g" ,".png"), height = 4, width = 5)
ggsave(plot = rhum_Pap_col_temp_g, filename = paste0(imsi_save_dir,"/", ITE_save_version, "/graph/", "rhum_Pap_col_temp_g" ,".png"), height = 4, width = 5)
ggsave(plot = rhum_col_temp_g, filename = paste0(imsi_save_dir,"/", ITE_save_version, "/graph/", "rhum_all_temp_g" ,".png"), height = 8, width = 16)


rhum_temp_col_Alt_g <- scatter_function_col_M(Data = df_analysis, X_var = rhum_flower_sampling_10days, Y_var = temp_flower_sampling_10days, Point_color = Alt_ITE,
                                                 Title = "Alternaria", X_name = "Rhum_flower", Y_name = "Temp_flower", Point_color_name = "Alt_ITE")
rhum_temp_col_Epi_g <- scatter_function_col_M(Data = df_analysis, X_var = rhum_flower_sampling_10days, Y_var = temp_flower_sampling_10days, Point_color = Epi_ITE,
                                              Title = "Epicoccum", X_name = "Rhum_flower", Y_name = "Temp_flower", Point_color_name = "Epi_ITE")
rhum_temp_col_Han_g <- scatter_function_col_M(Data = df_analysis, X_var = rhum_flower_sampling_10days, Y_var = temp_flower_sampling_10days, Point_color = Han_ITE,
                                              Title = "Hannaella", X_name = "Rhum_flower", Y_name = "Temp_flower", Point_color_name = "Han_ITE")
rhum_temp_col_Per_g <- scatter_function_col_M(Data = df_analysis, X_var = rhum_flower_sampling_10days, Y_var = temp_flower_sampling_10days, Point_color = Per_ITE,
                                              Title = "Periconia", X_name = "Rhum_flower", Y_name = "Temp_flower", Point_color_name = "Per_ITE")
rhum_temp_col_Cla_g <- scatter_function_col_M(Data = df_analysis, X_var = rhum_flower_sampling_10days, Y_var = temp_flower_sampling_10days, Point_color = Cla_ITE,
                                              Title = "Cladosporium", X_name = "Rhum_flower", Y_name = "Temp_flower", Point_color_name = "Cla_ITE")
rhum_temp_col_Pap_g <- scatter_function_col_M(Data = df_analysis, X_var = rhum_flower_sampling_10days, Y_var = temp_flower_sampling_10days, Point_color = Pap_ITE,
                                              Title = "Papiliotrema", X_name = "Rhum_flower", Y_name = "Temp_flower", Point_color_name = "Pap_ITE")
rhum_col_M_g <- plot_grid(
  rhum_temp_col_Alt_g,
  rhum_temp_col_Epi_g,
  rhum_temp_col_Han_g,
  rhum_temp_col_Per_g,
  rhum_temp_col_Cla_g,
  rhum_temp_col_Pap_g,
  ncol = 3,
  nrow = 2,
  align = "hv"
)
ggsave(plot = rhum_temp_col_Alt_g, filename = paste0(imsi_save_dir,"/", ITE_save_version, "/graph/", "rhum_temp_col_Alt_g" ,".png"), height = 4, width = 5)
ggsave(plot = rhum_temp_col_Epi_g, filename = paste0(imsi_save_dir,"/", ITE_save_version, "/graph/", "rhum_temp_col_Epi_g" ,".png"), height = 4, width = 5)
ggsave(plot = rhum_temp_col_Han_g, filename = paste0(imsi_save_dir,"/", ITE_save_version, "/graph/", "rhum_temp_col_Han_g" ,".png"), height = 4, width = 5)
ggsave(plot = rhum_temp_col_Per_g, filename = paste0(imsi_save_dir,"/", ITE_save_version, "/graph/", "rhum_temp_col_Per_g" ,".png"), height = 4, width = 5)
ggsave(plot = rhum_temp_col_Cla_g, filename = paste0(imsi_save_dir,"/", ITE_save_version, "/graph/", "rhum_temp_col_Cla_g" ,".png"), height = 4, width = 5)
ggsave(plot = rhum_temp_col_Pap_g, filename = paste0(imsi_save_dir,"/", ITE_save_version, "/graph/", "rhum_temp_col_Pap_g" ,".png"), height = 4, width = 5)
ggsave(plot = rhum_col_M_g, filename = paste0(imsi_save_dir,"/", ITE_save_version, "/graph/", "rhum_all_M_g" ,".png"), height = 8, width = 16)

df_analysis$incidence <- df_analysis$incidence*100
incidence_vs_Alt <- scatter_function_inc_vs_M(Data = df_analysis, X_var = incidence, Y_var = Alt_ITE, Title = "Alternaria", X_name = "Incidence", Y_name = "Alt_ITE")
incidence_vs_Epi <- scatter_function_inc_vs_M(Data = df_analysis, X_var = incidence, Y_var = Epi_ITE, Title = "Epicoccum", X_name = "Incidence", Y_name = "Epi_ITE")
incidence_vs_Han <- scatter_function_inc_vs_M(Data = df_analysis, X_var = incidence, Y_var = Han_ITE, Title = "Hannaella", X_name = "Incidence", Y_name = "Han_ITE")
incidence_vs_Per <- scatter_function_inc_vs_M(Data = df_analysis, X_var = incidence, Y_var = Per_ITE, Title = "Periconia", X_name = "Incidence", Y_name = "Per_ITE")
incidence_vs_Cla <- scatter_function_inc_vs_M(Data = df_analysis, X_var = incidence, Y_var = Cla_ITE, Title = "Cladosporium", X_name = "Incidence", Y_name = "Cla_ITE")
incidence_vs_Pap <- scatter_function_inc_vs_M(Data = df_analysis, X_var = incidence, Y_var = Pap_ITE, Title = "Papiliotrema", X_name = "Incidence", Y_name = "Pap_ITE")
incidence_vs_all_M <- plot_grid(
  incidence_vs_Alt,
  incidence_vs_Cla,
  incidence_vs_Epi,
  incidence_vs_Han,
  incidence_vs_Pap,
  incidence_vs_Per,
  ncol = 3,
  nrow = 2,
  align = "hv"
)
ggsave(plot = incidence_vs_Alt, filename = paste0(imsi_save_dir,"/", ITE_save_version, "/graph/", "incidence_vs_Alt" ,".png"), height = 4, width = 4)
ggsave(plot = incidence_vs_Epi, filename = paste0(imsi_save_dir,"/", ITE_save_version, "/graph/", "incidence_vs_Epi" ,".png"), height = 4, width = 4)
ggsave(plot = incidence_vs_Han, filename = paste0(imsi_save_dir,"/", ITE_save_version, "/graph/", "incidence_vs_Han" ,".png"), height = 4, width = 4)
ggsave(plot = incidence_vs_Per, filename = paste0(imsi_save_dir,"/", ITE_save_version, "/graph/", "incidence_vs_Per" ,".png"), height = 4, width = 4)
ggsave(plot = incidence_vs_Cla, filename = paste0(imsi_save_dir,"/", ITE_save_version, "/graph/", "incidence_vs_Cla" ,".png"), height = 4, width = 4)
ggsave(plot = incidence_vs_Pap, filename = paste0(imsi_save_dir,"/", ITE_save_version, "/graph/", "incidence_vs_Pap" ,".png"), height = 4, width = 4)
ggsave(plot = incidence_vs_all_M, filename = paste0(imsi_save_dir,"/", ITE_save_version, "/graph/", "incidence_vs_all_M" ,".png"), height = 8, width = 14)


##-------------------------------------
rhum_75 <- df_scaled %>% filter(rhum_flower_sampling_10days < 75) %>% filter(rhum_flower_sampling_10days > 70)
rhum_75 <- df_for_more_analysis %>% filter(rhum_flower_sampling_10days < 75) %>% filter(rhum_flower_sampling_10days > 70)
rhum_80 <- df_for_more_analysis %>% filter(rhum_flower_sampling_10days > 75)# %>% filter(year == 2023)
rhum_75$loc
rhum_80$loc

rhum_80 <- df_scaled %>% filter(rhum_flower_sampling_10days > 75) #%>% filter(rhum_flower_sampling_10days > 70)

boxplot(rhum_75$incidence, rhum_80$incidence)
t.test(rhum_75$incidence, rhum_80$incidence)

hist(df_scaled$Epicoccum)
hist(rhum_75_df$Epicoccum)

rhum_75_df_low_Epi <- rhum_75_df %>% filter(Epicoccum < 5) %>% filter(temp_flower_sampling_10days > 15)
rhum_75_df_high_Epi <- rhum_75_df %>% filter(Epicoccum > 5)%>% filter(temp_flower_sampling_10days > 15)
boxplot(rhum_75_df_low_Epi$incidence, rhum_75_df_high_Epi$incidence)
t.test(rhum_75_df_low_Epi$incidence, rhum_75_df_high_Epi$incidence)
plot(df_scaled$Epicoccum, df_scaled$rhum_)

library(corrplot)
cite_rhum <- results[["Epicoccum"]]$CITE_rhum_draws_mean
cite_rhum <- results[["Hannaella"]]$CITE_rhum_draws_mean
cite_rhum <- results[["Periconia"]]$CITE_rhum_draws_mean
cite_rhum <- results[["Cladosporium"]]$CITE_rhum_draws_mean
cite_rhum <- results[["Papiliotrema"]]$CITE_rhum_draws_mean

cor_matrix <- cor(
  cite_rhum,
  method = "pearson",
  use = "complete.obs"
)

corrplot(
  cor_matrix,
  method = "color",
  type = "upper",
  addCoef.col = "black",
  tl.col = "black",
  tl.srt = 0,
  number.cex = 1.0
)

pairs(
  cite_rhum,
  pch = 19,
  cex = 0.7
)

plot(results[["Hannaella"]]$CITE_rhum_draws_mean$rhum_60, results[["Hannaella"]]$CITE_rhum_draws_mean$rhum_80)
plot(results[["Periconia"]]$CITE_rhum_draws_mean$rhum_60, results[["Periconia"]]$CITE_rhum_draws_mean$rhum_80)

### Supplementary figure ---------------------------------------------------------------------------------------------------------
plot(df_for_more_analysis$rhum_flower_sampling_10days, df_for_more_analysis$incidence)
plot(df_for_more_analysis$temp_flower_sampling_10days, df_for_more_analysis$incidence)


cor.test(df_for_more_analysis$temp_flower_sampling_10days, df_for_more_analysis$incidence, method = "spearman")
cor.test(df_for_more_analysis$rhum_flower_sampling_10days, df_for_more_analysis$incidence, method = "spearman")
cor.test(df_for_more_analysis$rhum_flower_sampling_10days, df_for_more_analysis$incidence, method = "pearson")

temp_vs_inc_g <-  ggplot(df_for_more_analysis, aes(x = temp_flower_sampling_10days, y = incidence)) +
  geom_point(
    # aes(color = "black"),
    color = "#DA3B3B",
    size = 3
  ) +
  labs(
    title = "spearman = 0.74 / p < 0.001", 
    x = "\nTemp flower (℃)" ,#paste0(X_name), 
    y = "FHB incidence (%)\n"
  ) +
  theme_classic() + 
  theme(
    plot.title = element_text(size = 10, face = "bold"),
    axis.title.x = element_text(size = 15, hjust = 0.5, face = "bold"),
    axis.title.y = element_text(size = 15, hjust = 0.5, face = "bold"),
    axis.text.x = element_text(
      angle = 0, hjust = 0.5, vjust = 0.4,
      size = 15, face = "bold", color = "black"
    ),
    axis.text.y = element_text(
      size = 15, face = "bold", color = "black"
    )
  )
temp_vs_inc_g

rhum_vs_inc_g <-  ggplot(df_for_more_analysis, aes(x = rhum_flower_sampling_10days, y = incidence)) +
  geom_point(
    # aes(color = "black"),
    color = "#187711",
    size = 3
  ) +
  labs(
    title = "spearman = 0.79 / p < 0.001", 
    x = "\nRhum flower (%)" ,#paste0(X_name), 
    y = "FHB incidence (%)\n"
  ) +
  theme_classic() + 
  theme(
    plot.title = element_text(size = 10, face = "bold"),
    axis.title.x = element_text(size = 15, hjust = 0.5, face = "bold"),
    axis.title.y = element_text(size = 15, hjust = 0.5, face = "bold"),
    axis.text.x = element_text(
      angle = 0, hjust = 0.5, vjust = 0.4,
      size = 15, face = "bold", color = "black"
    ),
    axis.text.y = element_text(
      size = 15, face = "bold", color = "black"
    )
  )
temp_vs_inc_g
rhum_vs_inc_g


temp_hist_g <-  ggplot(df_for_more_analysis, aes(x = temp_flower_sampling_10days)) +
  geom_histogram(bins = 10, fill = "#DA3B3B", color = "white") +
  labs(
    title = "", 
    x = "\nTemp flower (℃)" ,
    y = "Count\n" 
  ) +
  theme_classic() + 
  theme(
    plot.title = element_text(size = 10, face = "bold"),
    axis.title.x = element_text(size = 15, hjust = 0.5, face = "bold"),
    axis.title.y = element_text(size = 15, hjust = 0.5, face = "bold"),
    axis.text.x = element_text(
      angle = 0, hjust = 0.5, vjust = 0.4,
      size = 15, face = "bold", color = "black"
    ),
    axis.text.y = element_text(
      size = 15, face = "bold", color = "black"
    )
  )


rhum_hist_g <-  ggplot(df_for_more_analysis, aes(x = rhum_flower_sampling_10days)) +
  geom_histogram(bins = 10, fill = "#187711", color = "white") +
  labs(
    title = "", 
    x = "\nRhum flower (%)" ,#paste0(X_name), 
    y = "Count\n" 
  ) +
  theme_classic() + 
  theme(
    plot.title = element_text(size = 10, face = "bold"),
    axis.title.x = element_text(size = 15, hjust = 0.5, face = "bold"),
    axis.title.y = element_text(size = 15, hjust = 0.5, face = "bold"),
    axis.text.x = element_text(
      angle = 0, hjust = 0.5, vjust = 0.4,
      size = 15, face = "bold", color = "black"
    ),
    axis.text.y = element_text(
      size = 15, face = "bold", color = "black"
    )
  )


# save
temp_vs_inc_g
rhum_vs_inc_g
temp_hist_g
rhum_hist_g

ggsave(temp_vs_inc_g , filename = "P:/000.codes/MB_causal_inference/output/v5.3/graph/Fig.S_temp_vs_inc_g.png", height = 5, width = 5)
ggsave(rhum_vs_inc_g , filename = "P:/000.codes/MB_causal_inference/output/v5.3/graph/Fig.S_rhum_vs_inc_g.png", height = 5, width = 5)
ggsave(temp_hist_g , filename = "P:/000.codes/MB_causal_inference/output/v5.3/graph/Fig.S_temp_hist_g.png", height = 5, width = 5)
ggsave(rhum_hist_g , filename = "P:/000.codes/MB_causal_inference/output/v5.3/graph/Fig.S_rhum_hist_g.png", height = 5, width = 5)
