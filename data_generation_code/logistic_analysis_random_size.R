library(tidyverse)
library(ggplot2)
library(lmtest)

matrix_input_dir = '~/Desktop/data/matrices/'
output_dir = '~/Desktop/cse599_f_2019/correlation_analysis_ms_fd_size/analysis_input/'

list_of_files <- dir(matrix_input_dir ,pattern = "*.csv")

subject_name <- unlist(lapply(list_of_files,function(x)strsplit(x,'-')[[1]][[1]]))

cor_ms_size <- rep(NA, length(list_of_files))
sig_diff_pval_lm <- rep(NA, length(list_of_files))
risk_diff_adjusted <- rep(NA, length(list_of_files))
risk_diff_unadjusted <- rep(NA, length(list_of_files))
odd_ratio_log <- rep(NA, length(list_of_files))
odd_ratio_lcb <- rep(NA, length(list_of_files))
odd_ratio_ucb <- rep(NA, length(list_of_files))
odd_ratio_p_val <- rep(NA, length(list_of_files))
odds_ratio_unadjust <- rep(NA, length(list_of_files))
odds_ratio_unadjust_p_val <- rep(NA, length(list_of_files))
lrt_fit_logit <- rep(NA, length(list_of_files))
lrt_fit_lm <- rep(NA, length(list_of_files))




cor_size_FD <- rep(NA, length(list_of_files))
partial_cor_size_FD <- rep(NA, length(list_of_files))
cor_MS_FD <- rep(NA, length(list_of_files))
partial_cor_MS_FD <- rep(NA, length(list_of_files))

glm_converged<- rep(NA, length(list_of_files))

for (i in seq_along(list_of_files)){
  
  curr_file<-list_of_files[i]
  
  cat('processing ',curr_file,'\n')
  load_data_input <- paste0("/Users/yiqunc/Desktop/cse599_f_2019/random_logistics/random_0_0_2_",
                            curr_file,'.RData')
  load(load_data_input)
  
  reg_df <- reg_df[!is.na(as.logical(reg_df$FD)),]
  reg_df$FD <- as.numeric(as.logical(reg_df$FD))
  
  logit_reg_ms <- glm(FD~MS,reg_df,family = binomial(logit))
  logit_reg_ms_size <- glm(FD~MS+size,reg_df,family = binomial(logit))
  
  logit_lrt_p_val <- lrtest(logit_reg_ms,logit_reg_ms_size)$`Pr(>Chisq)`[2]

  
  lm_ms_size <- lm(FD~MS+size,reg_df)
  
  lm_ms <- lm(FD~MS,reg_df)
  
  
  mm1 = lm(FD~size,data=reg_df)
  res1 = mm1$residuals
  mm2 = lm(MS~size,data=reg_df)
  res2 = mm2$residuals
  
  partial_cor_MS_FD[i] <- cor(res1,res2)
  cor_MS_FD[i] <- cor(reg_df$FD,reg_df$MS)
  
  
  mm1 = lm(FD~MS,data=reg_df)
  res1 = mm1$residuals
  mm2 = lm(size~MS,data=reg_df)
  res2 = mm2$residuals
  
  partial_cor_size_FD[i] <- cor(res1,res2)
  cor_size_FD[i] <- cor(reg_df$FD,reg_df$size)
  
  lm_lrt_p_val <- lrtest(lm_ms,lm_ms_size)$`Pr(>Chisq)`[2]
  
  lrt_fit_logit[i] <- logit_lrt_p_val
  lrt_fit_lm[i] <-  lm_lrt_p_val
    
  curr_cor_ms_size <- cor(reg_df$MS,reg_df$size)

  cor_ms_size[i] <- curr_cor_ms_size
   
  sig_diff_pval_lm[i] <- summary(lm_ms_size)$coefficients['size',4]
   
  risk_diff_adjusted[i] <- lm_ms_size$coefficients['MS']
  risk_diff_unadjusted[i] <- lm_ms$coefficients['MS']
   
  model_log_summary <- summary(logit_reg_ms_size)
  model_log_ms_summary <- summary(logit_reg_ms)
  
  odd_ratio_log[i] <- exp(model_log_summary$coefficients['MS',1])
  odd_ratio_lcb[i] <-   exp(model_log_summary$coefficients['MS',1]-
                               model_log_summary$coefficients['MS',2]*qnorm(0.975))

   odd_ratio_ucb[i] <-   exp(model_log_summary$coefficients['MS',1]+
                               model_log_summary$coefficients['MS',2]*qnorm(0.975))

   odd_ratio_p_val[i] <- model_log_summary$coefficients['MS',4]
   
   odds_ratio_unadjust[i] <-  exp(model_log_ms_summary$coefficients['MS',1])
   odds_ratio_unadjust_p_val[i] <- model_log_ms_summary$coefficients['MS',4]
   
   glm_converged[i] <- logit_reg_ms_size$converged
  
}




random_size_new_result <- data.frame(subject_name=subject_name,
                                     list_of_files = list_of_files,
                                     cor_ms_size=cor_ms_size,
                              sig_diff_pval_lm=sig_diff_pval_lm, 
                              risk_diff_adjusted = risk_diff_adjusted,
                              risk_diff_unadjusted= risk_diff_unadjusted,
                              lrt_fit_logit = lrt_fit_logit,
                              glm_converged = glm_converged,
                              lrt_fit_lm = lrt_fit_lm,
            odd_ratio_log=odd_ratio_log,
            odd_ratio_lcb=odd_ratio_lcb,odd_ratio_ucb=odd_ratio_ucb,
            odd_ratio_p_val=odd_ratio_p_val,
            odds_ratio_unadjust=odds_ratio_unadjust,
            odds_ratio_unadjust_p_val=odds_ratio_unadjust_p_val,
            cor_size_FD= cor_size_FD,
            partial_cor_size_FD = partial_cor_size_FD,
            cor_MS_FD= cor_MS_FD,
            partial_cor_MS_FD = partial_cor_MS_FD)


write.csv(random_size_new_result,file = paste0(output_dir,'random_size_new_result.csv'))

#### volcano plot 
#### look at how many datasets are improved












