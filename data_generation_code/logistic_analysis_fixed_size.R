library(tidyverse)
library(ggplot2)
library(lmtest)

matrix_input_dir = '~/Desktop/data/matrices/'
data_input_dir <- "/Users/yiqunc/Desktop/cse599_f_2019/fixed_logistics/"
output_dir = '~/Desktop/cse599_f_2019/correlation_analysis_ms_fd_size/analysis_input/'

list_of_files <- sort(dir(matrix_input_dir ,pattern = "*.csv"))

subject_name <- unlist(lapply(list_of_files,function(x)strsplit(x,'-')[[1]][[1]]))


sampling_seq <- seq(0.025,0.2,by = 0.025)

for (sampling_ratio in sampling_seq){
    
  risk_diff_unadjusted <- rep(NA, length(list_of_files))
  odd_ratio_log <- rep(NA, length(list_of_files))
  odd_ratio_lcb <- rep(NA, length(list_of_files))
  odd_ratio_ucb <- rep(NA, length(list_of_files))
  cor_ms_fd <- rep(NA, length(list_of_files))
  valid_samples <-rep(NA, length(list_of_files))
  glm_converge  <-rep(NA, length(list_of_files))
  
  
  for (i in seq_along(list_of_files)){
    
    curr_file<-list_of_files[i]
    
    cat('processing ',curr_file,'\n')
    load_data_input <- paste0(data_input_dir,"fixed_size_",sampling_ratio,"_",
                              curr_file,'.RData')
    load(load_data_input)
    
    reg_df <- reg_df[!is.na(as.logical(reg_df$FD)),]
    reg_df$FD <- as.numeric(as.logical(reg_df$FD))
    
    logit_reg_ms <- glm(FD~MS,reg_df,family = binomial(logit))
    #logit_reg_ms_size <- glm(FD~MS+size,reg_df,family = binomial(logit))
    
    #logit_lrt_p_val <- lrtest(logit_reg_ms,logit_reg_ms_size)$`Pr(>Chisq)`[2]
    cor_ms_fd[i] <- cor(reg_df$FD,reg_df$MS)
    
    #lm_ms_size <- lm(FD~MS+size,reg_df)
    
    lm_ms <- lm(FD~MS,reg_df)
    
    #lm_lrt_p_val <- lrtest(lm_ms,lm_ms_size)$`Pr(>Chisq)`[2]
    
    #lrt_fit_logit[i] <- logit_lrt_p_val
    #lrt_fit_lm[i] <-  lm_lrt_p_val
    glm_converge[i] <-logit_reg_ms$converged
    #curr_cor_ms_size <- cor(reg_df$MS,reg_df$size)
    
    #cor_ms_size[i] <- curr_cor_ms_size
    
    #sig_diff_pval_lm[i] <- summary(lm_ms_size)$coefficients['size',4]
    
    #risk_diff_adjusted[i] <- lm_ms_size$coefficients['MS']
    risk_diff_unadjusted[i] <- lm_ms$coefficients['MS']
    
    #model_log_summary <- summary(logit_reg_ms_size)
    model_log_summary <- summary(logit_reg_ms)
    
    if (logit_reg_ms$converged){
      odd_ratio_log[i] <- exp(model_log_summary$coefficients['MS',1])
      odd_ratio_lcb[i] <-   exp(model_log_summary$coefficients['MS',1]-
                                  model_log_summary$coefficients['MS',2]*qnorm(0.975))
      
      odd_ratio_ucb[i] <-   exp(model_log_summary$coefficients['MS',1]+
                                  model_log_summary$coefficients['MS',2]*qnorm(0.975))
    }
  
    
    
  }
  
  
  fixed_size_new_result <- data.frame(subject_name=subject_name,
                                      risk_diff_unadjusted=risk_diff_unadjusted,
                                      odd_ratio_log = odd_ratio_log,
                                      odd_ratio_lcb = odd_ratio_lcb,
                                      odd_ratio_ucb = odd_ratio_ucb,
                                      cor_ms_fd = cor_ms_fd,
                                      glm_converge = glm_converge,
                                      fixed_ratio = sampling_ratio,
                                      list_of_files = list_of_files
  )
  
  
  write.csv(fixed_size_new_result,file = paste0(output_dir,'fixed_size_',
                                                sampling_ratio,'_new_result.csv'))
  
}







