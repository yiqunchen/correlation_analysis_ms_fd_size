library(dplyr)
library(parallel)
#library(khb)

input_path = '/home/yiqun.chen/mutant_score/matrices/'
output_path = '/home/yiqun.chen/mutant_score/logistic_reg/'


# first part we are gonna look at the distributions of different dev

fault_detection_score <- function(input_df){
  # for a kill matrix, this function computes the 
  # number of mutants detected by at least one of 
  # the subsampled tests
  num_mutants_detected <- sum(input_df %>%  select(bug))
  # binarize it to get "at least one"
  return(num_mutants_detected>=1)
}

mutant_score <- function(input_df){
  # for a kill matrix, this function computes the 
  # number of mutants detected by at least one of 
  # the subsampled tests
  num_mutants_detected <- input_df %>%  select(.,starts_with('X')) %>%  colSums()
  # binarize it to get "at least one"
  MS <- sum(num_mutants_detected>0)
  return(MS)
}


list_of_files <- sort(dir(input_path, pattern = "*.csv"))

list_of_df <- vector('list',length = length(list_of_files))


list_of_df <- list_of_df[221:length(list_of_df)]
##### random sampling
set.seed(12345)
N_replica <- 10000

for (i in 1:length(list_of_files)){
  
  
  csv_file <- list_of_files[[i]]
  input_df <- read.csv(paste0(input_path,csv_file))
  ratio_sampling <-  runif(N_replica,min=0,max=0.1)
  sampled_size <- ceil(ratio_sampling*nrow(input_df))
  index_list <- c(1:nrow(input_df))
  sampled_index <- mclapply(sampled_size, function(x)sample(index_list, x,replace = FALSE),
                            mc.cores = 10)
  sampled_df_MS <- unlist(mclapply(sampled_index,function(x)mutant_score(input_df[x,]),
                            mc.cores = 10))
  sampled_df_FD <- unlist(mclapply(sampled_index,
                     function(x)fault_detection_score(input_df[x,]),  mc.cores = 10))
  
  reg_df <- data.frame(FD = sampled_df_FD, MS = sampled_df_MS, size = sampled_size)
  #log_reg <- glm(FD~MS+size,data=reg_df,family='binomial')
  #save(log_reg,csv_file , file =  paste0(output_path,'random_size_sample_',csv_file,'.RData'))
  log_reg_ms_size <- glm(FD~MS+size,data=reg_df,family='binomial')
  log_reg_ms <- glm(FD~MS,data=reg_df,family='binomial')
  
  log_comp_1 <- anova(log_reg_ms,log_reg_ms_size,test = 'Chisq')
  log_comp_p_val <- log_comp_1$`Pr(>Chi)`[2]
  
  lm_reg_ms_size <- lm(FD~MS+size,data=reg_df)
  lm_reg_ms <- lm(FD~MS,data=reg_df)
  
  lm_comp_1 <- anova(lm_reg_ms,lm_reg_ms_size,test = 'Chisq')
  lm_comp_p_val <- lm_comp_1$`Pr(>Chi)`[2]
  
  #khb_test_diff <- khb(log_reg_ms, log_reg_ms_size,med.sandwich=T,glm.sandwich=T)
  
  
  save(log_reg_ms_size,log_reg_ms,log_comp_p_val,lm_reg_ms_size,lm_reg_ms,lm_comp_p_val,csv_file ,
       file =  paste0(output_path,'random_size_sample_',csv_file,'.RData')) 

 
}



