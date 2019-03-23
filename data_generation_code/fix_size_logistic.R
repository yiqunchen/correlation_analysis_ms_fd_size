

library(dplyr)
library(parallel)

input_path = '/home/yiqun.chen/mutant_score/matrices/'
output_path = '/home/yiqun.chen/mutant_score/logistic_reg/fixed_size_sampling/'


# first part we are gonna look at the distributions of different dev

fault_detection_score <- function(input_df){
  # for a kill matrix, this function computes the 
  # number of mutants detected by at least one of 
  # the subsampled tests
  num_faults_detected <- sum(input_df %>%  select(bug))
  # binarize it to get "at least one"
  return(num_faults_detected>=1)
}


mutant_score <- function(input_df){
  # for a kill matrix, this function computes the 
  # number of mutants detected by at least one of 
  # the subsampled tests
  num_mutants_detected <- input_df %>%  select(.,starts_with('X')) %>%  colSums()
  # binarize it to get "at least one"
  MS <- sum(num_mutants_detected>=1)
  return(MS)
}


args <- commandArgs(trailingOnly = TRUE)
csv_file <- args[1]
ratio <- args[2]

#### 
input_df <- read.csv(paste0(input_path,csv_file))
cat('reading',csv_file,'\n')
# Chart-18f.csv

##### random sampling
set.seed(12345)
N_replica <- 10000

#seq(0.025,0.2,by=0.025)

ratio_sampling <- rep(ratio, times = N_replica) 
sampled_size <- ceil(ratio_sampling*nrow(input_df))

index_list <- c(1:nrow(input_df))



sampled_index <- mclapply(sampled_size,
                          function(x)sample(index_list, x,replace = FALSE),
                          mc.cores = 10)

sampled_df_MS <- unlist(mclapply(sampled_index,
                                 function(x)mutant_score(input_df[x,]),
                                 mc.cores = 10))

sampled_df_FD <- unlist(mclapply(sampled_index,
                                 function(x)fault_detection_score(input_df[x,]), 
                                 mc.cores = 10))

reg_df <- data.frame(FD = sampled_df_FD, MS = sampled_df_MS, size = sampled_size)

log_reg_ms <- glm(FD~MS,data=reg_df,family='binomial')

risk_reg_ms <- glm(FD~MS,data=reg_df,family='poisson')

lm_reg_ms <- lm(FD~MS,data=reg_df)



save(reg_df,log_reg_ms,risk_reg_ms,lm_reg_ms,csv_file ,ratio,
     file =  paste0(output_path,'fixed_size_',ratio,'_',csv_file,'.RData')) 





