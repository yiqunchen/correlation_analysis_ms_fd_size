
###### this script takes in  fault_coupling_prob.csv
###### and visualize/tabulate interesting and relevant statistics
##### regarding fault-mutant coupling

# input: run coupling_prob.R to get the csv file


library(tidyverse)
library(latex2exp)

input_dir <- '/Users/yiqunc/Desktop/cse599_f_2019/mutants-faults-revisited-data/'

input_dir = '~/Desktop/cse599_f_2019/correlation_analysis_ms_fd_size/analysis_input/'
plot_ouput_dir = '~/Desktop/cse599_f_2019/correlation_analysis_ms_fd_size/plot_output/'

random_size_plot_df <- read.csv(file = paste0(input_dir,'random_size_new_result.csv'))



coupling_prob_df <- read.csv( file = paste0(input_dir,"fault_coupling_prob.csv"))
coupling_prob_df <- coupling_prob_df %>% column_to_rownames('X')

# first we look at how many couplings are present in each test
# type
task_type <- sapply(strsplit(rownames(coupling_prob_df),'-'),function(x)x[1])

coupling_prob_only <- coupling_prob_df %>% select(starts_with('X'))%>% select_if(is.numeric)
num_of_perfect_couples <- rowSums(coupling_prob_only==1,na.rm = T)
ratio_of_perfect_couples <- num_of_perfect_couples/coupling_prob_df$n_mutants

num_of_perfect_decouples <- rowSums(coupling_prob_only==0,na.rm = T)
ratio_of_perfect_decouples <- num_of_perfect_decouples/coupling_prob_df$n_mutants

max_couples <- apply(coupling_prob_only, 1, function(x)max(x,na.rm = T))


num_good_test <- read.csv(file = paste0(input_dir,'num_good_test.csv'))


coupling_prob_plot <- coupling_prob_only%>%
  mutate(task =task_type )
coupling_prob_plot_gather <- coupling_prob_plot %>%
  gather(key = 'mutant',value = 'cond_prob',starts_with("X"))%>%
  na.omit()

couple_decouple_plot <- data.frame(task=task_type, 
                                   num_of_perfect_couples = num_of_perfect_couples,
                                   ratio_of_perfect_couples=ratio_of_perfect_couples,
                                   num_of_perfect_decouples=num_of_perfect_decouples,
                                   ratio_of_perfect_decouples = ratio_of_perfect_decouples,
                                   n_mutants = coupling_prob_df$n_mutants,
                                   n_good_tests = num_good_test$n_good_test,max_couple_prob = max_couples)

couple_decouple_plot%>% group_by(task) %>% 
  summarise(mean_pcp = mean(ratio_of_perfect_couples),
            med_pcp = median(ratio_of_perfect_couples),
            max_pcp = max(ratio_of_perfect_couples))

max_beta_proj <- c("Chart-25f.csv" ,"Closure-28f.csv", "Lang-13f.csv" ,
                  "Math-70f.csv"  ,  "Time-24f.csv" )
couple_decouple_plot[rownames(couple_decouple_plot)%in%max_beta_proj,] %>% 
  select(ratio_of_perfect_couples)


max_or_proj <- c("Chart-24f.csv" ,"Closure-101f.csv", "Lang-13f.csv" ,
                   "Math-65f.csv"  ,  "Time-21f.csv" )

couple_decouple_plot[rownames(couple_decouple_plot)%in%max_or_proj,] %>% 
  select(ratio_of_perfect_couples)

df_1 <- couple_decouple_plot %>% 
  mutate(list_of_files = rownames(couple_decouple_plot)) %>%
  select(ratio_of_perfect_couples,max_couple_prob,list_of_files,task)

df_2 <- random_size_plot_df %>% select(list_of_files,risk_diff_adjusted,partial_cor_MS_FD,
          glm_converged,odd_ratio_log)%>% left_join(df_1,by='list_of_files')

png(filename = paste0(plot_ouput_dir,"couple_reg_coef.png"),
    width = 8,height=4,res=300,units='in')


df_2 %>% ggplot(aes(x=ratio_of_perfect_couples, y=risk_diff_adjusted,
            color=task))+  geom_point(alpha=0.7)+
  facet_wrap(~task,scales = 'free')+
  scale_color_brewer(palette="Dark2")+
  theme_bw()+
  guides(color=guide_legend(title='Project name'))+
  ggtitle('Scatterplot of reg coef of MS and FD versus % of perfect coupling') +
  xlab('Percent of perfectly coupled mutants')+
  ylab(TeX('$\\hat{\\beta}_{MS}$'))+
  theme(plot.title = element_text(hjust=0.5,size=10))

dev.off()



png(filename = paste0(plot_ouput_dir,"couple_OR.png"),
    width = 8,height=4,res=300,units='in')

df_2 %>% filter(glm_converged) %>% ggplot(aes(x=ratio_of_perfect_couples, y=odd_ratio_log,
                    color=task))+  
  geom_point(alpha=0.7)+
  facet_wrap(~task,scales = 'free')+
  scale_color_brewer(palette="Dark2")+
  theme_bw()+
  guides(color=guide_legend(title='Project name'))+
  ggtitle('Scatterplot of OR versus % of perfect coupling') +
  xlab('Percent of perfectly coupled mutants')+
  ylab(TeX('$\\hat{OR}_{MS}$'))+
  theme(plot.title = element_text(hjust=0.5,size=10))
dev.off()


png(filename = paste0(plot_ouput_dir,"couple_partial_cor.png"),
    width = 8,height=4,res=300,units='in')

df_2 %>% ggplot(aes(x=ratio_of_perfect_couples, y=partial_cor_MS_FD,
                                              color=task))+  geom_point(alpha=0.7)+
  facet_wrap(~task,scales = 'free')+
  scale_color_brewer(palette="Dark2")+
  theme_bw()+
  guides(color=guide_legend(title='Project name'))+
  ggtitle('Scatterplot of partial correlation of MS and FD versus % of perfect coupling') +
  xlab('Percent of perfectly coupled mutants')+
  ylab('Partial correlation(MS,FD)')+
  theme(plot.title = element_text(hjust=0.5,size=10))

dev.off()



png(file = paste0(plot_ouput_dir,'perfect_coupling_prop_box.png'),
    width = 6,height=4,res=300,units='in')

couple_decouple_plot%>%
ggplot( aes(x=task, y=ratio_of_perfect_couples)) +
  geom_boxplot()+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  theme_bw()+
  ggtitle('Proportion of perfectly coupled mutants on the entire test suite')+
  ylab('Percent of perfectly coupled mutants')+
  xlab('Project name')+
  theme(plot.title = element_text(hjust=0.5,size=10))

dev.off()

png(file = paste0(plot_ouput_dir,'distri_coupling_prop.png'),
    width = 8,height=5,res=300,units='in')

coupling_prob_plot_gather%>%
  ggplot( aes( x=cond_prob,fill = task)) +
  geom_histogram(aes(y=..density..))+
  #geom_density(aes(y=..density..))+
  theme_bw()+
  ggtitle('Distribution of coupling probabilities for all mutants on the entire test suite')+
  ylab('Coupling probabilities')+
  facet_wrap(~task,scales = 'free')+
  scale_fill_brewer(palette="Dark2")+
  guides(fill=guide_legend(title='Project name'))+
  theme(plot.title = element_text(hjust=0.5,size=10))

dev.off()

# sample coupling pro
#  largest adjusted: 
# or adjusted Chart-24f.csv    Closure-101f.csv Lang-13f.csv     Math-65f.csv     Time-21f.csv   



