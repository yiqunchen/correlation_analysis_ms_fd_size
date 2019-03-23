
library(RColorBrewer)
library(tidyverse)
library(ggplot2)
library(lmtest)
library(latex2exp)


input_dir = '~/Desktop/cse599_f_2019/correlation_analysis_ms_fd_size/analysis_input/'
plot_ouput_dir = '~/Desktop/cse599_f_2019/correlation_analysis_ms_fd_size/plot_output/'

random_size_plot_df <- read.csv(file = paste0(input_dir,'random_size_new_result.csv'))

fixed_size_0.025_new_result <- read.csv(file = paste0(input_dir,'fixed_size_0.025_new_result.csv'))
fixed_size_0.05_new_result <- read.csv(file = paste0(input_dir,'fixed_size_0.05_new_result.csv'))
fixed_size_0.075_new_result <- read.csv(file = paste0(input_dir,'fixed_size_0.075_new_result.csv'))
fixed_size_0.1_new_result <- read.csv(file = paste0(input_dir,'fixed_size_0.1_new_result.csv'))
fixed_size_0.125_new_result <- read.csv(file = paste0(input_dir,'fixed_size_0.125_new_result.csv'))
fixed_size_0.15_new_result <- read.csv(file = paste0(input_dir,'fixed_size_0.15_new_result.csv'))
fixed_size_0.175_new_result <- read.csv(file = paste0(input_dir,'fixed_size_0.175_new_result.csv'))
fixed_size_0.2_new_result <- read.csv(file = paste0(input_dir,'fixed_size_0.2_new_result.csv'))


join_0.025 <- fixed_size_0.025_new_result %>% 
  select(subject_name,odd_ratio_log,list_of_files,fixed_ratio,glm_converge)
join_0.05 <- fixed_size_0.05_new_result %>% select(subject_name,
                                                   odd_ratio_log,list_of_files,fixed_ratio,glm_converge)
join_0.075 <- fixed_size_0.075_new_result %>% select(subject_name,
                                                     odd_ratio_log,list_of_files,fixed_ratio,glm_converge)
join_0.1 <- fixed_size_0.1_new_result %>% select(subject_name,
                                                 glm_converge,odd_ratio_log,list_of_files,fixed_ratio)
join_0.125 <- fixed_size_0.125_new_result %>% select(subject_name,
                  glm_converge,odd_ratio_log,list_of_files,fixed_ratio)
join_0.15 <- fixed_size_0.15_new_result %>% select(subject_name,
                                                   odd_ratio_log,list_of_files,fixed_ratio,glm_converge)
join_0.175 <- fixed_size_0.175_new_result %>% select(subject_name,
                                                     odd_ratio_log,list_of_files,fixed_ratio,glm_converge)
join_0.2 <- fixed_size_0.2_new_result %>% select(subject_name,
                                                 odd_ratio_log,list_of_files,fixed_ratio,glm_converge)

join_fixed_or <- bind_rows(join_0.025,join_0.05,join_0.075,join_0.1,
                             join_0.125,join_0.15,join_0.175,join_0.2) %>% 
  rename(odds_ratio = odd_ratio_log) %>%
  mutate(fixed_ratio = as.character(fixed_ratio))


join_random_adj <- random_size_plot_df %>% 
  select(subject_name,list_of_files, glm_converged,odd_ratio_log) %>% 
  mutate(fixed_ratio='random_adjusted') %>%
  rename(odds_ratio = odd_ratio_log,glm_converge = glm_converged)


join_random_unadj <- random_size_plot_df %>% 
  select(subject_name,list_of_files,glm_converged, odds_ratio_unadjust) %>%
  mutate(fixed_ratio='random_unadjusted') %>%
  rename(odds_ratio = odds_ratio_unadjust,glm_converge = glm_converged)

plot_or <- bind_rows(join_fixed_or,join_random_adj,join_random_unadj) 

plot_random_fixed <- plot_or$fixed_ratio
plot_random_fixed[!is.na(as.numeric(plot_random_fixed))] <- 'fixed'

plot_or <- plot_or %>% mutate(plot_random_fixed=plot_random_fixed)


png(filename = paste0(plot_ouput_dir,"beta_OR_size.png"),
    width = 8.5,height=4,res=300,units='in')


plot_or%>% filter(glm_converge) %>%
  mutate(plot_random_fixed = fct_recode(plot_random_fixed,
     'logit(FD) ~ MS, fixed size'= 'fixed', 'logit(FD) ~ MS, random size' = 'random_unadjusted',
           'logit(FD) ~ MS + Size, random size' = 'random_adjusted')) %>%
  ggplot( aes(x=subject_name, y=odds_ratio,
        fill = plot_random_fixed)) + 
  geom_boxplot(alpha=0.6)+
  ggtitle('Odds ratio of Mutation score and real fault detection with varying sizes') +
  #xlab('Project Name')+
  ylab(TeX('$\\hat{OR}_{MS}$'))+theme_bw()+scale_fill_brewer(palette="Dark2")+
  scale_color_brewer(palette="Dark2")+
  facet_wrap(subject_name~.,nrow = 2,scale='free')+
  theme(plot.title = element_text(hjust=0.5,size=10))+
  guides(fill=guide_legend(title='Regression models'))


dev.off()


###### boxplot replication by each category
png(filename = paste0(plot_ouput_dir,"stratified_OR_coef_size.png"),
    width = 8.5,height=4,res=300,units='in')


plot_or%>% filter(glm_converge) %>%
  filter(fixed_ratio%in%c(0.025,0.075,0.15,0.2,'random_unadjusted','random_adjusted'))%>%
  mutate(fixed_ratio = fct_recode(fixed_ratio,
                                  'logit(FD) ~ MS, 2.5% sampling'=  '0.025' ,    
                                  'logit(FD) ~ MS, 7.5% sampling'=  '0.075' , 
                                  'logit(FD) ~ MS, 15% sampling'=  '0.15' , 
                                  'logit(FD) ~ MS, 20% sampling'=  '0.2' , 
                                  'logit(FD) ~ MS, random size' = 'random_unadjusted',
                                  'logit(FD) ~ MS + Size, random size' = 'random_adjusted')) %>%
  ggplot( aes(x=fixed_ratio, y=odds_ratio,
              fill = fixed_ratio,colour=fixed_ratio)) + 
  geom_boxplot(alpha=0.6)+
  ggtitle('Estimated Odds ratio for MS across different models, stratified by projects') +
  #xlab('Project Name')+
  ylab(TeX('$\\hat{OR}_{MS}$'))+theme_bw()+
  facet_wrap(subject_name~.,nrow = 2,scale='free')+
  scale_fill_brewer(palette='Dark2')+
  scale_color_brewer(palette='Dark2')+
  theme(plot.title = element_text(hjust=0.5,size=10),axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


dev.off()


################ 
max_or <- plot_or %>% filter(!is.na(odds_ratio),plot_random_fixed=='random_adjusted') %>% 
  group_by(subject_name)%>% summarize(max(odds_ratio))

plot_or$list_of_files[plot_or$odds_ratio%in%max_or$`max(odds_ratio)`]

# Chart-24f.csv    Closure-101f.csv Lang-13f.csv     Math-65f.csv     Time-21f.csv   



