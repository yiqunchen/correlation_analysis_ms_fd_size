
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
  select(subject_name,risk_diff_unadjusted,list_of_files,fixed_ratio)
join_0.05 <- fixed_size_0.05_new_result %>% select(subject_name,
                     risk_diff_unadjusted,list_of_files,fixed_ratio)
join_0.075 <- fixed_size_0.075_new_result %>% select(subject_name,
                      risk_diff_unadjusted,list_of_files,fixed_ratio)
join_0.1 <- fixed_size_0.1_new_result %>% select(subject_name,
                                                 risk_diff_unadjusted,list_of_files,fixed_ratio)
join_0.125 <- fixed_size_0.125_new_result %>% select(subject_name,
                                                     risk_diff_unadjusted,list_of_files,fixed_ratio)
join_0.15 <- fixed_size_0.15_new_result %>% select(subject_name,
                                                   risk_diff_unadjusted,list_of_files,fixed_ratio)
join_0.175 <- fixed_size_0.175_new_result %>% select(subject_name,
                                                     risk_diff_unadjusted,list_of_files,fixed_ratio)
join_0.2 <- fixed_size_0.2_new_result %>% select(subject_name,
            risk_diff_unadjusted,list_of_files,fixed_ratio)

join_fixed_risk <- bind_rows(join_0.025,join_0.05,join_0.075,join_0.1,
                              join_0.125,join_0.15,join_0.175,join_0.2) %>% 
  rename(risk_diff = risk_diff_unadjusted) %>%
  mutate(fixed_ratio = as.character(fixed_ratio))


join_random_adj <- random_size_plot_df %>% 
  select(subject_name,list_of_files, risk_diff_adjusted) %>% 
  mutate(fixed_ratio='random_adjusted') %>%
  rename(risk_diff = risk_diff_adjusted)


join_random_unadj <- random_size_plot_df %>% 
  select(subject_name,list_of_files, risk_diff_unadjusted) %>%
  mutate(fixed_ratio='random_unadjusted') %>%
  rename(risk_diff = risk_diff_unadjusted)

plot_risk_diff <- bind_rows(join_fixed_risk,join_random_adj,join_random_unadj) 

plot_random_fixed <- plot_risk_diff$fixed_ratio
plot_random_fixed[!is.na(as.numeric(plot_random_fixed))] <- 'fixed'
  
plot_risk_diff <- plot_risk_diff %>% mutate(plot_random_fixed=plot_random_fixed)


png(filename = paste0(plot_ouput_dir,"beta_coef_size.png"),
    width = 7,height=4,res=300,units='in')

plot_risk_diff %>%
  mutate(plot_random_fixed = fct_recode(plot_random_fixed,
              'FD ~ MS, fixed size'= 'fixed', 'FD ~ MS, random size' = 'random_unadjusted',
              'FD ~ MS + Size, random size' = 'random_adjusted')) %>%
ggplot( aes(x=subject_name, y=risk_diff,
                           fill = plot_random_fixed)) + 
  geom_boxplot(alpha=0.6)+
  ggtitle('Difference between regression coefficients for different regression models') +
  xlab('Project Name')+
  ylab(TeX('$|\\hat{\\beta}_{MS}|$'))+
  theme_bw()+scale_fill_brewer(palette="Dark2")+
  scale_color_brewer(palette="Dark2")+
  theme(plot.title = element_text(hjust=0.5,size=10))+
  guides(fill=guide_legend(title='Regression models'))
  

dev.off()


###### boxplot replication by each category


png(filename = paste0(plot_ouput_dir,"stratified_beta_coef_size.png"),
    width = 7,height=4,res=300,units='in')

plot_risk_diff%>% 
  filter(fixed_ratio%in%c(0.025,0.075,0.15,0.2,'random_unadjusted','random_adjusted'))%>%
  mutate(fixed_ratio = fct_recode(fixed_ratio,
                            'FD ~ MS, 2.5% sampling'=  '0.025' ,    
                            'FD ~ MS, 7.5% sampling'=  '0.075' , 
                            'FD ~ MS, 15% sampling'=  '0.15' , 
                            'FD ~ MS, 20% sampling'=  '0.2' , 
                   'FD ~ MS, random size' = 'random_unadjusted',
              'FD ~ MS + Size, random size' = 'random_adjusted')) %>% 
  ggplot( aes(x=fixed_ratio, y=risk_diff,
              fill = fixed_ratio,colour=fixed_ratio)) + 
  geom_boxplot(alpha=0.6)+
  ggtitle('Regression coefficient of MS for different models, stratified by projects') +
  #xlab('Project Name')+
  ylab(TeX('$|\\hat{\\beta}_{MS}|$'))+theme_bw()+
  facet_wrap(subject_name~.,nrow = 2,scale='free')+
  scale_fill_brewer(palette='Dark2')+
  scale_color_brewer(palette='Dark2')+
  theme(plot.title = element_text(hjust=0.5,size=10),axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

dev.off()



################ look at highest 

max_coef <- plot_risk_diff %>% filter(!is.na(risk_diff),plot_random_fixed=='random_adjusted') %>% 
  group_by(subject_name)%>% summarize(max(risk_diff))

plot_risk_diff$list_of_files[plot_risk_diff$risk_diff%in%max_coef$`max(risk_diff)`]
 # largest adjusted: Chart-25f.csv   Closure-28f.csv Lang-13f.csv    Math-70f.csv    Time-24f.csv   

