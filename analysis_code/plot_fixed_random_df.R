

library(RColorBrewer)
library(tidyverse)
library(ggplot2)
library(lmtest)


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
      select(subject_name,cor_ms_fd,list_of_files,fixed_ratio)
join_0.05 <- fixed_size_0.05_new_result %>% select(subject_name,
                                                   cor_ms_fd,list_of_files,fixed_ratio)
join_0.075 <- fixed_size_0.075_new_result %>% select(subject_name,
                                                     cor_ms_fd,list_of_files,fixed_ratio)
join_0.1 <- fixed_size_0.1_new_result %>% select(subject_name,
                                                 cor_ms_fd,list_of_files,fixed_ratio)
join_0.125 <- fixed_size_0.125_new_result %>% select(subject_name,
                                                     cor_ms_fd,list_of_files,fixed_ratio)
join_0.15 <- fixed_size_0.15_new_result %>% select(subject_name,
                                                   cor_ms_fd,list_of_files,fixed_ratio)
join_0.175 <- fixed_size_0.175_new_result %>% select(subject_name,
                                                     cor_ms_fd,list_of_files,fixed_ratio)
join_0.2 <- fixed_size_0.2_new_result %>% select(subject_name,
                                                 cor_ms_fd,list_of_files,fixed_ratio)

join_random <- random_size_plot_df %>% 
  select(subject_name,list_of_files, cor_MS_FD) %>%
  mutate(fixed_ratio=NA)
colnames(join_random) <- c("subject_name", "list_of_files" ,"cor_ms_fd",'fixed_ratio')


fixed_random_cor <- bind_rows(join_0.025,join_0.05,join_0.075,join_0.1,
          join_0.125,join_0.15,join_0.175,join_0.2,join_random)

fixed_random_cor[is.na(fixed_random_cor$fixed_ratio),'fixed_ratio'] <- 'random'


fixed_random_cor$random <- (fixed_random_cor$fixed_ratio=='random')



###### boxplot replication
png(filename = paste0(plot_ouput_dir,"reproduce_fig_3.png"),
    width = 6,height=4,res=300,units='in')

fixed_random_cor %>% 
  mutate(random_plot = if_else(random,'Random','Fixed')) %>%
  ggplot( aes(x=subject_name, y=cor_ms_fd,fill = random_plot)) + 
  geom_boxplot(alpha=0.6)+
 # geom_point(position = position_jitterdodge())+
  ggtitle('Pearson Cor between mutation score and real fault detection') +
  xlab('Project Name')+
  ylab('Pearson Correlation')+theme_bw()+scale_fill_brewer(palette="Dark2")+
  scale_color_brewer(palette="Dark2")+
  theme(plot.title = element_text(hjust=0.5,size=12))+
  guides(fill=guide_legend(title='Size distribution'))

dev.off()
  

###### boxplot replication by each category

png(filename = paste0(plot_ouput_dir,"reproduce_fig_3_stratify_size.png"),
    width = 8,height=4,res=300,units='in')

fixed_random_cor%>% filter(fixed_ratio%in%c(0.025,0.075,0.125,0.175,0.2,'random'))%>%
ggplot( aes(x=fixed_ratio, y=cor_ms_fd,
                             fill = fixed_ratio)) + 
  geom_boxplot(alpha=0.6)+
  ggtitle('Pearson Cor between mutation score and real fault detection for varying sampling ratio') +
  xlab('Sampling ratio')+
  ylab('Pearson Correlation')+theme_bw()+
  facet_wrap(subject_name~.,nrow = 2)+
  scale_fill_brewer(palette='Dark2')+
  scale_color_brewer(palette='Dark2')+
  theme(plot.title = element_text(hjust=0.5,size=10))+
  guides(fill=guide_legend(title='Size distribution'))

dev.off()


#########violin plot

fixed_random_cor%>% filter(fixed_ratio%in%c(0.025,0.075,0.125,0.175,0.2,'random'))%>%
  ggplot( aes(x=fixed_ratio, y=cor_ms_fd,
              fill = fixed_ratio,colour=fixed_ratio)) + 
  geom_violin(alpha=0.6)+
  ggtitle('Pearson Correlation between unnormalized MS and FD with random sizes') +
  xlab('Project Name')+
  ylab('Pearson Correlation')+theme_bw()+
  facet_wrap(subject_name~.,nrow = 2)+
  scale_fill_brewer(palette='Dark2')+
  scale_color_brewer(palette='Dark2')+
  theme(plot.title = element_text(hjust=0.5))


random_size_plot_df %>% select(subject_name,
                               odd_ratio_log,odds_ratio_unadjust) %>%
  filter(odd_ratio_log<100) %>%
  ggplot( aes(x=odd_ratio_log, y=odds_ratio_unadjust,
              color=subject_name)) +
  geom_point()+ scale_color_brewer(palette="Dark2")+
  facet_wrap(~subject_name,scales = "free")+theme_bw()+
  geom_abline(slope=1,intercept = 0) 




########## Now let's look at a joint box plot for partial correlation and cor

png(filename = paste0(plot_ouput_dir,"partial_cor_plot.png"),
    width = 8,height=4,res=300,units='in')


random_size_plot_df %>% 
  select(subject_name,list_of_files, cor_MS_FD,partial_cor_MS_FD,cor_size_FD,
         partial_cor_size_FD)  %>% 
  gather('cor_MS_FD','partial_cor_MS_FD','cor_size_FD','partial_cor_size_FD',
          key='type',value ='correlation') %>% 
  mutate(type = factor(type, levels = c('cor_MS_FD','partial_cor_MS_FD',
           'cor_size_FD','partial_cor_size_FD'))) %>%
  mutate( type = fct_recode(type, "Cor(MS,FD)" = "cor_MS_FD"  
    , "PCor(MS,FD)" = "partial_cor_MS_FD" ,
    "Cor(Size,FD)"="cor_size_FD", "PCor(Size,FD)" = "partial_cor_size_FD"))%>%
  ggplot( aes(x=type, y=correlation,
             fill = type))+  geom_boxplot(alpha=0.6)+
  scale_fill_brewer(palette="Dark2")+
  facet_wrap(~subject_name)+theme_bw()+
  theme(plot.title = element_text(hjust=0.5,size=10),
        axis.text.x = element_text(angle = 60,hjust=1))+
  guides(fill=guide_legend(title='Type of correlation'))+
  ggtitle('(Partial) Pearson Cor for mutation score, test suite size, and real fault detection with random sizes') +
  xlab('Type of correlation')+
  ylab('Correlation value')

dev.off()

########## now try partial cor versus actual cor

png(filename = paste0(plot_ouput_dir,"partial_cor_ms_fault_scatter.png"),
    width = 7,height=4,res=300,units='in')

random_size_plot_df %>% 
  select(subject_name,list_of_files, cor_MS_FD,partial_cor_MS_FD,cor_size_FD,
         partial_cor_size_FD)  %>%
  ggplot( aes(x=partial_cor_MS_FD, y=cor_MS_FD,
              color=subject_name))+  geom_point(alpha=0.6)+
  scale_color_brewer(palette="Dark2")+
  facet_wrap(~subject_name)+theme_bw()+
  geom_abline(slope=1,intercept = 0)+
  guides(color=guide_legend(title='Project name'))+
  ggtitle('Change of Cor(MS,FD) after adjusting for size for random size distributions') +
  xlab('Partial correlation(MS,FD)')+
  ylab('Cor(MS,FD)')+
  theme(plot.title = element_text(hjust=0.5,size=10))


dev.off()

########

join_random_unadj <- random_size_plot_df %>% 
  select(subject_name,list_of_files, cor_MS_FD) %>%
  mutate(fixed_ratio='random_size')%>%
  rename(cor_ms_fd=cor_MS_FD)

join_random_partial_cor <- random_size_plot_df %>% 
  select(subject_name,list_of_files, partial_cor_MS_FD) %>%
  mutate(fixed_ratio='random_size_partial') %>%
  rename(cor_ms_fd=partial_cor_MS_FD)


fixed_random_cor_adj <- bind_rows(join_0.025,join_0.05,join_0.075,join_0.1,
           join_0.125,join_0.15,join_0.175,join_0.2) %>% 
  mutate(fixed_ratio= as.character(fixed_ratio))

fixed_random_cor_adj <- bind_rows(fixed_random_cor_adj,join_random_unadj,join_random_partial_cor)


fixed_random_cor_adj %>% filter(fixed_ratio%in%c(0.025,0.075,
         0.125,0.175,0.2,'random_size','random_size_partial'))%>%
  ggplot( aes(x=fixed_ratio, y=cor_ms_fd,
              fill = fixed_ratio)) +  
  geom_boxplot()+
  scale_fill_brewer(palette="Dark2")+
  facet_wrap(~subject_name)+theme_bw()


fixed_random_cor_adj %>% group_by(fixed_ratio) %>% 
  summarise(test_n = sum(is.na(cor_ms_fd)))


remove_na_files <- as.character(unique(fixed_random_cor_adj$list_of_files[is.na(fixed_random_cor_adj$cor_ms_fd)]))



###### let's plot this as well 


png(filename = paste0(plot_ouput_dir,"partial_cor_random_fixed.png"),
    width = 8,height=4,res=300,units='in')

fixed_random_cor_adj %>% filter((!list_of_files%in%remove_na_files),fixed_ratio%in%c(0.025,0.075,
      0.125,0.175,0.2,'random_size','random_size_partial'))%>%
  mutate(fixed_ratio =fct_recode(fixed_ratio,
          "Random"= "random_size","Random partial cor"='random_size_partial') )%>%
  ggplot( aes(x=fixed_ratio, y=cor_ms_fd,
              fill = fixed_ratio)) +  
  geom_boxplot(alpha=0.7)+
  scale_fill_brewer(palette="Dark2")+
  facet_wrap(~subject_name)+theme_bw()+
  guides(fill=guide_legend(title='Types of correlation'))+
  ggtitle('Distribution Cor(MS,FD) for different sizes and types of correlation') +
  xlab('Types of correlation')+
  ylab('Correlation value')+
  theme(plot.title = element_text(hjust=0.5,size=10),
        axis.text.x = element_text(angle = 60,hjust=1,size=5))

dev.off()


########## partial correlation
png(filename = paste0(plot_ouput_dir,"partial_cor_SIZE_fault_scatter.png"),
    width = 7,height=4,res=300,units='in')

random_size_plot_df %>% 
  select(subject_name,list_of_files, cor_MS_FD,partial_cor_MS_FD,cor_size_FD,
         partial_cor_size_FD)  %>%
  ggplot( aes(x=partial_cor_size_FD, y=cor_size_FD,
              color=subject_name))+  geom_point(alpha=0.7)+
  scale_color_brewer(palette="Dark2")+
  facet_wrap(~subject_name)+theme_bw()+
  geom_abline(slope=1,intercept = 0)+
  guides(color=guide_legend(title='Project name'))+
  ggtitle('Change of Cor(Size,FD) after adjusting for MS for random size distributions') +
  xlab('Partial correlation(Size,FD)')+
  ylab('Cor(Size,FD)')+
  theme(plot.title = element_text(hjust=0.5,size=10))


dev.off()



