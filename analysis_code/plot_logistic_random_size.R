

library(tidyverse)
library(ggplot2)
library(lmtest)

input_dir = '~/Desktop/cse599_f_2019/correlation_analysis_ms_fd_size/analysis_input/'
plot_ouput_dir = '~/Desktop/cse599_f_2019/correlation_analysis_ms_fd_size/plot_output/random_size/'

random_size_plot_df <- read.csv(file = paste0(input_dir,'random_size_new_result.csv'))

png(filename = paste0(plot_ouput_dir,"plot_cor_un_norm_ms_size.png"),
    width = 7,height=4,res=300,units='in')

ggplot(random_size_plot_df, aes(x=subject_name, y=cor_MS_FD)) + 
  geom_boxplot()+ geom_jitter(shape=16, position=position_jitter(0.2))+
 ggtitle('Pearson Correlation between unnormalized MS and FD with random sizes') +
  xlab('Project Name')+
  ylab('Pearson Correlation')+theme_bw()+
  theme(plot.title = element_text(hjust=0.5))

dev.off()

png(filename = "~/Desktop/cse599_f_2019/plot_final_week/plot_cor_ms_size.png",
    width = 7,height=4,res=300,units='in')
ggplot(random_size_plot_df, aes(x=subject_name, y=cor_ms_size)) + 
  geom_boxplot()+ geom_jitter(shape=16, position=position_jitter(0.2))+
  theme_bw()
dev.off()


### correlation between size and ms
png(filename = "~/Desktop/cse599_f_2019/plot_final_week/plot_cor_ms_size.png",
    width = 7,height=4,res=300,units='in')
ggplot(random_size_plot_df, aes(x=subject_name, y=cor_ms_size)) + 
  geom_boxplot()+ geom_jitter(shape=16, position=position_jitter(0.2))+
  theme_bw()
dev.off()
### log p value
# 
alpha_cutoff = log(0.05)

png(filename = "~/Desktop/cse599_f_2019/plot_final_week/plot_sign_diff.png",
    width = 7,height=4,res=300,units='in')

ggplot(random_size_plot_df, aes(x=subject_name, y=log(sig_diff_pval_lm))) + 
  geom_boxplot()+ geom_jitter(shape=16, position=position_jitter(0.2))+
  geom_hline(aes(yintercept=alpha_cutoff), colour="#990000")+theme_bw()

dev.off()
### some multiple testing correction
mean(random_size_plot_df$sig_diff_pval_lm <= 0.05)

mean(random_size_plot_df$lrt_fit_logit <= 0.05)


mean(random_size_plot_df$lrt_fit_lm <= 0.05)

mean(p.adjust(random_size_plot_df$lrt_fit_lm,method = 'BH')<= 0.05)

mean(p.adjust(random_size_plot_df$lrt_fit_logit,method = 'BH')<= 0.05)

mean(p.adjust(random_size_plot_df$lrt_fit_lm,method = 'bonferroni')<= 0.05)


plot_risk_diff <- random_size_plot_df %>% select(subject_name,
      risk_diff_adjusted,risk_diff_unadjusted) %>%
gather('risk_diff_adjusted','risk_diff_unadjusted',
                                key='adjust',value ='RD')
  
plot_risk_diff$adjust[plot_risk_diff$adjust=='risk_diff_adjusted'] <- 'adjusted'
plot_risk_diff$adjust[plot_risk_diff$adjust=='risk_diff_unadjusted'] <-'unadjusted'

png(filename = "~/Desktop/cse599_f_2019/plot_final_week/plot_rd_diff.png",
    width = 7,height=4,res=300,units='in')

ggplot(plot_risk_diff, aes(x=adjust, y=RD,
              fill=adjust)) + 
  geom_boxplot()+ geom_jitter(shape=16, position=position_jitter(0.2))+
  scale_fill_brewer(palette="Dark2")+facet_wrap(~subject_name)+theme_bw()
dev.off()


png(filename = "~/Desktop/cse599_f_2019/plot_final_week/plot_rd_diff_scatter.png",
    width = 7,height=4,res=300,units='in')
 random_size_plot_df %>% select(subject_name,
         risk_diff_adjusted,risk_diff_unadjusted) %>%
  ggplot( aes(x=risk_diff_adjusted, y=risk_diff_unadjusted,
                             color=subject_name)) +
  geom_point()+ scale_color_brewer(palette="Dark2")+
  facet_wrap(~subject_name,scales = 'free')+theme_bw()+
 geom_abline(slope=1,intercept = 0) 
dev.off()


############ OR

plot_odds_ratio <- random_size_plot_df %>% select(subject_name,
                          odd_ratio_log,odds_ratio_unadjust) %>%
  filter(odd_ratio_log<100) %>%
  gather('odd_ratio_log','odds_ratio_unadjust',
         key='adjust',value ='OR')  

plot_odds_ratio$adjust[plot_odds_ratio$adjust=='odd_ratio_log'] <- 'adjusted'
plot_odds_ratio$adjust[plot_odds_ratio$adjust=='odds_ratio_unadjust'] <-'unadjusted'

png(filename = "~/Desktop/cse599_f_2019/plot_final_week/plot_or_diff_scatter.png",
    width = 7,height=4,res=300,units='in')
 random_size_plot_df %>% select(subject_name,
         odd_ratio_log,odds_ratio_unadjust) %>%
  filter(odd_ratio_log<100) %>%
  ggplot( aes(x=odd_ratio_log, y=odds_ratio_unadjust,
              color=subject_name)) +
  geom_point()+ scale_color_brewer(palette="Dark2")+
  facet_wrap(~subject_name,scales = "free")+theme_bw()+
  geom_abline(slope=1,intercept = 0) 
dev.off()
#############

plot_ci_odds <- random_size_plot_df %>% select(odd_ratio_log,
                    odd_ratio_lcb,odd_ratio_ucb,odds_ratio_unadjust) %>%
                    filter(odd_ratio_log<5)%>%
                mutate(contains_unadjusted= (odds_ratio_unadjust<odd_ratio_ucb)&(odd_ratio_lcb<odds_ratio_unadjust))
mean(plot_ci_odds$contains_unadjusted)

plot_ci_odds$csv<-rownames(plot_ci_odds)


png(filename = "~/Desktop/cse599_f_2019/plot_final_week/plot_ci_adjust.png",
    width = 7,height=4,res=300,units='in')

ggplot(plot_ci_odds, aes(x=csv, y=odd_ratio_log, color = contains_unadjusted)) + 
  geom_errorbar(aes(ymin=odd_ratio_lcb, ymax=odd_ratio_ucb), width=.1) +
  geom_point() +
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  scale_color_manual(values=c('grey','red'))+
  labs(x = 'task/bug', y = 'adjusted OR',color = 'Contains unadjusted  OR')+
  theme(plot.title = element_text(hjust = 0.5,size=10))

dev.off()

##############

volcano_df <- random_size_plot_df%>% 
  select(subject_name,odd_ratio_log,odd_ratio_p_val) %>%
  filter(odd_ratio_log<100) %>%
  mutate(minus_log_p_val = -log(odd_ratio_p_val))

png(filename = "~/Desktop/cse599_f_2019/plot_final_week/volcano_stratify.png",
    width = 7,height=4,res=300,units='in')
volcano_df %>%filter(odd_ratio_log<5) %>%
  ggplot( aes(x=odd_ratio_log, y=minus_log_p_val,
              color=subject_name)) +
  geom_point()+ scale_color_brewer(palette="Dark2")+
  facet_wrap(~subject_name,scales = "free")+theme_bw()+
  geom_abline(slope=0,intercept =-log(0.05) ) 
dev.off()

##############




