# Yaniv Brandvain 
# Altruistic load
# Plotting results from SLiMulations
# March 15th 2022

library(tidyverse)
library(janitor)
library(conflicted)
library(patchwork)
conflict_prefer_all("dplyr",quiet = TRUE)


### Figure 3
all_gen <- read_csv("allGenerations.csv") %>% 
  clean_names() %>%
  mutate(selfing_rate = factor(str_remove(as.character(self_r),"0")),
         mutation_rate = factor(mut_rate, labels = c(bquote(mu==5~x~10^-8), bquote(mu==5~x~10^-7), bquote(mu==5~x~10^-6))),
         inbreeding_depression = factor(in_d, labels = c(bquote(delta==0.5), bquote(delta==0.75), bquote(delta==1))),
         selection_coeff = factor(selec_r, labels = c(bquote(s==1), bquote(s==0.8), bquote(s==0.6), bquote(s==0.4), bquote(s==0.2)))) %>%
  mutate(mutation_rate = fct_rev(mutation_rate),
         selection_coeff = fct_rev(selection_coeff))



plotForDeltab <- function(delta){
  all_gen %>%
    filter( mutation_rate %in% c("mu == 5 ~ x ~ 10^-8","mu == 5 ~ x ~ 10^-6"),in_d==delta,selec_r%in% c(-.2,-.6,-1))%>%
    group_by(replicate, pair_type, mut_rate, selec_r , in_d, self_r) %>%
    filter( generation == max(generation )) %>%
    ungroup()%>%
    mutate( mutation_rate  = fct_rev( mutation_rate ))%>%
    ggplot(aes(x = selfing_rate, y =1000 * avg_mut_freq, color = pair_type, lty = pair_type, shape = pair_type))+
    geom_point()+
    facet_grid( mutation_rate ~  selection_coeff,labeller = "label_parsed")+
    stat_summary(fun = "mean",geom = "line", aes(group = pair_type))+ 
    theme(legend.position = "bottom" )+
    scale_y_continuous(breaks =   c(0,2,4,6),limits = c(0,7))+
    scale_linetype_manual(values = c(2,1))+
    guides(alpha = "none")
}

delta1 <- plotForDeltab(1) 
delta1 <-delta1 +labs(title = expression(Inbreeding~depression:~delta==1),
                      y = expression(Mean~mutation~freq~x~10^3),
                      x = expression(Selfing~rate~(sigma)))

delta.75 <- plotForDeltab(.75)
delta.75 <-delta.75 +labs(title = expression(Inbreeding~depression:~delta==.75),
                          y = expression(Mean~mutation~freq~x~10^3),
                          x = expression(Selfing~rate~(sigma)))

z   <- delta1  +  delta.75 & theme(legend.position = "right") 
z + plot_layout(guides = "collect") + 
  plot_annotation(tag_levels = 'A')

ggsave("SLiM_mut_freq.png",width = 8.5, height = 3.5)



### Figure S1



plotForDeltabSupp <- function(delta){
  all_gen %>%
    filter( mutation_rate %in% c("mu == 5 ~ x ~ 10^-8","mu == 5 ~ x ~ 10^-6"),in_d==delta,selec_r%in% c(-.2,-.6,-1))%>%
    group_by(replicate, pair_type, mut_rate, selec_r , in_d, self_r) %>%
    filter( generation == max(generation )) %>%
    ungroup()%>%
    mutate( mutation_rate  = fct_rev( mutation_rate ))%>%
    ggplot(aes(x = selfing_rate, y =avg_muts_per_ind, color = pair_type, lty = pair_type, shape = pair_type))+
    geom_point()+
    facet_grid( mutation_rate ~  selection_coeff,labeller = "label_parsed", scales = "free_y")+
    stat_summary(fun = "mean",geom = "line", aes(group = pair_type))+ 
    theme(legend.position = "bottom" )+
    #scale_y_continuous(breaks =   c(0,2,4,6),limits = c(0,7))+
    scale_linetype_manual(values = c(2,1))+
    guides(alpha = "none")
}

delta1Supp <- plotForDeltabSupp(1) 
delta1Supp <-delta1Supp +labs(title = expression(Inbreeding~depression:~delta==1),
                      y = "# mutations per individual",
                      x = expression(Selfing~rate~(sigma)))

delta.75Supp <- plotForDeltabSupp(.75)
delta.75Supp <-delta.75Supp +labs(title = expression(Inbreeding~depression:~delta==.75),
                                  y = "# mutations per individual",
                                  x = expression(Selfing~rate~(sigma)))

z   <- delta1Supp  +  delta.75Supp & theme(legend.position = "right") 
z + plot_layout(guides = "collect") + 
  plot_annotation(tag_levels = 'A')

ggsave("SLiM_mut_count.png",width = 8.5, height = 3.5)



