# Yaniv Brandvain 
# Altruistic load
# Finding equilibrium allele frequencies by numerical iteration [Figure 2]
# March 15th 2022

library(tidyverse)
library(viridis)


#Calculate mean fitness
findW <- function(self, h, s, paa, pAa, pAA, p, replaceaa, replaceAa, replaceAA){
 (1/2)*(self - 1)*(2*paa*(replaceaa + 1)*(h*p*s - 1) + pAa*(replaceAa + 1)*(s*(h + p) - 2) - 2*pAA*(replaceAA + 1)*((h - 1)*p*s - h*s + 1))
}

# Find replacement probabilities in each maternal family
replaces <- function(p,h,s,self){
  c(aa = (1 - self) *(h*s)* p ,
    Aa = (s/2) * (self* (1/ 2 - p) +  h + p),
    AA = s *(self + (1 - self) * (p + h - h * p)) )
}

# Find genotype frequencies every generation
newFreqs <- function(geno_freqs, h, s, self){
   paa       <- geno_freqs["aa"]  %>% as.numeric()
   pAa       <- geno_freqs["Aa"]  %>% as.numeric()
   pAA       <- geno_freqs["AA"]  %>% as.numeric()
   p         <- sum(geno_freqs * c(0,.5,1))
  this_replace <- replaces(p = p,  h = h, s = s, self = self)
   replaceaa <- this_replace["aa"]  %>% as.numeric()
   replaceAa <- this_replace["Aa"]  %>% as.numeric()
   replaceAA <- this_replace["AA"]  %>% as.numeric()
  WBAR <- findW(self = self, h = h, s = s, 
        paa = paa, pAa = pAa , pAA = pAA, p = p, 
        replaceaa = replaceaa, replaceAa = replaceAa, replaceAA = replaceAA) %>%
    as.numeric()
  new_freqs <- c(aa =  ((1 - p)*(1 - self)*((1 - p) + paa*replaceaa + pAa*(replaceAa/2)))/WBAR,
                 Aa = ((1 - h*s)*(p*paa*replaceaa*(1 - self) + 
                          p*paa*(1 - self) + (1 - p)*pAA*replaceAA*(1 - self) + 
                          (1 - p)*pAA*(1 - self) + (1/2)*pAa*
                          replaceAa*(1 - self) + (1/2)*pAa*(1 - self)))/WBAR,
                 AA = (1 - s)*((pAa*(1 - self)*p*(1/2) + pAa*replaceAa*(1 - self)*p*(1/2) + 
                       pAA*(1 - self)*p + pAA*replaceAA*(1 - self)*p)/WBAR) )
    return(new_freqs)
}

# Find equilbrium allele byt iteration of recursion equations
findEq <- function(initial_aa, initial_Aa, initial_AA, self, h, s){
  geno_freqs <- c(aa = initial_aa , Aa = initial_Aa, AA = initial_AA)
  g   <- 0 
  diff_genos <- 0
  while(g < 100 | diff_genos > 1e-13){
    g <- g+1
    new_freqs  <- newFreqs(geno_freqs = geno_freqs, h = h, s = s, self = self)
    diff_genos <- sum(abs(new_freqs - geno_freqs))
    geno_freqs <- new_freqs  / sum(new_freqs)
  }
  return(c(g,geno_freqs))
}

# Make a table of parameter combinations to loop over
params <- crossing(s = seq(0.04,1,length = 25),
                   selfing_rate = c(.01,.02,.04,.08,.16,.32),
                   h = seq(0,.075,length = 16))

# Looping over paramter combinations
# CAUTION.This takes a long time
eqfreqs <- apply(params, ,1, function(X){
        tmp <- findEq(initial_aa = .9999, initial_Aa = .0001, initial_AA = 0, 
               self = X[["selfing_rate"]], h = X[["h"]], s = X[["s"]])
        data.frame(self = X[["selfing_rate"]],
                   h = X[["h"]],  gen = tmp[1], s = X[["s"]],
                   paa = tmp[["aa"]],  pAa = tmp[["Aa"]],  pAA = tmp[["AA"]], p   = tmp[["Aa"]]/2 + tmp[["AA"]] )
      }) %>% bind_rows()

# Saving the output
#write_csv(x = eqfreqs, file = "eq_freqs.csv")


read_csv("eq_freqs.csv") %>%
  rename(`selfing rate` = self) %>%
  ggplot(aes(x = s, y = p, color = h))+
  geom_point(size = .5)+
  facet_wrap(~`selfing rate`, labeller = "label_both")+
  scale_color_viridis_c(option = "inferno")+
  scale_y_continuous(limits = c(0,.07))+
  labs(x = "Selection coefficient (s)",
       y = "Equilibrium allele frequency (p*)",
       color = "Dominance\ncoefficient (h)  ")+
  scale_x_continuous(breaks = seq(0,1,.25), labels = c("0","1/4","1/2","3/4","1")) + 
  theme(#legend.position = c(0.85, 0.82), legend.direction = "horizontal",
        legend.position = "bottom",
        strip.text = element_text(size=11),
        axis.title =  element_text(size=12)
        #legend.background  = element_rect(fill = "transparent", colour = "transparent") 
        )
+
  guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5))


