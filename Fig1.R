# Yaniv Brandvain 
# Altruistic load
# Plotting invasion criteria  [Figure 1]
# March 15th 2022

library(tidyverse)
library(viridis)
crossing(h = seq(0.0005,.1,length = 200),
         s = seq(0.005 ,1,length = 200))%>% 
  mutate(crit.self = 2*h *(3 + h *s) / (1-h*s)) %>%
  ggplot(aes(x = s, y = h, fill=crit.self))+
  geom_raster()+
  scale_fill_viridis_c()+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  labs(x = expression(Selection~coefficient~(s)),
       y = expression(Dominance~coefficient~(h)),
       title = "Selfing rate required for invasion")+
  theme(        axis.title =  element_text(size=12))

