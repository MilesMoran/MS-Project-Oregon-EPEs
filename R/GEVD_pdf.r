################################################################################

library(ggplot2)
library(dplyr)
library(tidyr)
library(SpatialExtremes)

# set WD to the location of this file
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

################################################################################

y <- 
    data.frame(X = seq(-4,4,0.01)) %>% 
    mutate(
        sn0.5 = dgev(X, 0, 1, shape=-0.50),
        s0.5  = dgev(X, 0, 1, shape=+0.50),
        s0  = dgev(X, 0, 1, shape=+0.00)
    ) %>%  
    pivot_longer(cols=-c("X"), names_to="Shape", values_to="Density") 
    

Shape.labs <- c("ξ = 0", "ξ = +0.5", "ξ = -0.5")
names(Shape.labs) <- c("s0", "s0.5", "sn0.5")

ggplot(y, aes(x=X, y=Density, color=Shape)) + 
    geom_line(show.legend=FALSE, lwd=2) + 
    xlim(-3,3) + ylim(0,0.5) + 
    facet_wrap(~Shape, labeller=labeller(Shape=Shape.labs)) + 
    theme(strip.text.x = element_text(size = 20))

ggsave("../Plots/GEVD_pdfs.png", width=1200, height=300, dpi=100, units="px")

################################################################################