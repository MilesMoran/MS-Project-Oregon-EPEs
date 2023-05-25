## =============================================================================
## load libraries and functions ================================================
## =============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(magrittr)
library(lubridate)

library(RColorBrewer)

library(climextRemes)
library(SpatialExtremes)

# set WD to the location of this file
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## =============================================================================
## import and clean data =======================================================
## =============================================================================

# (1) import contents of csv file
# (2) convert precipitation to mm
# (3) remove observations that are extremely large (>300mm)
# (4) flag observations that are suspiciously large (>200mm and not removed)
# (5) create temporal blocks as in Risser (2019)

startYear <- 1979

dat <- 
    readRDS("../Data/rain-daily.rds") %>% 
    rename(Precip = Precipitation) %>% 
    mutate(
        Precip = Precip/10,
        Precip = ifelse(Precip > 300, NA, Precip),
        Errant = ifelse(Precip > 200 & !is.na(Precip), 1, 0),
        Blocks = case_when(
            (Months %in% c(12, 1, 2)) ~ 1,  # DJF (Winter)
            (Months %in% c( 3, 4, 5)) ~ 2,  # MAM (Spring)
            (Months %in% c( 6, 7, 8)) ~ 3,  # JJA (Summer)
            (Months %in% c( 9,10,11)) ~ 4   # SON (Autumn)
        ),
        BlocksId = Blocks + 4*(Years-startYear) + 4*(Months%%12==0),
        YMD = ymd(paste(Years, Months, Days, sep="-"))
    ) 

## =============================================================================
## summarize the data for each temporal block ==================================
## =============================================================================

# (1) for eack block, get the block-maximum
# (2) for each block, identify how many observations were NA or Errant
# (3) only keep blocks in which the # of OK observations was >= 60
# (4) re-add the "Blocks" variable which denotes season
# (5) re-add lat/lon data from original dataset

stations <- 
    dat %>% 
    select(Id, lat, long) %>% 
    distinct()

dat.blockSummary <- 
    dat %>% 
    group_by(Id, BlocksId) %>% 
    summarize(
        maxPrecip = max(Precip, na.rm=TRUE),
        BlockSizeTotal = length(Precip),
        BlockSizeNA = sum(is.na(Precip)),
        BlockSizeErrant = sum(Errant),
        BlockSizeOK = (BlockSizeTotal - BlockSizeNA - BlockSizeErrant),
        .groups="drop"
    ) %>% 
    filter(BlockSizeOK >= 60) %>% 
    mutate(Blocks = ((BlocksId-1) %% 4)+1) %>% 
    filter(Blocks == 1) %>%             ### this is block-specific
    mutate(BlocksId = (BlocksId-1)/4)   ### this is block-specific

## =============================================================================
## Estimate GEV Parameters for Each Station (Assuming μ Varies With Time) ======
## =============================================================================

# also give 20 year return periods

dat.MLEs <- 
    dat.blockSummary %>% 
    filter(Id != "USC00352493") %>%  ### why does this one not converge?
    use_series(Id) %>% 
    unique() %>% 
    purrr::map(function(group) {
        filter(dat.blockSummary, Id==group) %>% 
            { fit_gev(y=.$maxPrecip, x=data.frame(t=unique(.$BlocksId)),
                      locationFun=~t, getParams=T, returnPeriod=20) 
            } %>%  
            { c(list("Id"=group),.) }
    }) %>% 
    purrr::map_df(unlist) %>% 
    select(-c("nllh", "info.convergence", "info.counts.function", 
              "info.failure", "info.counts.gradient",)) %>% 
    mutate_at(vars(-c("Id")), as.numeric) %>% 
    left_join(stations, by="Id")

## =============================================================================
## Plot GEVs Densities for Select Locations (Based on Estimated Params) ========
## =============================================================================

mm_to_inch <- function(x) { round(x/25.4, 2) }

example.stations <- 
    c("USW00024229",  # PDX
      "USW00024130",  # Baker City (Far North-West OR)
      "USW00024225",  # Medford 
      "USW00024230")  # Redmond (close to Bend)

x <- seq(0, 100, 0.01)
y <- data.frame(X=x)

for(stationId in example.stations) {
    ### Get PDF for Parameter Estimates in 1980
    IdYearCombo <- paste(stationId, "1980", sep="-")
    y[,IdYearCombo] <-
        filter(dat.MLEs, Id==stationId) %>% 
        { dgev(x, .$mle.mu0+1*.$mle.mu1, .$mle.scale, .$mle.shape) }

    ### Get PDF for Parameter Estimates in 2010
    IdYearCombo <- paste(stationId, "2010", sep="-")
    y[,IdYearCombo] <-
        filter(dat.MLEs, Id==stationId) %>% 
        { dgev(x, .$mle.mu0+40*.$mle.mu1, .$mle.scale, .$mle.shape) }
}

y %>% 
    set_colnames(
        c("X", "(1) Portland Airport-1980", "(1) Portland Airport-2010",
          "(2) Baker City-1980", "(2) Baker City-2010",
          "(3) Medford-1980", "(3) Medford-2010",
          "(4) Redmond-1980", "(4) Redmond-2010")
    ) %>% 
    pivot_longer(
        cols = -c("X"),
        names_to = c("Id", "Year"),
        names_pattern = "(.+)-(.+)",
        values_to = "Density"
    ) %>% 
    ggplot(aes(x=X, y=Density, color=Year)) + geom_line() + 
        scale_x_continuous(limits=c(0,100)) + 
        scale_y_continuous(limits=c(0,0.14)) + 
        labs(x="Precipitation (mm)") + 
        theme(axis.line=element_line(), legend.position="bottom") + 
        facet_wrap(~Id)

ggsave("../Plots/WinterMaxPcp-PDFs.png", width=1200, height=800, dpi=120, units="px")

## =============================================================================
## Map GEV Parameter Estimates =================================================
## =============================================================================

library(latex2exp)

borders <- 
    map_data("state") %>% 
    filter(region == "oregon") %>% 
    mutate(long=-long)
    
map.oregon <- 
    ggplot(mapping=aes(x=-long, y=lat)) + 
    geom_polygon(aes(group=group), data=borders, fill="white", color="black") + 
    coord_fixed(1.3) + 
    theme_void()

plot.mu0 <- 
    map.oregon + 
    geom_point(aes(color=mle.mu0), data=dat.MLEs) + 
    scale_color_viridis_c() + 
    labs(color=TeX("$\\hat{\\mu}_0$"))

plot.mu1 <- 
    map.oregon + 
    geom_point(aes(color=mle.mu1), data=dat.MLEs) + 
    scale_color_distiller(palette="RdBu", limits=c(-0.65, +0.65), oob = scales::squish) + 
    labs(color=TeX("$\\hat{\\mu}_1$"))

plot.sigma <- 
    map.oregon + 
    geom_point(aes(color=mle.scale), data=dat.MLEs) + 
    scale_color_viridis_c() + 
    labs(color=TeX("$\\hat{\\sigma}$"))

plot.xi <-
    map.oregon + 
    geom_point(aes(color=mle.shape), data=dat.MLEs) + 
    scale_color_distiller(palette="RdBu", limits=c(-0.65, +0.65)) + 
    labs(color=TeX("$\\hat{\\xi}$"))

p <- gridExtra::grid.arrange(plot.mu0, plot.sigma, plot.mu1, plot.xi, ncol=2, nrow=2)

ggsave("../Plots/WinterMaxPcp-Estimates-Nonstationary.png", plot=p,
       width=2000, height=1545, dpi=300, units="px")

## =============================================================================
## Map Estimated Return Values =================================================
## =============================================================================

# these are the return-values for the parameters as they stand in t ≈ 2010
# (since μ_t = μ_0 + [t × μ_1] )

plot.returnvals <- 
    map.oregon + 
    geom_point(aes(color=returnValue20), data=dat.MLEs) + 
    scale_color_viridis_c() + 
    labs(color=TeX("$\\hat{r}_{20}$"))

ggsave("../Plots/WinterMaxPcp-20YrReturn-2000-Nonstationary.png", plot=plot.returnvals,
       width=1200, height=800, dpi=300, units="px", bg="white")

## =============================================================================
## =============================================================================
