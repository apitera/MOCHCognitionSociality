# Code to make Figures 2 & 3 in Pitera et al., in preparation: "Spatiotemporal variation in cognitive phenotype, social network position, and distribution of social associations in a food-caching bird"

# These figures show snow depth data (cm) from high (ca., 2400 m) and low elevation (ca., 1900 m) sites at Sagehen Experimental Forest over 5 winter seasons (referring to non-breeding seasons for the study species of interest: the mountain chickadee, Poecile gambeli)

# Libraries etc...------
source("https://gist.githubusercontent.com/benmarwick/2a1bb0133ff568cbe28d/raw/fb53bd97121f7f9ce947837ef1a4c65a73bffb3f/geom_flat_violin.R") 

library(cowplot)
library(ggplot2)
library(lubridate)

# Load in data ------
dat.path <- "./icloud/LabThings/MOCHData/AnimalSocialNetworks/Assort/Cleaned/"
weat.dat <- "AssortCog_Fig2_Fig3_Data.RData"

load(paste0(dat.path, weat.dat))


# ================ Fig. 2a ==================
#  Snow depth over time for all seasons + elevations: --

# ad <-
ggplot(data = all.weather, aes(x = Day.of.Winter, y = snow.depth, color = season, linetype = location, shape = location)) +
  geom_line(size = 0.5) +
  xlab("Day of season") +
  ylab("Snow depth (cm)") +
  scale_color_manual(values = c("#230038" , # 2015-16
                                "#4E53B3" , # 2016-17
                                "#C58FFF" , # 2017-18
                                "#8ADBD7" , # 2018-19
                                "#5BA696"), 
      labels = c("2015-16", "2016-17", "2017-18", "2018-19", "2019-20"), name = "")+
  scale_x_continuous(limits = c(0, 200), expand = c(0,0), breaks = seq(0, 200, by = 15)) +
  guides(linetype = "none")+
  theme_classic() +
  theme(plot.margin = margin(.25,.25,.25,.25,"cm"), 
      legend.position = c(.5,.99),
      legend.direction = "horizontal",
      legend.spacing.x = unit(.001, "cm"),
      legend.text = element_text(size = 12, margin = margin(r = 15)),
      legend.background = element_rect(fill = "transparent", color = NA),
      axis.ticks = element_line(color = "black"),
      axis.title.y = element_text(vjust = 2.5, size = 15),
      axis.text.y = element_text(size = 13, color = "black", margin = margin(r = 2)),
      axis.title.x = element_text(size = 15),
      axis.text.x = element_text(size = 13, color = "black", margin = margin(t = 2)),
      plot.background = element_rect(fill = "transparent", color = NA),
      panel.background = element_rect(fill = "transparent")) 


# ========================= Fig. 2b =====================================
# Plot snow depth within each season + elevation *for network data days ONLY*--------


new.labs <- c("   2015-16", "    2016-17", "2017-18     ", "2018-19", "   2019-20")
names(new.labs) <- c("2015_16", "2016_17", "2017_18", "2018_19", "2019_20")

# rc <-
ggplot(data = net.dates, aes(y = snow.depth, x = elevation, color = season, fill= paste(season, elevation), shape = elevation)) +
  ylab("Snow depth (cm)") +
  facet_grid(cols = vars(season), switch = 'x', labeller = labeller(season = new.labs)) +
  geom_flat_violin(position = position_nudge(x = .3, y = 0)) +
  geom_boxplot(position = position_nudge(x = 0.2, y = 0),
               width = .1, outlier.shape = NA) +
  geom_point(aes(y = snow.depth, size = elevation, shape = elevation),
             position = position_jitter(width = 0.1),
              alpha = 0.3) +
  scale_shape_manual(values = c(18, 5)) +
  scale_size_manual(values = c(2.5, 1.5)) +
  scale_fill_manual(values = alpha(c("#230038", "transparent",
                                     "#4E53B3", "transparent",
                                     "#C58FFF",
                                     "#8ADBD7",
                                     "#5BA696",  "transparent"), 0.5)) +
  scale_color_manual(values = c("#230038",
                                "#4E53B3",
                                "#C58FFF",
                                "#8ADBD7",
                                "#5BA696")) +
  theme_classic() +
  geom_hline(yintercept = 100, color = "snow3", linetype = 3) + 
  theme(plot.margin = margin(.25,.25,.25,.25,"cm"), 
        legend.position = "none",
        axis.title.y = element_text(vjust = 2.5, size = 15, color = "black"),
        axis.text.y = element_text(size = 13, color = "black"),
        axis.ticks.y = element_line(color = "black"),
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text.x = element_text(vjust = 2.5, size = 13, color = "black"),
        panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.spacing.x = unit(0, "null"),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent")) 

# fig_2 <- plot_grid(ad, rc, nrow = 2, ncol = 1, labels = "auto")
# ggsave(fig_2, filename = paste0(dat.path,"AssortCog_Fig2.tiff"), width = 174, height = 150, units = "mm")



# =====================Fig. 3 - Appendix ================ 
# Supplemental figure: snow depth over time for network days ONLY 
nd <-
ggplot(data = net.dates, aes(x = Day.of.Winter, y = snow.depth, color = season, linetype = elevation, fill = season, size = elevation, shape = elevation)) +
  geom_line(linewidth = 0.4, alpha = 0.7, aes(linetype = elevation), show.legend = F) +
  geom_point(aes(size = elevation), stroke = .5) +
  scale_shape_manual(values = c(23, 5))+
  scale_size_manual(values = c(1.05, 1)) +
  scale_fill_manual(values = alpha(c("#230038" , # 2015-16
                                     "#4E53B3" , # 2016-17
                                     "#C58FFF" , # 2017-18
                                     "#8ADBD7" , # 2018-19
                                     "#5BA696"), 
                                   0.5), 
                    labels = c("2015-16", "2016-17", "2017-18", "2018-19", "2019-20"), 
                    name = "")+ #2019-20)
  scale_color_manual(values = c("#230038", # 2015-16
                                "#4E53B3", # 2016-17
                                "#C58FFF",  #2017-18
                                "#8ADBD7",# 2018-19
                                "#5BA696"),
                     labels = c("2015-16", "2016-17", "2017-18", "2018-19", "2019-20"), 
                     name = "")+ 
  ylab("Snow depth (cm)") +
  xlab("Day of season") +
  geom_hline(yintercept = 100, linetype = 3, color = "snow3") + 
  scale_x_continuous(limits = c(0, 185), expand = c(0,0), breaks = seq(0, 185, by = 15)) +
  guides(shape = "none", size = "none", linetype = "none", 
         color = guide_legend(override.aes = list(shape = 23, size = 2))) +
  theme_classic() +
  theme(plot.margin = margin(.25,0,.25,.25,"cm"), 
        legend.position = c(.5, 1),
        legend.direction = "horizontal",
        legend.spacing.x = unit(.001, "cm"),
        legend.text = element_text(size = 12, margin = margin(r = 15)),
        axis.title.y = element_text(vjust = 2.5, size = 15),
        axis.text.y = element_text(size = 13, color = "black", margin = margin(r = 2)),
        axis.title.x = element_text(size = 15),
        axis.text.x = element_text(size = 13, color = "black", margin = margin(t = 2)),
        legend.background = element_rect(fill = "transparent", color = NA),
        axis.ticks = element_line(color = "black"),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent")) 

