# Testing script

#### Define parameters for data pull ##########################################
NEONsites <- c("KING")
startdate <- "2019-01"
enddate <- "2019-12"

#### Pull raw data and K from NEON API ########################################
#data <- request_NEON(NEONsites, startdate, enddate)

# Save raw data
#saveRDS(data, file = "data/Raw_KING_2019-01_2019-12.rds")
# Load
data <- readRDS(file = "data/Raw_KING_2019-01_2019-12.rds")

#### Visualize raw data #######################################################
library(ggplot2)
# Discharge v. K600
ggplot(data = data$k600_clean, aes(x = meanQ, y = k600)) +
  geom_smooth(method = "lm", color = "black") +
  ggpubr::stat_regline_equation(label.y = max(data$k600_clean$k600)) +
  ggpubr::stat_cor(label.y = max(data$k600_clean$k600)*0.9) +
  facet_wrap(vars(siteID)) +
  geom_point(shape = 21, fill = "white", size = 2) +
  # Theme adjustments
  theme_classic() +
  theme(axis.text = element_text(color = "black"),
        panel.background = element_rect(color = "black"),
        legend.position = "bottom")
# Datetime v. O2
ggplot(data = data$data, aes(x = solarTime, y = DO_mgL,
                             color = horizontalPosition)) +
  facet_grid(rows = vars(horizontalPosition)) +
  geom_line(alpha = 0.8) +
  geom_point(shape = 21, fill = "white", size = 0.2) +
  # Theme adjustments
  scale_x_datetime(date_breaks = "month", date_labels = "%b %Y") +
  scale_color_brewer(palette = 6, type = "qual") +
  theme_classic() +
  theme(axis.text = element_text(color = "black"),
        panel.background = element_rect(color = "black"),
        legend.position = "bottom")
# Time period for two station modeling
ggplot(data = data$data, aes(x = solarTime, y = modelingStrategy,
                           color = modelingStrategy)) +
  #facet_grid(rows = vars(horizontalPosition)) +
  geom_line(alpha = 0.8, color = "grey") +
  geom_point(shape = 21, fill = "white") +
  # Theme adjustments
  scale_x_datetime(date_breaks = "month", date_labels = "%b %Y") +
  scale_color_brewer(palette = 6, type = "qual") +
  theme_classic() +
  theme(axis.text = element_text(color = "black"),
        panel.background = element_rect(color = "black"),
        legend.position = "bottom")

#### Fit 2 station metabolism model to raw data ###############################


