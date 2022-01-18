# Testing script
library (ggplot2)

#### Define parameters for data pull ##########################################
NEONsites <- c("KING")
startdate <- "2019-01"
enddate <- "2019-12"

#### Pull raw data and K from NEON API ########################################
data <- request_NEON(NEONsites, startdate, enddate)

# Save raw data
#saveRDS(data, file = "data/Raw_KING_2019-01_2019-12.rds")
# Load

#### Visualize raw data #######################################################
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

#### Fit 2 station metabolism model to raw data ###############################
results <- fit_twostation(data, nbatch = 1e4)

ggplot(data = results$results, aes(x = date)) +
  geom_ribbon(aes(ymin = GPP.lower, ymax= GPP.upper, fill = "GPP"), alpha = 0.5) +
  geom_line(aes(y = GPP, color = "GPP"))+
  geom_point(aes(y = GPP, color = "GPP"), shape = 21, fill = "white") +
  # Theme adjustments
  scale_x_date(date_breaks = "month", date_labels = "%b %Y") +
  scale_color_brewer(palette = 6, type = "qual") +
  theme_classic() +
  theme(axis.text = element_text(color = "black"),
        panel.background = element_rect(color = "black"),
        legend.position = "bottom")

ggplot(data = results$results, aes(x = date)) +
  geom_ribbon(aes(ymin = K.lower, ymax= K.upper, fill = "K"), alpha = 0.5) +
  geom_line(aes(y = K, color = "K"))+
  geom_point(aes(y = K, color = "K"), shape = 21, fill = "white") +
  # Theme adjustments
  scale_x_date(date_breaks = "month", date_labels = "%b %Y") +
  scale_color_brewer(palette = 6, type = "qual") +
  theme_classic() +
  theme(axis.text = element_text(color = "black"),
        panel.background = element_rect(color = "black"),
        legend.position = "bottom")
