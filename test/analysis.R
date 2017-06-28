library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)

setwd("/opt/utils/adaptSpotEnergies/test/")
##--------
##
data <- as.data.table(read.table(file = "directions_check.dat", header = FALSE, sep = ",",
                      colClasses = "numeric",
                      col.names  = c("tid",
                      "x", "y", "z",
                      "v.x", "v.y", "v.z", "cond")))
data <- data %>% mutate(cond = ifelse(cond == 1, "end", "start"))

data_long <- data %>% gather(type, val, x:v.z, factor_key = TRUE) %>% unite(cond, cond, type)
data_wide <- data_long %>% spread(cond, val) %>% select(-end_v.x, -end_v.y, -end_v.z)

data_wide <- data_wide %>% mutate(r.x = end_x - start_x, r.y = end_y - start_y, r.z = end_z - start_z)
data_wide <- data_wide %>% mutate(r = sqrt(r.x^2 + r.y^2 + r.z^2), v = sqrt(start_v.x^2 + start_v.y^2 + start_v.z^2))
data_wide <- data_wide %>% mutate(cos = (r.x*start_v.x + r.y*start_v.y + r.z*start_v.z)/(r*v))
data_wide <- data_wide %>% mutate(beamid = ifelse(tid > 102, 1, 0))
data_wide <- data_wide %>% mutate(spotid = ifelse(tid > 102, tid - 102, tid))
data_wide$beamid <- as.factor(data_wide$beamid)

ggplot(data = data_wide, aes(x = tid, y = cos)) +
    geom_line()

ggplot(data = data_wide, aes(x = spotid, y = r, color = beamid)) +
    geom_line()

data_wide2 <- data_wide %>% filter(tid >= 38 & tid <= 82)

##--------
##
data <- as.data.table(read.table(file = "energy_shifts.dat", header = FALSE, sep = ",",
                                 colClasses = "numeric",
                                 col.names  = c("energy")))
data <- data %>% mutate(energy = energy / 1000000)
ggplot(data, aes(energy)) +
    geom_histogram(binwidth = 0.5)
