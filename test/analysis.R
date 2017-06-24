library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)

setwd("/opt/utils/adaptSpotEnergies/test/")
data <- as.data.table(read.table(file = "data.dat", header = FALSE, sep = " ",
                   colClasses = "numeric",
                   col.names  = c("start.x", "start.y", "start.z",
                                  "start2.x", "start2.y", "start2.z",
                                  "end.x", "end.y", "end.z", "wepl",
                                  "dCos.x", "dCos.y", "dCos.z",
                                  "energy", "beamid", "spotid")))
data <- data %>% mutate(spotid = ifelse(beamid == 1, spotid - 102, spotid))
data <- data %>% mutate(r.x = end.x - start2.x) %>% mutate(r.y = end.y - start2.y) %>% mutate(r.z = end.z - start2.z)
data <- data %>% mutate(r = sqrt(r.x^2 + r.y^2 + r.z^2))
data <- data %>% mutate(dCos = sqrt(dCos.x^2 + dCos.y^2 + dCos.z^2))
data <- data %>% mutate(trace = sqrt((start.x - start2.x)^2 + (start.y - start2.y)^2 + (start.z - start2.z)^2))
data <- data %>% mutate(ok = (r.x*dCos.x + r.y*dCos.y + r.z*dCos.z)/(dCos*r))

data$meta.x <- as.factor(data$meta.x)

data <- data %>% mutate(dCos = sqrt(dCos.x^2 + dCos.y^2 + dCos.z^2))

hist(data$dCos, breaks = 10)
hist(data$trace, breaks = 10)
plot(data$trace)

ggplot(data = data, aes(x = meta.y, y = trace, color = meta.x)) +
    geom_line()

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
