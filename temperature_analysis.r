#############################################################################################################################
############################################## 4- EXPERIMENTAL WARMING ANALYSES #############################################
#############################################################################################################################

data <- read.csv(file=".../temperature_data.csv", header = TRUE)

data$plot = factor (data$plot)
data$date = factor (data$date)
data$country = factor(data$country)
data$habitat = factor(data$habitat)
data$treatment = factor(data$treatment)
data$site = factor(data$site)
data$year = factor(data$year)
data$month = factor(data$month)
data$day = factor(data$day)

data[sapply(data, is.infinite)] <- NA # replace inf and -inf values by NA

summary(data)


library(lme4)
library(lmerTest)
library(pbkrtest)

mod_mean = lmer(mean ~ habitat*treatment + (1|year/month/day) + (1|country/site), data=data)
mod_min = lmer(min ~ habitat*treatment + (1|year/month/day) + (1|country/site), data=data)
mod_max = lmer(max ~ habitat*treatment + (1|year/month/day) + (1|country/site), data=data)

summary(mod_mean)
summary(mod_min) #non significative interaction --> remove from the model
summary(mod_max)

mod_min2 = lmer(min ~ habitat + treatment + (1|year/month/day) + (1|country/site), data=data)
summary(mod_min2)


# See differences
library (effects)
plot (effect ("habitat:treatment", mod_mean))
effect ("habitat:treatment", mod_mean)

plot (effect ("habitat:treatment", mod_max))
effect ("habitat:treatment", mod_max)

plot(effect("habitat", mod_min2))
effect("habitat", mod_min2)

plot(effect("treatment", mod_min2))
effect("treatment", mod_min2)


############ Analysing forests and hedgerows separately
data <- read.csv(file="C:/Users/pablo/Desktop/submitt!!/figshare/final docs/4 - temperature analyses/temperature_data.csv", header = TRUE)

data$plot = factor (data$plot)
data$date = factor (data$date)
data$country = factor(data$country)
data$habitat = factor(data$habitat)
data$treatment = factor(data$treatment)
data$site = factor(data$site)
data$year = factor(data$year)
data$month = factor(data$month)
data$day = factor(data$day)

data[sapply(data, is.infinite)] <- NA # replace inf and -inf values by NA

forest = data[data$habitat == "F",]
hedgerow = data[data$habitat == "H",]

## Run models
library(lme4)
library(lmerTest)
library(pbkrtest)
library (effects)

# 1 - Forests
mod_min_f = lmer(min ~ treatment + (1|year/month/day) + (1|country/site), data=forest)
mod_mean_f = lmer(mean ~ treatment + (1|year/month/day) + (1|country/site), data=forest)
mod_max_f = lmer(max ~ treatment + (1|year/month/day) + (1|country/site), data=forest)

summary(mod_min_f)
summary(mod_mean_f)
summary(mod_max_f)

# See differences
plot (effect ("treatment", mod_min_f))
effect ("treatment", mod_min_f)
mean(forest$min[forest$treatment == "C"], na.rm = TRUE)
mean(forest$min[forest$treatment == "T"], na.rm = TRUE)

plot (effect ("treatment", mod_mean_f))
effect ("treatment", mod_mean_f)
mean(forest$mean[forest$treatment == "C"], na.rm = TRUE)
mean(forest$mean[forest$treatment == "T"], na.rm = TRUE)

plot (effect ("treatment", mod_max_f))
effect ("treatment", mod_max_f)
mean(forest$max[forest$treatment == "C"], na.rm = TRUE)
mean(forest$max[forest$treatment == "T"], na.rm = TRUE)

# 2 - Hedgerows
mod_min_h = lmer(min ~ treatment + (1|year/month/day) + (1|country/site), data=hedgerow)
mod_mean_h = lmer(mean ~ treatment + (1|year/month/day) + (1|country/site), data=hedgerow)
mod_max_h = lmer(max ~ treatment + (1|year/month/day) + (1|country/site), data=hedgerow)

summary(mod_min_h)
summary(mod_mean_h) # no significative differences
summary(mod_max_h)

# See differences
plot (effect ("treatment", mod_min_h))
effect ("treatment", mod_min_h)
mean(hedgerow$min[hedgerow$treatment == "C"], na.rm = TRUE)
mean(hedgerow$min[hedgerow$treatment == "T"], na.rm = TRUE)

plot (effect ("treatment", mod_mean_h))
effect ("treatment", mod_mean_h)
mean(hedgerow$mean[hedgerow$treatment == "C"], na.rm = TRUE)
mean(hedgerow$mean[hedgerow$treatment == "T"], na.rm = TRUE)
mean(hedgerow$mean, na.rm= TRUE) # as no significative differences found

plot (effect ("treatment", mod_max_h))
effect ("treatment", mod_max_h)
mean(hedgerow$max[hedgerow$treatment == "C"], na.rm = TRUE)
mean(hedgerow$max[hedgerow$treatment == "T"], na.rm = TRUE)

## Average heating

av_heat_f = data.frame(
  habitat = "forest",
  min = (mean(forest$min[forest$treatment == "T"], na.rm = TRUE) - mean(forest$min[forest$treatment == "C"], na.rm = TRUE)),
  mean = (mean(forest$mean[forest$treatment == "T"], na.rm = TRUE)- mean(forest$mean[forest$treatment == "C"], na.rm = TRUE)),
  max = (mean(forest$max[forest$treatment == "T"], na.rm = TRUE)- mean(forest$max[forest$treatment == "C"], na.rm = TRUE))
)

av_heat_h = data.frame(
  habitat = "hedgerow",
  min = (mean(hedgerow$min[hedgerow$treatment == "T"], na.rm = TRUE) - mean(hedgerow$min[hedgerow$treatment == "C"], na.rm = TRUE)),
  mean = (mean(hedgerow$mean[hedgerow$treatment == "T"], na.rm = TRUE)- mean(hedgerow$mean[hedgerow$treatment == "C"], na.rm = TRUE)),
  max = (mean(hedgerow$max[hedgerow$treatment == "T"], na.rm = TRUE)- mean(hedgerow$max[hedgerow$treatment == "C"], na.rm = TRUE))
)

av_heat = rbind(av_heat_f, av_heat_h) ## Average heating


