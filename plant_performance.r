#############################################################################################################################
################################################ VITAL RATE REGRESSIONS #####################################################
#############################################################################################################################

library(lme4)
library(lmerTest)
library(pbkrtest)
library(readr)
data <- read_csv(".../data.csv")

data$habitat = factor(data$habitat)
data$treatment = factor(data$treatment)
data$transplant <- factor(data$transplant)
data$site <- factor(data$site)
data$survival = factor (data$survival)
data$rep_prob = factor (data$rep_prob)
data$seed_prod_2019 = factor (data$seed_prod_2019)
summary(data)

## Plant height
height_model_glmm <- lmer(z ~ habitat*treatment + (1|transplant/site), data= data)
summary(height_model_glmm)
anova(height_model_glmm)

height_model_glmm2 <- lmer(z ~ habitat + treatment + (1|transplant/site), data= data)
summary(height_model_glmm2)
anova(height_model_glmm2)
# In this case not necessary but it is with poisson regressions
height_noT <- lmer(z ~ habitat + (1|transplant/site), data= data)
height_noH <- lmer(z ~ treatment + (1|transplant/site), data= data)
anova(height_model_glmm2, height_noT, test="LRT") # treatment p value
anova(height_model_glmm2, height_noH, test="LRT") # habitat p value

# See differences
library (effects)
plot (effect ("habitat", height_model_glmm2)) # Forest plants have higher height
plot (effect ("treatment", height_model_glmm2)) # Treatment plants have higher height


## Flowers number
flower_model_glmm <- glmer(rep_nr ~ habitat*treatment + (1|transplant/site), data= data, family= poisson)
summary(flower_model_glmm)
anova(flower_model_glmm)

flower_model_glmm2 <- glmer(rep_nr ~ habitat + treatment + (1|transplant/site), data= data, family= poisson)
summary(flower_model_glmm2)
anova(flower_model_glmm2) # No p values

flower_noT <- glmer(rep_nr ~ habitat + (1|transplant/site), data= data, family= poisson)
flower_noH <- glmer(rep_nr ~ treatment + (1|transplant/site), data= data, family= poisson)
anova(flower_model_glmm2, flower_noT, test="LRT") # treatment p value
anova(flower_model_glmm2, flower_noH, test="LRT") # habitat p value

# See differences
plot (effect ("habitat", flower_model_glmm2)) # Forest plants have higher flower nr
plot (effect ("treatment", flower_model_glmm2)) # Treatment plants have higher flower nr


## seed number
seed_model_glmm <- glmer(seed_nr_2019 ~ habitat*treatment + (1|transplant/site), data= data, family= poisson)
summary(seed_model_glmm)
anova(seed_model_glmm) #no p value

seed_noInt <- glmer(seed_nr_2019 ~ habitat + treatment + (1|transplant/site), data= data, family = poisson) # no interaction model
anova(seed_model_glmm, seed_noInt, test = "LRT")


# See differences
plot (effect ("habitat:treatment", seed_model_glmm)) # + to - number of seeds: FC > Hc > FT > HT

## Now, analyses of the probabilities (survival and reproduction probability)

## Survival
surv_mod1 <- glmer (survival ~ habitat*treatment + (1|transplant/site), data = data, family = binomial, na.action = na.omit)
summary(surv_mod1)#no significative differences
anova(surv_mod1)

surv_mod2 <- glmer (survival ~ habitat + treatment + (1|transplant/site), data = data, family = binomial, na.action = na.omit)
summary(surv_mod2)
anova(surv_mod2)
########################NO SIGNIFICATIVE DIFFERENCES AMONG GROUPS FOR SURVIVAL


## Reprorduction
rep_mod1 <- glmer (rep_prob ~ habitat*treatment + (1|transplant/site), data = data, family = binomial, na.action = na.omit)
summary(rep_mod1)#no significative differences
anova(rep_mod1)

rep_mod2 <- glmer (rep_prob ~ habitat + treatment + (1|transplant/site), data = data, family = binomial, na.action = na.omit)
summary(rep_mod2)
anova(rep_mod2)

########################NO SIGNIFICATIVE DIFFERENCES AMONG GROUPS FOR REPRODUCTION


### Calculate the p values of these models
## For models with interactions
anova(surv_mod1, surv_mod2, test = "LRT") #p = 0.8959 - survival
anova(rep_mod1, rep_mod2, test = "LRT") #p = 0.2964 - reproduction


## For models without interactions
#Survival
surv_mod_noT <- glmer (survival ~ habitat + 1 + (1|transplant/site), data = data, family = binomial, na.action = na.omit)
surv_mod_noH <- glmer (survival ~ 1 + treatment + (1|transplant/site), data = data, family = binomial, na.action = na.omit)

anova(surv_mod2, surv_mod_noT)# p = 0.1926 for treatment
anova(surv_mod2, surv_mod_noH)# p = 0.3927 for habitat

#Reproduction
rep_mod_noT <- glmer (rep_prob ~ habitat + 1 + (1|transplant/site), data = data, family = binomial, na.action = na.omit)
rep_mod_noH <- glmer (rep_prob ~ 1 + treatment + (1|transplant/site), data = data, family = binomial, na.action = na.omit)

anova(rep_mod2, rep_mod_noT)# p = 0.1991 for treatment
anova(rep_mod2, rep_mod_noH)# p = 0.6026 for habitat


### Average values
surv_1 = length(data$survival [data$survival == 1 & !is.na(data$survival)]) # Get all 1 from dataset, without na
surv_0 = length(data$survival [data$survival == 0 & !is.na(data$survival)]) # Get all 0 from dataset, without na
surv_rate = surv_1/(surv_0 + surv_1)
surv_rate #0.8244681

rep_1 = length(data$rep_prob [data$rep_prob == 1 & !is.na(data$rep_prob)]) # Get all 1 from dataset, without na
rep_0 = length(data$rep_prob [data$rep_prob == 0 & !is.na(data$rep_prob)]) # Get all 0 from dataset, without na
rep_rate = rep_1/(rep_0 + rep_1)
rep_rate #0.9251337
