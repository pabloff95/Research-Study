############################################################################################################################################
#################################################### 2 - MODEL SELECTION ###################################################################
############################################################################################################################################
rm(list=ls()) #clean global environment

library(lme4)

data <- read.csv(file= ".../data.csv")
data$site <- factor(data$site)
data$individual<- factor (data$individual)
data$survival<- factor (data$survival)
data$rep_prob<- factor (data$rep_prob)
data$seed_prod_2019<- factor (data$seed_prod_2019)
summary(data)


#create data subsets
forest_control = data[which(data$habitat == "F" & data$treatment == "C"),]
forest_warming = data[which(data$habitat == "F" & data$treatment == "T"),]
hedgerow_control = data[which(data$habitat == "H" & data$treatment == "C"),]
hedgerow_warming = data[which(data$habitat == "H" & data$treatment == "T"),]

#create dataframes for analysing seed_nr (needs dataframe with only individuals that produced seeds: seed_prod =1) 
#and rep_nr (needs dataframe with only individuals that reproducted: rep_prob =1)  for each factor combination 
df_rep_nr_fc = subset (forest_control, forest_control$rep_prob == 1)
df_rep_nr_fw = subset (forest_warming, forest_warming$rep_prob == 1)
df_rep_nr_hc = subset (hedgerow_control, hedgerow_control$rep_prob == 1)
df_rep_nr_hw = subset (hedgerow_warming, hedgerow_warming$rep_prob == 1)

df_seed_nr_fc = subset (forest_control, forest_control$seed_prod_2019 == 1)
df_seed_nr_fw = subset (forest_warming, forest_warming$seed_prod_2019 == 1)
df_seed_nr_hc = subset (hedgerow_control, hedgerow_control$seed_prod_2019 == 1)
df_seed_nr_hw = subset (hedgerow_warming, hedgerow_warming$seed_prod_2019 == 1)

#create data frames to introduce the AICs of the models
model_type = c("1", "z", "z + I(z^2)")
forest_control_aic  = data.frame (factor= "forests_control", model_type = model_type)
forest_warming_aic  = data.frame (factor= "forests_warming", model_type = model_type)
hedgerow_control_aic  = data.frame (factor= "hedgerow_control", model_type = model_type)
hedgerow_warming_aic  = data.frame (factor= "hedgerow_warming", model_type = model_type)


### 1- Survival

#habitat = forest;  treatment = control
surv_forest_control_1 = glm (survival ~ 1 , family= binomial, data = forest_control, na.action = na.omit)
surv_forest_control_2 = glm (survival ~ z , family= binomial, data = forest_control, na.action = na.omit)
surv_forest_control_3 = glm (survival ~ z + I(z^2) , family= binomial, data = forest_control, na.action = na.omit)

#habitat = forest;  treatment = warming
surv_forest_warming_1 = glm (survival ~ 1 , family= binomial, data = forest_warming, na.action = na.omit)
surv_forest_warming_2 = glm (survival ~ z , family= binomial, data = forest_warming, na.action = na.omit)
surv_forest_warming_3 = glm (survival ~ z + I(z^2) , family= binomial, data = forest_warming, na.action = na.omit)

#habitat = hedgerow;  treatment = control
surv_hedgerow_control_1 = glm (survival ~ 1 , family= binomial, data = hedgerow_control, na.action = na.omit)
surv_hedgerow_control_2 = glm (survival ~ z , family= binomial, data = hedgerow_control, na.action = na.omit)
surv_hedgerow_control_3 = glm (survival ~ z + I(z^2) , family= binomial, data = hedgerow_control, na.action = na.omit)

#habitat = hedgerow;  treatment = warming
surv_hedgerow_warming_1 = glm (survival ~ 1 , family= binomial, data = hedgerow_warming, na.action = na.omit)
surv_hedgerow_warming_2 = glm (survival ~ z , family= binomial, data = hedgerow_warming, na.action = na.omit)
surv_hedgerow_warming_3 = glm (survival ~ z + I(z^2) , family= binomial, data = hedgerow_warming, na.action = na.omit)

#store results (AIC)
forest_control_aic <- cbind(forest_control_aic, data.frame(survival = (c( AIC(surv_forest_control_1), AIC(surv_forest_control_2), AIC(surv_forest_control_3) ))))
forest_warming_aic <- cbind(forest_warming_aic, data.frame(survival = (c( AIC(surv_forest_warming_1), AIC(surv_forest_warming_2), AIC(surv_forest_warming_3) ))))
hedgerow_control_aic <- cbind(hedgerow_control_aic, data.frame(survival = (c( AIC(surv_hedgerow_control_1), AIC(surv_hedgerow_control_2), AIC(surv_hedgerow_control_3) ))))
hedgerow_warming_aic <- cbind(hedgerow_warming_aic, data.frame(survival = (c( AIC(surv_hedgerow_warming_1), AIC(surv_hedgerow_warming_2), AIC(surv_hedgerow_warming_3) ))))


### 2- Growth

#habitat = forest;  treatment = control
growth_forest_control_1 = lmer (z1 ~ 1 + (1|transplant/site),  data = forest_control, na.action = na.omit)
growth_forest_control_2 = lmer (z1 ~ z + (1|transplant/site),  data = forest_control, na.action = na.omit)
growth_forest_control_3 = lmer (z1 ~ z + I(z^2) + (1|transplant/site),  data = forest_control, na.action = na.omit)

#habitat = forest;  treatment = warming
growth_forest_warming_1 = lmer (z1 ~ 1 + (1|transplant/site),  data = forest_warming, na.action = na.omit)
growth_forest_warming_2 = lmer (z1 ~ z + (1|transplant/site),  data = forest_warming, na.action = na.omit)
growth_forest_warming_3 = lmer (z1 ~ z + I(z^2) + (1|transplant/site),  data = forest_warming, na.action = na.omit)

#habitat = hedgerow;  treatment = control
growth_hedgerow_control_1 = lmer (z1 ~ 1 + (1|transplant/site),  data = hedgerow_control, na.action = na.omit)
growth_hedgerow_control_2 = lmer (z1 ~ z + (1|transplant/site),  data = hedgerow_control, na.action = na.omit)
growth_hedgerow_control_3 = lmer (z1 ~ z + I(z^2) + (1|transplant/site),  data = hedgerow_control, na.action = na.omit)

#habitat = hedgerow;  treatment = warming
growth_hedgerow_warming_1 = lmer (z1 ~ 1 + (1|transplant/site),  data = hedgerow_warming, na.action = na.omit)
growth_hedgerow_warming_2 = lmer (z1 ~ z + (1|transplant/site),  data = hedgerow_warming, na.action = na.omit)
growth_hedgerow_warming_3 = lmer (z1 ~ z + I(z^2) + (1|transplant/site),  data = hedgerow_warming, na.action = na.omit)

#store results (AIC): heare the AIC has been calculated with anova() as the result from AIC was significantly different
forest_control_aic <- cbind(forest_control_aic, data.frame(growth = (c( AIC(growth_forest_control_1), AIC(growth_forest_control_2), AIC(growth_forest_control_3) ))))
forest_warming_aic <- cbind(forest_warming_aic, data.frame(growth = (c( AIC(growth_forest_warming_1), AIC(growth_forest_warming_2), AIC(growth_forest_warming_3) ))))
hedgerow_control_aic <- cbind(hedgerow_control_aic, data.frame(growth = (c( AIC(growth_hedgerow_control_1), AIC(growth_hedgerow_control_2), AIC(growth_hedgerow_control_3) ))))
hedgerow_warming_aic <- cbind(hedgerow_warming_aic, data.frame(growth = (c( AIC(growth_hedgerow_warming_1), AIC(growth_hedgerow_warming_2), AIC(growth_hedgerow_warming_3) ))))

### 3 - rep probability

#habitat = forest;  treatment = control
rep_prob_forest_control_1 = glm (rep_prob ~ 1 , family= binomial, data = forest_control, na.action = na.omit)
rep_prob_forest_control_2 = glm (rep_prob ~ z , family= binomial, data = forest_control, na.action = na.omit)
rep_prob_forest_control_3 = glm (rep_prob ~ z + I(z^2)  , family= binomial, data = forest_control, na.action = na.omit)

#habitat = forest;  treatment = warming
rep_prob_forest_warming_1 = glm (rep_prob ~ 1 , family= binomial, data = forest_warming, na.action = na.omit)
rep_prob_forest_warming_2 = glm (rep_prob ~ z , family= binomial, data = forest_warming, na.action = na.omit)
rep_prob_forest_warming_3 = glm (rep_prob ~ z + I(z^2) , family= binomial, data = forest_warming, na.action = na.omit)

#habitat = hedgerow;  treatment = control
rep_prob_hedgerow_control_1 = glm (rep_prob ~ 1 , family= binomial, data = hedgerow_control, na.action = na.omit)
rep_prob_hedgerow_control_2 = glm (rep_prob ~ z , family= binomial, data = hedgerow_control, na.action = na.omit)
rep_prob_hedgerow_control_3 = glm (rep_prob ~ z + I(z^2) , family= binomial, data = hedgerow_control, na.action = na.omit)

#habitat = hedgerow;  treatment = warming
rep_prob_hedgerow_warming_1 = glm (rep_prob ~ 1 , family= binomial, data = hedgerow_warming, na.action = na.omit)
rep_prob_hedgerow_warming_2 = glm (rep_prob ~ z , family= binomial, data = hedgerow_warming, na.action = na.omit)
rep_prob_hedgerow_warming_3 = glm (rep_prob ~ z + I(z^2) , family= binomial, data = hedgerow_warming, na.action = na.omit)

#store results (AIC)
forest_control_aic  <- cbind(forest_control_aic,  data.frame(rep_prob = (c( AIC(rep_prob_forest_control_1), AIC(rep_prob_forest_control_2), AIC(rep_prob_forest_control_3)))))
forest_warming_aic <- cbind(forest_warming_aic, data.frame(rep_prob = (c( AIC(rep_prob_forest_warming_1), AIC(rep_prob_forest_warming_2), AIC(rep_prob_forest_warming_3) ))))
hedgerow_control_aic  <- cbind(hedgerow_control_aic,  data.frame(rep_prob = (c( AIC(rep_prob_hedgerow_control_1), AIC(rep_prob_hedgerow_control_2), AIC(rep_prob_hedgerow_control_3) ))))
hedgerow_warming_aic  <- cbind(hedgerow_warming_aic,  data.frame(rep_prob = (c( AIC(rep_prob_hedgerow_warming_1), AIC(rep_prob_hedgerow_warming_2), AIC(rep_prob_hedgerow_warming_3) ))))



### 4 - rep number

#habitat = forest;  treatment = control
rep_nr_forest_control_1 = glmer (rep_nr ~ 1 + (1|transplant/site)  , family= poisson, data = df_rep_nr_fc, na.action = na.omit)
rep_nr_forest_control_2 = glmer (rep_nr~ z + (1|transplant/site)  , family= poisson, data = df_rep_nr_fc, na.action = na.omit)
rep_nr_forest_control_3 = glmer (rep_nr~ z + I(z^2) + (1|transplant/site)  , family= poisson, data = df_rep_nr_fc, na.action = na.omit)

#habitat = forest;  treatment = warming
rep_nr_forest_warming_1 = glmer (rep_nr~ 1 + (1|transplant/site)  , family= poisson, data = df_rep_nr_fw, na.action = na.omit)
rep_nr_forest_warming_2 = glmer (rep_nr~ z + (1|transplant/site)  , family= poisson, data = df_rep_nr_fw, na.action = na.omit)
rep_nr_forest_warming_3 = glmer (rep_nr~ z + I(z^2) + (1|transplant/site)  , family= poisson, data = df_rep_nr_fw, na.action = na.omit)

#habitat = hedgerow;  treatment = control
rep_nr_hedgerow_control_1 = glmer (rep_nr~ 1 + (1|transplant/site)  , family= poisson, data = df_rep_nr_hc, na.action = na.omit)
rep_nr_hedgerow_control_2 = glmer (rep_nr~ z + (1|transplant/site)  , family= poisson, data = df_rep_nr_hc, na.action = na.omit)
rep_nr_hedgerow_control_3 = glmer (rep_nr~ z + I(z^2) + (1|transplant/site)  , family= poisson, data = df_rep_nr_hc, na.action = na.omit)

#habitat = hedgerow;  treatment = warming
rep_nr_hedgerow_warming_1 = glmer (rep_nr~ 1 + (1|transplant/site)  , family= poisson, data = df_rep_nr_hw, na.action = na.omit)
rep_nr_hedgerow_warming_2 = glmer (rep_nr~ z + (1|transplant/site)  , family= poisson, data = df_rep_nr_hw, na.action = na.omit)
rep_nr_hedgerow_warming_3 = glmer (rep_nr~ z + I(z^2) + (1|transplant/site)  , family= poisson, data = df_rep_nr_hw, na.action = na.omit)

#store results (AIC)
forest_control_aic  <- cbind(forest_control_aic,  data.frame(rep_nr = (c( AIC(rep_nr_forest_control_1), AIC(rep_nr_forest_control_2), AIC(rep_nr_forest_control_3) ))))
forest_warming_aic <- cbind(forest_warming_aic, data.frame(rep_nr = (c( AIC(rep_nr_forest_warming_1), AIC(rep_nr_forest_warming_2), AIC(rep_nr_forest_warming_3) ))))
hedgerow_control_aic  <- cbind(hedgerow_control_aic,  data.frame(rep_nr = (c( AIC(rep_nr_hedgerow_control_1), AIC(rep_nr_hedgerow_control_2), AIC(rep_nr_hedgerow_control_3) ))))
hedgerow_warming_aic  <- cbind(hedgerow_warming_aic,  data.frame(rep_nr = (c( AIC(rep_nr_hedgerow_warming_1), AIC(rep_nr_hedgerow_warming_2), AIC(rep_nr_hedgerow_warming_3) ))))


### 5 - seed number

#habitat = forest;  treatment = control
seed_nr_forest_control_1 = glmer (seed_nr_2019 ~ 1 + (1|transplant/site), family= poisson, data = df_seed_nr_fc, na.action = na.omit)
seed_nr_forest_control_2 = glmer (seed_nr_2019 ~ z1 + (1|transplant/site), family= poisson, data = df_seed_nr_fc, na.action = na.omit)
seed_nr_forest_control_3 = glmer (seed_nr_2019 ~ z1 + I(z1^2) + (1|transplant/site), family= poisson, data = df_seed_nr_fc, na.action = na.omit)

#habitat = forest;  treatment = warming
seed_nr_forest_warming_1 = glmer (seed_nr_2019 ~ 1 + (1|transplant/site), family= poisson, data = df_seed_nr_fw, na.action = na.omit)
seed_nr_forest_warming_2 = glmer (seed_nr_2019 ~ z1 + (1|transplant/site), family= poisson, data = df_seed_nr_fw, na.action = na.omit)
seed_nr_forest_warming_3 = glmer (seed_nr_2019 ~ z1 + I(z1^2) + (1|transplant/site), family= poisson, data = df_seed_nr_fw, na.action = na.omit)

#habitat = hedgerow;  treatment = control
seed_nr_hedgerow_control_1 = glmer (seed_nr_2019 ~ 1 + (1|transplant/site), family= poisson, data = df_seed_nr_hc, na.action = na.omit)
seed_nr_hedgerow_control_2 = glmer (seed_nr_2019 ~ z1 + (1|transplant/site), family= poisson, data = df_seed_nr_hc, na.action = na.omit)
seed_nr_hedgerow_control_3 = glmer (seed_nr_2019 ~ z1 + I(z1^2) + (1|transplant/site), family= poisson, data = df_seed_nr_hc, na.action = na.omit)

#habitat = hedgerow;  treatment = warming
seed_nr_hedgerow_warming_1 = glmer (seed_nr_2019 ~ 1 + (1|transplant/site), family= poisson, data = df_seed_nr_hw, na.action = na.omit)
seed_nr_hedgerow_warming_2 = glmer (seed_nr_2019 ~ z1 + (1|transplant/site), family= poisson, data = df_seed_nr_hw, na.action = na.omit)
seed_nr_hedgerow_warming_3 = glmer (seed_nr_2019 ~ z1 + I(z1^2) + (1|transplant/site), family= poisson, data = df_seed_nr_hw, na.action = na.omit)

#store results (AIC)
forest_control_aic  <- cbind(forest_control_aic,  data.frame(seed_nr = (c( AIC(seed_nr_forest_control_1), AIC(seed_nr_forest_control_2), AIC(seed_nr_forest_control_3) ))))
forest_warming_aic <- cbind(forest_warming_aic, data.frame(seed_nr = (c( AIC(seed_nr_forest_warming_1), AIC(seed_nr_forest_warming_2), AIC(seed_nr_forest_warming_3) ))))
hedgerow_control_aic  <- cbind(hedgerow_control_aic,  data.frame(seed_nr = (c( AIC(seed_nr_hedgerow_control_1), AIC(seed_nr_hedgerow_control_2), AIC(seed_nr_hedgerow_control_3) ))))
hedgerow_warming_aic  <- cbind(hedgerow_warming_aic,  data.frame(seed_nr = (c( AIC(seed_nr_hedgerow_warming_1), AIC(seed_nr_hedgerow_warming_2), AIC(seed_nr_hedgerow_warming_3) ))))

##### Summarize AIC results in one dataframe and get the models with the lowest AIC in "final_models" dataframe #####
#put together all the dataframes with the AIC values
results <- rbind(forest_control_aic, forest_warming_aic, hedgerow_control_aic, hedgerow_warming_aic)

#create the final data frame with the best models in it
min_data     = aggregate(results[, 3:7], data.frame(results$factor), min)   #this creates a data frame with the minimum AIC of the 3 models of each factor
final_models = data.frame (matrix(NA, ncol = 6, nrow =4))                   #create an empty dataframe of same dimensions as the above df
names(final_models) = names (min_data)                                      #name it as min_data
#get the best models for each factor 
for (j in 2: 6) {
  for (i in 1:4) {
    value = min_data [,j] [i]                   #for "i" and for "j" --> to locate the value in "min_data"
    for (z in 3:7) {
      column = results [,z]                     #for "z" to select the column (the population parameter model) in "results"
      for (k in 1:12) {                         #for "k" to go through all the values of the column above
        if (value == column [k]) {
          final_models [,j] [i] = as.character(results [,2] [k])      #if the ij value of "min_data" is equal to any value of column [k], print it to the final data frame "final_ models"
          next #reset iteration
        }
      }
    }
  }
}
#finish first column of the data frame
final_models$results.factor = min_data$results.factor
colnames(final_models)[1] = "factor"

#replace the z in seed_nr model by "z1"
final_models$seed_nr <- gsub("z", "z1", final_models$seed_nr)


#Save the best models data frame
write.csv(final_models, file= "C:/Users/pablo/Desktop/submitt!!/figshare/final docs/2 - model selection/best_models.csv", row.names = FALSE) # Get best fitting models
write.csv(results, "C:/Users/pablo/Desktop/submitt!!/figshare/final docs/2 - model selection/all_AIC.csv") # Get all AICs

