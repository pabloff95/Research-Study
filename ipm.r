############################################################################################################################################
####################################################### IPM FUNCTIONS  #####################################################################
############################################################################################################################################
rm(list=ls()) #clean global environment

###### 1 - Demographic functions

#Kernel: K(z',z) = P(z', z) + F(z', z) = s(z) � G(z', z) + p rep(z) � b rep (z) � p seed (z) � b seed (z) � p germination � recruit_size (z')

#survival and growth kernel: P(z', z)
P_z1z <- function (z1, z, par) {
  return (s_z(z, par) * G_z1z(z1, z, par))
}

#components:
#survival probability: S(z)
s_z <- function (z, par) {
  linear.p <- as.numeric(par["survival_int"]) + as.numeric(par["survival_x"])*z + as.numeric(par["survival_x_2"])*(z^2)   #linear predictor: calculates the predicted y, in logit scale, for the values of size (z)
  p <- 1/(1 + exp (-linear.p))                          #inverse logit transformation to probability: transformation of Y from logit scale to probability
  return (p)                                            #returns survival probability
}

#probability density function for size at the next census, conditional on current size z: G(z',z)
G_z1z <- function (z1, z, par) {
  mean <- as.numeric(par["growth_int"]) + as.numeric(par["growth_x"])*z + as.numeric(par["growth_x_2"])*(z^2) #this calculates the mean size for next year
  sd <- as.numeric(par ["growth_sd"])                                       #sd 
  pd_growth <- dnorm(z1, mean= mean, sd= sd)                                #probability density function of new size z1, for current size z
  return (pd_growth)
}

#fecundity kernel:
F_z1z <- function (z1, z, par, p_germin = p_germ) {    
  return (prep_z (z, par) * brep_z (z, par) * bseed_z (z, par) * p_germin * recruit_size_pb (z1))
}

#components
#rep probability: p rep(z)
prep_z <- function (z, par) {
  lin.p <- as.numeric(par["rep_prob_int"]) + as.numeric(par["rep_prob_x"])*z + as.numeric(par["rep_prob_x_2"])*(z^2)
  p <- 1/(1+exp(-lin.p))
  return (p)                                #returns rep probability
}

#rep mean number
brep_z <- function (z, par) {
  x <- exp (as.numeric(par["rep_nr_int"]) + as.numeric(par["rep_nr_x"])*z + as.numeric(par["rep_nr_x_2"])*(z^2))
  return(x)
}

#seed mean number
bseed_z <- function (z, par) {
  x <- exp (as.numeric(par["seed_nr_int"]) + as.numeric(par["seed_nr_x"])*z + as.numeric(par["seed_nr_x_2"])*(z^2))
  return(x)
}

#recruitment size
recruit_size_pb <- function (z1, recruit_mean = rec_mean, recruit_sd= rec_sd) {
  pd_rec <- dnorm(z1, mean= recruit_mean, sd= recruit_sd)
  return (pd_rec)
}

##### 2 - Matrix function

#Function for the interation matrix (applies the midllepoint rule)
mk_K <- function(m, par, L, U, U1=U, L1=L) {
  h <- (U - L)/m                                                          #heith to apply middlepoint rule
  meshpts <- L + ((1:m) - 1/2) * h                                        #sets the pint where the heith is measured
  P <- h * (outer(meshpts, pmin(pmax(meshpts,L1),U1), P_z1z, par = par))  #outer help to compute the matrix
  F <- h * (outer(meshpts, pmin(pmax(meshpts,L1),U1), F_z1z, par = par))  #pmin and pmax to set a floor and ceiling for the integration limits
  K <- P + F
  return(list(K = K, meshpts = meshpts, P = P, F = F))
}

###### 3 - Confidence interval function 

confidence_interval <- function(vector, interval) {
  # Standard deviation of sample
  vec_sd <- sd(vector)
  # Sample size
  n <- length(vector)
  # Mean of sample
  vec_mean <- mean(vector)
  # Error according to t distribution
  error <- qt((interval + 1)/2, df = n - 1) * vec_sd / sqrt(n)
  # Confidence interval as a vector
  result <- c("lower" = vec_mean - error, "upper" = vec_mean + error)
  return(result)
}


############################################################################################################################################
######################################################## IPM CODE ##########################################################################
############################################################################################################################################

#Packages
library(lme4)
library(ggplot2)
library(gridExtra)

#set seed
set.seed(321) 

#load data file
data <- read.csv(file= ".../data.csv")
data$site <- factor(data$site)
data$individual<- factor (data$individual)
data$survival<- factor (data$survival)
data$rep_prob<- factor (data$rep_prob)
data$seed_prod_2019<- factor (data$seed_prod_2019)


#constant terms
p_germ = 0.14                          #germination probability (Verheyen and Henry, 2004)
rec_mean = 7                           #Bieberich et al 2018: recruitment mean size
rec_sd = 2                             #Bieberich et al 2018: recruitment sd size

#IPM computation
n.rep = 1000                                                       #number of iterations of the boostrapp loop
nBigMatrix = 100                                                   #number of mesh points 


results = data.frame()  #used to store the results

#create the subsets for each factor combination
forest_control = data[which(data$habitat == "F" & data$treatment == "C"),]
forest_warming = data[which(data$habitat == "F" & data$treatment == "T"),]
hedgerow_control = data[which(data$habitat == "H" & data$treatment == "C"),]
hedgerow_warming = data[which(data$habitat == "H" & data$treatment == "T"),]

#set limitis of integration
lower_fc =  min(min(forest_control$z, na.rm= TRUE), min(forest_control$z1, na.rm= TRUE))
upper_fc =  max(max(forest_control$z, na.rm= TRUE), max(forest_control$z1, na.rm= TRUE))
lower_fw =  min(min(forest_warming$z, na.rm= TRUE), min(forest_warming$z1, na.rm= TRUE))
upper_fw =  max(max(forest_warming$z, na.rm= TRUE), max(forest_warming$z1, na.rm= TRUE))
lower_hc =  min(min(hedgerow_control$z, na.rm= TRUE), min(hedgerow_control$z1, na.rm= TRUE))
upper_hc =  max(max(hedgerow_control$z, na.rm= TRUE), max(hedgerow_control$z1, na.rm= TRUE))
lower_hw =  min(min(hedgerow_warming$z, na.rm= TRUE), min(hedgerow_warming$z1, na.rm= TRUE))
upper_hw =  max(max(hedgerow_warming$z, na.rm= TRUE), max(hedgerow_warming$z1, na.rm= TRUE))
int_limits <- list (
  fc_lim = c(lower_fc, upper_fc),
  fw_lim = c(lower_fw, upper_fw),
  hc_lim = c(lower_hc, upper_hc),
  hw_lim = c(lower_hw, upper_hw)
)


for (bt in 1:n.rep) {
  #sample individuals from the dataset for bootstrap
  data.s_fc <- sample(1:nrow(forest_control), size = dim(forest_control)[1],replace = TRUE) #takes the row number of 48 (dim (data)) individuals of the original dataset
  #replace= TRUE: the same individual can be picked more than once -> bootstrap
  data.s_fw <- sample(1:nrow(forest_warming), size = dim(forest_warming)[1],replace = TRUE)
  data.s_hc <- sample(1:nrow(hedgerow_control), size = dim(hedgerow_control)[1],replace = TRUE)
  data.s_hw <- sample(1:nrow(hedgerow_warming), size = dim(hedgerow_warming)[1],replace = TRUE)
  
  data_boot_fc <- forest_control [data.s_fc,]                                        #construct the re-sampled data in a new dataframe
  data_boot_fw <- forest_warming [data.s_fw,]        
  data_boot_hc <- hedgerow_control [data.s_hc,]        
  data_boot_hw <- hedgerow_warming [data.s_hw,]        
  
  data_list = list(data_boot_fc, data_boot_fw, data_boot_hc, data_boot_hw)
  
  ### 1- create subsets for each factor combination that is going to be analysed
  
  #create data frames for modelling the seed_nr and rep_nr for each factor combination (FC, FW, HC, HW)
  df_rep_nr_fc = subset (data_boot_fc, data_boot_fc$rep_prob == 1)
  df_rep_nr_fw = subset (data_boot_fw, data_boot_fw$rep_prob == 1)
  df_rep_nr_hc = subset (data_boot_hc, data_boot_hc$rep_prob == 1)
  df_rep_nr_hw = subset (data_boot_hw, data_boot_hw$rep_prob == 1)
  
  df_seed_nr_fc = subset (data_boot_fc, data_boot_fc$seed_prod_2019 == 1)
  df_seed_nr_fw = subset (data_boot_fw, data_boot_fw$seed_prod_2019 == 1)
  df_seed_nr_hc = subset (data_boot_hc, data_boot_hc$seed_prod_2019 == 1)
  df_seed_nr_hw = subset (data_boot_hw, data_boot_hw$seed_prod_2019 == 1)
  
  
  ### 2- make models (with the best models already implemented in the code)
  
  # Forest - Control
  fc_surv      = glm   (survival       ~ z                                 , family= binomial, data = data_boot_fc , na.action = na.omit)
  fc_growth    = lmer  (z1             ~ 1  + (1|transplant/site)          ,                   data = data_boot_fc , na.action = na.omit)
  fc_rep_prob  = glm   (rep_prob       ~ z                                 , family= binomial, data = data_boot_fc , na.action = na.omit)
  
  fc_rep_nr    = glmer (rep_nr         ~ z  + (1|transplant/site)          , family= poisson,  data = df_rep_nr_fc , na.action = na.omit)
  fc_seed_nr   = glmer (seed_nr_2019   ~ z1 + (1|transplant/site)          , family= poisson,  data = df_seed_nr_fc, na.action = na.omit)
  
  # Forest - Warming
  fw_surv      = glm   (survival       ~ z                                 , family= binomial, data = data_boot_fw , na.action = na.omit)
  fw_growth    = lmer  (z1             ~ z  + (1|transplant/site)          ,                   data = data_boot_fw , na.action = na.omit)
  fw_rep_prob  = glm   (rep_prob       ~ z                                 , family= binomial, data = data_boot_fw , na.action = na.omit)
  
  fw_rep_nr    = glmer (rep_nr         ~ z  + (1|transplant/site)          , family= poisson,  data = df_rep_nr_fw , na.action = na.omit)
  fw_seed_nr   = glmer (seed_nr_2019   ~ 1  + (1|transplant/site)          , family= poisson,  data = df_seed_nr_fw, na.action = na.omit)
  
  # Hedgerow - Control
  hc_surv      = glm   (survival       ~ z                                 , family= binomial, data = data_boot_hc , na.action = na.omit)
  hc_growth    = lmer  (z1             ~ 1  + (1|transplant/site)          ,                   data = data_boot_hc , na.action = na.omit)
  hc_rep_prob  = glm   (rep_prob       ~ z                                 , family= binomial, data = data_boot_hc , na.action = na.omit)
  
  hc_rep_nr    = glmer (rep_nr         ~ 1  + (1|transplant/site)          , family= poisson,  data = df_rep_nr_hc , na.action = na.omit)
  hc_seed_nr   = glmer (seed_nr_2019   ~ z1 + I(z1^2) + (1|transplant/site), family= poisson,  data = df_seed_nr_hc, na.action = na.omit)
  
  # Hedgerow - Warming
  hw_surv      = glm   (survival       ~ 1                                 , family= binomial, data = data_boot_hw , na.action = na.omit)
  hw_growth    = lmer  (z1             ~ z  + (1|transplant/site)          ,                   data = data_boot_hw , na.action = na.omit)
  hw_rep_prob  = glm   (rep_prob       ~ z                                 , family= binomial, data = data_boot_hw , na.action = na.omit)
  
  hw_rep_nr    = glmer (rep_nr         ~ z  + (1|transplant/site)          , family= poisson,  data = df_rep_nr_hw , na.action = na.omit)
  hw_seed_nr   = glmer (seed_nr_2019   ~ z1 + (1|transplant/site)          , family= poisson,  data = df_seed_nr_hw, na.action = na.omit)
  
  #create list with all the models
  lista = list (fc_surv,fc_growth, fc_rep_prob, fc_rep_nr, fc_seed_nr,
                fw_surv,fw_growth, fw_rep_prob, fw_rep_nr, fw_seed_nr,
                hc_surv,hc_growth, hc_rep_prob, hc_rep_nr, hc_seed_nr,
                hw_surv,hw_growth, hw_rep_prob, hw_rep_nr, hw_seed_nr)
  
  ### 3 - Store the parameters in a dataframe
  params = data.frame(int= numeric(), x= numeric(), x_2= numeric(), sd= numeric())
  
  for (i in 1: (5*4)) {  # 5 = total number of vital rates regressed ; 4 = total number of factors analized
    x = data.frame(summary(lista [[i]])$coef) [,1]
    sd = summary(lista [[i]])$sigma
    if (is.null(sd)) {sd = 0}
    new.row = data.frame (int = x [1], x= x [2], x_2= x [3], sd = sd)
    params = rbind (params, new.row)
  }
  #add the informative variables
  params$factor= c(rep("forest_control", 5), rep("forest_warming", 5), rep("hedgerow_control", 5), rep("hedgerow_warming", 5))
  params$parameter = rep (c("survival", "growth", "rep_prob", "rep_nr", "seed_nr"), 4)
  #replace NAs by 0
  params [is.na(params)] = 0
  
  
  ### 4 - Create the parameter variables (to be stored in a list: param_list)
  param_list = list()
  names_param_list = character()
  
  for (j in 1:length(unique(params$factor))) {
    
    #iniciate the variables to store data and select the factor (FC, FW, HC and HW)
    all_factors = unique(params$factor)                                  #combination that is going to be selected
    factor_df = subset(params, params$factor == all_factors [j])         #get the data from that factor
    factor_parameters = numeric()                                        #vector where the the parameters of all the  kernel functions of the factor are going to be stored
    
    #collect all the parameters of 1 factor
    for (i in 1:length(factor_df$parameter)) {
      int = paste0(factor_df$parameter[i],"_",names(factor_df)[1])       #these functions create the names for each parameter (int, x, x2, sd) (ex.: rep_prob_x_2)
      x   = paste0(factor_df$parameter[i],"_",names(factor_df)[2])
      x2  = paste0(factor_df$parameter[i],"_",names(factor_df)[3])
      sd  = paste0(factor_df$parameter[i],"_",names(factor_df)[4])
      factor_df.names = c(int, x, x2, sd)
      
      values = c(as.numeric(factor_df [i,] [c(1:4)]))                   #here all the values (int, x, x2, sd) of each parameter are collected
      names (values) = factor_df.names
      
      factor_parameters = c(factor_parameters, values)                  #store all the parameters (int, x, x2, sd) of each kernel component (surv, growth, ...) of one unique factor (FC, FW, HC, HW) in a variable
    }
    
    #include the kernel function parameters of each factor in a list
    variable_name = paste0(all_factors[j],"_parameters") #create variable a name for each factor, ex: hedgerow_warming_parameters
    param_list [[j]] =  factor_parameters                #include the vector
    names_param_list [j] = as.character(variable_name)   #names of eahc vector to be added when the loop finishes 
  }
  names(param_list) = names_param_list
  
  ### 5 - Calculate the IPM matrix
  for (w in 1: length(param_list)){                  #creates the integration matrices
    limit = int_limits [[w]]
    factor_combination= as.character(unique(params$factor) [w])
    name = paste0("IPM_",factor_combination)
    assign(name, mk_K(m= nBigMatrix, par= param_list [[w]], L= limit [1] , U= limit [2] ))
  }
  IPM_list = list (IPM_forest_control, IPM_forest_warming, IPM_hedgerow_control, IPM_hedgerow_warming) #store IPM matrices in list
  
  ### 6 - Calculate growth rate
  for (i in 1:length(IPM_list)) {
    factor_comb = c("forest-control", "forest-warming", "hedgerow-control", "hedgerow-warming")
    ipm = IPM_list [[i]]
    eigenval = Re(eigen(ipm$K)$values [1])
    new_gr = data.frame(Factor = factor_comb [i], Growth_rate = eigenval)
    results = rbind(results, new_gr)
  }
  cat(bt, "\n")
}



### 7 - Summarize results in "summary_results" dataframe
summary_results = rbind(
  data.frame(
    factor = "forest-control",
    mean = mean(results$Growth_rate[results$Factor == "forest-control"]),
    sd = sd(results$Growth_rate[results$Factor == "forest-control"]),
    lower_CI = confidence_interval(vector = results$Growth_rate[results$Factor == "forest-control"], 0.95) [1],
    upper_CI = confidence_interval(vector = results$Growth_rate[results$Factor == "forest-control"], 0.95) [2]
  ),
  data.frame(
    factor = "forest-warming",
    mean = mean(results$Growth_rate[results$Factor == "forest-warming"]),
    sd = sd(results$Growth_rate[results$Factor == "forest-warming"]),
    lower_CI = confidence_interval(vector = results$Growth_rate[results$Factor == "forest-warming"], 0.95) [1],
    upper_CI = confidence_interval(vector = results$Growth_rate[results$Factor == "forest-warming"], 0.95) [2]
  ),
  data.frame(
    factor = "hedgerow-control",
    mean = mean(results$Growth_rate[results$Factor == "hedgerow-control"]),
    sd = sd(results$Growth_rate[results$Factor == "hedgerow-control"]),
    lower_CI = confidence_interval(vector = results$Growth_rate[results$Factor == "hedgerow-control"], 0.95) [1],
    upper_CI = confidence_interval(vector = results$Growth_rate[results$Factor == "hedgerow-control"], 0.95) [2]
  ),
  data.frame(
    factor = "hedgerow-warming",
    mean = mean(results$Growth_rate[results$Factor == "hedgerow-warming"]),
    sd = sd(results$Growth_rate[results$Factor == "hedgerow-warming"]),
    lower_CI = confidence_interval(vector = results$Growth_rate[results$Factor == "hedgerow-warming"], 0.95) [1],
    upper_CI = confidence_interval(vector = results$Growth_rate[results$Factor == "hedgerow-warming"], 0.95) [2]
  )
)


#save the growth rates (+ all the growth rates generated during the bootstrap)
write.csv(summary_results, file= "C:/Users/pablo/Desktop/submitt!!/figshare/final docs/6 - IPM/growth_rates.csv", row.names = FALSE)
write.csv(results, file= "C:/Users/pablo/Desktop/submitt!!/figshare/final docs/6 - IPM/growth_rates_all_bootstrap.csv", row.names = FALSE)


### 8- Plot growth rates
library(ggplot2)
grw_rate <- read.csv(file= "C:/Users/pablo/Desktop/submitt!!/figshare/final docs/6 - IPM/growth_rates_all_bootstrap.csv") ### = "results" dataframe

grw_rate$habitat = NA
for (i in 1:length(grw_rate$Factor)) {
  if (grw_rate$Factor [i] == "forest-control" || grw_rate$Factor [i] == "forest-warming"  ) {grw_rate$habitat [i] = "Forest"}
  if (grw_rate$Factor [i] == "hedgerow-control" || grw_rate$Factor [i] == "hedgerow-warming"  ) {grw_rate$habitat [i] = "Hedgerow"} 
}

grw_rate$treatment = NA
for (i in 1:length(grw_rate$Factor)) {
  if (grw_rate$Factor [i] == "forest-control" || grw_rate$Factor [i] == "hedgerow-control"  ) {grw_rate$treatment [i] = "Control"}
  if (grw_rate$Factor [i] == "forest-warming" || grw_rate$Factor [i] == "hedgerow-warming"  ) {grw_rate$treatment [i] = "Warming"} 
}

growth_rate_plot = ggplot(grw_rate, aes(x= treatment, y= Growth_rate, fill= treatment ))+ geom_violin(trim = FALSE, fill= "grey85", colour= "grey85") +
  geom_boxplot(width = 0.2, size= 1) +
  ylab("Growth rate") + xlab ("") +
  facet_grid(.~habitat) + 
  theme(legend.position="none",
        panel.background = element_rect(colour="black", fill= "snow1"),
        panel.grid.major = element_line (colour= "gray90"),
        panel.grid.minor = element_line (colour= "gray90"),
        axis.title = element_text(size=13),
        strip.text.y = element_text(size = 15, color = "black", face = "bold.italic"),
        strip.text.x = element_text(size = 15, color = "black", face = "bold.italic"),
        strip.background = element_rect(color="black", size=1.2, linetype="solid"),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.25)),
        axis.title.y = element_text(size = rel(1.05))) +
  coord_cartesian(ylim = c(0, 30)) +                             # ylim = 12 for better visualization
  geom_hline(yintercept=1, color = "red", size= 0.8)

growth_rate_plot


### Find significative differences on growth rates: Tukey test
mod = lm(Growth_rate ~  Factor ,data = grw_rate)
anova(mod)

library(multcomp)
library(abind)
adhoc <- glht(mod, linfct = mcp(Factor = "Tukey"))
summary(adhoc)

