############################################################################################################################################
#################################################### 4 - IPM FUNCTIONS  ####################################################################
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
###################################################  SENSIVITY AND ELASTICITY ##############################################################
############################################################################################################################################

library(lme4)
library(fields)
library(viridis)
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
nBigMatrix = 100                                                   #number of mesh points 

##### 1 - Run the IPM without bootstrap (only the original data) 
#[the IPM code is the same, just modified to remove the boostrap part]

{
  #create the subsets for each factor combination
  forest_control = data[which(data$habitat == "F" & data$treatment == "C"),]
  forest_warming = data[which(data$habitat == "F" & data$treatment == "T"),]
  hedgerow_control = data[which(data$habitat == "H" & data$treatment == "C"),]
  hedgerow_warming = data[which(data$habitat == "H" & data$treatment == "T"),]
  
  data_list = list (forest_control, forest_warming, hedgerow_control, hedgerow_warming)
  
  #integration limits
  lower_fc = min(min(forest_control$z, na.rm= TRUE), min(forest_control$z1, na.rm= TRUE))
  upper_fc = max(max(forest_control$z, na.rm= TRUE), max(forest_control$z1, na.rm= TRUE))
  lower_fw = min(min(forest_warming$z, na.rm= TRUE), min(forest_warming$z1, na.rm= TRUE))
  upper_fw = max(max(forest_warming$z, na.rm= TRUE), max(forest_warming$z1, na.rm= TRUE))
  lower_hc = min(min(hedgerow_control$z, na.rm= TRUE), min(hedgerow_control$z1, na.rm= TRUE))
  upper_hc = max(max(hedgerow_control$z, na.rm= TRUE), max(hedgerow_control$z1, na.rm= TRUE))
  lower_hw = min(min(hedgerow_warming$z, na.rm= TRUE), min(hedgerow_warming$z1, na.rm= TRUE))
  upper_hw = max(max(hedgerow_warming$z, na.rm= TRUE), max(hedgerow_warming$z1, na.rm= TRUE))
  int_limits <- list (
    fc_lim = c(lower_fc, upper_fc),
    fw_lim = c(lower_fw, upper_fw),
    hc_lim = c(lower_hc, upper_hc),
    hw_lim = c(lower_hw, upper_hw)
  )
  
  
  ### 1- create data frames for modelling the seed_nr and rep_nr for each factor combination (FC, FW, HC, HW)
  df_rep_nr_fc = subset (forest_control, forest_control$rep_prob == 1)
  df_rep_nr_fw = subset (forest_warming, forest_warming$rep_prob == 1)
  df_rep_nr_hc = subset (hedgerow_control, hedgerow_control$rep_prob == 1)
  df_rep_nr_hw = subset (hedgerow_warming, hedgerow_warming$rep_prob == 1)
  
  df_seed_nr_fc = subset (forest_control, forest_control$seed_prod_2019 == 1)
  df_seed_nr_fw = subset (forest_warming, forest_warming$seed_prod_2019 == 1)
  df_seed_nr_hc = subset (hedgerow_control, hedgerow_control$seed_prod_2019 == 1)
  df_seed_nr_hw = subset (hedgerow_warming, hedgerow_warming$seed_prod_2019 == 1)
  
  
  ### 2- make models (with the best models already implemented in the code)
  
  # Forest - Control
  fc_surv      = glm  (survival       ~ z                                 , family= binomial, data = forest_control, na.action = na.omit)
  fc_growth    = lmer (z1             ~ 1            + (1|transplant/site),                   data = forest_control, na.action = na.omit)
  fc_rep_prob  = glm  (rep_prob       ~ z                                 , family= binomial, data = forest_control, na.action = na.omit)
  
  fc_rep_nr    = glmer(rep_nr         ~ z            + (1|transplant/site), family= poisson,  data = df_rep_nr_fc  , na.action = na.omit)
  fc_seed_nr   = glmer(seed_nr_2019   ~ z1           + (1|transplant/site), family= poisson,  data = df_seed_nr_fc , na.action = na.omit)
  
  # Forest - Warming
  fw_surv      = glm  (survival       ~ z                                 , family= binomial, data = forest_warming, na.action = na.omit)
  fw_growth    = lmer (z1             ~ z            + (1|transplant/site),                   data = forest_warming, na.action = na.omit)
  fw_rep_prob  = glm  (rep_prob       ~ z                                 , family= binomial, data = forest_warming, na.action = na.omit)
  
  fw_rep_nr    = glmer(rep_nr         ~ z            + (1|transplant/site), family= poisson,  data = df_rep_nr_fw  , na.action = na.omit)
  fw_seed_nr   = glmer(seed_nr_2019   ~ 1            + (1|transplant/site), family= poisson,  data = df_seed_nr_fw , na.action = na.omit)
  
  # Hedgerow - Control
  hc_surv      = glm  (survival       ~ z                                 , family= binomial, data = hedgerow_control, na.action = na.omit)
  hc_growth    = lmer (z1             ~ 1            + (1|transplant/site),                   data = hedgerow_control, na.action = na.omit)
  hc_rep_prob  = glm  (rep_prob       ~ z                                 , family= binomial, data = hedgerow_control, na.action = na.omit)
  
  hc_rep_nr    = glmer(rep_nr         ~ 1            + (1|transplant/site), family= poisson , data = df_rep_nr_hc    , na.action = na.omit)
  hc_seed_nr   = glmer(seed_nr_2019   ~ z1 + I(z1^2) + (1|transplant/site), family= poisson , data = df_seed_nr_hc   , na.action = na.omit)
  
  # Hedgerow - Warming
  hw_surv      = glm  (survival       ~ 1                                 , family= binomial, data = hedgerow_warming, na.action = na.omit)
  hw_growth    = lmer (z1             ~ z            + (1|transplant/site),                   data = hedgerow_warming, na.action = na.omit)
  hw_rep_prob  = glm  (rep_prob       ~ z                                 , family= binomial, data = hedgerow_warming, na.action = na.omit)
  
  hw_rep_nr    = glmer(rep_nr         ~ z            + (1|transplant/site), family= poisson , data = df_rep_nr_hw    , na.action = na.omit)
  hw_seed_nr   = glmer(seed_nr_2019   ~ z1           + (1|transplant/site), family= poisson , data = df_seed_nr_hw   , na.action = na.omit)
  
  #create list with all the models
  lista = list (fc_surv,fc_growth, fc_rep_prob, fc_rep_nr, fc_seed_nr,
                fw_surv,fw_growth, fw_rep_prob, fw_rep_nr, fw_seed_nr,
                hc_surv,hc_growth, hc_rep_prob, hc_rep_nr, hc_seed_nr,
                hw_surv,hw_growth, hw_rep_prob, hw_rep_nr, hw_seed_nr)
  
  ### 3 - Store the parameters in a dataframe
  params = data.frame(int= numeric(), x= numeric(), x_2= numeric(), sd= numeric())
  
  for (i in 1: (5*4)) {  # 5 = total number of parameters estimated ; 4 = total number of factors analized
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
    all_factors = unique(params$factor)                                  #factor that is going to be picked
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
    #here L and U are selected from the min / max value of the columns z and z1
  }
  IPM_list = list (IPM_forest_control, IPM_forest_warming, IPM_hedgerow_control, IPM_hedgerow_warming)
  
}



##### 2 - Elasticity

#growth rate (dominant eigenvalue)
growth_fc = Re(eigen(IPM_forest_control$K)$values[1])
growth_fw = Re(eigen(IPM_forest_warming$K)$values[1])
growth_hc = Re(eigen(IPM_hedgerow_control$K)$values[1])
growth_hw = Re(eigen(IPM_hedgerow_warming$K)$values[1])

#dominant right eigen vector
w_fc = Re(eigen(IPM_forest_control$K)$vectors [,1])
w_fw = Re(eigen(IPM_forest_warming$K)$vectors [,1])
w_hc = Re(eigen(IPM_hedgerow_control$K)$vectors [,1])
w_hw = Re(eigen(IPM_hedgerow_warming$K)$vectors [,1])

#normalize to get the stable stage distribution
stable_dis_fc = w_fc / sum (w_fc)
stable_dis_fw = w_fw / sum (w_fw)
stable_dis_hc = w_hc / sum (w_hc)
stable_dis_hw = w_hw / sum (w_hw)

#dominant left eigenvector
v_fc = Re(eigen(t(IPM_forest_control$K))$vectors [,1])
v_fw = Re(eigen(t(IPM_forest_warming$K))$vectors [,1])
v_hc = Re(eigen(t(IPM_hedgerow_control$K))$vectors [,1])
v_hw = Re(eigen(t(IPM_hedgerow_warming$K))$vectors [,1])

#normalize to get the reproductive value
rep_val_fc = v_fc / v_fc [1]
rep_val_fw = v_fw / v_fw [1]
rep_val_hc = v_hc / v_hc [1]
rep_val_hw = v_hw / v_hw [1]


#get width of mesh points cells
h_fc = diff (IPM_forest_control$meshpts [1:2])
h_fw = diff (IPM_forest_warming$meshpts [1:2])
h_hc = diff (IPM_hedgerow_control$meshpts [1:2])
h_hw = diff (IPM_hedgerow_warming$meshpts [1:2])

x_fc = sum(stable_dis_fc*rep_val_fc)*h_fc
x_fw = sum(stable_dis_fw*rep_val_fw)*h_fw
x_hc = sum(stable_dis_hc*rep_val_hc)*h_hc
x_hw = sum(stable_dis_hw*rep_val_hw)*h_hw

#sensivity
sensivity_fc = outer (rep_val_fc, stable_dis_fc)/x_fc
sensivity_fw = outer (rep_val_fw, stable_dis_fw)/x_fw
sensivity_hc = outer (rep_val_hc, stable_dis_hc)/x_hc
sensivity_hw = outer (rep_val_hw, stable_dis_hw)/x_hw

#elasticity
elasticity_fc = matrix(as.vector(sensivity_fc)*as.vector(IPM_forest_control$K)/growth_fc, nrow = nBigMatrix)
elasticity_fw = matrix(as.vector(sensivity_fw)*as.vector(IPM_forest_warming$K)/growth_fw, nrow = nBigMatrix)
elasticity_hc = matrix(as.vector(sensivity_hc)*as.vector(IPM_hedgerow_control$K)/growth_hc, nrow = nBigMatrix)
elasticity_hw = matrix(as.vector(sensivity_hw)*as.vector(IPM_hedgerow_control$K)/growth_hw, nrow = nBigMatrix)


#compute the separate elasticity surfaces associated with P and F
# Survival / growth component
pvals_fc <- IPM_forest_control$P / h_fc
pelasticity_fc <- pvals_fc * sensivity_fc / growth_fc

pvals_fw <- IPM_forest_warming$P / h_fw
pelasticity_fw <- pvals_fw * sensivity_fw / growth_fw

pvals_hc <- IPM_hedgerow_control$P / h_hc
pelasticity_hc <- pvals_hc * sensivity_hc / growth_hc

pvals_hw <- IPM_hedgerow_warming$P / h_hw
pelasticity_hw <- pvals_hw * sensivity_hw / growth_hw

# Fecundity component
fvals_fc <- IPM_forest_control$F / h_fc
felasticity_fc <- fvals_fc * sensivity_fc / growth_fc

fvals_fw <- IPM_forest_warming$F / h_fw
felasticity_fw <- fvals_fw * sensivity_fw / growth_fw

fvals_hc <- IPM_hedgerow_control$F / h_hc
felasticity_hc <- fvals_hc * sensivity_hc / growth_hc

fvals_hw <- IPM_hedgerow_warming$F / h_hw
felasticity_hw <- fvals_hw * sensivity_hw / growth_hw



##### 3 - Plots

### 1 - All together
dev.new(width=5, height=5)
{
  par(mfrow= c(4,2))
  #forest control
  image.plot(col = rev(plasma(12)), IPM_forest_control$meshpts, IPM_forest_control$meshpts, t(pvals_fc), xlab= "Size (t)", ylab= "Size (t+1)", main= "Elasticity Survival/Growth - Forest control")
  image.plot(col = rev(plasma(12)), IPM_forest_control$meshpts, IPM_forest_control$meshpts, t(fvals_fc), xlab= "Size (t)", ylab= "Size (t+1)", main= "Elasticity Fecundity - Forest control")

  #forest warming
  image.plot(col = rev(plasma(12)), IPM_forest_warming$meshpts, IPM_forest_warming$meshpts, t(pvals_fw), xlab= "Size (t)", ylab= "Size (t+1)", main= "Elasticity Survival/Growth - Forest Warming")
  image.plot(col = rev(plasma(12)), IPM_forest_warming$meshpts, IPM_forest_warming$meshpts, t(fvals_fw), xlab= "Size (t)", ylab= "Size (t+1)", main= "Elasticity Fecundity - Forest warming")

  #hedgerow control
  image.plot(col = rev(plasma(12)), IPM_hedgerow_control$meshpts, IPM_hedgerow_control$meshpts, t(pvals_hc), xlab= "Size (t)", ylab= "Size (t+1)", main= "Elasticity Survival/Growth - Hedgerow control")
  image.plot(col = rev(plasma(12)), IPM_hedgerow_control$meshpts, IPM_hedgerow_control$meshpts, t(fvals_hc), xlab= "Size (t)", ylab= "Size (t+1)", main= "Elasticity Fecundity - Hedgerow control")

  #hedgerow warming
  image.plot(col = rev(plasma(12)), IPM_hedgerow_warming$meshpts, IPM_hedgerow_warming$meshpts, t(pvals_hw), xlab= "Size (t)", ylab= "Size (t+1)", main= "Elasticity Survival/Growth - Hedgerow Warming")
  image.plot(col = rev(plasma(12)), IPM_hedgerow_warming$meshpts, IPM_hedgerow_warming$meshpts, t(fvals_hw), xlab= "Size (t)", ylab= "Size (t+1)", main= "Elasticity Fecundity - Hedgerow warming")
  
}

