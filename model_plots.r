############################################################################################################################################
###################################################### 3 - MODEL PLOTS (Fig 3) #####################################################################
############################################################################################################################################
rm(list=ls()) #clean global environment

#Packages
library(lme4)
library(ggplot2)
library(effects)
library(gridExtra)
library(grid)

#read data
data <- read.csv(file= ".../data.csv")
data$site <- factor(data$site)
data$individual<- factor (data$individual)
data$survival<- factor (data$survival)
data$rep_prob<- factor (data$rep_prob)
data$seed_prod_2019<- factor (data$seed_prod_2019)

forest = data[which(data$habitat == "F"),]
hedgerow = data[which(data$habitat == "H"),]

### Survival
survival_f = ggplot(forest[!is.na(forest$z),], aes(x=z, y=as.numeric(as.character(survival)), colour= treatment)) + 
  geom_smooth(aes(group=treatment, fill=treatment),method = "glm", se=TRUE, method.args= list(family="binomial"), size=1.25, formula= y ~ x,show.legend = FALSE) + #creates the model line
  geom_point(show.legend = FALSE)+
  scale_color_manual(values = c("black", "red")) + #selects the colours of the dots and line
  scale_fill_manual(values = c("grey", "rosybrown2"))+   #selects the colours of the errror bands
  ylab("Survival rate")+
  theme(panel.background = element_rect(colour="black", fill= "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank (),
        
        axis.title = element_text(size=13),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank()) +
  coord_cartesian(xlim = c(min(forest$z, na.rm=T), max(forest$z, na.rm= T)), ylim= c(0,1)) 

survival_h = ggplot(hedgerow[!is.na(hedgerow$z),], aes(x=z, y=as.numeric(as.character(survival)), colour= treatment)) + 
  geom_smooth(data= hedgerow[hedgerow$treatment == "C",],method = "glm", se=TRUE, method.args= list(family="binomial"), fill = "grey", size=1.25, formula= y ~ x,show.legend = FALSE) + 
  geom_smooth(data= hedgerow[hedgerow$treatment == "T",],method = "glm", se=TRUE, method.args= list(family="binomial"), fill = "rosybrown2",size=1.25, formula= y ~ 1,show.legend = FALSE) + # 2 geom_smooth because control and treatment had different best fitting models
  geom_point(show.legend = FALSE) +
  scale_color_manual(values = c("black", "red")) + 
  theme(panel.background = element_rect(colour="black", fill= "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank (),
        axis.title = element_text(size=13),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank()) +
  coord_cartesian(xlim = c(min(hedgerow$z, na.rm=T), max(hedgerow$z, na.rm= T)), ylim= c(0,1)) 


### Rep_prob
rep_prob_f = ggplot(forest[!is.na(forest$z),], aes(x= z, y= as.numeric(as.character(rep_prob)), colour=treatment))+
  geom_smooth(aes(group=treatment, fill=treatment),method= "glm", se=TRUE, method.args= list(family="binomial"), size=1.25, formula= y ~ x,show.legend = FALSE) +
  geom_point(show.legend = FALSE)+
  scale_color_manual(values = c("black", "red"))+
  scale_fill_manual(values = c("grey", "rosybrown2"))+
  ylab("Reproduction rate") + 
  theme(panel.background = element_rect(colour="black", fill= "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank (),
        
        axis.title = element_text(size=13),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank()) +
  coord_cartesian(xlim = c(min(forest$z, na.rm=T), max(forest$z, na.rm= T)), ylim= c(0,1)) 

rep_prob_h = ggplot(hedgerow[!is.na(hedgerow$z),], aes(x= z, y= as.numeric(as.character(rep_prob)), colour=treatment))+
  geom_smooth(aes(group=treatment, fill=treatment),method= "glm", se=TRUE, method.args= list(family="binomial"), size=1.25, formula= y ~ x,show.legend = FALSE) +
  geom_point(show.legend = FALSE)+
  scale_color_manual(values = c("black", "red"))+
  scale_fill_manual(values = c("grey", "rosybrown2"))+
  ylab("Reproduction rate") + 
  theme(panel.background = element_rect(colour="black", fill= "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank (),
        
        axis.title = element_text(size=13),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank()) +
  coord_cartesian(xlim = c(min(hedgerow$z, na.rm=T), max(hedgerow$z, na.rm= T)), ylim= c(0,1)) 


### Seed_prob
seed_prob_f = ggplot(forest[!is.na(forest$z),], aes(x= z1, y= as.numeric(as.character(seed_prod_2019)), colour=treatment))+
  geom_smooth(aes(group=treatment, fill=treatment),method= "glm", se=TRUE, method.args= list(family="binomial"), size=1.25, formula= y ~ x,show.legend = FALSE) +
  geom_point(show.legend = FALSE)+
  scale_color_manual(values = c("black", "red"))+
  scale_fill_manual(values = c("grey", "rosybrown2"))+
  ylab("Seed rate") + 
  theme(panel.background = element_rect(colour="black", fill= "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank (),
        
        axis.title = element_text(size=13),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank()) +
  coord_cartesian(xlim = c(min(forest$z, na.rm=T), max(forest$z, na.rm= T)), ylim= c(0,1)) 

seed_prob_h = ggplot(hedgerow[!is.na(hedgerow$z),], aes(x= z1, y= as.numeric(as.character(seed_prod_2019)), colour=treatment))+
  geom_smooth(aes(group=treatment, fill=treatment),method= "glm", se=TRUE, method.args= list(family="binomial"), size=1.25, formula= y ~ x,show.legend = FALSE) +
  geom_point(show.legend = FALSE)+
  scale_color_manual(values = c("black", "red"))+
  scale_fill_manual(values = c("grey", "rosybrown2"))+
  ylab("Seed rate") + 
  theme(panel.background = element_rect(colour="black", fill= "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank (),
        
        axis.title = element_text(size=13),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank()) +
  coord_cartesian(xlim = c(min(hedgerow$z, na.rm=T), max(hedgerow$z, na.rm= T)), ylim= c(0,1)) 

### Growth
#forests
fc_data = forest[which(forest$treatment == "C"),]
fw_data = forest[which(forest$treatment == "T"),]

fc_growth    = lmer  (z1 ~ 1 + (1|transplant/site),  data = fc_data, na.action = na.omit)
fw_growth    = lmer  (z1 ~ z + (1|transplant/site),  data = fw_data, na.action = na.omit)
#get fitted values of not null models
fitted_val_fw_growth = as.data.frame(Effect("z", fw_growth, xlevels = length(fw_data$z)))
fitted_val_fw_growth$z1 = fw_data$z1
fitted_val_fw_growth$treatment = "T"

growth_f = ggplot(forest[!is.na(forest$z),], aes(z, z1, colour = treatment)) + 
  geom_point(show.legend = FALSE) +
  geom_ribbon(data = fitted_val_fw_growth, aes(ymin=lower, ymax=upper), colour = NA, fill = "rosybrown2", alpha=0.5) + #errror bands of GLMM
  geom_line(y=summary(fc_growth)$coef [1], size = 1.25, colour = "black",show.legend = FALSE) + #line of null model (treatment = C)
  geom_line(data = fitted_val_fw_growth, aes(z, fit), size = 1.25,show.legend = FALSE) + #line of glmm model (treatment = T)
  scale_color_manual(values = c("black", "red"))+
  ylab("Size t+1 (cm)") + 
  theme(panel.background = element_rect(colour="black", fill= "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank (),
        axis.title = element_text(size=13),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank()) +
  coord_cartesian(xlim = c(min(forest$z, na.rm=T), max(forest$z, na.rm= T)), ylim= c(min(data$z1, na.rm = T),max(data$z1, na.rm = T))) 

#hedgerows
hc_data = hedgerow[which(forest$treatment == "C"),]
hw_data = hedgerow[which(forest$treatment == "T"),]
hc_growth    = lmer  (z1 ~ 1 + (1|transplant/site),  data = hc_data, na.action = na.omit)
hw_growth    = lmer  (z1 ~ z + (1|transplant/site),  data = hw_data, na.action = na.omit)
#get fitted values of not null models
fitted_val_hw_growth = as.data.frame(Effect("z", hw_growth, xlevels = length(hw_data$z)))
fitted_val_hw_growth$z1 = hw_data$z1
fitted_val_hw_growth$treatment = "T"

growth_h = ggplot(hedgerow[!is.na(hedgerow$z),], aes(z, z1, colour = treatment)) + 
  geom_point(show.legend = FALSE) +
  geom_ribbon(data = fitted_val_hw_growth, aes(ymin=lower, ymax=upper), colour = NA, fill = "rosybrown2", alpha=0.5) + #errror bands of GLMM
  geom_line(y=summary(hc_growth)$coef [1], size = 1.25, colour = "black",show.legend = FALSE) + #line of null model (treatment = C)
  geom_line(data = fitted_val_hw_growth, aes(z, fit), size = 1.25,show.legend = FALSE) + #line of glmm model (treatment = T)
  scale_color_manual(values = c("black", "red"))+
  theme(panel.background = element_rect(colour="black", fill= "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank (),
        axis.title = element_text(size=13),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank()) +
  coord_cartesian(xlim = c(min(hedgerow$z, na.rm=T), max(hedgerow$z, na.rm= T)), ylim= c(min(data$z1, na.rm = T),max(data$z1, na.rm = T))) 

### Rep number
#forests
df_rep_nr_fc = subset (fc_data, fc_data$rep_prob == 1)
df_rep_nr_fw = subset (fw_data, fw_data$rep_prob == 1)
fc_rep_nr    = glmer (rep_nr~ z + (1|transplant/site)  , family= poisson, data = df_rep_nr_fc, na.action = na.omit)
fw_rep_nr    = glmer (rep_nr~ z + (1|transplant/site)  , family= poisson, data = df_rep_nr_fw, na.action = na.omit)
#get fitted values
fitted_val_fc_rep_nr = as.data.frame(Effect("z", fc_rep_nr, xlevels = length(df_rep_nr_fc$z)))
fitted_val_fc_rep_nr$rep_nr = df_rep_nr_fc$rep_nr
fitted_val_fc_rep_nr$treatment = "C"

fitted_val_fw_rep_nr = as.data.frame(Effect("z", fw_rep_nr, xlevels = length(df_rep_nr_fw$z)))
fitted_val_fw_rep_nr$rep_nr = df_rep_nr_fw$rep_nr
fitted_val_fw_rep_nr$treatment = "T"

fitted_val_rep_nr = rbind(fitted_val_fc_rep_nr, fitted_val_fw_rep_nr)

forest_rep = rbind(df_rep_nr_fc, df_rep_nr_fw)

rep_nr_f = ggplot(forest_rep[!is.na(forest_rep$z),], aes(z, rep_nr, colour = treatment)) + 
  geom_point(show.legend = FALSE) +
  geom_ribbon(data = fitted_val_rep_nr[which(fitted_val_rep_nr$treatment == "C"),], aes(ymin=lower, ymax=upper), colour= NA, fill = "grey", alpha=0.5) + 
  geom_ribbon(data = fitted_val_rep_nr[which(fitted_val_rep_nr$treatment == "T"),], aes(ymin=lower, ymax=upper), colour = NA, fill = "rosybrown2", alpha=0.5) +
  geom_line(data = fitted_val_rep_nr[which(fitted_val_rep_nr$treatment == "C"),], aes(z, fit), size = 1.25,show.legend = FALSE) +
  geom_line(data = fitted_val_rep_nr[which(fitted_val_rep_nr$treatment == "T"),], aes(z, fit), size = 1.25,show.legend = FALSE) +
  scale_color_manual(values = c("black", "red"))+
  ylab("Number of flower heads") + 
  theme(panel.background = element_rect(colour="black", fill= "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank (),
        axis.title = element_text(size=13),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank()) +
  coord_cartesian(xlim = c(min(forest$z, na.rm=T), max(forest$z, na.rm= T)), ylim= c(min(data$rep_nr, na.rm = T),max(data$rep_nr, na.rm = T))) 

#hedgerows
df_rep_nr_hc = subset (hc_data, fc_data$rep_prob == 1)
df_rep_nr_hw = subset (hw_data, fw_data$rep_prob == 1)
hc_rep_nr    = glmer (rep_nr~ 1 + (1|transplant/site)  , family= poisson, data = df_rep_nr_hc, na.action = na.omit)
hw_rep_nr    = glmer (rep_nr~ z + (1|transplant/site)  , family= poisson, data = df_rep_nr_hw, na.action = na.omit)
#get fitted values
fitted_val_hw_rep_nr = as.data.frame(Effect("z", hw_rep_nr, xlevels = length(df_rep_nr_hw$z)))
fitted_val_hw_rep_nr$rep_nr = df_rep_nr_hw$rep_nr
fitted_val_hw_rep_nr$treatment = "T"

hedgerow_rep = rbind(df_rep_nr_hc, df_rep_nr_hw)

rep_nr_h = ggplot(hedgerow_rep[!is.na(hedgerow_rep$z),], aes(z, rep_nr, colour = treatment)) + 
  geom_point(show.legend = FALSE) +
  geom_ribbon(data = fitted_val_hw_rep_nr, aes(ymin=lower, ymax=upper), colour = NA, fill = "rosybrown2", alpha=0.5) +
  geom_line(y=summary(hc_rep_nr)$coef [1], size = 1.25, colour = "black") + #line of null model (treatment = C)
  geom_line(data = fitted_val_hw_rep_nr, aes(z, fit), size = 1.25,show.legend = FALSE) +
  scale_color_manual(values = c("black", "red"))+
  theme(panel.background = element_rect(colour="black", fill= "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank (),
        axis.title = element_text(size=13),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank()) +
  coord_cartesian(xlim = c(min(hedgerow$z, na.rm=T), max(hedgerow$z, na.rm= T)), ylim= c(min(data$rep_nr, na.rm = T),max(data$rep_nr, na.rm = T))) 

### Seed number
#forests
df_seed_nr_fc = subset (fc_data, fc_data$seed_prod_2019 == 1)
df_seed_nr_fw = subset (fw_data, fw_data$seed_prod_2019 == 1)
fc_seed_nr   = glmer (seed_nr_2019 ~ z1 + (1|transplant/site), family= poisson, data = df_seed_nr_fc, na.action = na.omit)
fw_seed_nr   = glmer (seed_nr_2019 ~ 1 + (1|transplant/site), family= poisson, data = df_seed_nr_fw, na.action = na.omit)

#get fitted values
fitted_val_fc_seed = as.data.frame(Effect("z1", fc_seed_nr, xlevels = length(df_seed_nr_fc$z1)))
fitted_val_fc_seed$seed_nr_2019 = df_seed_nr_fc$seed_nr_2019
fitted_val_fc_seed$treatment = "C"

forest_seed = rbind(df_seed_nr_fc, df_seed_nr_fw)

seed_nr_f = ggplot(forest_seed[!is.na(forest_seed$z1),], aes(z1, seed_nr_2019, colour = treatment)) + 
  geom_point(show.legend = FALSE) +
  geom_ribbon(data = fitted_val_fc_seed, aes(ymin=lower, ymax=upper), colour= NA, fill = "grey", alpha=0.5) + #fc
  geom_line(y=summary(fw_seed_nr)$coef [1], size = 1.25, colour = "red",show.legend = FALSE) + #line of null model (treatment = T)
  geom_line(data = fitted_val_fc_seed, aes(z1, fit), size = 1.25,show.legend = FALSE) +
  scale_color_manual(values = c("black", "red"))+
  ylab("Seed number") + 
  xlab("Size (cm)") + 
  theme(panel.background = element_rect(colour="black", fill= "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank (),
        axis.title = element_text(size=13)) +
  coord_cartesian(xlim = c(min(forest$z, na.rm=T), max(forest$z, na.rm= T)), ylim= c(min(data$seed_nr_2019, na.rm = T),max(data$seed_nr_2019, na.rm = T))) 

#hedgerows
df_seed_nr_hc = subset (hc_data, hc_data$seed_prod_2019 == 1)
df_seed_nr_hw = subset (hw_data, hw_data$seed_prod_2019 == 1)
hc_seed_nr   = glmer (seed_nr_2019 ~ z1 + I(z1^2) + (1|transplant/site), family= poisson, data = df_seed_nr_hc, na.action = na.omit)
hw_seed_nr   = glmer (seed_nr_2019 ~ z1 + (1|transplant/site), family= poisson, data = df_seed_nr_hw, na.action = na.omit)

#get fitted values
fitted_val_hc_seed = as.data.frame(Effect("z1", hc_seed_nr, xlevels = length(df_seed_nr_hc$z1)))
fitted_val_hc_seed$seed_nr_2019 = df_seed_nr_hc$seed_nr_2019
fitted_val_hc_seed$treatment = "C"

fitted_val_hw_seed = as.data.frame(Effect("z1", hw_seed_nr, xlevels = length(df_seed_nr_hw$z1)))
fitted_val_hw_seed$seed_nr_2019 = df_seed_nr_hw$seed_nr_2019
fitted_val_hw_seed$treatment = "T"

fitted_val_h_seed = rbind(fitted_val_hc_seed, fitted_val_hw_seed)

hedgerow_seed = rbind(df_seed_nr_hc, df_seed_nr_hw)

seed_nr_h = ggplot(hedgerow_seed[!is.na(hedgerow_seed$z1),], aes(z1, seed_nr_2019, colour = treatment)) + 
  geom_point(show.legend = FALSE) +
  geom_ribbon(data = fitted_val_hc_seed, aes(ymin=lower, ymax=upper), colour= NA, fill = "grey", alpha=0.5) + #hc
  geom_ribbon(data = fitted_val_hw_seed, aes(ymin=lower, ymax=upper), colour= NA, fill = "rosybrown2", alpha=0.5) + #hw
  geom_line(data = fitted_val_hc_seed, aes(z1, fit), size = 1.25,show.legend = FALSE) + #hc line
  geom_line(data = fitted_val_hw_seed, aes(z1, fit), size = 1.25,show.legend = FALSE) + #hw line
  scale_color_manual(values = c("black", "red"))+
  xlab("Size (cm)") + 
  theme(panel.background = element_rect(colour="black", fill= "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank (),
        axis.title = element_text(size=13),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank()) +
  coord_cartesian(xlim = c(min(hedgerow$z, na.rm=T), max(hedgerow$z, na.rm= T)), ylim= c(min(data$seed_nr_2019, na.rm = T),max(data$seed_nr_2019, na.rm = T))) 



##grid.arrange(survival_f, survival_h,rep_prob_f, rep_prob_h, seed_prob_f, seed_prob_h,growth_f, growth_h,rep_nr_f, rep_nr_h, seed_nr_f, seed_nr_h, ncol = 2, nrow = 6 )

## Align plots and show them
plots <- list(survival_f, survival_h, 
              rep_prob_f, rep_prob_h, 
              #seed_prob_f, seed_prob_h,
              growth_f, growth_h,
              rep_nr_f, rep_nr_h,
              seed_nr_f, seed_nr_h)
grobs <- list()
widths <- list()

#collect the widths for each grob of each plot
for (i in 1:length(plots)){
  grobs[[i]] <- ggplotGrob(plots[[i]])
  widths[[i]] <- grobs[[i]]$widths[2:5]
}
#get the max width
maxwidth <- do.call(grid::unit.pmax, widths)
#assign the max width to each plot
for (i in 1:length(grobs)){
  grobs[[i]]$widths[2:5] <- as.list(maxwidth)
}
#Generate the plots
##do.call("grid.arrange", c(grobs, ncol = 2))

surv_f = grobs[[1]]
surv_h = grobs[[2]]
rep_pr_f = grobs[[3]]
rep_pr_h = grobs[[4]]
growth_rate_f = grobs[[5]]
growth_rate_h = grobs[[6]]
rep_num_f = grobs [[7]]
rep_num_h = grobs [[8]]
seed_num_f = grobs [[9]]
seed_num_h = grobs [[10]]

grid.arrange(
  arrangeGrob(surv_f, left= textGrob(expression(bold("a)")), gp = gpar(fontsize = 13, col= "black"), rot= 360, y=0.9), top = textGrob(expression(bold("Forests")), gp = gpar(fontsize = 13, col= "black"))),
  arrangeGrob(surv_h, top= textGrob(expression(bold("Hedeerows")), gp = gpar(fontsize = 13, col= "black"))),
  
  arrangeGrob(rep_pr_f, left= textGrob(expression(bold("b)")), gp = gpar(fontsize = 13, col= "black"), rot= 360, y=0.9)),
  arrangeGrob(rep_pr_h),
  
  #arrangeGrob(seed_rate_f, left = textGrob(expression(bold("c)")), gp = gpar(fontsize = 13, col= "black"), rot= 360, y=0.9)),
  #arrangeGrob(seed_rate_h),
  
  arrangeGrob(growth_rate_f, left= textGrob(expression(bold("c)")), gp = gpar(fontsize = 13, col= "black"), rot= 360, y=0.9)),
  arrangeGrob(growth_rate_h),
  
  arrangeGrob(rep_num_f, left= textGrob(expression(bold("d)")), gp = gpar(fontsize = 13, col= "black"), rot= 360, y=0.9)),
  arrangeGrob(rep_num_h),
  
  arrangeGrob(seed_num_f, left= textGrob(expression(bold("e)")), gp = gpar(fontsize = 13, col= "black"), rot= 360, y=0.9)),
  arrangeGrob(seed_num_h),
  
  ncol= 2 )




