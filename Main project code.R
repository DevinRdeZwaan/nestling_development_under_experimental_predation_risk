##########################################################################################
### R code for:																			
###	"Hierarchical fear: parental behaviour and corticosterone mediate nestling growth 
###  in response to predation risk"

###	DR de Zwaan and K Martin
###	April 2020	

### Run with R version 3.6.3
##########################################################################################  
  
  
  
### Required packages
library(piecewiseSEM)
library(lme4)
library(lmerTest)
library(multcompView)
library(car)
library(ggplot2)
library(MuMIn)
library(doBy)
library(plyr)
library(simsem)
library(lavaan)
library(ggplot2)
library(visreg)
library(dplyr)

### Set working directory

setwd("C:/Users/Devin/Dropbox/Working manuscripts/Chapter 2")

### Read in file

nestlings <- read.csv(file.choose(), stringsAsFactors = FALSE,
                      strip.white = TRUE, na.strings = "") # load nest summary combined

### Check variables
head(nestlings)
str(nestlings)

#####################################
### Create 2 datasets without NAs: A) only nestlings with CORT data, B) all nestlings
### with complete data: (A) used for path analysis, (B) used for all other analysis 

nestlings_cort <- subset(nestlings, CORT_data > 0)

nestlings_complete <- subset(nestlings, complete_data >0)


#################################################################################################
################## Summary calculations for first paragraph of Results ##########################
#################################################################################################

### Sample sizes

table(nestlings_complete$complete_data) # sample size = 188
table(nestlings_cort$CORT_data)    # sample size = 112

length(table(nestlings_complete$nest_id)) # 54 nests
length(table(nestlings_cort$nest_id)) # 56 nests


### Calculate avg growth at 7 days relative to adult size
### Average adult wing = 105, tarsus = 22.35, mass = 34.09 (de Zwaan, unpublished data)

mean(nestlings_complete$wing_7d)/105 # 37.4%

mean(nestlings_complete$tarsus_7d)/22.35 # 88.6%

mean(nestlings_complete$mass_7d)/34.09 # 60.0%


### Calculate average growth between 5 and 7 days post-hatch

nestlings_complete$wing_growth <- nestlings_complete$wing_7d - nestlings_complete$wing_5d

nestlings_complete$tarsus_growth <- nestlings_complete$tarsus_7d - nestlings_complete$tarsus_5d

nestlings_complete$mass_growth <- nestlings_complete$mass_7d - nestlings_complete$mass_5d

### Average and SE of growth for all traits

mean(nestlings_complete$wing_growth) #13.8
sd(nestlings_complete$wing_growth)/sqrt(length(nestlings_complete$wing_growth)) # 0.2

mean(nestlings_complete$tarsus_growth) # 3.31
sd(nestlings_complete$tarsus_growth)/sqrt(length(nestlings_complete$tarsus_growth)) # 0.08

mean(nestlings_complete$mass_growth) # 3.99
sd(nestlings_complete$mass_growth)/sqrt(length(nestlings_complete$mass_growth)) # 0.19


### Average growth between 5-7 days relative to growth over all 7 days

mean(nestlings_complete$wing_growth)/mean(nestlings_complete$wing_7d) # 35.3%

mean(nestlings_complete$tarsus_growth)/mean(nestlings_complete$tarsus_7d) # 16.7%

mean(nestlings_complete$mass_growth)/mean(nestlings_complete$mass_7d) # 19.5%


#################################################################################################
### Relationship between nestling size and age at fledge - second Results paragraph #############
#################################################################################################

### Turn age at fledge into a factor

nestlings_complete$age_at_fledge_f <- factor(nestlings_complete$age_at_fledge, levels = c("8", "9", "10", "11"))

### Wing model

AAF_wing_lmer <- lmer(wing_7d ~ age_at_fledge_f + brood_size + clutch_initiation 
                      + (1|nest_id), data=nestlings_complete)

Anova(AAF_wing_lmer, type=3, test="F") # F = 3.23, P = 0.03

### Tarsus model 

AAF_tarsus_lmer <- lmer(tarsus_7d ~ age_at_fledge_f + brood_size + clutch_initiation 
                        + (1|nest_id), data=nestlings_complete)

Anova(AAF_tarsus_lmer, type=3, test="F") # F = 2.31, P = 0.09

### Mass model

AAF_mass_lmer <- lmer(mass_7d ~ age_at_fledge_f + brood_size + clutch_initiation 
                      + (1|nest_id), data=nestlings_complete)

Anova(AAF_mass_lmer, type=3, test="F") # F = 1.31, P = 0.28


### Wing load

nestlings_complete$wing_load <- summary(lm(wing_7d~ mass_7d, data=nestlings_complete))$residuals

AAF_wing_load_lmer <- lmer(wing_load ~ age_at_fledge_f + brood_size + clutch_initiation
                           + (1|nest_id), data=nestlings_complete)

Anova(AAF_wing_load_lmer, type=3, test="F") # F = 2.88, P = 0.04


### Correlation between wing and tarsus

cor.test(nestlings_complete$wing_7d, nestlings_complete$tarsus_7d) # r = 0.86



###################################
### Average, min, max, and SE of age at fledge
###################################

### Average within nest

within_nest_summary <-ddply(nestlings_complete, .(nest_id), summarize, 
                            age_at_fledge = mean(age_at_fledge))


### Remove NAs
within_nest_subset <- subset(within_nest_summary, !is.na(age_at_fledge))

aaf_summary <- summaryBy(age_at_fledge ~ 1, data = within_nest_subset, FUN = c(mean,min,max,sd,length))
aaf_summary$se <- aaf_summary$age_at_fledge.sd / sqrt(aaf_summary$age_at_fledge.length)
aaf_summary

### mean = 9.2, min = 8, max = 11, SE = 0.2


#################################################################################################
######################################## Path analysis ##########################################
#################################################################################################

#####################
### Calculate variables for path analysis
#####################

### Change in provisioning / 10 min
nestlings_cort$prov_dif <- (nestlings_cort$prov_rate_control - nestlings_cort$prov_rate_treatment)*10

### Make treatment into a factor & set SAVS as the baseline
nestlings_cort$treatment_f <- relevel(factor(nestlings_cort$treatment, ordered=FALSE), ref="SAVS")


##########################
### Important note - prior to running path models, you need to make a small update to 
### an internal function in piecewiseSEM to avoid a 'missing column' error.
### To do this, load and run the entire R script called 'piecewiseSEM code fix'
#########################


##############################
#### Two hypothesized paths - A) causal, and B) correlational
##############################


### A) causal model

wing_causal <- psem(lmer(wing_5d ~ brood_size + clutch_initiation + temp_sub_10_pre + CORT_pg_per_mm_pre + (1|nest_id), data=nestlings_cort),
                  lmer(wing_7d ~ wing_5d + CORT_pg_per_mm_post*clutch_initiation + prov_dif + (1|nest_id), data=nestlings_cort),
                  lm(prov_dif ~ brood_size + treatment_f, data=nestlings_cort),
                  lmer(CORT_pg_per_mm_post ~ CORT_pg_per_mm_pre + prov_dif+ treatment_f + clutch_initiation + (1|nest_id), data=nestlings_cort))


summary(wing_causal) # Fishers' C = 13.28 on 24 df; AIC = 69.28


# correlational model

wing_cor <- psem(lmer(wing_5d ~ brood_size + clutch_initiation + temp_sub_10_pre + CORT_pg_per_mm_pre +(1|nest_id), data=nestlings_cort),
                  lmer(wing_7d ~ wing_5d + prov_dif + treatment_f + (1|nest_id), data=nestlings_cort),
                  lm(prov_dif ~ brood_size + treatment_f, data=nestlings_cort),
                  lmer(CORT_pg_per_mm_post ~ CORT_pg_per_mm_pre + prov_dif + treatment_f + clutch_initiation + (1|nest_id), data=nestlings_cort),
                  wing_7d %~~%  CORT_pg_per_mm_post) 

summary(wing_cor) # Fisher's C = 20.11 on 24 DF; AIC = 74.11


#####
### Conclusion - retain causal model based on AIC and Fisher's C

### See Supplementary Appendix for the path structure of the causal and correlational model
#####


#####
### Model selection procedure to pare down global causal model to final model
### can be found in separate R code named 'Path model selection'.
### The final model after selection is included below.
#####


#####
### Final model used in manuscript
#####

wing_final <- psem(lmer(wing_5d ~ clutch_initiation + temp_sub_10_pre + CORT_pg_per_mm_pre + (1|nest_id), data=nestlings_cort),
                  lmer(wing_7d ~ wing_5d + CORT_pg_per_mm_post*clutch_initiation + (1|nest_id), data=nestlings_cort),
                  lm(prov_dif ~ treatment_f, data=nestlings_cort),
                  lmer(CORT_pg_per_mm_post ~ CORT_pg_per_mm_pre + prov_dif + treatment_f + clutch_initiation + (1|nest_id), data=nestlings_cort))

summary(wing_final) # Fisher's C = 13.65; AIC = 63.65


###################
### Path model power validation
###################

# remove NAs from dataset

nestlings_subset <- subset(nestlings_cort, !is.na(prov_dif))


# create theoretical population of 10,000 individuals based on existing distributions

expanded_pop_cort <- data.frame(wing_5d = rnorm(10000, mean(nestlings_cort$wing_5d), sd(nestlings_cort$wing_5d)),
                                wing_7d = rnorm(10000, mean(nestlings_cort$wing_7d), sd(nestlings_cort$wing_7d)),
                                CORT_pg_per_mm_pre = rnorm(10000, mean(nestlings_cort$CORT_pg_per_mm_pre), sd(nestlings_cort$CORT_pg_per_mm_pre)),
                                CORT_pg_per_mm_post = rnorm(10000, mean(nestlings_cort$CORT_pg_per_mm_post), sd(nestlings_cort$CORT_pg_per_mm_post)),
                                clutch_initiation = rnorm(10000, mean(nestlings_cort$clutch_initiation), sd(nestlings_cort$clutch_initiation)),
                                temp_sub_10_pre = rnorm(10000, mean(nestlings_cort$temp_sub_10_pre), sd(nestlings_cort$temp_sub_10_pre)),
                                treatment_dummy = sample(c(0, 1, 2), size = 10000, replace = TRUE, prob = c(0.39, 0.32, 0.29)),
                                # create dummy variable (SAVS =0, Fox =1, Raven = 2) - prob based on proportion in dataset
                                prov_dif = rnorm(10000, mean(nestlings_subset$prov_dif), sd(nestlings_subset$prov_dif)))

### Check
head(expanded_pop_cort)


### Create model formula to match final path model.
### Note - nest.id as a random effect is not necessary because this is a theoretical population
### where individual nestlings are not associated with each other in the same way as a natural
### population (i.e., nestlings are not grouped within nests)

nestling.mcmc <- "
wing_5d ~ temp_sub_10_pre + CORT_pg_per_mm_pre + clutch_initiation
wing_7d ~ wing_5d +  CORT_pg_per_mm_post*clutch_initiation
CORT_pg_per_mm_post ~ CORT_pg_per_mm_pre + treatment_dummy + prov_dif + clutch_initiation
prov_dif ~ treatment_dummy
"

### Run model 1000 times and restrict sample size to observed data (n = 112)
output.f <- sim(1000, model=nestling.mcmc, n=112, rawData=expanded_pop_cort, lavaanfun = "cfa")

plotCutoff(output.f, 0.05) # graphically shows how many iterations were outside the norm
# Most iterations before the cut-off (red line)


### Create the same CFA model drawing from original data
### First make dummy variable in original data set corresponding to treatment levels

nestlings_cort$treatment_dummy <- ifelse(nestlings_cort$treatment_f == "SAVS", 0,
                                  ifelse(nestlings_cort$treatment_f == "Fox", 1,2))

### CFA model
nestling_fit <- cfa(nestling.mcmc, data = nestlings_cort)

### Compare simulated model with observed model to evaluate sample size adequacy
pValue(nestling_fit, output.f) 
# AIC = 1, Chisq = 0.36, RMSEA = 0.35
# Note: Because distribution is randomly generated, these test statistics will vary slightly
# among attempts.

### Conclusion: the two models are not significantly different (> 0.05) looking at AIC, Chisq, 
### and RMSEA, indicating a robust fit based on our sample size


#####################################
### Extract standardized coefficients and SEs from final path model
######################################

### A) for numerical predictors

coefs <-coefs(wing_final, standardize="scale") # extract coefficients

coefs$Std.Estimate <- as.numeric(coefs$Std.Estimate) # need to make numerical
coefs$Std.Error <- as.numeric(coefs$Std.Error)
coefs$Estimate <- as.numeric(coefs$Estimate)

### Note: Warning about "NAs introduced by coercion" refers to the factor levels.
### Next step will remove these.

### Remove the rows that do not have standard estimates
coefs
coefs.updated <- coefs[-c(8:11,15:18),]
str(coefs.updated)

### SE of coefs are not scaled, so do so now

coefs.updated$scaled.se <- (coefs.updated$Std.Error * coefs.updated$Std.Estimate)/coefs.updated$Estimate # calculate scaled standard errors

### coefs and SE for numerical predictors
coefs.table <- coefs.updated[c("Response", "Predictor", "Std.Estimate", "scaled.se", "P.Value")]
coefs.table

### B) for categorical predictors (treatment)

### Note - path analysis output gives marginal means for each treatment level.
### To calculate standardized coefficients so comparable to numerical predictors:
### 1) Subtract marginal mean from reference mean (i.e., Fox - SAVS) to get coefficient
### 2) Calculate pair-wise SD using Hedge's G formula
### 3) Divide calculated coefficient by pair-wise SD for standradized coefficient
### 4) Divide standard error by pair-wise SD for standardized SE


### Provisioning rate ~ predation risk

nestlings_subset <- subset(nestlings_cort, !is.na(prov_dif)) # remove NAs from dataset

summaryBy(prov_dif ~ treatment_f, data=nestlings_subset, FUN=c(length,sd)) # calculate SD & sample size

# extract marginal means and SE from path analysis output
marginal_means <- coefs[c(8:11,15:18), c('Response', 'Predictor', 'Estimate', 'Std.Error')]
marginal_means

## Fox - Hedge's G formula
sqrt((33*(0.8417366)^2 + 33*(0.8216543)^2)/(34+34-2)) #0.83

## Fox - standardized B coefficient
(0.5735 - -0.3941) / 0.83 # 1.17

## Fox - standardized SE
0.1433/0.83 #0.17

## Raven - Hedge's G formula
sqrt((27*(0.8440266)^2 + 33*(0.8216543)^2)/(28+34-2)) # 0.83

## Raven - standardized B coefficient
(0.8357 - -0.3941) / 0.83 # 1.48

## Raven - standardized SE
0.1579/0.83 # 0.19



### Corticosterone ~ predation risk

summaryBy(CORT_pg_per_mm_post ~ treatment_f, data=nestlings_cort, FUN=c(length,sd)) # calculate SD & sample size

## Fox - Hedge's G formula
sqrt((35*(5.155177)^2 + 43*(4.670002)^2)/(36+44-2)) #4.89

## Fox - standardized B coefficient
(11.8715 - 12.9681) / 4.89 # -0.22

## Fox - standardized SE
0.8154/4.89 #0.17

## Raven - Hedge's G formula
sqrt((31*(3.801787)^2 + 43*(4.670002)^2)/(32+44-2)) # 4.33

## Raven - standardized B coefficient
(10.3826 - 12.9681) / 4.33 # -0.60

## Raven - standardized SE
0.9165/4.33 # 0.21

### The above coefficients and SE can be found in Fig 3 in the manuscript




#################################################################################################
######################################### Figures ###############################################
#################################################################################################

#################
#### Figure 2 - boxplots showing relationship between age at fledge and size at 7 d
#################

### Clean dataset to use for plotting 

# Remove age at fledge NAs
nestlings_complete_plot <- subset(nestlings_complete, !is.na(age_at_fledge))

# Convert to factor
nestlings_complete_plot$age_at_fledge_f <- factor(nestlings_complete_plot$age_at_fledge, 
                                                  levels = c("8", "9", "10", "11"))

### Plot

### Panel A - wing length

wing_plot <- ggplot(nestlings_complete_plot, aes(x=age_at_fledge_f, y=wing_7d, group=age_at_fledge_f)) + 
  geom_boxplot(fill="lightskyblue3", color="black", lwd=0.5, outlier.shape = NA)  +
  geom_jitter(size = 1.5, shape=19, color = "black", position=position_jitter(0.2)) + # because of jitter, points will not be in exactly the same position
  scale_y_continuous(limits=c(18,60), expand=c(0.05,0.05))+ 
  labs(y="Wing at 7 days (mm)", x = "Age at fledge") +
  theme_classic()+  
  theme(axis.line = element_line(colour="black", size=1.3),   
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.border = element_blank(), 
        axis.ticks.length = unit(0.2, "cm"),
        axis.title.x = element_text(size=22, colour="black", vjust= 1.5), 
        axis.text.x  = element_text(size=18, colour="black"), 
        axis.title.y=element_text(size=22, colour="black", vjust=1.8), 
        axis.text.y=element_text(size=18, colour="black")) 

wing_plot

# Export
png(file="fig2_panelA_wing_length.png",width=3000,height=2500, res=600)
wing_plot
dev.off()



### Panel B - wing load

wing_load <- ggplot(nestlings_complete_plot, aes(x=age_at_fledge_f, y=wing_load, group=age_at_fledge_f)) + 
  geom_boxplot(fill="lightskyblue3", color="black", lwd=0.5, outlier.shape = NA)  +
  geom_jitter(size = 1.5, shape=19, position=position_jitter(0.2)) +
  scale_y_continuous(limits=c(-15,13), breaks= c(-15,-10,-5,0,5,10), expand=c(0.05,0.05))+
  labs(y="Wing load at 7 days", x = "Age at fledge (days)") +
  theme_classic()+  
  theme(axis.line = element_line(colour="black", size=1.3),  
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.border = element_blank(), 
        axis.ticks.length = unit(0.2, "cm"), 
        axis.title.x = element_text(size=22, colour="black"),
        axis.text.x  = element_text(size=18, colour="black"), 
        axis.title.y=element_text(size=22, colour="black", vjust=1.5),
        axis.text.y=element_text(size=18, colour="black")) 

wing_load

# Export
png(file="fig2_panelB_wing_load.png",width=3000,height=2500, res=600)
wing_load
dev.off()


### Panel C - tarsus length

tarsus_plot <- ggplot(nestlings_complete_plot, aes(x=age_at_fledge_f, y=tarsus_7d, group=age_at_fledge_f)) + 
  geom_boxplot(fill="lightskyblue3", color="black", lwd=0.5, outlier.shape = NA)  +
  geom_jitter(size = 1.5, shape=19, position=position_jitter(0.2)) +
  scale_y_continuous(limits=c(10,30), breaks=c(10,15,20,25,30), expand=c(0.05,0.05))+
  labs(y="Tarsus at 7 days (mm)", x = "Age at fledge") +
  theme_classic()+  
  theme(axis.line = element_line(colour="black", size=1.3),  
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        axis.ticks.length = unit(0.2, "cm"),
        axis.title.x = element_text(size=22, colour="black"), 
        axis.text.x  = element_text(size=18, colour="black"), 
        axis.title.y=element_text(size=22, colour="black", vjust=1.8),
        axis.text.y=element_text(size=18, colour="black"))

tarsus_plot

png(file="fig2_panelC_tarsus_length.png",width=3000,height=2500, res=600)
tarsus_plot
dev.off()


### Panel D - mass

mass_plot <- ggplot(nestlings_complete_plot, aes(x=age_at_fledge_f, y=mass_7d, group=age_at_fledge_f)) + 
  geom_boxplot(fill="lightskyblue3", color="black", lwd=0.5, outlier.shape = NA)  +
  geom_jitter(size = 1.5, shape=19, position=position_jitter(0.2)) +
  scale_y_continuous(limits=c(8,30), expand=c(0.05,0.05))+
  labs(y="Mass at 7 days (g)", x = "Age at fledge") +
  theme_classic()+  
  theme(axis.line = element_line(colour="black", size=1.3),  
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), 
        panel.border = element_blank(), 
        axis.ticks.length = unit(0.2, "cm"), 
        axis.title.x = element_text(size=22, colour="black"), 
        axis.text.x  = element_text(size=18, colour="black"), 
        axis.title.y=element_text(size=22, colour="black", vjust=1.8), 
        axis.text.y=element_text(size=18, colour="black")) 

mass_plot

png(file="fig2_panelD_mass.png",width=3000,height=2500, res=600)
mass_plot
dev.off()


##############################
## Figure 4. Relationship between DFE*CORT and wing lengh figure
##############################

### Look at spread of CORT data

hist(nestlings_cort$CORT_pg_per_mm_post, 20)

### Calculate 1 sd below and 2 sd above median to capture spread.
### Representative of nestlings with low versus high CORT concentration.
### Equivalent to 6 pg/mm and 20 pg/mm.

sd_below <- round(median(nestlings_cort$CORT_pg_per_mm_post) - sd(nestlings_cort$CORT_pg_per_mm_post))
sd_above <- round(median(nestlings_cort$CORT_pg_per_mm_post) + (2 * sd(nestlings_cort$CORT_pg_per_mm_post)))

### Run sub-model from path analysis
wing_cort <- lmer(wing_7d ~ wing_5d + CORT_pg_per_mm_post*clutch_initiation + (1|nest_id), data=nestlings_cort)


### Create colours
colour.1 <- rgb(216/255,179/255,101/255, alpha=0.5)
colour.2 <- rgb(90/255,180/255,172/255, alpha=0.5)


png(file="Fig4.png", width=3000,height=2500, res=600)
visreg(wing_cort, "clutch_initiation", by="CORT_pg_per_mm_post", type="conditional", overlay=TRUE, 
       breaks=c(sd_below, sd_above), rug=FALSE, yaxt="n", xaxt="n",
       line.par = list(col = c("black", "black"), lty=c("solid", "dashed")),
       fill = list(col = c(colour.1, colour.2)),
       points=list(cex=1.2, pch=c(16,1), col="black"),
       xlab = list("Julian date", cex=1.5), ylab = list("Wing length (mm)", cex=1.5))
axis(1, cex.axis=1.18)
axis(2, cex.axis=1.18)
dev.off()



##########################
## Figure 5 - influence of predation risk on nestling wing and DTF
##########################

# Note:
# Values were calculated by hand for each of Fox and Raven using the SD of the prov_dif variable
# and multiplying it through to calculate the SD change in wing length (indirect).
# We kept the influence of predator on CORT fixed but varied the influence of CORT on wing length 
# by time of year based on the interaction (three times - early (-1sd), average clutch initiation,
# and late (+1sd) (direct effects).
# The predicted points are displayed in increments of 1 sd change in provisioning rate based
# on the observed range of values for each treatment (~ -1 visit/10 min increments)

####
# First, create table for nestling growth in response to fox
####

prov_response_f <- factor(c('0','-1','-2','-3','0','-1','-2','-3','0','-1','-2','-3'), 
                          ordered =TRUE, levels = c('0', '-1', '-2', '-3'))

wing_growth_f <-c(2.74,-0.25,-3.24,-6.23,1.31,-0.12,-1.56,-2.93,0.12,0,-0.12,-0.25) # early, avg, late

clutch_initiation_f <- factor(c("early", "early", "early", "early", "avg", "avg","avg","avg", 
                         "late","late","late","late"), ordered = TRUE, levels = c('early', 'avg', 'late'))

fox_table <- bind_cols(list(prov_response_f, wing_growth_f, clutch_initiation_f))
fox_table <- as.data.frame(fox_table)

colnames(fox_table) <- c('prov_response', 'wing_growth', 'clutch_initiation')


### Second, create raven table

prov_response_r <- factor(c('0','-1','-2','-3','0','-1','-2','-3','0','-1','-2','-3'), 
                          ordered =TRUE, levels = c('0', '-1', '-2', '-3'))

wing_growth_r <- as.vector(c(8.12,4.87,1.62,-1.62,3.86,2.30,0.74,-0.74,0.34,0.20,0.07,-0.07)) # early, avg, late

clutch_initiation_r <- factor(c("early", "early", "early", "early", "avg", "avg","avg","avg", 
                          "late","late","late","late"), ordered = TRUE, levels = c('early', 'avg', 'late'))

raven_table <- bind_cols(list(prov_response_r, wing_growth_r, clutch_initiation_r))
raven_table <- as.data.frame(raven_table)

colnames(raven_table) <- c('prov_response', 'wing_growth', 'clutch_initiation')


#########
## Fig 5A - Fox
#########

dodge <- position_dodge(0.08) #jitter points slightly

fox_wing_plot <- ggplot(fox_table, aes(prov_response, wing_growth, group=clutch_initiation)) +
  geom_hline(yintercept=0, linetype="dashed", color = "grey65", size=1.5) +
  geom_line(position=dodge, colour="black", size=0.8) +
  geom_point(aes(fill=clutch_initiation), shape=21, color = "black", stroke = 1, size=6, position=dodge)+
  scale_fill_manual(values=c('#e5f5f9','#99d8c9', '#2ca25f')) +
  scale_y_continuous(limits=c(-7,3), expand=c(0.05,0.05))+
  labs(y="Wing growth (mm)", x = "Relative provisioning response (nest visits/10 min)",
       fill = "") + 
  theme_classic()+ 
  theme(axis.line = element_line(colour="black", size=1.3),   
        panel.border = element_blank(), 
        legend.text=element_text(size=16),
        axis.ticks.length = unit(0.2, "cm"),
        axis.title.x = element_text(size=22, colour="black"), 
        axis.text.x  = element_text(size=18, colour="black"), 
        axis.title.y=element_text(size=22, colour="black", vjust=2), 
        axis.text.y=element_text(size=18, colour="black")) 

fox_wing_plot

# Export
png(file="Fig5A_fox.png",width=3000,height=2500, res=600)
fox_wing_plot
dev.off()


###################
## Fig 5B - raven plot 
####################

dodge <- position_dodge(0.08) #jitters points slightly

raven_wing_plot <- ggplot(raven_table, aes(prov_response, wing_growth, group=clutch_initiation)) +
  geom_hline(yintercept=0, linetype="dashed", color = "grey65", size=1.5) +
  geom_line(position=dodge, colour="black", size=0.8) +
  geom_point(aes(fill=clutch_initiation), shape=23, color = "black", stroke = 1, size=6, position=dodge)+
  scale_fill_manual(values=c('#ffeda0','#feb24c', '#f03b20')) +
  scale_y_continuous(limits=c(-2,9), expand=c(0.05,0.05))+
  labs(y="Wing growth (mm)", x = "Relative provisioning response (nest visits/10 min)",
       fill = "") + 
  theme_classic()+ 
  theme(axis.line = element_line(colour="black", size=1.3),   
        panel.border = element_blank(), 
        legend.text=element_text(size=16),
        axis.ticks.length = unit(0.2, "cm"),
        axis.title.x = element_text(size=22, colour="black"), 
        axis.text.x  = element_text(size=18, colour="black"), 
        axis.title.y=element_text(size=22, colour="black", vjust=2), 
        axis.text.y=element_text(size=18, colour="black"))

raven_wing_plot


#Export
png(file="Fig5B_raven.png",width=3000,height=2500, res=600)
raven_wing_plot
dev.off()


