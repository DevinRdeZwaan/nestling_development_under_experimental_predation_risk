#################################################################################################
################### Model selection procedure to identify final best model ######################
#################################################################################################

#####
## Note - full details and model outputs can be found in the Supplementary Appendix.
## Covariates such as clutch initiation and temperature were retained as controls based on 
## previous research.
## The main predictor 'treatment' was also retained in all cases because we are testing
## the influence of predation risk on wing length and feather CORT

### Required packages

library(MuMIn)
library(piecewiseSEM)

#############
### 1) Start by removing terms from wing (5d) submodel
#############

### Determine order of removal

wing_5d_model <- lmer(wing_5d ~ brood_size + clutch_initiation + temp_sub_10_pre + CORT_pg_per_mm_pre + (1|nest_id), data=nestlings_cort)
wing_5d_wo_nestlings <- lmer(wing_5d ~ clutch_initiation + temp_sub_10_pre + CORT_pg_per_mm_pre + (1|nest_id), data=nestlings_cort)
wing_5d_wo_CORT <- lmer(wing_5d ~ brood_size + clutch_initiation + temp_sub_10_pre + (1|nest_id), data=nestlings_cort)

AICc(wing_5d_model, wing_5d_wo_nestlings, wing_5d_wo_CORT)
### Order of removal - nestlings, CORT


### Global path model
wing_causal <- psem(lmer(wing_5d ~ brood_size + clutch_initiation + temp_sub_10_pre + CORT_pg_per_mm_pre + (1|nest_id), data=nestlings_cort),
                    lmer(wing_7d ~ wing_5d + CORT_pg_per_mm_post*clutch_initiation + prov_dif + (1|nest_id), data=nestlings_cort),
                    lm(prov_dif ~ brood_size + treatment_f, data=nestlings_cort),
                    lmer(CORT_pg_per_mm_post ~ CORT_pg_per_mm_pre + prov_dif+ treatment_f + clutch_initiation + (1|nest_id), data=nestlings_cort))


summary(wing_causal) # Fishers' C = 13.28 on 24 df; AIC = 69.28


wing_causal_wo_brood_size <- psem(lmer(wing_5d ~ #brood_size + 
                            clutch_initiation + temp_sub_10_pre + CORT_pg_per_mm_pre + (1|nest_id), data=nestlings_cort),
                    lmer(wing_7d ~ wing_5d + CORT_pg_per_mm_post*clutch_initiation + prov_dif + (1|nest_id), data=nestlings_cort),
                    lm(prov_dif ~ brood_size + treatment_f, data=nestlings_cort),
                    lmer(CORT_pg_per_mm_post ~ CORT_pg_per_mm_pre + prov_dif+ treatment_f + clutch_initiation + (1|nest_id), data=nestlings_cort))


summary(wing_causal_wo_brood_size) # Fishers' C = 16.34; AIC = 70.34
### Remove nestlings


wing_causal_wo_cort <- psem(lmer(wing_5d ~ #brood_size + CORT_pg_per_mm_pre + 
                                clutch_initiation + temp_sub_10_pre + (1|nest_id), data=nestlings_cort),
                    lmer(wing_7d ~ wing_5d + CORT_pg_per_mm_post*clutch_initiation + prov_dif + (1|nest_id), data=nestlings_cort),
                    lm(prov_dif ~ brood_size + treatment_f, data=nestlings_cort),
                    lmer(CORT_pg_per_mm_post ~ CORT_pg_per_mm_pre + prov_dif+ treatment_f + clutch_initiation + (1|nest_id), data=nestlings_cort))


summary(wing_causal_wo_cort) # Fishers' C = 24.98; AIC = 76.98
### Keep CORT



###########
### 2) Now remove terms from provisioning rate submodel
###########

### No need to determine order of term removal because only testing one variable

### Best model so far
wing_causal_2 <- psem(lmer(wing_5d ~ clutch_initiation + temp_sub_10_pre + CORT_pg_per_mm_pre + (1|nest_id), data=nestlings_cort),
                                 lmer(wing_7d ~ wing_5d + CORT_pg_per_mm_post*clutch_initiation + prov_dif + (1|nest_id), data=nestlings_cort),
                                 lm(prov_dif ~ brood_size + treatment_f, data=nestlings_cort),
                                 lmer(CORT_pg_per_mm_post ~ CORT_pg_per_mm_pre + prov_dif+ treatment_f + clutch_initiation + (1|nest_id), data=nestlings_cort))


summary(wing_causal_2) # Fishers' C = 16.34; AIC = 70.34


wing_causal_2_wo_brood_size <- psem(lmer(wing_5d ~ clutch_initiation + temp_sub_10_pre + CORT_pg_per_mm_pre + (1|nest_id), data=nestlings_cort),
                      lmer(wing_7d ~ wing_5d + CORT_pg_per_mm_post*clutch_initiation + prov_dif + (1|nest_id), data=nestlings_cort),
                      lm(prov_dif ~ #brood_size + 
                           treatment_f, data=nestlings_cort),
                      lmer(CORT_pg_per_mm_post ~ CORT_pg_per_mm_pre + prov_dif+ treatment_f + clutch_initiation + (1|nest_id), data=nestlings_cort))


summary(wing_causal_2_wo_brood_size) # Fisher's C = 11.54; AIC = 63.54
### Remove nestlings from provisioning rate submodel


###############
### 3) Now remove terms from CORT (7 d) submodel
###############

### Determine order of term removal
CORT_7d_model <- lmer(CORT_pg_per_mm_post~ CORT_pg_per_mm_pre + prov_dif+ treatment_f + clutch_initiation + (1|nest_id), data=nestlings_cort)
CORT_7d_wo_CORT5 <- lmer(CORT_pg_per_mm_post~ prov_dif+ treatment_f + clutch_initiation + (1|nest_id), data=nestlings_cort)
CORT_7d_wo_prov <- lmer(CORT_pg_per_mm_post~ CORT_pg_per_mm_pre + treatment_f + clutch_initiation + (1|nest_id), data=nestlings_cort)


AICc(CORT_7d_model, CORT_7d_wo_CORT5, CORT_7d_wo_prov)
### Order of removal: CORT 5d, prov


### Best model so far
wing_causal_3 <- psem(lmer(wing_5d ~ clutch_initiation + temp_sub_10_pre + CORT_pg_per_mm_pre + (1|nest_id), data=nestlings_cort),
                      lmer(wing_7d ~ wing_5d + CORT_pg_per_mm_post*clutch_initiation + prov_dif + (1|nest_id), data=nestlings_cort),
                      lm(prov_dif ~ treatment_f, data=nestlings_cort),
                      lmer(CORT_pg_per_mm_post ~ CORT_pg_per_mm_pre + prov_dif+ treatment_f + clutch_initiation + (1|nest_id), data=nestlings_cort))

summary(wing_causal_3) # Fisher's C = 11.54; AIC = 63.54


wing_causal_3_wo_cort <- psem(lmer(wing_5d ~ clutch_initiation + temp_sub_10_pre + CORT_pg_per_mm_pre + (1|nest_id), data=nestlings_cort),
                      lmer(wing_7d ~ wing_5d + CORT_pg_per_mm_post*clutch_initiation + prov_dif + (1|nest_id), data=nestlings_cort),
                      lm(prov_dif ~ treatment_f, data=nestlings_cort),
                      lmer(CORT_pg_per_mm_post ~ #CORT_pg_per_mm_pre + 
                            prov_dif+ treatment_f + clutch_initiation + (1|nest_id), data=nestlings_cort))

summary(wing_causal_3_wo_cort) # Fisher's C = 15.76; AIC = 65.76
### Retain CORT


wing_causal_3_wo_prov <- psem(lmer(wing_5d ~ clutch_initiation + temp_sub_10_pre + CORT_pg_per_mm_pre + (1|nest_id), data=nestlings_cort),
                      lmer(wing_7d ~ wing_5d + CORT_pg_per_mm_post*clutch_initiation + prov_dif + (1|nest_id), data=nestlings_cort),
                      lm(prov_dif ~ treatment_f, data=nestlings_cort),
                      lmer(CORT_pg_per_mm_post ~ CORT_pg_per_mm_pre + #prov_dif + 
                            treatment_f + clutch_initiation + (1|nest_id), data=nestlings_cort))

summary(wing_causal_3_wo_prov) # Fisher's C = 23.96; AIC = 73.96
### Retain provisioning rate


###############
### 4) Now remove terms from wing (7 d) submodel
##############

### Determine order of term removal
wing_model <- lmer(wing_7d ~ wing_5d + CORT_pg_per_mm_post*clutch_initiation + prov_dif + (1|nest_id), data=nestlings_cort)
wing_model_wo_5d <- lmer(wing_7d ~ CORT_pg_per_mm_post*clutch_initiation + prov_dif + (1|nest_id), data=nestlings_cort)
wing_model_wo_CORT_int <- lmer(wing_7d ~ wing_5d + CORT_pg_per_mm_post + prov_dif + (1|nest_id), data=nestlings_cort)
wing_model_wo_CORT <- lmer(wing_7d ~ wing_5d + prov_dif + (1|nest_id), data=nestlings_cort)
wing_model_wo_prov <- lmer(wing_7d ~ wing_5d + CORT_pg_per_mm_post*clutch_initiation + prov_dif + (1|nest_id), data=nestlings_cort)

AICc(wing_model, wing_model_wo_5d, wing_model_wo_CORT_int, wing_model_wo_CORT, wing_model_wo_prov)
### Order of removal: int, cort, prov, wing 5d

### Best model so far
wing_causal_4 <- psem(lmer(wing_5d ~ clutch_initiation + temp_sub_10_pre + CORT_pg_per_mm_pre + (1|nest_id), data=nestlings_cort),
                      lmer(wing_7d ~ wing_5d + CORT_pg_per_mm_post*clutch_initiation + prov_dif + (1|nest_id), data=nestlings_cort),
                      lm(prov_dif ~ treatment_f, data=nestlings_cort),
                      lmer(CORT_pg_per_mm_post ~ CORT_pg_per_mm_pre + prov_dif+ treatment_f + clutch_initiation + (1|nest_id), data=nestlings_cort))

summary(wing_causal_4) # Fisher's C = 11.54; AIC = 63.54


wing_causal_4_wo_int <- psem(lmer(wing_5d ~ clutch_initiation + temp_sub_10_pre + CORT_pg_per_mm_pre + (1|nest_id), data=nestlings_cort),
                      lmer(wing_7d ~ wing_5d + CORT_pg_per_mm_post + prov_dif + (1|nest_id), data=nestlings_cort),
                      lm(prov_dif ~ treatment_f, data=nestlings_cort),
                      lmer(CORT_pg_per_mm_post ~ CORT_pg_per_mm_pre + prov_dif+ treatment_f + clutch_initiation + (1|nest_id), data=nestlings_cort))

summary(wing_causal_4_wo_int) # Fisher's C = 20.65; AIC = 68.65
### Retain interaction with clutch initiation date


wing_causal_4_wo_cort <- psem(lmer(wing_5d ~ clutch_initiation + temp_sub_10_pre + CORT_pg_per_mm_pre + (1|nest_id), data=nestlings_cort),
                             lmer(wing_7d ~ wing_5d + #CORT_pg_per_mm_post 
                                    prov_dif + (1|nest_id), data=nestlings_cort),
                             lm(prov_dif ~ treatment_f, data=nestlings_cort),
                             lmer(CORT_pg_per_mm_post ~ CORT_pg_per_mm_pre + prov_dif+ treatment_f + clutch_initiation + (1|nest_id), data=nestlings_cort))

summary(wing_causal_4_wo_cort) # Fisher's C = 26.64; AIC = 72.64
### Retain CORT


wing_causal_4_wo_prov <- psem(lmer(wing_5d ~ clutch_initiation + temp_sub_10_pre + CORT_pg_per_mm_pre + (1|nest_id), data=nestlings_cort),
                      lmer(wing_7d ~ wing_5d + CORT_pg_per_mm_post*clutch_initiation + #prov_dif +
                             (1|nest_id), data=nestlings_cort),
                      lm(prov_dif ~ treatment_f, data=nestlings_cort),
                      lmer(CORT_pg_per_mm_post ~ CORT_pg_per_mm_pre + prov_dif+ treatment_f + clutch_initiation + (1|nest_id), data=nestlings_cort))

summary(wing_causal_4_wo_prov) # Fisher's C = 13.65; AIC = 63.65
### Remove provisioning rate


wing_causal_4_wo_wing5 <- psem(lmer(wing_5d ~ clutch_initiation + temp_sub_10_pre + CORT_pg_per_mm_pre + (1|nest_id), data=nestlings_cort),
                              lmer(wing_7d ~ #wing_5d + 
                                     CORT_pg_per_mm_post*clutch_initiation + (1|nest_id), data=nestlings_cort),
                              lm(prov_dif ~ treatment_f, data=nestlings_cort),
                              lmer(CORT_pg_per_mm_post ~ CORT_pg_per_mm_pre + prov_dif+ treatment_f + clutch_initiation + (1|nest_id), data=nestlings_cort))

summary(wing_causal_4_wo_wing5) # Fisher's C = 168.51; AIC = 216.51
### Retain wing (5 d)


############
### Final selected model
############

wing_final <- psem(lmer(wing_5d ~ clutch_initiation + temp_sub_10_pre + CORT_pg_per_mm_pre + (1|nest_id), data=nestlings_cort),
                   lmer(wing_7d ~ wing_5d + CORT_pg_per_mm_post*clutch_initiation + (1|nest_id), data=nestlings_cort),
                   lm(prov_dif ~ treatment_f, data=nestlings_cort),
                   lmer(CORT_pg_per_mm_post ~ CORT_pg_per_mm_pre + prov_dif + treatment_f + clutch_initiation + (1|nest_id), data=nestlings_cort))

summary(wing_final) # Fisher's C = 13.65; AIC = 63.65


#####
### Model selection results can be found in the Supplementary Appendix
#####
