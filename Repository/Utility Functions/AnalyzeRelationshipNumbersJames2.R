par(mfrow=c(3,2))
re_load <- FALSE
jrep <- 3
nbreaks <- 20
graph_hists <- FALSE
cum_only<- TRUE
lower_age <- 50
upper_age <- 52
lower_age_overlay <- 30
upper_age_overlay <- 32
if (graph_hists == TRUE) par(mfrow=c(4,4))

if (re_load == TRUE) {
  load("H:/TasP2/Exp5_Equal_Rel_Durs_Lower_Mean_Duration/txt_opto_Exp5_2paramsEqualRelDursN10000__under25_under30_random_prop1.RData")
  pl_low <- evomodel
  load("H:/TasP2/Exp1_Base10000_second_try/txt_opto_Exp1_3_paramsBase10000_more_reps__CD4_nadir_under500_random_prop1.RData")
  pl_default <- evomodel
  load("H:/TasP2/Exp51_Group1_Super_High_Risk/txt_opto_Exp51_paramsGroup1_SuperHighRisk__random_prop1.RData")
  pl_very_high <- evomodel
  load("H:/TasP2/Exp51_Group1_Super_High_Risk/txt_opto_Exp51_3_paramsGroup1_SuperHighRisk_v3__random_prop1.RData")
  pl_high <- evomodel
  load("H:/TasP2/Exp51_Group1_Super_High_Risk/txt_opto_Exp51_2_paramsGroup1_SuperHighRisk_v2__random_prop1.RData")
  pl_ultra <- evomodel
}

main_title <- "Baseline: D1 = 730, D2 = 3650 Days"
nbreaks = 15
curr_model <- pl_default
source('H:/TasP2/Exp51_Group1_Super_High_Risk/JamesPatner_list_Script2.R')
nbreaks = 8
source('H:/TasP2/Exp51_Group1_Super_High_Risk/JamesPatner_list_Script2_overlay.R')

main_title <- "Perturbation 2: D1 = D2 = 3*365 Days"
nbreaks = 20
curr_model <- pl_low 
source('H:/TasP2/Exp51_Group1_Super_High_Risk/JamesPatner_list_Script2.R')
nbreaks = 15
source('H:/TasP2/Exp51_Group1_Super_High_Risk/JamesPatner_list_Script2_overlay.R')

main_title <- "Perturbation 3: D1 = 10, D2= 6*365 Days,\n No transitions "
nbreaks = 20
curr_model <- pl_very_high
source('H:/TasP2/Exp51_Group1_Super_High_Risk/JamesPatner_list_Script2.R')
nbreaks = 15
source('H:/TasP2/Exp51_Group1_Super_High_Risk/JamesPatner_list_Script2_overlay.R')


main_title <- "Perturbation 4: D1 = 100, D2=3650 Days,\n42% enter into Group 1"
nbreaks = 15
curr_model <- pl_high
source('H:/TasP2/Exp51_Group1_Super_High_Risk/JamesPatner_list_Script2.R')
source('H:/TasP2/Exp51_Group1_Super_High_Risk/JamesPatner_list_Script2_overlay.R')


main_title <- "Perturbation 5: D1 = 10, D2=3650 Days,\n42% enter into Group 1"
nbreaks = 20
curr_model <- pl_ultra
source('H:/TasP2/Exp51_Group1_Super_High_Risk/JamesPatner_list_Script2.R')
source('H:/TasP2/Exp51_Group1_Super_High_Risk/JamesPatner_list_Script2_overlay.R')




