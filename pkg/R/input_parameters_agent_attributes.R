#' @title Title
#'
#' @description Description
#'
#' @param x A number.
#' @param y A number.
#' @return return value here.
#' @details
#' Additional details here
#' @examples
#' example function call here

#' @export
input_parameters_agent_attributes <-function(){
  
  #Description:
  # Vector of the individual agent attributes that populate the “pop” list
  
  list( popVariables= c(
    
#--  Virulence model parameters ------------------------------ #

    "s",                               "Status",           
    "NumRecipients",                    
    "ViralContribToLogSP0" ,
    "EnvirContribToLogSP0",            "LogSetPoint",           
    "SetPoint",                        "d_acute",  
    "Generation",                      "RandomTimeToAIDS",
    "Time_Inf",                        "Time_Inf_Adj",      "V",   
     "r0",                             "age_infection",
    "vl_phase2_trans",                 "rate_phase2",
    "PPP",                              "vl_peak_agent",
    "vl_at_test",                      "cd4_at_test",
    "vl_expected",
    "Donors_V",                        "Donors_treated",  
    "Donors_treated_2nd_line",         "Donors_CD4",
    "Donors_ViralContribToLogSP0",     "Donors_EnvirContribToLogSP0",
    "Donors_Total_Time_Inf_At_Trans",  "Donors_Generation",
    "Donors_Index",                    "Donors_age",                   
    "Donors_LogSetPoint",              "Donors_SetPoint",
    "Donors_d_acute",                  "Time_Death",
    "Donors_diag_status",              "Donors_age", 
    
    "V_vec",                           "I_vec",
    "M_vec",                           "L_vec",
    "K",                               "Imm_Trig", 
    "CD4count",                        "CD4tot",
    "ChronPhase",                      "OnDrug",
    "Adherence1","Adherence2",         "Adherence3","Adherence4",
    "Drug1", "Drug2",                  "Drug3", "Drug4", 
    "aim3_no_muts",                    "adherence_start",
    "adherence_type",                  "CD4",                                             
    "CD4_time",                        "CD4_initial_value",
    "CD4_treatment_delay_index",       "spvl_cat",
    "CD4_time_death",                   "CD4_nadir",
    "start_aids_cd4",                   "start_max_aids",
    "virus_sens_vacc",                  "virus_sens_drug",
    "virus_part_res_drug",              "virus_3_plus_drug_muts",
    "virus_3_plus_drug_muts",             "Aim3RoundingErrors",               
     "aim3_mutations_long",             "CYP_6_slow",
    
# -- testing for and treatment of drug resistant viruses (aim 3) --- #    
    "eligible_2nd_line_ART",
    "treated_2nd_line",  
    "diag_resist_status",
    "diag_resist_time", 
    "last_neg_resist_test",
    "time_init_2nd_line",
    
#--  Vital dynamics /social/treatment ------------------------------ #
    "vaccinated",                      "vacc_init_time",           
    "age","sqrt_age",                   "arrival_time",
    "last_neg_test",                   "diag_status",
    "diag_time",                       "disclosure_status",
    "id",                              "eligible_care",
    "att1",                            "ai_prob",
    "treated",                         "tx_init_time",
    "circum",                           "tx_schedule",               
    "sti_status",                       "sex",
    "total_acts",                        "role",
    "min_time_tx",                      "enhanced_testing",
    "time_hiv_sex_act",                 "time_hiv_sex",
     "last_disc_sex",                   "tx_stop_time",              
    "rand_prob_test",                   "ever_enhanced_testing",
    "rand_prob_test_init",              "tx_dropout",
     "pos_partner_duration",             "no_partners_past_prep",
    "no_partners_now_prep",              "have_diag_partner",
    "have_disc_partner",                 "on_prep",                      
    "prep_decrease",                     "eligible_for_prep",
    "prep_list",                         "known_pos_partner_duration",
    "condom_user")
    ##############################
  ) #end of  list
  
}