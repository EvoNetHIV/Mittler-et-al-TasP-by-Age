
#specify output location and experiment number
local_path <- getwd()
exp_number <- "51"

#--------------------------------------------------------------
#load evonet
library(evonet)
#--------------------------------------------------------------
#Load default parameters
primary_parameters  <- input_parameters_primary()
cd4_data            <- input_parameters_cd4_data()

#--- combine individual parameters into single list
evoparams <- c(primary_parameters, cd4_data)

#--------------------------------------------------------------
# User specific parameters 
evoparams$nsims=16
evoparams$QA_QC_pause_time = 0
evoparams$fast_edgelist = TRUE
evoparams$save_partner_list = TRUE
evoparams$n_steps          = 40*365

# Treatment campaign / limits on the number who can be treated
evoparams$proportion_treated = 1 # 0.22 # used with tx_limit="percentage"
evoparams$prob_care                = 0.95
evoparams$prob_care_after_campaign = 1.0  # This applies only for social_bring_new_groups_into_care_module
evoparams$prob_eligible_ART        = 1.0
evoparams$prob_eligible_2nd_line_ART = 1.0
evoparams$start_treat_before_big_campaign = 5*365 # % treated at "start_treat_before_big_campaign"
evoparams$start_treatment_campaign      = 10*365
evoparams$start_scale_up_campaign      = 5*365
evoparams$proportion_treated_begin = 0.50*evoparams$proportion_treated
evoparams$tx_type = c("all_diag") # Treat everyone who gets diagnosed.  ** THis gets over-written below
evoparams$tx_limit = "absolute_num" # Choices: "absolute_num" or percentage" 
evoparams$tx_schedule_props  =c("F"=1.0,"V"=0.0,"N"=0,"P"=0)
evoparams$max_num_treated = 200  # used with tx_limit="absolute_num"
evoparams$yearly_incr_tx     = 0.02 # 0.05 means 5% more people get treated each year once campaign starts
evoparams$attr_treatment_threshold = 2
evoparams$min_age_recruit_for_care = 25
evoparams$max_age_recruit_for_care = 40
evoparams$vl_full_supp           = 1e-3 # Assume highly suppressiver therapy
evoparams$prob_tx_droput = 0.05    # Probability of "V"-type agent discounting therapy over the course of a year

# Testing for HIV
evoparams$testing_model            = "memoryless"
evoparams$mean_test_interval_male  =  365  #  Default value for people who aren't already requesting testing
evoparams$mean_test_interval_female = 365 #  50*365 would mean ~2% of nontesters come in for testing per year.
evoparams$mean_test_interval_under25 = 365 #  50*365 would mean ~2% of nontesters come in for testing per year.
evoparams$test_result_delay = 30
evoparams$mean_trtmnt_delay = 15
evoparams$prob_enhanced_testing_before_campaign = 1.0 # Maximum percent of population under care (treatment and/or regular testing with follow-through)

# used with "baserand"
evoparams$prob_enhanced_testing_after_campaign = 1.0 # Maximum percent of population under care (treatment and/or regular testing with follow-through)
evoparams$reduction_test_interval_enhanced = 1 # 0.0002*180000 = every 36 days

# Population / demographic parameters
evoparams$initial_pop      = 2000
evoparams$initial_infected = 200  
evoparams$model_sex <- "hetero"
evoparams$asmr_data_male= "south_africa_male"
evoparams$asmr_data_female= "south_africa_female"
evoparams$initial_agedata_male = "stable_age_no_hiv_dist"
evoparams$initial_agedata_female = "stable_age_no_hiv_dist"
evoparams$birth_model = "exponential_growth"
evoparams$baseline_input_exp_growth = 0.007*(evoparams$initial_pop/100)
evoparams$pop_growth_rate_annual   = 0.01
evoparams$min_age = 16
evoparams$max_age = 100

# Factors affecting transmission probabilities (note that many of these are overwritten below)
evoparams$trans_lambda  =  0.000247 * 10 # 20X Hughes et al. (b/c they may be biased) 8/17/17: 20% lower with update_cd4_daily_with_age for exp 36
evoparams$condom_use_rel_dur = FALSE
evoparams$condom_use_age = TRUE
evoparams$age_condom_use_halves = 35
evoparams$condom_prob <- 1.00

evoparams$prob_sex_by_age  = TRUE
evoparams$prob_sex_age_19  = 0.40 # used when prob_sex_by_age == TRUE
#evoparams$mean_sex_acts_day <- 1.0 # NEED TO CHANGE BACK TO 0.36
evoparams$circum_prob <- 0.4
evoparams$RR_cond_male_concurrent <- 1.0
evoparams$RR_cond_fem_concurrent <- 1.0
evoparams$min_age <- 16
evoparams$prob_elig_vl_test <- 0.45

# Printing parameters
evoparams$QA_QC_pause_time = 0
evoparams$print_frequency = 365
evoparams$popsumm_frequency = 365
evoparams$network_print_frequency = 500000
evoparams$plot_nw=F
# Specify output files to create
evoparams$save_vl_list    = FALSE        # Set to true to see each agents VL history
evoparams$save_coital_acts = FALSE      # Set to true see coital acts
evoparams$save_network = FALSE          # Set to true see relationship histories (small runs only)
evoparams$save_infection_matrix = FALSE # Set this AND "save_coital_acts" to TRUE to see the phylogenies
evoparams$save_RData_file = TRUE       # Set to true to save Rdata file
evoparams$save_summary_figs = TRUE      # Set to true to save summary figures file

##########################
### Set up risk groups ###
##########################
evoparams$generic_nodal_att_no_categories = 3
evoparams$generic_nodal_att_values        = 1:3  # short term, medium-term, and long-term relatioships
# Define proportions of people in different risk groups at "birth"
f3_0 = 0.1             # proportion of newly entering agents with a preference for short-term relationships
f2_0 = 0.8             # proportion of newly entering agents with a preference for medium-term relationships
f1_0 = 1 - f2_0 - f3_0 # proportion of newly entering agents with a preference for long-term relationships
evoparams$generic_nodal_att_values_props_births = c(f1_0, f2_0, f3_0)
# Probabilities of changing risk group status over time 
p_long <- 0.00015  # Per day probability of an agent switching from short- or medium-term relationship preferences to a preference for long-term relationships
p_medium  <- 0.0003   # Per day probability of an agent switching from short-term to medium-term relationship preference
# 0.00015 corresponds to a ~5% chance per year of transitioning out
# 0.0003 corresponds to ~11% chance per year of transitining out
# Old setup
#evoparams$generic_nodal_att_trans_mat   <- matrix(
#  c(1 - p_medium - p_long,  p_medium,  p_long,   # 1--> 1, 1-->2, 1-->3
#    0,                      1-p_long,  p_long,   # 2--> 1, 2-->2, 2-->3
#    0,                      0,         1     ),  # 3--> 1, 3-->2, 3-->3
#  nrow=3,byrow=T,dimnames=list(c("G1","G2","G3"))
#)
# New formulation with group 1 now defined to be the long-duration group
evoparams$generic_nodal_att_trans_mat   <- matrix(
  c(1,        0,         0,                        
    p_long,   1-p_long,  0,                        
    p_long,   p_medium,  1 - p_medium - p_long  ),
  nrow=3,byrow=T,dimnames=list(c("G1","G2","G3"))
)
# Define proportions of people in different risk groups at t = 0
# Empirically derived values that match the steady-states that result from the parameters (f1_0,f2_0,f3_0,p_long,p_medium) above
f3 = 0.015       # proportion of agents with a preference for short-term relationships
f2 = 0.35         # proportion of agents with a preference for medium-term relationships
f1 = 1 - f2 - f3  # proportion of agents with a preference for long-term relationships
# Eventually, I would would like to replace these empircally derived values for f1, f2, and f3 with some formulas
# Those formulas would look something like this (under development)...
# f1 =  birth_rate*f1_0                          /  ( 2 - p_medium - p_long - ave_death_short )
# f2 = (birth_rate*f2_0 + p_medium*f1)           /  ( 2 - p_long - ave_death_medium )
# f3 = (birth_rate*f3_0 + p_long*f1 + p_long*f2) /  ( 2 - ave_death_long )
evoparams$generic_nodal_att_values_props  = c(f1, f2, f3)
# Set up formation terms for people in the different risk groups  ###
#evoparams$nw_form_terms <- "~edges + nodemix('att1', base=1) + absdiff('age',pow=1.5) +concurrent +offset(nodematch('sex', diff=FALSE))"
evoparams$nw_form_terms <- "~edges + nodemix('att1', base=1) + absdiff('age',pow=1.0) +offset(nodematch('sex', diff=FALSE))"
md= 0.70 # mean degree
abs_diff_age = 4   # Difference in age
# Matrix formula for determining edge counts assuming all risk groups have the same instaneous mean degree (where E = total edges) 
#         G1           G2         G3
#  G1   E*f1*f1    2*E*f1*f2   2*E*f1*f3
#  G2                E*f2*f2   2*E*f2*f3
#  G3                            E*f3*f3    # Note this formula assumes random mixing (does that matter?)
tot_edges = md * evoparams$initial_pop / 2   # Would add up to 150 for n=1000 and md=0.35
abs_diff_age_param <- tot_edges * abs_diff_age
# Previous verson had some Formulas for correcting mean degrees between the threee groups (which I expect wont be neeed anymore)
evoparams$target_stats <- c(tot_edges, 2*tot_edges*f1*f2, tot_edges*f2*f2, 2*tot_edges*f1*f3, 2*tot_edges*f2*f3, tot_edges*f3*f3,
                            abs_diff_age_param)#0.06*tot_edges/2)  # last one (blanked out here) was a concurrency term
evoparams$nw_coef_form <- -Inf
# Define the dissolution terms 
evoparams$dissolution <- "~offset(edges) + offset(nodemix('att1', base=1))"
D1 = 10*365  # Duration preference for people in group 1-1 pairs (long-term relationships)
D2 = 365    # Duration preference for people in group 2-2 pairs (medium-term relationships)
D3 = 10      # Duration preference for people in group 3-3 pairs (short-term relationships)
# Matrix of durations assuming a geometric mean for relationship duration
#         G1       G2             G3
#  G1     D1    (D1*D2)^0.5    (D1*D3)^0.5 
#  G2              D2          (D2*D3)^0.5    
#  G3                             D3 
evoparams$relation_dur=c(D1, sqrt(D1*D2), D2, sqrt(D1*D3), sqrt(D2*D3), D3)
#add parameters that are functions of other input parameters

evoparams  <- input_parameters_derived(evoparams)

#convert raw parameter list into EpiModel object
evoparams <- do.call(EpiModel::param.net,evoparams)

#params to loop through----------------------------

criterion_vec = list(
  c("random"),
  c("under25","random"),
  c("over30","random"),
  c("CD4_nadir_under500","random"),      # Immunological
  c("S6.0","S5.5","S5.0","S4.5","S4.0","random"),   # Virological
  c("under25","under30","random"),
  c("under25","under30","under35","random")   # Virological
)
criterion_names <- c(1:length(criterion_vec))
for (i in 1:length(criterion_vec)) {
  criterion_names[i] = paste(criterion_vec[i][[1]],collapse="_")
}


#make sure param vecs same length
param_vec = list(0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0)
param_names <- c(1:length(param_vec))
for (i in 1:length(param_vec)) {
  param_names[i] = paste(param_vec[i][[1]],collapse="_")
}


nw = setup_initialize_network(evoparams)

  #estimate initial network (create argument list, then call fxn)
  netest_arg_list <- list(
    nw            =  nw,
    formation     =  as.formula(evoparams$nw_form_terms),
    target.stats  =  evoparams$target_stats,
    coef.form     =  evoparams$nw_coef_form,
    constraints   =  as.formula(evoparams$nw_constraints),
    verbose       =  FALSE,
    coef.diss     =  dissolution_coefs( dissolution =  as.formula(evoparams$dissolution),
                                        duration    =  evoparams$relation_dur,
                                        d.rate      =  3e-05) )
  
  estimated_nw <- do.call(EpiModel::netest, netest_arg_list)

for(ii in 1:1)
{  
  evoparams$tx_type <- criterion_vec[ii]
  
  for(jj in 1:length(param_vec))
  {
    model_name = paste("txt_opto_",exp_number,"_",criterion_names[ii],"_prop",param_names[jj],sep="")
    evoparams$proportion_treated  = param_vec[jj][[1]]
    evoparams$proportion_treated_begin = 0.5 * evoparams$proportion_treated 
    evoparams$yearly_incr_tx     = 0.02 # 0.05 means 5% more people get treated each year once campaign starts

    #caclulate derived parameters (values based on other input parameters)
    evoparams  <- input_parameters_derived(evoparams)
    #convert raw parameter list into EpiModel object
    evoparams <- do.call(EpiModel::param.net,evoparams)

  
  #--------------------------------------------------------------
  #-- create initial vector of infection status as an epimodel object
  infected_list <- EpiModel::init.net(i.num=evoparams$initial_infected,
                                      status.rand = FALSE)
  
  #--------------------------------------------------------------
  #---  Create list of modules to run for input into epimodel_control_fxn() below
  # ***   Note: initialize fxn must always be first and verbose fxn last AND death fxn
  # ***   must precede birth fxn (these conditions may change in future)
  # ***   treatment_fxn must be before update_vl and update_cd4
  
  evo_module_list<- list(
    "initialize.FUN"     = initialize_module,
    "restart.FUN"        = restart_module,
    #"plot_nw.FUN"        = plot_network_fxn, #turn off plotting to save time
    "aging.FUN"          = vital_aging_module,
    "testing.FUN"        = social_testing_diagnosis_module,
    "adherence.FUN"      = social_adherence,
    "resist.testing.FUN" = social_resistance_testing_module,
    "dropout.FUN"        = social_treatment_dropout_john,
    "treatment.FUN"      = social_treatment_module_multiple_criteria_v2, # social_treatment_module, social_treatment_module_john, social_treatment_module_multiple_criteria
    "treat.2nd.line.FUN" = social_2nd_line_treatment_module, # social_treatment_module, social_treatment_module_john
    "update_vl.FUN"      = viral_update_gamma_john, # viral_update_gamma, viral_update_aim3
    "update_cd4_FUN"     = viral_update_cd4_daily, # viral_update_cd4_diff_eqn, viral_update_cd4_simple_diff_eqn, 
    "coital_acts.FUN"    = social_coital_acts_module,
    "trans.FUN"          = transmission_main_module,
    "trans_book.FUN"     = transmission_bookkeeping_module,
    "trans_cd4.FUN"      = transmission_cd4_module,
    "deaths.FUN"         = vital_deaths_module_john,
    "births.FUN"         = vital_births_module,
    "social_trans.FUN"   = social_attribute_transition_module,
    "summary.FUN"        = summary_module,
    "resim_nets.FUN"     = EpiModel::resim_nets,
    "verbose.FUN"        = NULL)
  
  #--- call epimodel's control fxn (load evonet modules into epimodel)
  evocontrol <- setup_epimodel_control_object(evonet_params = evoparams,
                                              module_list   = evo_module_list)
  
  #--------------------------------------------------------------

    evomodel  <- EpiModel::netsim(x = estimated_nw,
                                  param = evoparams,
                                  init = infected_list,
                                  control = evocontrol)

  plots_popsumm(evomodel,outpath=outpath,
                name=model_name,nw_stats=TRUE,max_points_rep=100,
                evoparams$popsumm_frequency)

  assign(model_name,evomodel)
  file_name <- paste(model_name,".RData",sep="")
  save(list=model_name,
       file = file.path(outpath,file_name) )
  remove(evomodel)
  remove(list=model_name)
}#end of loop
}  
#--------------------------------------------------------------
