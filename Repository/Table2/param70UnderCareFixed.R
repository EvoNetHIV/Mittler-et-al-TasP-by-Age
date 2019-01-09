# Deviations for default parameters used in revised TasP-by-Age simulations
# Note: Some of these will be modified by the calling program

# Load evonet
library(evonet)

# Create list
evoparams=list()

# Top-level simulation parameters

evoparams$initial_pop      = 2000
evoparams$initial_infected = 150
evoparams$n_steps          = 45*365
evoparams$nsims            = 16
evoparams$ncores            = 16


# Treatment campaign / limits on the number who can be treated
evoparams$start_treatment_campaign      = 20*365   # Start time of big TasP campaign
evoparams$start_treat_before_big_campaign = 19*365 # Start time of gradual ramp-up prior to the start of the big TasP campaign
evoparams$proportion_treated = 1 # After the start of the TasP campaign *** Typically changed by calling script ***
evoparams$proportion_treated_begin = 0 # Treated before gradual ramp-up 
evoparams$prob_care                = 0.70 # Percent population that could get treated given "all out" campaign
evoparams$prob_eligible_ART        = 1.0
evoparams$tx_limit = "absolute_num" # Choices: "absolute_num" or percentage"
evoparams$yearly_incr_tx  = 0.02 # 0.05 means 5% more people get treated each year once campaign starts
evoparams$vl_full_supp    = 1e-3 # Assume highly suppressiver therapy
evoparams$tx_schedule_props  = c("F"=0.5,"V"=0.5,"N"=0,"P"=0) # Percentage of people who always (F), sometimes (V), or never (N) take therapy
evoparams$prob_tx_droput  = 0.00   # Probability of "V"-type agent discontinuing therapy over the course of a year
evoparams$mean_trtmnt_delay = 15    # Delay betweent diagnosis and start of treatment

# Criteria for targeted TasP campaign (requires evonet function: "targeted_treatment2")
TasP_criteria = list(
  c("random"),
  c("under25","random"),
  c("CD4_nadir_under500","random"),      # Immunological
  c("S6.0","S5.5","S5.0","S4.5","S4.0","random"),   # Virological
  c("under25","under30","random"),
  c("under25","under30","under35","random")   # Virological
)

# Testing for HIV
evoparams$testing_model            = "memoryless"
evoparams$mean_test_interval_male  =  365  #  Default value for people who aren't already requesting testing
evoparams$mean_test_interval_female = 365 #  50*365 would mean ~2% of nontesters come in for testing per year.
evoparams$mean_test_interval_under25 = 365 #  50*365 would mean ~2% of nontesters come in for testing per year.
evoparams$test_result_delay = 30    # Amount of time that it takes to get test results.

# Demographic parameters
evoparams$model_sex <- "hetero"
evoparams$asmr_data_male= "south_africa_male"
evoparams$asmr_data_female= "south_africa_female"
evoparams$initial_agedata_male = "stable_age_no_hiv_dist"
evoparams$initial_agedata_female = "stable_age_no_hiv_dist"
evoparams$birth_model = "exponential_growth"
evoparams$baseline_input_exp_growth = 0.007*(evoparams$initial_pop/100) # Calibrated by hand (assumes that initial_agedata_male = "stable_age_no_hiv_dist" )
evoparams$pop_growth_rate_annual   = 0.01
evoparams$min_age = 16
evoparams$max_age = 99

# Factors affecting transmission probabilities (note that some may be overwritten by calling program)
evoparams$trans_lambda  =  0.000247 *15 # some value times Hughes et al. (b/c they may be biased)  
evoparams$trans_RR_receptive_vaginal = 1.5
evoparams$condom_use_rel_dur = FALSE
evoparams$condom_use_age = TRUE
evoparams$age_condom_use_halves = 35
evoparams$condom_prob <- 0.80 # When condom_use_age == TRUE this gives probability that 16-year old uses condom, otherwise this is the probability for everyone

evoparams$prob_sex_by_age  = TRUE
evoparams$prob_sex_age_19  = 0.20 # used when prob_sex_by_age == TRUE
evoparams$circum_prob <- 0.4
evoparams$RR_cond_male_concurrent <- 1.0
evoparams$RR_cond_fem_concurrent <- 1.0
evoparams$min_age <- 16

evoparams$prop_AI                  = 0.0
evoparams$mean_prop_acts_AI        = 0.0

######################################################
############ Set up risk groups ######################
######################################################

md= 0.75 # mean degree
male_female_age_diff = 4 # 4 means that women partner with men who are, on average, 4 years older than them
abs_diff_age = 4   # Difference in age AFTER accounting for the male-female age difference

# Assume two risk groups hat differ in their relationship-duration tendancies
evoparams$generic_nodal_att_no_categories = 2
evoparams$generic_nodal_att_values        = 1:2  # 1 = short-term, 2 = long-term relatioship
D1 = 2*365     # Duration preference for people in group 1-1 pairs (short-term relationships)
D2 = 10*365  # Duration preference for people in group 2-2 pairs (long-term relationships)
f1_0 = 0.9         # proportion of newly entering agents that have a preference/tendancy for short-term relationships
f2_0 = 0.1         # proportion of newly entering agents that have a preference/tendancy for long-term relationships
p_long <- 0.00011  # Per day probability of an switching from short-to long-term preference (0.00011 amounts to ~4% probability/year)
f1 = 0.42          # Steady-state proportion of agents with a tendancy to form short-term relationships.  This value was determined empirically!
f2 = 1 - f1        # proportion of agents with a preference/tendancy to form long-term relationships

concur_women = 0.07 * evoparams$initial_pop / 2 #  Approx midpoint of studies collected by Kathryn
concur_men = 0.17 * evoparams$initial_pop / 2  # Approx midpoint of studies collected by Kathryn rounded to nearest tenth

# Functional equations
evoparams$generic_nodal_att_values_props_births = c(f1_0, f2_0)
evoparams$generic_nodal_att_trans_mat   <- matrix(
  c( 1- p_long, p_long,
     0,          1),
  nrow=2,byrow=T,dimnames=list(c("G1","G2"))
)
evoparams$generic_nodal_att_values_props  = c(f1, f2)
evoparams$nw_form_terms <- "~edges + nodemix('att1', base=1) + absdiffby('age','sex',4,c('f','m')) + concurrent(by='sex')  +offset(nodematch('sex', diff=FALSE))"
#evoparams$nw_form_terms  <- "~edges + nodemix('att1', base=1) + absdiff('age',pow=1.0) +  offset(nodematch('sex', diff=FALSE))"
#evoparams$nw_form_terms <- "~edges + nodemix('att1', base=1) + absdiff('age',pow=1.0) + concurrent(by='sex')+ offset(nodematch('sex', diff=FALSE))"

tot_edges = md * evoparams$initial_pop / 2   # Would add up to 150 for n=1000 and md=0.35
abs_diff_age_param <- tot_edges * abs_diff_age
# Matrix formula for determining edge counts assuming all risk groups have the same instaneous mean degree (where E = total edges)
#         G1           G2
#  G1   E*f1*f1    2*E*f1*f2
#  G2                E*f2*f2
#evoparams$target_stats <- c(tot_edges, 2*tot_edges*f1*f2, tot_edges*f2*f2,   abs_diff_age_param)  # Original idea
evoparams$target_stats <- c(tot_edges, 2*tot_edges*f1*f2, tot_edges*f2*f2,   abs_diff_age_param, concur_women, concur_men)  # Original idea
evoparams$nw_coef_form <- -Inf
evoparams$dissolution <- "~offset(edges) + offset(nodemix('att1', base=1))"

# Matrix of durations assuming a geometric mean for relationship duration
#         G1       G2
#  G1     D1    (D1*D2)^0.5
#  G2              D2
evoparams$relation_dur=c(D1, sqrt(D1*D2), D2)

# Specify evonet functions
modules <- c(
  "aging",
  "social_testing_diagnosis_eligibles_only_module",
  "adherence",
  "treatment_dropout",
  "targeted_treatment2",
  "viral_update_delayed_rebound",
  "cd4_update2", # viral_update_cd4_daily
  "coital_acts",
  "transmission",
  "deaths",
  "births",
  "risk_group_changes",
  "summary_module")

# List of names of TasP criteria for loop script
criterion_names <- c(1:length(TasP_criteria))
for (i in 1:length(TasP_criteria)) {
  criterion_names[i] = paste(TasP_criteria[i][[1]],collapse="_")
}

# Specify output files to create
evoparams$save_vl_list    = FALSE        # Set to true to see each agents VL history
evoparams$save_coital_acts = FALSE      # Set to true see coital acts
evoparams$save_network = FALSE          # Set to true see relationship histories (small runs only)
evoparams$save_infection_matrix = FALSE # Set this AND "save_coital_acts" to TRUE to see the phylogenies
evoparams$save_RData_file = TRUE       # Set to true to save Rdata file
evoparams$save_summary_figs = TRUE      # Set to true to save summary figures file
evoparams$save_partner_list = TRUE

# Printing parameters
evoparams$QA_QC_pause_time = 0
evoparams$print_frequency = 365
evoparams$popsumm_frequency = 365
evoparams$network_print_frequency = 500000
evoparams$plot_nw=F

# System options
options(error=NULL)
options(warn=-1)

# Network simulation method
evoparams$fast_edgelist = TRUE


# Create consolidated "evoparams" list 
param_list=evoparams
evoparams <- do.call(evonet_setup,param_list) # Default parameters for underlying evonet programs






