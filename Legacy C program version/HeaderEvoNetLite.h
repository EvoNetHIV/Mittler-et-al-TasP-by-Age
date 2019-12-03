/*
 *  NetworkHeader.h
 *  
 *
 *  Created by John Mittler on 12/6/12.
 *  Copyright 2012 University of Washington. All rights reserved.
 *
 */
// New parameters added 3/25/18

long NetStatsCount;
double PercentNonTargetedYear1;
double ActualNonTargetedYear1;
double baseline_treat_prob, total_prob;
double ActualNonTargetedYear1;
double treatment_goal_start; 
double Start_Spontaneous_Treatment;
long repl;
long RandomIndex[2000001];
long old_random_number_seed;
int first_run = 1; // used to ensure that big outputs are made only on first of many runs

long replicates = 1;
long N0, Infected0;

// Run options
int randomized_parameters = 0; // if 1 multiple parameters get randomized before the simulation
int gradual_tx_incr = 1;  // If 1 the number of newly treated agents increases at a constant rate per year.  If 0, tx increases suddenly
                          // though the model allows for a yearly incr afterwards and a pre_TasP ramp-up.

// New code for testing
void HIVTesting(void);
int Diagnosed[2000001];
long TimeDiag[2000001];
double daily_prob_diagnosis = 1.0/365.0;
long delay_diag_tx = 14.0;
long delay_test_antibodies = 35.0;

// Treatment strategies
long under25 = 0, under30 = 0, under35 = 0, cd4_low = 0, cd4_super_low = 0, vl_high = 0, vl_very_high = 0, vl_super_high = 0, men = 0, women = 0, men27women23 = 0;
long men23women27 = 0, men33women27 = 0, men27women33 = 0;
double lowest_tx = 0.0, highest_tx = 1.0, tx_incr = 0.1;
double fold_incr_prob_with_targeting = 10.0;
double yearly_incr_tx = 0.01;
int TasP_Strategy = 1;

// Tracking variables for statistics
long susc_days = 0, susc_days_u50 = 0, susc_days_u25 = 0, susc_days_Mid = 0, susc_days_o50 = 0;
long new_infections = 0, new_infections_u50 = 0, new_infections_u25 = 0, new_infections_Mid = 0, new_infections_o50 = 0;;
long Time_All_Treated = 100000000, Last_Time_Not_All_Treated = 0;
long susc_days_last_5_years = 0, susc_days_last_5_years_u50 = 0, new_infections_last_5_years = 0, new_infections_last_5_years_u50 = 0;
long InfectedTreated = 0;
long pill_count, TotalTreated, NonTargetedTx, NonTargetedTxYear1, TotalTreatedYear1, TotalDiedAIDSAfterTasP=0;
long NumPartners[2000001], NumPartnersLastTime[2000001], Time_AIDS_Death[2000001];
double ActTxStartTasP = 0.0; // Percent actually treated at the start of the TasP campaign (mainly useful for gradual treatment increase)
long TrackedAgentsGroup1start = 10000000, TrackedAgentsGroup1end = 10000000, TrackedAgentsGroup2start = 10000000, TrackedAgentsGroup2end = 10000000; 
long tx_person_days = 0, cd4_cat1_person_days = 0, cd4_cat2_person_days = 0, cd4_cat3_person_days = 0, cd4_cat4_person_days = 0, cd4_cat5_person_days = 0; // For DALY calcs

// DALY estimates
double DALY00, DALY01, DALY02, DALY03, DALY04, DALY05, DALY06, DALY07, DALY08, DALY10, DALY15, DALY20; // DALY03 means DALYs with 3% discount rate
double life_expectancy = 60.0;
double cost_hiv_treated = 0.003;
double cost_cd4_gt_350 = 0.02;
double cost_cd4_200_350 = 0.05;
double cost_cd4_lt_200 = 0.5;
double cost_died_AIDS = 1.0;

// Plotting options
int plt_GS=1;
int plt_NS= 1; 
int plt_out=0;
int plt_H=0;
int plt_R=0;
int plt_AH=0;
int plt_AM=1;
int plt_AD=1;
int plt_VL=0;
int plt_spVL=0;
int plt_warn=1;
long print_frequency = 365;

// Infection risk factors
double MaxAgeSex = 75.0;
double circum_prob = 0.4;
double RR_female_recipient = 1.5;
double RR_circum = 0.5;
double RR_youth = 2.0;
double RR_condom = 0.2;
double condom_prob16 = 0.69;
double age_condom_use_halves = 34.0; // Approx fit 2012 SA Nat Prev, Incid & Behav Surveyjj
double number_tx_spontaneous = 0.0;
double prob_sex_in_AIDS = 0.1;
double male_concurrency = 0.40;
double female_concurrency = 0.04;
double prob_conc_more_partners = 0.0;
double LongevityFudgeFactor = 1.0;
double increased_duration_per_year = 60.0; // In days.  If 50, duration will increase 500 days in 10 years, 2500 days (=6 years) in 50 years
double rate_stops_being_concurrent = 0.04 * 1.0/365; // 2.5% chance per year ==> ~50% drop over 25 years (as in Maugh-Brown et al. (2015 or 2015)  
double perc_tx_start = 0.80; 
double percent_high_risk = 0.0;
double reduced_linkage_low_risk = 1.0;

// Standard C functions
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

//  Key evonet_Lite functions
void GetParameters(void), PrintStats(void);
void OpenFiles(void);
void GetParameters(void);
void InitializePopulation(void);
void InitialAgeDistribution(void);
void CalculateLifeExpectancies(void);
void Births(void);
void Aging(void);
void UpdateDALYs(void);
void AddInfecteds(void);
void DefineSexualContactNetwork(void);
void PrintEdgeList(void);
void PrintHeaders(void);
void PrintStats(void);
void PrintHeritabilityStats(void);
void PrintPartnershipDistributions(void);
void InitiateTreatment(void);
void InitiateTreatmentGradual(void);
void TreatmentDropout(void);
void CollectStatistics(void);
void InitiateGradualTreatment(void);
void GetTreatmentProbs(double baseline_treat_prob);
void GetTreatmentProbsMultiplicative(double baseline_treat_prob);
void TreatAgents(int targeted_treatments, long num_newly_treated);
void GetRandomTreatmentProbs(double baseline_treat_prob);
void UpdateCD4Counts(void);
void UpdateCD4CountsStochastic(void);
void UpdateViralLoads(void);
void UpdateViralLoadsAim3(void);
void SimulateTransmission(void);
void BirthOfNewSusceptibles(void);
void NaturalDeath(void);
void DeathOfAIDSPatients(void);
void RemoveLinksToDeadPeople(void);
void SimulatePartnerShipDissolution(void);
void SimulateNewPartnershipFormation(void);
void CheckEdgeListForImpossibleValues(void);  // Does some elementary error checking
void RecordStatusOfLivingHIVPatientsForRegistryFile(void); 
void PrintFirstAndLastContacts(void);
void PrintAgeDistributions(void);
void PrintTrackedAgents(void);
void ConstantInputNewSusceptibles(void);
void VaccinatePopulation(void);
void SimulateVaccineDecay(void);
void UsualInitialValues(long i);
void PrintVLdistribution(void);
void AddNewLinks(long LinksToAdd);
void AddNewLinks2(long LinksToAdd);
void AddLinkForSinglePerson(long person1);
void PrintAgeMatches(void);
void PrintRelationshipRegistry(void);

// General use functions
double norm_rand(double rn_mean, double rn_variance);
double jmpow(double a, double b);
double GetGammaDelay(double k, double theta);
long GetPoisson(long num, double pois_prob);
long max(long a, long b);
long rpois(double lamdba);

// Special case parameters
long Flat_Viral_Load = 0; // Setting this to 1 forces VL to be the same value during the entire primary infection period
double PrimaryInfectionLevel = 1.0e7; // Average viral load during primary infection
double expected_time_treated;
long Progression_Model;

// Printing parameters
long VL_print_time; // Instructs program to print out VL distributions every VL_print_time days.
double VL_print_lower, VL_print_upper, VL_print_interval;  // Lower and upper end of the distribution together with bin size
long VL_print_count = 0;
long printcount;
int printOutput;

// Network Parameters
double OrigAverageDegree;
long BreakUpList1[3000001], BreakUpList2[3000001], BreakUpList3[3000001], BreakUpList4[3000001]; // These arrays track indices of partnerships to be broken up
long RevisedNumLinks[3000001], TimeLastBreakUp[2000001];
long num_cases_no_link_found = 0;
int time_varying_mean_degree = 1;
double d_rate = 0.00025;

// New lines for vaccination
int Vaccinated[3000001], ResistantToVaccine[3000001];
long LengthVaccinated[3000001];
double PercentVaccinated, VaccineDuration, PercentResistantToVaccine, RR_Vaccinated;

// Condom-campaign stuff
double MutationVariance, H, h;
long OrigN;
long Start_Condom_Campaign;
double Percent_using_condums_after_condom_campaign;
double orig_prob_sex;


/* Files used by this program */
FILE *ParFile; // Contains parameters
FILE *ParOut;  // Output file that contains parameters only (used to verify that parameters were read in correctly) 
FILE *OutputFile; // Records average and standard deviation of viral loads (and other parameters) at each time step
FILE *NumPartnersFile; // Histogram num partners.  1st col = 0 partners, 2nd = 1 part, 21st col = 20+ partners. rows: total, women, men, concur-women, concur-men 
FILE *HeritOutputFile; // Records relationship between donors and recipients at a single timepoint (tprintHerit)
FILE *PatientRegistryOutputFile; // Records characteristics of each patient who died and state of living patients at tfinal
FILE *VLOutputFile; // Gives a histogram of log viral loads at different times
FILE *spVLOutputFile; // Gives a histogram of (log-transformed) set point viral loads at different times
FILE *NetworkStatsOutputFile; // Gives summaries of mean degree, etc. 
FILE *AgeDistFile; // Gives age distributions at various times during the simulation
FILE *RelLengthFile; // Gives average relationship lengths (uncensored) for people in different age groups
FILE *AgeMatchFile; // Gives list ages of all couples to test age-related homophily terms
FILE *RelationshipRegistryFile; // Gives list ages of all couples to test age-related homophily terms
FILE *GrandSummaryFile; // Gives final statistics from each run
FILE *AgentHistoryFile; // Detailed time course for selected agents
FILE *WarningsFile; // List of warning messages from each run

/*Epidemiological parameters */
double Inf_Prob; // Probability of infection (calculated from quantities below)
double MaxInfRate, VHalfMaxInfRate, HillCoeffInfRate, beta; // Maximum probability of infection,  viral load at which infection prob is 0.5 of the max
double BaselineDeathRate, MaxAIDSDeathRate, VHalfMaxAIDSDeathRate, HillCoeffAIDSDeathRate, alpha; // Rate of inputs and deaths of people
double PopGrowth; // Rate that population increases each year
double yearly_rate, daily_rate; // rates of progression to AIDS using John's new formula
double prog_rate; // prog_rate is rate of increase of viral load per year
long Inf_Rate_Function; // 1 = power function from Lingappa et al, 0 = hill function  
double InfRateBaseline, InfRateExponent; // Prob(Infmission) = InfRateBaseline * [log(VL)] ^ InfRateExponent
double BirthRate; // Rate of entry of new susceptibles into the population.
double NaturalDeathRate; // Rate at which people normally die (of causes other than AIDS)
long MaxN = 2000000; // NOTE: This is the extent of the remaining arrays minus 1
long MaxNReached = 0; // Flag to indicate if population was about to exceed the maximum array size
double tMaxNReached = 0.0; // Time that births or reincarnations were halted
long NewAIDSDeaths = 0; // Number of deaths in each time period (only used when reincarnation is allowed)
long Widowed_partner;
long newbirths; // Internal counter to count number of births each day
double Active; // Internal -- number of sexually active people
double asmr[101]; // Age-specific mortality rate.
 
/* Virological parameters */
double V[2000001], I[2000001], M[2000001], L[2000001];
double Donors_V[2000001], s[2000001];
double V0, r0, V_peak, Max_VL_AIDS, VL_Incr_AIDS, t_peak, t_acute, Vss;
double d_acute[2000001], Donors_d_acute[2000001];
double death_rate_constant, death_rate_exponent; // Death rate parameters for daily probability death function (used only when Gamma_Death !=1)
double VL_Start_AIDS[2000001]; 
long   Time_Start_AIDS[20000001];
long   TimeInAIDS[20000001];

double ExpectedDelayTime; // Assume gamma distributed delays
double Dmax, Dk, D50; // Parameters governing Fraser's gamma distributed times to AIDS
double shape_parameter,theta; // These parameters control the mean and variance of the gamma delay

/* Infection status parameters */
long N, Infected, Susceptible; // Initial number of people, Initial number infected, Number who could be infected
long currN;  // Value of N used to calculate the number of births that will occur
long Dead, DiedAIDS, DiedNat; // Cumulative Number Died of AIDS or natural causes
double prob_sex; // Probability that partners will have sex each day
int age_dependent_sex = 1; // Probablity of sex drops linearly until age 75 if set to 1
int age_dependent_mixing = 1; // People more likely to partner with someone their own age if set to 1
double Time_Inf[2000001], Donors_Total_Time_Inf_At_Trans[2000001]; // Time infected
int Status[2000001]; // Status of each person (Infected=1, Suscp = 0, Died = -1, DiedAIDS = -2)
long Generation[2000001], Donors_Generation[2000001], Donors_Index[2000001], NumRecipients[2000001];
long HoneymoonDays;
double AgeAIDS[2000001];
double IncreasedProbSexHoneymoonDays;

/* Risk group status */
int HadSexAlready[2000001], HighRisk[2000001];
double percent_high_risk;
long Sex[2000001], Concurrent[2000001], FSW[2000001];
double Age[2000001];
int Circumcised[2000001];
 
/* Parameters governing heritability of VL setpoint */
double AverageLogSP0, VarianceLogSP0, Heritability; // Primary inputs
double SetPoint[20000010], Donors_SetPoint[20000010]; // Patient specific values (drawn from random number generator)
double LogSetPoint[20000010], Donors_LogSetPoint[20000010]; // Patient specific values (drawn from random number generator)
double ViralContribToLogSP0[2000001],EnvirContribToLogSP0[2000001]; // Viral and environmental deviations that contribute to final setpoint
double Donors_ViralContribToLogSP0[2000001], Donors_EnvirContribToLogSP0[2000001]; // Viral and environmental deviations that contribute to final setpoint

/* CD4 Dynamics */
int CD4[2000001];
int CD4_nadir[2000001];
int CD4_initial_value[2000001];
double CD4time[2000001];
int SpvlCat[2000001];
int CD4_Exp_Flag = 1;
double CD4_TimeToAIDS[2000001];
double CD4_TimeToAIDS_exp[2000001][5];
int enhanced_progression_age = 0;
long CD4_treatment_delay_index[2000001];
double prob_cd4_rebound, prob_cd4_rebound_further;

/* Look up table for CD4 counts */
double CD4_lookup[10][5];//cd4 as function of SPVL and time since infection
void CD4_Category();
int SPVL_Category(double SPVL_param);
double cd4_initProbs[10][4];
void cd4InitProbs();
void CreateCD4Table();
int initialCD4(int spvl_level);
int CD4_After_Treatment; // Viral load after treatment starts
int CD4_Determines_Treatment; //if 0, vl determines treament, if 1, then cd4

/* Statistical results (and internal params needed to get those statistics) */
double Vsum, Vave, Vstd, Vstdsum; // Statistics of viral load
double set_ave, set_sum, set_std, set_stdsum; // Statistics on logSetPoint
double d_ave, d_sum, d_std, d_stdsum; // Statistics on decay rate of infected cells after peak (our way of setting setpoints)
double G_ave, G_sum, G_std, G_stdsum; // Statistics on viral generation
double vcount; // Needed for averaging above quantities
double AveLinks,SumLinks;  // Average number of links per person

/* Treatment Parameters */
int first_tx_strategy = 1, last_tx_strategy = 1;  // Strategy 1 = random, strategy 2 = under age 25, ... stategy 9 = under age 25 and CD4 > 500
long UnderCare[2000001]; // Pre-defined list of patients who would get treated after infection (e.g., they live near a clinic)
double ProbDropout[2000001]; // Pre-defined list of patients who would get treated after infection (e.g., they live near a clinic)
long Treated[2000001]; // Patients currently being treated
long Time_Treated[2000001]; // Day that patient started therapy
long TargetedTreated[2000001]; // Patients currently being treated specifically b/c of a targeted TasP campaign
double treatment_prob[2000001], rel_treatment_prob[20000001], cum_treatment_prob[2000001]; // Used in Initiate Treatment routine
long tx_days[2000001]; // Number of days that each agent has been treated. (Useful for adjusting VL after a treatment interruption)
long tx_days_aids[2000001]; // Number of days agent treated after onset of AIDS  (Useful for adjusting VL after a treatment interruption)
int ever_aids[2000001]; // Agent ever progressed to AIDS (even with therapy)
double Start_TasP_Campaign; // Time (day of epidemic) when public health authorities start new treatment campaign
double prob_dropout; //  Per year probability of a treated person discontinuing therapy
double percent_under_care; // Percent of people who get treated after being infected
double Increment_Prob_Under_Care_High_CD4; // Increase in probability of going under care for people with low CD4 counts (= high CD4 categories)
double VL_After_Treatment; // Viral load after treatment starts
double Start_SexReduction_Campaign; // # Day of epidemic when public health authorities start treatment campaign
double Reduction_Mean_Degree; //  Reduction in mean degree resulting from this campaign
double Start_Faithfulness_Campaign; // Day of epidemic when public health authorities start campaign to reduce partnership turnover
double Increased_Duration; // Increase in partnership durations resulting from this campaign
double dropout_prob; // Per per probability of discontinuing therapy
/* Time-keeping variables */
long time, tfinal, tprintHerit;

/* Stopping paramters */
int StopEarly = 0;   // This variable is set to 1 when the program is unable to find a new link

/* Internal variables */
long i, j, k, count, daycount, deadcount;
double epsilon = 1e-10;
char junk[12],descript[30];

/* Variables required for random number generation*/
double nrand(void); // Normal random numbers with mean 0 and stdev 1
double randnum, normalrand; // internal variables
long random_number_seed; // This is the randum number seed
void init_genrand(unsigned long SSS); // All of the rest are existing MT functions
long genrand_int31(void);
unsigned long genrand_int32(void);
double genrand_real1(void);
double genrand_real2(void);
double genrand_real3(void);
double genrand_res53(void);

/* Sexual network parameters */
long NumLinks[2000001], Links[2000001][10]; // Number of sexual links for each person, Identities of those links (up to 5)
long TimeLinked[2000001][10];
long TotalEdges; // Total number of sexual links in the entire population (note that one link links two individuals
double AverageDegree; // Average number of sexual partners person
long MaxLinks; // Maximum number of links per person (MaxLinks = 1 implies monogamy)
long NewLinks; // Temporary variable giving the number of links for each new person entering the population
long Link_back; // Internal variable to find links back from widow's of the dead person back to the dead person
double dp_TotalEdges; // double precision version of TotalEdges (for reading in only, immediately converted to long)
double Duration[2000001]; // Length that person i tends to stay in a relationship
double MinDuration, MaxDuration; // Mins and Maxs used to draw a uniform random number for each person
double PartnershipDuration;  // Duration of each partnership (mean of the underlygin partnership duration of each partner)
double ProbPartnersBreakUp;  // Prob of partnerships breaking up 
long PartnersOfPerson1; // Internal parameter: index of each person's partners 
long LinksToRemove; // Temporary variable indicating number of links per dead person that should be purged
long AlreadyLinked; // Bookkeeping variable to prevent same persons from linking up more than once (but still counting this as two links)
long failed_to_find_a_new_link = 0; // Set to one if algorithm is unable to find a new pair (indicates some kind of programming problem)
double tr; // Turnover rate = Links_Broken_Per_Day/TotalEdges;
long escape_counter; // Used to prevent endless loops (in case of programming error)
long random_partner; // Which of person i's partners (if more than one) does person i have sex with today
long index_to_person2, found_link, index_to_person1; // Used internally to add and break links
long link_break; // Used internally to count up links that need to broken
long breakups; // Internal parameter: how many people broke up each day
double NLinks; // Number of people linked after deaths, births, and breakups
double dpLinksToAdd; // Number of links to add after accounting for deaths, births, and breakups
long LinksToAdd; // Number of links to add after accounting for deaths, births, and breakups
double ExpectedLinks; // Expected number of links that you would have had had their been now births, deaths, or breakups 
long HadDifficutlyFindingLinks = 0; // Flag set when the program had difficulty finding new links

