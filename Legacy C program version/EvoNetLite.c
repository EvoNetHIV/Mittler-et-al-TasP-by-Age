/************************************************************************/
/*********************    HIV NETWORK SIMULATOR    **********************/
/************************************************************************/
/*                                                                      */
/*  Goal: Simulate effect of targeted treatment on HIV prevalence       */
/*                                                                      */
/*  Dependencies: Variables declared in a header file listed below      */
/*                                                                      */
/*  Input: Reads parameters from a file called "parameters.txt"         */
/*                                                                      */
/*  Outputs: To console & ~13 "__.txt" files defined in the header file */ 
/*                                                                      */
/************************************************************************/

#include "HeaderEvoNetLite.h"  // Includes names for variables and subroutines

double OrigAverageDegree; 
double LifeExpectancy[101];

int main()
{
 printf("Starting!!!\n");
 /* Initialize the population */
 OpenFiles();
 GetParameters();
 OrigAverageDegree = AverageDegree;
 BirthRate = BirthRate*N0/2000;
 old_random_number_seed = random_number_seed;
 init_genrand(random_number_seed);
 PrintHeaders();
 // Experiment 2 specific code -- randomize a bunch of epidemiological parameters */
 if (randomized_parameters == 1) {
   MaxInfRate = MaxInfRate * (1 + 0.5*genrand_real2() - 0.25);
   MinDuration = MinDuration * (1 + 0.5*genrand_real2() - 0.25); 
   MaxDuration = MaxDuration * (1 + 0.5*genrand_real2() - 0.25);
   male_concurrency = male_concurrency * (1 + 0.5*genrand_real2() - 0.25);
   female_concurrency = female_concurrency * (1 + 0.5*genrand_real2() - 0.25);
   LongevityFudgeFactor = LongevityFudgeFactor * (1 + genrand_real2());
   Infected = Infected * (1 + 0.5*genrand_real2() - 0.25);
   circum_prob = circum_prob * (1 + 0.5*genrand_real2() - 0.25);
   increased_duration_per_year = increased_duration_per_year * (1 + 0.5*genrand_real2() - 0.25);
   prob_sex = prob_sex * (1 + 0.5*genrand_real2() - 0.25);
 }
 //for (yearly_incr_tx = 0.0; yearly_incr_tx <= 10.01; yearly_incr_tx = yearly_incr_tx + 1.0) {
 for (TasP_Strategy = first_tx_strategy ; TasP_Strategy <= last_tx_strategy; TasP_Strategy++) {
 
   init_genrand(old_random_number_seed); // All strategies start with the same random number seed 
   printf("For Strategy %d, seed was %ld: first random number is %lf\n",TasP_Strategy,old_random_number_seed,genrand_real2()); 
   if (first_run == 0) {plt_NS = 0; plt_AD = 0; plt_AM = 0;}
   cd4_low = 0; under25 = 0; under30 = 0; under35 = 0; vl_high = 0; vl_very_high = 0; vl_super_high = 0; men27women23 = 0; men = 0; women = 0; cd4_super_low = 0;
   // Strategy 1 -- no targeting
   if (TasP_Strategy == 2) under25 = 1;
   if (TasP_Strategy == 3) {under25 = 1; under30 = 1;}
   if (TasP_Strategy == 4) {under25 = 1; under30 = 1; under35 = 1;}
   if (TasP_Strategy == 5) cd4_low = 1;
   if (TasP_Strategy == 6) {cd4_low = 1; under25 = 1;}
   if (TasP_Strategy == 7) vl_super_high = 1;
   if (TasP_Strategy == 8) {vl_super_high = 1; vl_very_high = 1;}
   if (TasP_Strategy == 9) {vl_super_high = 1; vl_very_high = 1; vl_high = 1;}
   if (TasP_Strategy == 10) {men27women23 = 1;}
   if (TasP_Strategy == 11) {men = 1;}
   if (TasP_Strategy == 12) {women = 1;}
   if (TasP_Strategy == 13) {cd4_super_low = 1;}
 
   for (perc_tx_start = lowest_tx; perc_tx_start <= highest_tx; perc_tx_start = perc_tx_start + tx_incr) {
     first_run = 0;   
     for (repl=1;repl<=replicates;repl++) {
       time = 0;
       printf("\nStarting new simulation: Percent treated = %lf, yearly_incr_tx = %lf, Under25 = %ld, Under30 = %ld, Under35 = %ld, cd4_low = %ld, cd4_aids = %ld, vl_high = %ld, vl_very_high = %ld, vl_super_high = %ld, Replicate = %ld\n",
             perc_tx_start,yearly_incr_tx,under25,under30,under35,cd4_low,cd4_super_low,vl_high,vl_very_high,vl_super_high,repl);   
       //printf("Repl\tTime\tN\tLinks\tUnInf\tInf\tDiag\tTx\tNoTx\tDead\tdAIDS\tdNat\tAge\tLogV\t\tVstd\t\tSet_Ave\t\tSet_std\t\td_ave\t\td_std\t\tMD\tUnLink\tNewInfs\tDiscCouples\n");
       printf("Repl\tTime\tN\tLinks\tUnInf\tInf\tDiag\tTx\tNoTx\tDead\tdAIDS\tdNat\tAge\tLogV\tSPVL\tMD\tUnLink\tNewInfs\tDiscCouples\n");
       orig_prob_sex = prob_sex;
       H = Heritability;
       h = sqrt(Heritability);
       Dead = 0; DiedAIDS = 0, TotalDiedAIDSAfterTasP = 0; DiedNat = 0; Time_All_Treated = 3*tfinal, Last_Time_Not_All_Treated = 0;
       new_infections=0; new_infections_u50 = 0; new_infections_u25=0; new_infections_Mid = 0; new_infections_o50 = 0; 
       susc_days = 0; susc_days_u50 = 0; susc_days_u25 = 0; susc_days_Mid = 0; susc_days_o50 = 0;
       tx_person_days = 0, cd4_cat1_person_days = 0, cd4_cat2_person_days = 0, cd4_cat3_person_days = 0, cd4_cat4_person_days = 0, cd4_cat5_person_days = 0;; 
       cd4InitProbs(); 
       CreateCD4Table(); 
       InitialAgeDistribution();
       InitializePopulation();  
       CalculateLifeExpectancies();
       DefineSexualContactNetwork();
       AddInfecteds();
       UpdateViralLoads();
       //PrintEdgeList(); PrintFirstAndLastContacts();
       pill_count = 0; NonTargetedTxYear1 = 0; TotalTreatedYear1 = 0; NonTargetedTx = 0; TotalTreated = 0; 
       PrintVLdistribution();
       PrintAgeDistributions();
       PrintAgeMatches();
       //PrintStats();
       //PrintTrackedAgents();
       AverageDegree = OrigAverageDegree;
       /* Time Loop */
       do {
         if (VL_print_count >= VL_print_time) {
            PrintVLdistribution();
            VL_print_count = 0;
         }
         if (time == Start_TasP_Campaign) PrintPartnershipDistributions();
         if (time == Start_SexReduction_Campaign) AverageDegree = AverageDegree/Reduction_Mean_Degree;
         if ((time == 365) || (time == 10*365) || (time == 20*365) || (time == 30*365) || (time==40*365)) {
           PrintAgeDistributions();
           PrintAgeMatches();
         }
         if (time <= tfinal - 5*365) {
           susc_days_last_5_years = 0;
           susc_days_last_5_years_u50 = 0;
           new_infections_last_5_years = 0;
           new_infections_last_5_years_u50 = 0;
         }
         VL_print_count++;
         CollectStatistics();
	 HIVTesting();
         if (gradual_tx_incr == 1) InitiateTreatmentGradual();
                              else InitiateTreatment();
         TreatmentDropout();
         UpdateCD4CountsStochastic();
         UpdateViralLoads();
         SimulateTransmission(); 
         Births(); //ConstantInputNewSusceptibles();
         Aging();
         NaturalDeath();
         DeathOfAIDSPatients();
         UpdateDALYs();
         RemoveLinksToDeadPeople();
         SimulatePartnerShipDissolution(); 
         SimulateNewPartnershipFormation();
         SimulateVaccineDecay();
         CheckEdgeListForImpossibleValues();  // Does some elementary error checking
         if ( time % print_frequency == 0) PrintStats();
         PrintHeritabilityStats();
         PrintTrackedAgents();
         time = time + 1;
         
         if (time == Start_Condom_Campaign) prob_sex = orig_prob_sex * (1.0-Percent_using_condums_after_condom_campaign);
         fflush(stdout);
       } while ( (time <= tfinal) && (StopEarly == 0)); // End time loop
       if (StopEarly == 1) {
         PrintStats();
       }
       if (Start_TasP_Campaign > tfinal) PrintPartnershipDistributions();
      
       //RecordStatusOfLivingHIVPatientsForRegistryFile(); PrintFirstAndLastContacts();
       if (StopEarly == 1) {
         printf("Note: Program halted early -- most likely because of difficulties in finding partners using current search algorithm.  Consider increasing MaxLinks.\n");
       }
       if (num_cases_no_link_found > 0) if (plt_warn==1) fprintf(WarningsFile,"Warning at time %ld: Number of times no link was found = %ld\n",time, num_cases_no_link_found);
       if (TotalTreatedYear1 > 0) ActualNonTargetedYear1 = ((double) NonTargetedTxYear1) / ((double) TotalTreatedYear1); 
                             else ActualNonTargetedYear1 = 0.0;
       if (plt_GS==1) fprintf(GrandSummaryFile,"%d\t%ld\t%4.3lf\t%4.3lf\t%4.3lf\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%lf\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
                              TasP_Strategy,under25,perc_tx_start,ActTxStartTasP,yearly_incr_tx,repl,TotalDiedAIDSAfterTasP,Infected,Time_All_Treated,Last_Time_Not_All_Treated,
                              Infected - InfectedTreated,ActualNonTargetedYear1,pill_count,susc_days_last_5_years,new_infections_last_5_years,
                              cd4_cat1_person_days, cd4_cat2_person_days, cd4_cat3_person_days, cd4_cat4_person_days, cd4_cat5_person_days, tx_person_days,
                              DALY00, DALY01, DALY02, DALY03, DALY05, DALY07, DALY08, DALY10, DALY15, DALY20); 
     } // Replicates
   } // tx percent
 } // Strategies

 printf("DONE (New) !!! %c%c",7,7); // print bell sounds (ascii character 7)
 
} // Main


void OpenFiles()
{
 if ( (ParFile = fopen("Parameters.txt","r")) == NULL) { printf("Cannot open file 'Parameters.txt'.  Hit ctrl-C to quit : "); scanf("%s",junk); }
 if ( (ParOut = fopen("ParameterReguritation.txt","w")) == NULL) { printf("Cannot open file 'ParameterReguritation.txt'.  Hit ctrl-C to quit : "); scanf("%s",junk); }
 if ( (OutputFile = fopen("ViralLoadOutput.txt","w")) == NULL) { printf("Cannot open file 'ViralLoadOutput.txt'.  Hit ctrl-C to quit : "); scanf("%s",junk); }
 if ( (NumPartnersFile = fopen("NumPartners.txt","w")) == NULL) { printf("Cannot open file 'NumPartners.txt'.  Hit ctrl-C to quit : "); scanf("%s",junk); }
 if ( (HeritOutputFile = fopen("HeritabilityOutput.txt","w")) == NULL) { printf("Cannot open file 'HeritabilityOutput.txt'.  Hit ctrl-C to quit : "); scanf("%s",junk); }
 if ( (PatientRegistryOutputFile = fopen("PatientRegistryOutput.txt","w")) == NULL) { printf("Cannot open file 'PatientRegistryOutput.txt'.  Hit ctrl-C to quit : "); scanf("%s",junk); }
 if ( (VLOutputFile = fopen("VLHistogram.txt","w")) == NULL) { printf("Cannot open file 'VLOutput.txt'.  Hit ctrl-C to quit : "); scanf("%s",junk); }
 if ( (spVLOutputFile = fopen("spVLHistogram.txt","w")) == NULL) { printf("Cannot open file 'spVLOutput.txt'.  Hit ctrl-C to quit : "); scanf("%s",junk); }
 if ( (NetworkStatsOutputFile = fopen("NetworkStats.txt","w")) == NULL) { printf("Cannot open file 'NetworkStats.txt'.  Hit ctrl-C to quit : "); scanf("%s",junk); }
 if ( (AgeDistFile = fopen("AgeDistribution.txt","w")) == NULL) { printf("Cannot open file 'AgeDistribution.txt'.  Hit ctrl-C to quit : "); scanf("%s",junk); }
 if ( (RelLengthFile = fopen("RelationshipLengths.txt","w")) == NULL) { printf("Cannot open file 'RelationshipLenghts.txt'.  Hit ctrl-C to quit : "); scanf("%s",junk); }
 if ( (AgeMatchFile = fopen("AgeMatches.txt","w")) == NULL) { printf("Cannot open file 'AgeMatches.txt'.  Hit ctrl-C to quit : "); scanf("%s",junk); }
 if ( (RelationshipRegistryFile = fopen("RelationshipRegistry.txt","w")) == NULL) { printf("Cannot open file 'RelationshipRegistry.txt'.  Hit ctrl-C to quit : "); scanf("%s",junk); }
 if ( (GrandSummaryFile = fopen("GrandSummary.txt","w")) == NULL) { printf("Cannot open file 'GrandSummary.txt'.  Hit ctrl-C to quit : "); scanf("%s",junk); }
 if ( (AgentHistoryFile = fopen("AgentHistory.txt","w")) == NULL) { printf("Cannot open file 'AgentHistory.txt'.  Hit ctrl-C to quit : "); scanf("%s",junk); }
 if ( (WarningsFile = fopen("Warnings.txt","w")) == NULL) { printf("Cannot open file 'Warnings.txt'.  Hit ctrl-C to quit : "); scanf("%s",junk); }
}

void PrintTrackedAgents()
{ 
  long i;
  if (plt_AH == 1) {
    for (i=1; i<=N; i++) {
      if (Status[i] == 1) {
        if (i >= TrackedAgentsGroup1start && i <= TrackedAgentsGroup1end)
           fprintf(AgentHistoryFile,"%ld\t%ld\t%lf\t%d\t%lf\t%lf\t%ld\t%ld\t%ld\t%lf\t%lf\t%lf\t%d\n",time,i, SetPoint[i], Status[i], Time_Inf[i], s[i], Treated[i], tx_days[i], tx_days_aids[i], Age[i], ProbDropout[i], V[i],CD4[i]);
        if (i >= TrackedAgentsGroup2start && i <= TrackedAgentsGroup2end)
           fprintf(AgentHistoryFile,"%ld\t%ld\t%lf\t%d\t%lf\t%lf\t%ld\t%ld\t%ld\t%lf\t%lf\t%lf\t%d\n",time,i, SetPoint[i], Status[i], Time_Inf[i], s[i], Treated[i], tx_days[i], tx_days_aids[i], Age[i], ProbDropout[i], V[i],CD4[i]);
      }
    }
  }
}

void PrintAgeDistributions()
{
  long i, age, approx_age;
  long AgeDist[102]; 
  //printf("Entering PrintAgeDistributions\n");
  for (age=1; age<=100; age++) AgeDist[age] = 0; 
  for (i=1; i<=N; i++) {
    if (Status[i] >= 0) {
      approx_age = (long) Age[i];
      if ((approx_age < 0) || (approx_age > 100)) {
        printf("Problem in PrintAgeDistributions: approx_age = %ld, i = %ld\n",approx_age,i);
      }
      AgeDist[approx_age] = AgeDist[approx_age] + 1; 
    }
  }
  if (time <= 1) {
    for (age=16; age<=100; age++) if (plt_AD==1) fprintf(AgeDistFile,"\t%ld",age);
    if (plt_AD==1) fprintf(AgeDistFile,"\n");
  }
  if (plt_AD==1) fprintf(AgeDistFile,"%ld",time/365);
  for (age=16; age<=100; age++) if (plt_AD==1) fprintf(AgeDistFile,"\t%ld",AgeDist[age]);
  if (plt_AD==1) fprintf(AgeDistFile,"\n");
  //printf("Exiting PrintAgeDistributions\n");
}

void PrintPartnershipDistributions() 
{
  long i, j, nLinks;
  long num_partners[21];
  long num_partners_women[21];
  long num_partners_men[21];
  long num_partners_concurrent_women[21];
  long num_partners_concurrent_men[21];
 
  for (i = 0; i<= 20; i++) {
    num_partners[i] = 0;
    num_partners_women[i] = 0;
    num_partners_men[i] = 0;
    num_partners_concurrent_women[i] = 0;
    num_partners_concurrent_men[i] = 0;
  }

  for (i = 1; i<= N; i++) {
    if (Status[i] >= 0) {
      nLinks = NumLinks[i];
      if (nLinks > 20) nLinks = 20;
      num_partners[nLinks] = num_partners[nLinks] + 1;
      if (Sex[i] == 1) num_partners_women[nLinks] = num_partners_women[nLinks] + 1;
      if (Sex[i] == 2) num_partners_men[nLinks] = num_partners_men[nLinks] + 1;
      if (Sex[i] == 1 && Concurrent[i] == 1) num_partners_concurrent_women[nLinks] = num_partners_concurrent_women[nLinks] + 1;
      if (Sex[i] == 2 && Concurrent[i] == 1) num_partners_concurrent_men[nLinks] = num_partners_concurrent_men[nLinks] + 1;
    }
  }

  for (i = 0; i<=20; i++) fprintf(NumPartnersFile,"\t%ld",num_partners[i]); fprintf(NumPartnersFile,"\n");
  for (i = 0; i<=20; i++) fprintf(NumPartnersFile,"\t%ld",num_partners_women[i]); fprintf(NumPartnersFile,"\n");
  for (i = 0; i<=20; i++) fprintf(NumPartnersFile,"\t%ld",num_partners_men[i]); fprintf(NumPartnersFile,"\n");
  for (i = 0; i<=20; i++) fprintf(NumPartnersFile,"\t%ld",num_partners_concurrent_women[i]); fprintf(NumPartnersFile,"\n");
  for (i = 0; i<=20; i++) fprintf(NumPartnersFile,"\t%ld",num_partners_concurrent_men[i]); fprintf(NumPartnersFile,"\n");
} 

void PrintAgeMatches()
{
  long i, j;
  for (i=1; i<=N; i++) {
    if (Status[i] >= 0) {
      for (j=1; j<= NumLinks[i]; j++) {
        if (Status[i] >= 0 && Status[Links[i][j]] >= 0) {
          if (Sex[i] == 1) {
            if (plt_AM==1) fprintf(AgeMatchFile,"%ld\t%lf\t%lf\t%ld\t%ld\n",time/365,Age[i],Age[Links[i][j]],Sex[i],Sex[Links[i][j]]);
            if (Sex[Links[i][j]] == 1) if (plt_warn==1) fprintf(WarningsFile,"Warning at time %ld: Sex[%ld] = %ld linked to Sex[%ld] = %ld\n",time,i,Sex[i],Links[i][j],Sex[Links[i][j]]);
          } else {
            if (plt_AM==1) fprintf(AgeMatchFile,"%ld\t%lf\t%lf\t%ld\t%ld\n",time/365,Age[Links[i][j]],Age[i],Sex[Links[i][j]],Sex[i]);
            if (Sex[Links[i][j]] == 2) if (plt_warn==1) fprintf(WarningsFile,"Warning at time %ld: Sex[%ld] = %ld linked to Sex[%ld] = %ld\n",time,i,Sex[i],Links[i][j],Sex[Links[i][j]]);
          }
        }
      }
    }
  }
  //printf("Exiting PrintAgeMatches\n");
}

void PrintVLdistribution()
{
 double current_VL;
 long i;
 long Print_Array_VL[1000];
 long Print_Array_spVL[1000];
 long bin_count = 0;

 /* Print header giving the midpoint VL (on log scale) for each element in the histogram at time 0 */
 if (time == 0) {
   if (plt_VL==1) fprintf(VLOutputFile,"VL\t");
   if (plt_spVL==1) fprintf(spVLOutputFile,"spVL\t");
   for (current_VL = VL_print_lower; current_VL <= VL_print_upper; current_VL = current_VL + VL_print_interval) {
       if (plt_VL==1) fprintf(VLOutputFile,"%3.2lf\t", current_VL+VL_print_interval/2.0);
       if (plt_spVL==1) fprintf(spVLOutputFile,"%3.2lf\t", current_VL+VL_print_interval/2.0);
   }
   if (plt_VL==1) fprintf(VLOutputFile,"\n");
   if (plt_spVL==1) fprintf(spVLOutputFile,"\n");
 }


 /* Zero out elements of the histogram */
 bin_count = 0;
 for (current_VL = VL_print_lower; current_VL <= VL_print_upper; current_VL = current_VL + VL_print_interval) {
   bin_count++;
   Print_Array_VL[bin_count] = 0;
   Print_Array_spVL[bin_count] = 0;
 }

 /* For each infected person, increment histogram cell spanning that person's viral load */
 for (i=1; i<=N; i++) {
   if (Status[i] == 1) {
     bin_count = 0;
     for (current_VL = VL_print_lower; current_VL <= VL_print_upper; current_VL = current_VL + VL_print_interval) {
       bin_count++;
       if ( (log10(V[i]) >= current_VL) && (log10(V[i]) < current_VL+VL_print_interval)) {
         //if (time - Time_Inf[i] > 50.0) {
           Print_Array_VL[bin_count] = Print_Array_VL[bin_count] + 1;
         //}
       }
       if ( (LogSetPoint[i] >= current_VL) && (LogSetPoint[i] < current_VL+VL_print_interval)) {
         Print_Array_spVL[bin_count] = Print_Array_spVL[bin_count] + 1;
       }
     }
   }
 }

 /* Print histogram to a file */
 bin_count = 0;
 if (plt_VL==1) fprintf(VLOutputFile,"%ld\t", time);
 if (plt_spVL==1) fprintf(spVLOutputFile,"%ld\t", time);
 for (current_VL = VL_print_lower; current_VL <= VL_print_upper; current_VL = current_VL + VL_print_interval) {
   bin_count++;
   if (plt_VL==1) fprintf(VLOutputFile,"%ld\t", Print_Array_VL[bin_count]);
   if (plt_spVL==1) fprintf(spVLOutputFile,"%ld\t", Print_Array_spVL[bin_count]);
 }
 if (plt_VL==1) fprintf(VLOutputFile,"\n");
 if (plt_spVL==1) fprintf(spVLOutputFile,"\n");

}

void HIVTesting() {
  long i;
  for (i = 1; i <=N; i++) {
    if (Status[i] == 1 && (time-Time_Inf[i] > delay_test_antibodies) && Diagnosed[i] == 0) { // Assume ~35 days until sufficient antibodies for testing 
       if (genrand_real2() < daily_prob_diagnosis) {
         Diagnosed[i] = 1;
		 TimeDiag[i] = time;
       }
    }
  }

}
void VaccinatePopulation() {
  long i;
  for (i = 1; i <=N; i++) {
    if ((Status[i] == 0) && (Vaccinated[i]==0)) { // Questions about second condition (not vaccinating those already vaccinated)
                                                  // Could get weird if vaccine compaigns occure before vaccine wears off.
                                                  // Current algorithm will work if vaccine campaigns are longer than vaccine duration
       if (genrand_real2() < PercentVaccinated) {
         Vaccinated[i] = 1;
         LengthVaccinated[i] = 0;
       }
    }
  }

}

void SimulateVaccineDecay() {
  long i;
  for (i = 1; i <=N; i++) {
    if (Vaccinated[i]==1) {
      LengthVaccinated[i] = LengthVaccinated[i] + 1;
      if (LengthVaccinated[i] > VaccineDuration) {
        Vaccinated[i] = 0;
      }
    }
  }

}
void InitializePopulation() {
  /* Initalizes the s[], Status[], Numrecepients[], Duration[], UnderCare[], Sex[], Age[], and Treated[] arrays to their default values */
  long i;
  N = N0;
  for (i=1;i<=200000;i++) {
     NumLinks[i] = 0;
     s[i] = prog_rate; // Same for all patients in this version.  Could be made to vary from patient to patient in future studies
     Status[i] = 0;
     CD4[i] = 0;
     CD4_nadir[i] = 0;
     V[i] = 0.0; I[i] = 0.0; M[i] = 0.0; L[i] = 0.0; 
     CD4_TimeToAIDS[i]= 0;
     NumRecipients[i] = 0;
	 Diagnosed[i] = 0;
     Vaccinated[i] = 0;
     Circumcised[i] = 0;
     NumPartners[i] = 0;
     NumPartnersLastTime[i] = 0;
     ResistantToVaccine[i] = 0;
     LengthVaccinated[i] = 0;
     VL_Start_AIDS[i] = 0.0;
     Time_Start_AIDS[i] = 100000;
     Time_AIDS_Death[i] = 0;
     TimeLastBreakUp[i] = -100000;
     TimeInAIDS[i] = 0;
     AgeAIDS[i] = 0;
     NumLinks[i] = 0;
     DALY00 = 0.0, DALY01 = 0.0, DALY02 = 0.0, DALY03 = 0.0, DALY04 = 0.0, DALY05 = 0.0, DALY06 = 0.0, DALY07 = 0.0, DALY08 = 0.0, DALY10 = 0.0, DALY15 = 0.0, DALY20 = 0.0;
     for (j=1;j<=MaxLinks;j++) {
       Links[i][j] = 0;
       TimeLinked[i][j] = 0;
     }
  }
  for (i=1;i<=N;i++) {
     UsualInitialValues(i); // This does not set ages.:
  }
}

//double increased_duration_per_year = 50.0; // In days.  If 50, duration will increase 500 days in 10 years, 2500 days (=6 years) in 50 years

void UsualInitialValues (long i) {
  double Time_Concur_halves;
  if (rate_stops_being_concurrent > 0) Time_Concur_halves = 1.0/rate_stops_being_concurrent;
                                  else Time_Concur_halves = 1000000.0;
  if (genrand_real2() < 0.5) { 
    Sex[i] = 1; // Female
  } else {
    Sex[i] = 2; // Male
  }
  if (genrand_real2() < percent_high_risk) HighRisk[i] = 1;
                                      else HighRisk[i] = 0;
  if (Sex[i] == 1) {
    if (genrand_real2() < female_concurrency) Concurrent[i] = 1;
                                         else Concurrent[i] = 0;
  }
  if (Sex[i] == 2) {
    if (genrand_real2() < circum_prob) Circumcised[i] = 1;
    if (genrand_real2() < male_concurrency) Concurrent[i] = 1;
                                       else Concurrent[i] = 0;
  }
  if (Concurrent[i] == 1 && Age[i] > 45.0) {
    if (genrand_real2() < ((Age[i] - 45.0)*365.0)/((Age[i] - 45.0)*365 + Time_Concur_halves)) Concurrent[i] = 0;
  }
  ProbDropout[i] = 2.0*genrand_real2()*prob_dropout/365.0;
  NumLinks[i] = 0;
  Status[i] = 0;
  Diagnosed[i] = 0;
  TimeDiag[i] = 10000000;
  Duration[i] = genrand_real2()*2.0*MinDuration + (Age[i]-16.0)*increased_duration_per_year; // Distribution between 0 and 2*Min
  if (Duration[i] < 1.0) Duration[i] = 1.0;
  UnderCare[i] = 0; // Set to zero to start.  This will be updated in AddInfecteds() below
  Treated[i] = 0; // Assume no treatment for initial population
  Time_Treated[i] = 0; // Assume no treatment for initial population
  tx_days[i] = 0; // Assume no treatment for initial population
  tx_days_aids[i] = 0; // Assume no treatment for initial population
  ever_aids[i] = 0; // Assume no treatment for initial population
  TargetedTreated[i] = 0; // Assume no treatment for initial population
}

void Births()
{
   double per_day_birth_rate;
   long nBirths;
   double NewSusceptiblesPerDayDouble;
   long NewSusceptiblesPerDay;
   long i;
   newbirths = 0;
   if (MaxNReached ==0) {
     per_day_birth_rate = BirthRate*exp(PopGrowth*time/365);
     nBirths = rpois(per_day_birth_rate);
     NewSusceptiblesPerDay = nBirths;
     for (i=1;i<=NewSusceptiblesPerDay;i++) {
       if (N < MaxN) {
          N++;
          newbirths++;
          Age[N] = 16.0;
          UsualInitialValues(N);
       } else {
         if (MaxNReached ==0) {
            if (plt_warn==1) fprintf(WarningsFile,"Warning at time %ld: MaxN reached.  No more births or reincarnation events allowed\n",time);
            MaxNReached = 1;
            tMaxNReached = time;
         } // MaxNReached
       }  // else part of N < MaxN
     } // for i = 1 to N
   } // MaxNReached == 0
}

void Aging()
{
  long i;
  for (i=1; i<=N; i++) {
    if (Status[i] >= 0) {
      Age[i] = Age[i] + 1.0/365.0;
      Duration[i] = Duration[i] +  increased_duration_per_year * (1/365.0);
      if (Duration[i] >  MaxDuration) Duration[i] = MaxDuration;
      if (Concurrent[i] == 1) {
        if (genrand_real2() < rate_stops_being_concurrent) {
          Concurrent[i] = 0;
        }
      }
    }
    if (Age[i] > 100) Age[i] = 100;
  } 
}

void CalculateLifeExpectancies()
{
  long i, cur_age;
  double prob_live_to[101][101];
  printf("asmr[16] = %lf\n",asmr[16]);
  for (cur_age = 15; cur_age <= 100; cur_age++) {
     prob_live_to[cur_age][cur_age] = (1-0.5*asmr[cur_age]*365.0);
     // printf("%ld : ",cur_age);
     for (i=cur_age+1; i <= 99; i++) {
       prob_live_to[cur_age][i] = prob_live_to[cur_age][i-1]*(1-asmr[i]*365.0);
       // printf("%3.2lf ",prob_live_to[cur_age][i]);
     }
     prob_live_to[cur_age][100] = 0.0;
     //printf("\n");
  }
  for (cur_age = 15; cur_age <= 100; cur_age++) {
     LifeExpectancy[cur_age] = cur_age + 0.5; // Assuming prob_live_to[cur_age][cur_age+1] == 0, agent would live on average 0.5 years, hence the 0.5 term 
     for (i= cur_age+1; i<= 100; i++) {
       LifeExpectancy[cur_age] = LifeExpectancy[cur_age] + prob_live_to[cur_age][i];
     }
     //printf("Life Expectancy for someone age %ld is %lf\n",cur_age,LifeExpectancy[cur_age]);
  }
}

void UpdateDALYs()
{
  long i;
  double relevant_cost;
  for (i=1; i<=N; i++) {
    if (time > Start_TasP_Campaign) {
      if (Status[i] == 1) {
        relevant_cost =0.0;
        if (Treated[i]==1) relevant_cost = cost_hiv_treated/365.0;
        if (CD4[i] == 2 && Treated[i]==0) relevant_cost = cost_cd4_gt_350/365.0;
        if (CD4[i] == 3 && Treated[i]==0) relevant_cost = cost_cd4_200_350/365.0;
        if (CD4[i] == 4 && Treated[i]==0) relevant_cost = cost_cd4_lt_200/365.0;
        DALY00 = DALY00 + relevant_cost;
        DALY01 = DALY01 + relevant_cost * pow(0.99, (time-Start_TasP_Campaign)/365.0);
        DALY02 = DALY02 + relevant_cost * pow(0.98, (time-Start_TasP_Campaign)/365.0);
        DALY03 = DALY03 + relevant_cost * pow(0.97, (time-Start_TasP_Campaign)/365.0);
        DALY05 = DALY05 + relevant_cost * pow(0.95, (time-Start_TasP_Campaign)/365.0);
        DALY07 = DALY07 + relevant_cost * pow(0.93, (time-Start_TasP_Campaign)/365.0);
        DALY08 = DALY08 + relevant_cost * pow(0.92, (time-Start_TasP_Campaign)/365.0);
        DALY10 = DALY10 + relevant_cost * pow(0.90, (time-Start_TasP_Campaign)/365.0);
        DALY15 = DALY15 + relevant_cost * pow(0.85, (time-Start_TasP_Campaign)/365.0);
        DALY20 = DALY20 + relevant_cost * pow(0.80, (time-Start_TasP_Campaign)/365.0);
      }
      if (Status[i] == -2 & Time_AIDS_Death[i] == time) {
        relevant_cost = cost_died_AIDS * (LifeExpectancy[(long) Age[i]] - Age[i]);
        DALY00 = DALY00 + relevant_cost;
        DALY01 = DALY01 + relevant_cost * pow(0.99, (time-Start_TasP_Campaign)/365.0);
        DALY02 = DALY02 + relevant_cost * pow(0.98, (time-Start_TasP_Campaign)/365.0);
        DALY03 = DALY03 + relevant_cost * pow(0.97, (time-Start_TasP_Campaign)/365.0);
        DALY05 = DALY05 + relevant_cost * pow(0.95, (time-Start_TasP_Campaign)/365.0);
        DALY07 = DALY07 + relevant_cost * pow(0.93, (time-Start_TasP_Campaign)/365.0);
        DALY08 = DALY08 + relevant_cost * pow(0.92, (time-Start_TasP_Campaign)/365.0);
        DALY10 = DALY10 + relevant_cost * pow(0.90, (time-Start_TasP_Campaign)/365.0);
        DALY15 = DALY15 + relevant_cost * pow(0.85, (time-Start_TasP_Campaign)/365.0);
        DALY20 = DALY20 + relevant_cost * pow(0.80, (time-Start_TasP_Campaign)/365.0);
      } 
    }  
  }
}

long rpois(double lambda)
{
  double L, p, u;
  long k;
  L = exp(-lambda);
  k = 0;
  p = 1.0;
  do {
    k = k + 1.0;
    u = genrand_real2();
    p = p * u;
  } while (p > L);
  return(k-1);
}


void InitialAgeDistribution() {
  long jj, low_limit, high_limit, found;
  double a, b, r, sum_age_vec;
  double age_vec[101], rel_age_dist[101], cum_age_dist[101];
  double randnum;
  b = BirthRate; // Initial birth rate
  a = 1.0/365.0; //Per day aging rate
  r = PopGrowth/365.0; // pop_growth_rate_timestep;
 // From WHO http://apps.who.int/gho/data/?theme=main&vid=61540 (first number is male. Downloaded: 7/2/18) 
  for (i=15; i<=19;i++) asmr[i] = 0.5 * (0.002 + 0.001) / 365.0;
  for (i=20; i<=24;i++) asmr[i] = 0.5 * (0.002 + 0.002) / 365.0;
  for (i=25; i<=29;i++) asmr[i] = 0.5 * (0.004 + 0.003) / 365.0;
  for (i=30; i<=34;i++) asmr[i] = 0.5 * (0.006 + 0.005) / 365.0;
  for (i=35; i<=39;i++) asmr[i] = 0.5 * (0.009 + 0.008) / 365.0;
  for (i=40; i<=44;i++) asmr[i] = 0.5 * (0.011 + 0.007) / 365.0;
  for (i=45; i<=49;i++) asmr[i] = 0.5 * (0.013 + 0.008) / 365.0;
  for (i=50; i<=54;i++) asmr[i] = 0.5 * (0.018 + 0.009) / 365.0;
  for (i=55; i<=59;i++) asmr[i] = 0.5 * (0.024 + 0.013) / 365.0;
  for (i=60; i<=64;i++) asmr[i] = 0.5 * (0.037 + 0.019) / 365.0;
  for (i=65; i<=69;i++) asmr[i] = 0.5 * (0.054 + 0.028) / 365.0;
  for (i=70; i<=74;i++) asmr[i] = 0.5 * (0.077 + 0.042) / 365.0;
  for (i=75; i<=80;i++) asmr[i] = 0.5 * (0.108 + 0.062) / 365.0;
  for (i=80; i<=85;i++) asmr[i] = 0.5 * (0.150 + 0.101) / 365.0;
  for (i=85; i<=99;i++) asmr[i] = 0.5 * (0.237 + 0.184) / 365.0;
  asmr[100] =  1.0;

  /* female + male
  for (i=16; i<=19;i++) asmr[i] = 0.5 * (0.0013 + 0.0018) / 365.0;
  for (i=20; i<=24;i++) asmr[i] = 0.5 * (0.0035 + 0.0039) / 365.0;
  for (i=25; i<=29;i++) asmr[i] = 0.5 * (0.0072 + 0.0071) / 365.0;
  for (i=30; i<=34;i++) asmr[i] = 0.5 * (0.0113 + 0.0117) / 365.0;
  for (i=35; i<=39;i++) asmr[i] = 0.5 * (0.0130 + 0.0157) / 365.0;
  for (i=40; i<=44;i++) asmr[i] = 0.5 * (0.0129 + 0.0185) / 365.0;
  for (i=45; i<=49;i++) asmr[i] = 0.5 * (0.0127 + 0.0197) / 365.0;
  for (i=50; i<=54;i++) asmr[i] = 0.5 * (0.0125 + 0.0197) / 365.0;
  for (i=55; i<=59;i++) asmr[i] = 0.5 * (0.0127 + 0.0219) / 365.0;
  for (i=60; i<=64;i++) asmr[i] = 0.5 * (0.0194 + 0.0325) / 365.0;
  for (i=65; i<=69;i++) asmr[i] = 0.5 * (0.0269 + 0.0441) / 365.0;
  for (i=70; i<=74;i++) asmr[i] = 0.5 * (0.0379 + 0.0582) / 365.0;
  for (i=75; i<=80;i++) asmr[i] = 0.5 * (0.0563 + 0.0815) / 365.0;
  for (i=80; i<=99;i++) asmr[i] = 0.5 * (0.1403 + 0.1629) / 365.0;
  asmr[100] =  1.0;
  */

/*  From Evonet
 asmr =c(0.0013, 0.0013, 0.0013, 0.0013,      16-19 c
    0.0035, 0.0035, 0.0035, 0.0035, 0.0035,   20-24 c 
    0.0072, 0.0072, 0.0072, 0.0072, 0.0072,   25-29 c 
    0.0113, 0.0113, 0.0113, 0.0113, 0.0113,   30-34 c 
    0.0130, 0.0130, 0.0130, 0.0130, 0.0130,   35-39 c
    0.0129, 0.0129, 0.0129, 0.0129, 0.0129,   40-44 c 
    0.0127, 0.0127, 0.0127, 0.0127, 0.0127,   45-49 c 
    0.0125, 0.0125, 0.0125, 0.0125, 0.0125,   50-54 c 
    0.0127, 0.0127, 0.0127, 0.0127, 0.0127,   55-59 c 
    0.0194, 0.0194, 0.0194, 0.0194, 0.0194,   60-64 c 
    0.0269, 0.0269, 0.0269, 0.0269, 0.0269,   65-69 c
    0.0379, 0.0379, 0.0379, 0.0379, 0.0379,   70-74 c 
    0.0563, 0.0563, 0.0563, 0.0563, 0.0563,   75-79 c 
    0.1403, 0.1403, 0.1403, 0.1403, 0.1403,   80-84 c 
    0.1403, 0.1403, 0.1403, 0.1403, 0.1403,   85-89 c
    0.1403, 0.1403, 0.1403, 0.1403, 0.1403,   90-94 c 
    0.1403, 0.1403, 0.1403, 0.1403, 0.1403,   95-99 c
    0.1403                                    100

 asmr =c(0.0018, 0.0018, 0.0018, 0.0018,
    0.0039, 0.0039, 0.0039, 0.0039, 0.0039,
    0.0071, 0.0071, 0.0071, 0.0071, 0.0071, 
    0.0117, 0.0117, 0.0117, 0.0117, 0.0117, 
    0.0157, 0.0157, 0.0157, 0.0157, 0.0157, 
    0.0185, 0.0185, 0.0185, 0.0185, 0.0185, 
    0.0197, 0.0197, 0.0197, 0.0197, 0.0197, 
    0.0197, 0.0197, 0.0197, 0.0197, 0.0197, 
    0.0219, 0.0219, 0.0219, 0.0219, 0.0219, 
    0.0325, 0.0325, 0.0325, 0.0325, 0.0325, 
    0.0441, 0.0441, 0.0441, 0.0441, 0.0441, 
    0.0582, 0.0582, 0.0582, 0.0582, 0.0582, 
    0.0815, 0.0815, 0.0815, 0.0815, 0.0815, 
    0.1629, 0.1629, 0.1629, 0.1629, 0.1629, 
    0.1629, 0.1629, 0.1629, 0.1629, 0.1629, 
    0.1629, 0.1629, 0.1629, 0.1629, 0.1629, 
    0.1629, 0.1629, 0.1629, 0.1629, 0.1629, 
    0.1629),age_range=c(16,100))
*/
 
  for (jj = 1; jj <= 15; jj++) {
    age_vec[jj] = 0.0;
  }
  age_vec[16] = b/(a + r + asmr[16]);
  sum_age_vec = age_vec[16];
  age_vec[100] = 0.0; // Hard code 100% death at age 100
  for (jj = 17; jj<= 99; jj++) {
     age_vec[jj] = age_vec[jj-1] * a/(a + r + asmr[jj]);
     sum_age_vec = sum_age_vec + age_vec[jj];
  }
  for (jj = 1; jj<= 100; jj++) {
    rel_age_dist[jj] = age_vec[jj]/sum_age_vec;
  }
  //printf("age \t asmr \t age_vec \t Rel_prob \t CumProb \t(sum_age_vec = %lf, a = %lf, r = %lf)\n",sum_age_vec,a,r);
  cum_age_dist[1] = 0.0;
  for (jj = 2; jj<= 100; jj++) {
    cum_age_dist[jj] = rel_age_dist[jj] + cum_age_dist[jj-1];
    //printf("%ld \t %lf \t %lf \t %lf \t %lf (sum_age_vec = %lf)\n",jj,asmr[jj],age_vec[jj],rel_age_dist[jj],cum_age_dist[jj],sum_age_vec);
  }
  
  for (i=1;i<=N0;i++) {
    randnum = genrand_real2();
    found = 0;
    //printf("Determining age for agent %ld (of %ld), randnum = %lf\n",i,N0,randnum);
    for (jj=16; jj<= 100;jj++) {
      if ( (randnum <= cum_age_dist[jj]) && (found==0)) {
        found = 1;
        Age[i] = ((double) jj) +  genrand_real2();
        //printf("     agent %ld assigned to age %lf\n",i,Age[i]);
      }
    }
    if (found == 0) {
      printf("Problem in InitialAgeDistribution: No age assigned to agent %ld\n.",i);
      printf("Type ctrl-c to stop: ");
      scanf("%s",junk);
    }
  }
    
}

/* 
Code that worked in stand-alone script 
#   b = (a + r + d[1])*N[1]
#   for (age in 2: 83) {
#     N[age] <- N[age-1]*a/(a+r+d[age])
#   }
#  where N1 = lowest age class (=18 in most of our models)
#  a = aging rate
#  di = death rate of persons of age i
#  r = rate of growth in the birth rate
#  b = initial birth rate (rate of entry of persons of "age 1") 
#  Equations
#   dN1/dt = b*exp(r*t) - a*N1 - d1*N1
#   dN2/dt = a*N1 - a*N2 - d2*N2
#   dNi/dt = a*Ni-1 - a*Ni - di*Ni
# Equation only works when population is at steady state, meaning that b = (a+r+d[1])*N[1]
# In other words "baseline_input_exp_growth" (=b) needs to be tuned so that rate of initial input into the youngest
# age class equals the initial death and ageing rates of the youngest age class

b <- baseline_input_exp_growth  # Initial birth rate
  a <- 1/365  # Per day aging rate
  r <- pop_growth_rate_timestep
  d_f <- mort_per_timestep_female # Per day death rate (assume females and males the same)
  age_vec <- c(min_age:(max_age-1))  # Set up the vector (arbitrary numbers for now)
  age_vec[1] <- b/(a + r + d_f[min_age])
  low_limit <- min_age +1
upper_limit <- max_age -1
for (age in low_limit:upper_limit) {
  age_index <- age - min_age +1 # Youngest age class (e.g., 16) is represented as 1 in our model
  age_vec[age_index] <- age_vec[age_index-1]*a/(a +r + d_f[age_index])
}
age_vec[1] <- age_vec[1]/2  # Reduce youngest class by 50% to reflect fact that the average age of newly entering agents is min_age + 0.5.  
final_age_dist <- age_vec/sum(age_vec)
*/


void AddInfecteds()
{
 long i,ii, Reject_i, escape_counter = 0;
 double rand_num;

 Infected = Infected0; 
 for (ii=1;ii<=Infected0;ii++) {
   
   Reject_i = 0; escape_counter = 0;
   do {
     Reject_i = 0; // Assume to be okay until shown otherwise
     i = N * genrand_real2() + 1;
     if (Status[i] != 0) Reject_i = 1;
     rand_num = genrand_real2();
     escape_counter++;
     if (escape_counter > 100) {
       printf("Trouble in AddInfecteds: unable to seed population\n");
       printf("ii = %ld, i = %ld, N = %ld, Status[%ld] = %d, Age[%ld] = %lf\n",ii,i,N,i,Status[i],i,Age[i]);
       printf("Reject_i = %ld, escape_counter = %ld, rand_num=%lf\n",Reject_i,escape_counter,rand_num);
       printf("Type ctrl-C to stop:"); scanf("%s",junk);
     }     
   } while (Reject_i == 1);

   if (genrand_real2() < PercentResistantToVaccine) { 
     ResistantToVaccine[i] = 1;
   } else {
     ResistantToVaccine[i] = 0;
   }

   ViralContribToLogSP0[i] = norm_rand(AverageLogSP0, sqrt(H*VarianceLogSP0));
   
   EnvirContribToLogSP0[i] = norm_rand(0, sqrt((1-H)*VarianceLogSP0));
   
   LogSetPoint[i] = ViralContribToLogSP0[i] + EnvirContribToLogSP0[i];
   
   SpvlCat[i]=SPVL_Category(LogSetPoint[i] );
   CD4[i]=initialCD4(SpvlCat[i]);
   CD4_initial_value[i]=CD4[i];
      
   //adding exponential distribution for time in each cd4 category
   //can tighten this up once it gets working....
   CD4_TimeToAIDS_exp[i][1]=log(genrand_real2()) / (-1.0/CD4_lookup[SpvlCat[i]][1]) ;
   CD4_TimeToAIDS_exp[i][2]=log(genrand_real2()) / (-1.0/CD4_lookup[SpvlCat[i]][2]) ;
   CD4_TimeToAIDS_exp[i][3]=log(genrand_real2()) / (-1.0/CD4_lookup[SpvlCat[i]][3]) ;
   CD4_TimeToAIDS_exp[i][4]=log(genrand_real2()) / (-1.0/CD4_lookup[SpvlCat[i]][4]) ;

   if(CD4_Exp_Flag==1) {
     double ll=0.0;
     int kk;
     for(kk=CD4[i];kk<=3;kk++) {
       ll =  ll + CD4_TimeToAIDS_exp[i][kk];
     }
     CD4_TimeToAIDS[i]=ll;
   }

   //printf("%d----%f\n",i,CD4_TimeToAIDS_exp[i][1]);
   //printf("%d -- %d -- %d\n",i,SpvlCat[i], CD4[i]);
   //printf("-- %d --", CD4[i]);

   SetPoint[i] = pow(10.0,LogSetPoint[i]);
   d_acute[i] = 0.4; // Placeholder until exact value is calculated in UpdateViralLoads()
   Status[i] = 1;
   Generation[i] = 1;
   
   if (genrand_real2() < percent_under_care) {
      UnderCare[i] = 1;
   } else {
      UnderCare[i] = 0;
   }

   Time_Inf[i] = 0;
   CD4time[i]=-Time_Inf[i] ; // What does this mean???
   TimeInAIDS[i]= 0;
   V[i] = V0;
   
   Donors_V[i] = -1.0; // Flag absence of information with -1's or zeros.
   Donors_ViralContribToLogSP0[i] = -1.0; // Flag absence of information with -1's or zeros
   Donors_EnvirContribToLogSP0[i] = -1.0; // Flag absence of information with -1's or zeros
   Donors_LogSetPoint[i] = -1.0; // Flag absence of information with -1's or zeros
   Donors_SetPoint[i] = 0.0; // Flag absence of information with -1's or zeros
   Donors_d_acute[i] = -1.0; // Flag absence of information with -1's or zeros
   Donors_Total_Time_Inf_At_Trans[i] = -999.0; // Set unknown donors for the initial infected population to -1 year
   Donors_Generation[i] = 0;
   Donors_Index[i] = 0;
   //printf("%ld\t%ld\t%ld\t%6.5e\t%6.5e\t%6.5e\t%6.5e\t%6.5e\t%6.5e\n",
   //        time,i,Generation[i],V[i],LogSetPoint[i],ViralContribToLogSP0[i],EnvirContribToLogSP0[i],d_acute[i],time-Time_Inf[i]);
 }
}


double norm_rand(double rn_mean, double rn_variance)
{
  double rand1, rand2;
  rand1 = genrand_real2(); rand2 = genrand_real2();
  return(rn_mean + rn_variance*sqrt(-2.0*log(rand1))*cos(2.0*3.14159265358979*rand2));
}

void DefineSexualContactNetwork()
{
 long i,j;
 long Active =0;
 for (i=1;i<=N;i++) {
   NumLinks[i] = 0;
   for (j=1;j<=5;j++) {
     Links[i][j] = 0;
   }
   if (Status[i] >= 0) Active++;
 }
 TotalEdges = (long) ( N * AverageDegree / 2.0 + 0.5);
 //Links_Broken_Per_Day = N/(30);
 //Links_Broken_Per_Day = N/(3*365);
  AddNewLinks(TotalEdges);
}


/*void PrintEdgeList()
{
 long i,j;
 for (i=1;i<=N;i++) {
   if (Status[i] == 0) printf("Person %ld is linked to : ",i);
   if (Status[i] == 1) printf("Person %ld* is linked to : ",i);
   for (j=1;j<=NumLinks[i];j++) {
     if (j==1) printf("%ld",Links[i][j]);
     if (j>1) printf(", %ld",Links[i][j]);
     if (Status[Links[i][j]] == 1) printf("*");
   }
   printf("\n");
 }
 fflush(stdout);
//printf("Type 'c' to continue : ");
//scanf("%s",junk);
} */



void PrintFirstAndLastContacts()
{
if(printOutput==1)
{
 printf("At time %ld: the first (6) and last (6) contacts are: \n",time);
 if (N > 20) {
  for (i=1;i<=6;i++) {
   if (Status[i] == 0) printf("Person %ld is linked to : ",i);
   if (Status[i] == 1) printf("Person %ld* is linked to : ",i);
   if (Status[i] <0 ) printf("Person %ld : dead ",i);
   for (j=1;j<=NumLinks[i];j++) {
     if (j==1) printf("%ld",Links[i][j]);
     if (j>1) printf(", %ld",Links[i][j]);
     if (Status[Links[i][j]] == 1) printf("*");
   }
   printf("\n");
  }
  for (i=N/2-3;i<=N/2+3;i++) {
   if (Status[i] == 0) printf("Person %ld is linked to : ",i);
   if (Status[i] == 1) printf("Person %ld* is linked to : ",i);
   if (Status[i] <0 ) printf("Person %ld : dead ",i);
   for (j=1;j<=NumLinks[i];j++) {
     if (j==1) printf("%ld",Links[i][j]);
     if (j>1) printf(", %ld",Links[i][j]);
     if (Status[Links[i][j]] == 1) printf("*");
   }
   printf("\n");
  }
  for (i=N-6;i<=N;i++) {
   if (Status[i] == 0) printf("Person %ld is linked to : ",i);
   if (Status[i] == 1) printf("Person %ld* is linked to : ",i);
   if (Status[i] <0 ) printf("Person %ld : dead ",i);
   for (j=1;j<=NumLinks[i];j++) {
     if (j==1) printf("%ld",Links[i][j]);
     if (j>1) printf(", %ld",Links[i][j]);
     if (Status[Links[i][j]] == 1) printf("*");
   }
   printf("\n");
  }
 } else {
   for (i=1;i<=N;i++) {
   if (Status[i] == 0) printf("Person %ld is linked to : ",i);
   if (Status[i] == 1) printf("Person %ld* is linked to : ",i);
   if (Status[i] <0 ) printf("Person %ld : dead ",i);
   for (j=1;j<=NumLinks[i];j++) {
     if (j==1) printf("%ld",Links[i][j]);
     if (j>1) printf(", %ld",Links[i][j]);
     if (Status[Links[i][j]] == 1) printf("*");
   }
   printf("\n");
  }

 }
}
}

void CreateCD4Table()
{
   /* Look up table for CD4 counts from Cori,Fraser,Pickles TableS13 Supplement; years spent in cat.*/
   // 3/21/18 -- Updated to match evonet (john)
//logSPVL <= 3
   CD4_lookup[1][1]=6.08;
   CD4_lookup[1][2]=5.01;
   CD4_lookup[1][3]=3.60;
   CD4_lookup[1][4]=4.67;	//15.36 vs 14.51 (old notes not deleted on 3/21/18)
 //3 < logSPVL <= 3.5
   CD4_lookup[2][1]=4.69;
   CD4_lookup[2][2]=2.52;
   CD4_lookup[2][3]=3.68;
   CD4_lookup[2][4]=4.11;	//11.71 vs 11.99  Corrected from [1][4] to [2][4] on 4/8//18
 //3.5 < logSPVL <= 4
   CD4_lookup[3][1]=3.94;
   CD4_lookup[3][2]=4.07;
   CD4_lookup[3][3]=2.38;
   CD4_lookup[3][4]=3.54;	//11.56 vs 10.4
 //4 < logSPVL <= 4.5
   CD4_lookup[4][1]=2.96;
   CD4_lookup[4][2]=3.09;
   CD4_lookup[4][3]=3.81;
   CD4_lookup[4][4]=2.98;	//10.21 vs 9.91
 //4.5 < logSPVL <= 5
   CD4_lookup[5][1]=2.25;
   CD4_lookup[5][2]=2.32;
   CD4_lookup[5][3]=3.21;
   CD4_lookup[5][4]=2.42;	//8.41 vs 8.09
 //5 < logSPVL <= 5.5
   CD4_lookup[6][1]=1.47;
   CD4_lookup[6][2]=1.55;
   CD4_lookup[6][3]=2.27;
   CD4_lookup[6][4]=1.86;	//5.74 vs 5.68
 //5.5 < logSPVL <= 6
   CD4_lookup[7][1]=0.95;
   CD4_lookup[7][2]=1.19;
   CD4_lookup[7][3]=1.00;
   CD4_lookup[7][4]=1.29;	//5.58 vs 3.62
 //6 < logSPVL < 6.5
   CD4_lookup[8][1]=0.32;
   CD4_lookup[8][2]=0.59;
   CD4_lookup[8][3]=0.68;
   CD4_lookup[8][4]=0.73;	//1.51 vs 1.69
  //logSPVL > 6.5
   CD4_lookup[9][1]=0.30;
   CD4_lookup[9][2]=0.46;
   CD4_lookup[9][3]=0.37;
   CD4_lookup[9][4]=0.17;	//1.69 vs 1.27

   // Special test to see if increasing progression times (i.e., reducing death rates) can result in a higher prev/incid ratio (which seems off in sims to date)
   long i, j;
   for (i=1;i<=9;i++) {
     for (j=1;j<=4;j++) {
       CD4_lookup[i][j] = LongevityFudgeFactor*CD4_lookup[i][j];
     }
   }

 }

int SPVL_Category(double SPVL_param)
{
   int SPVL_cat;
   if(SPVL_param<=3.0){SPVL_cat = 1;}
   else{if(SPVL_param<=3.5){SPVL_cat = 2;}
   else{if(SPVL_param<=4.0){SPVL_cat = 3;}
   else{if(SPVL_param<=4.5){SPVL_cat = 4;}
   else{if(SPVL_param<=5.0){SPVL_cat = 5;}
   else{if(SPVL_param<=5.5){SPVL_cat = 6;}
   else{if(SPVL_param<=6.0){SPVL_cat = 7;}
   else{if(SPVL_param<=6.5){SPVL_cat = 8;}
   else{SPVL_cat = 9;}}}}}}}}
   return(SPVL_cat);

}


void cd4InitProbs(){

 //From Pickles CD4 report supplement, table3
   cd4_initProbs[1][1]=0.0;
   cd4_initProbs[1][2]=0.12;
   cd4_initProbs[1][3]=0.88;

   cd4_initProbs[2][1]=0.01;
   cd4_initProbs[2][2]=0.12;
   cd4_initProbs[2][3]=0.87;

   cd4_initProbs[3][1]=0.03;
   cd4_initProbs[3][2]=0.12;
   cd4_initProbs[3][3]=0.85;

   cd4_initProbs[4][1]=0.03;
   cd4_initProbs[4][2]=0.19;
   cd4_initProbs[4][3]=0.78;

   cd4_initProbs[5][1]=0.05;
   cd4_initProbs[5][2]=0.21;
   cd4_initProbs[5][3]=0.73;

   cd4_initProbs[6][1]=0.04;
   cd4_initProbs[6][2]=0.25;
   cd4_initProbs[6][3]=0.71;

   cd4_initProbs[7][1]=0.09;
   cd4_initProbs[7][2]=0.27;
   cd4_initProbs[7][3]=0.64;

   cd4_initProbs[8][1]=1.0;//having all 3 probs =1, forces into box 3
   cd4_initProbs[8][2]=1.0;
   cd4_initProbs[8][3]=1.0;

   cd4_initProbs[9][1]=1.0;
   cd4_initProbs[9][2]=1.0;
   cd4_initProbs[9][3]=1.0;

}


int initialCD4(int spvl_level){
  //can only be 1,2,3 initially
  int tempCat;
  double tempProb;
  tempProb=genrand_real2();

   if(tempProb<cd4_initProbs[spvl_level][1]){tempCat=3;}
   else{if(tempProb<cd4_initProbs[spvl_level][2]){tempCat=2;}
   else{tempCat=1;}}
  // printf("%d",tempCat);
   return(tempCat);
}

void UpdateCD4CountsStochastic()
{
   long i;
   double prog_rate, age_adjusted_prog_time;
   for (i=1;i<=N;i++) {
     if (Status[i] == 1) {
       if (CD4[i] > CD4_nadir[i]) CD4_nadir[i] = CD4[i];
       if( Treated[i] == 0) {
         if (CD4[i] <= 3) {
           age_adjusted_prog_time = CD4_lookup[SpvlCat[i]][CD4[i]];
	   if (Age[i] < 30.0) age_adjusted_prog_time = 1.23*age_adjusted_prog_time;
           if (Age[i] >= 30.0 && Age[i] < 35.0) age_adjusted_prog_time = 1.08*age_adjusted_prog_time;
           if (Age[i] >= 35.0 && Age[i] < 40.0) age_adjusted_prog_time = 0.92*age_adjusted_prog_time;
           if (Age[i] >= 40.0) age_adjusted_prog_time = 0.82*age_adjusted_prog_time;
           prog_rate = 1.0 / age_adjusted_prog_time;
	   prog_rate = prog_rate / 365.0;
	   if (genrand_real2() < prog_rate) {
             CD4[i]=CD4[i]+1;
           }
           TimeInAIDS[i] = 0;
           if (CD4[i] == 4 && TimeInAIDS[i] == 0) AgeAIDS[i] = Age[i];
          } else { // Deterministic dynamics for SPVL category 4 (AIDS)
	    age_adjusted_prog_time = 365*CD4_lookup[SpvlCat[i]][4];
	    if (AgeAIDS[i] < 30.0) age_adjusted_prog_time = 1.83*age_adjusted_prog_time;
            if (AgeAIDS[i] >= 30.0 && AgeAIDS[i] < 35.0) age_adjusted_prog_time = 1.85*age_adjusted_prog_time;
            if (AgeAIDS[i] >= 35.0 && AgeAIDS[i] < 40.0) age_adjusted_prog_time = 0.64*age_adjusted_prog_time;
            if (AgeAIDS[i] >= 40.0) age_adjusted_prog_time = 0.51*age_adjusted_prog_time;
            if( TimeInAIDS[i] > age_adjusted_prog_time) {
              CD4[i] = 5;
            } else {
              TimeInAIDS[i]=TimeInAIDS[i]+1;
            }
         } // CD4 >= 4 
       } //  Treated == 0. 
       if (time > 7300000 && i == 100) {printf("Inside UpdateCD4Stochastic: time = %ld, CD4[%ld] = %d, tx_days = %ld, Treated=%ld\n(type 'c' to continue)\n",time,i, CD4[i], tx_days[i], Treated[i]); scanf("%s",junk); }
       if (Treated[i] == 1) {
         if (CD4[i] >= 2 && CD4[i] == CD4_nadir[i]) { 
           if (genrand_real2() < prob_cd4_rebound) CD4[i] = CD4[i] - 1; //CD4 categories may decrease (i.e., CD4 counts increase) in patients on therapy
         }
         if (CD4[i] >= 2 && CD4[i] == CD4_nadir[i] - 1)  {
          if (genrand_real2() < prob_cd4_rebound_further) CD4[i] = CD4[i] - 1;
         }
       } // Treated
     } // Status == 1

     // Statistics used for DALYS
     if (Status[i] == 1 && Treated[i] == 0 & CD4[i] == 1 && time > Start_TasP_Campaign) cd4_cat1_person_days++;
     if (Status[i] == 1 && Treated[i] == 0 & CD4[i] == 2 && time > Start_TasP_Campaign) cd4_cat2_person_days++;
     if (Status[i] == 1 && Treated[i] == 0 & CD4[i] == 3 && time > Start_TasP_Campaign) cd4_cat3_person_days++;
     if (Status[i] == 1 && Treated[i] == 0 & CD4[i] == 4 && time > Start_TasP_Campaign) cd4_cat4_person_days++;
     if (Status[i] == 1 && Treated[i] == 1 && time > Start_TasP_Campaign) tx_person_days++;
     if (Status[i] == -2 && time > Start_TasP_Campaign) cd4_cat5_person_days++;

  } // i = 1 to N
}

void UpdateCD4Counts()
{
   long i;
   double tempTime,TimeParam;
   for (i=1;i<=N;i++) {
     if (Status[i] == 1) {
      if (CD4[i] > CD4_nadir[i]) CD4_nadir[i] = CD4[i];
      if (CD4[i] < 5) {
        TimeParam = (CD4time[i])/365.0;
        if (CD4_Exp_Flag==0) {
           tempTime = CD4_lookup[SpvlCat[i]][CD4[i]];
        } else {
           tempTime = CD4_TimeToAIDS_exp[i][CD4[i]];
          //printf("tempTime---%f",tempTime);
        }
        if (enhanced_progression_age == 1) {
           /* Notes from Evonet
             Data in Cori-Pickles suggests relative to their averages...
               1.23-fold more time in CD4 categories 1, 2, and 3 for people <30			
               1.08-fold more time in CD4 categories 1, 2, and 3 for people 30-35			
               0.92-fold less time in CD4 categories 1, 2, and 3 for people 35-40			
               0.82-fold less time in CD4 categories 1, 2, and 3 for people >40		
             Corresponding ratios for CD4 < 200 to death
                1.83-fold more time in CD4 category 4 people <30	
                1.85-fold more time in CD4 category 4 people 30-35	
                0.64-fold less time in CD4 category 4 people 35-40		
                0.51-fold less time in CD4 category 4 people >40
           */
           if (CD4[i] <= 3) {
             if (Age[i] < 30.0) tempTime = 1.23*tempTime;
             if (Age[i] >= 30.0 && Age[i] < 35.0) tempTime = 1.08*tempTime;
             if (Age[i] >= 35.0 && Age[i] < 40.0) tempTime = 0.92*tempTime;
             if (Age[i] >= 40.0) tempTime = 0.82*tempTime;
           } else {
             if (Age[i] < 30.0) tempTime = 1.83*tempTime;
             if (Age[i] >= 30.0 && Age[i] < 35.0) tempTime = 1.85*tempTime;
             if (Age[i] >= 35.0 && Age[i] < 40.0) tempTime = 0.64*tempTime;
             if (Age[i] >= 40.0) tempTime = 0.51*tempTime;
           } 
        }
        if( TimeParam > tempTime && Treated[i] == 0) {
          CD4[i]=CD4[i]+1;
          CD4time[i]=0;
          if (CD4[i]==4 && CD4time[i]==0 && CD4_Exp_Flag==0) CD4_TimeToAIDS[i]=(long)(time-Time_Inf[i]);
        } else {
          CD4time[i]=CD4time[i]+1;
        }
      }
      if (Treated[i] == 1) {
        if (CD4[i] >= 2 && CD4[i] == CD4_nadir[i]) {
           if (genrand_real2() < prob_cd4_rebound) CD4[i] = CD4[i] - 1;
        }
        if (CD4[i] >= 2 && CD4[i] == CD4_nadir[i] - 1)  {
           if (genrand_real2() < prob_cd4_rebound_further) CD4[i] = CD4[i] - 1;
        }
      }
    }

    // Statistics used for DALYS
    if (Status[i] == 1 && Treated[i] == 0 & CD4[i] == 1 && time > Start_TasP_Campaign) cd4_cat1_person_days++;
    if (Status[i] == 1 && Treated[i] == 0 & CD4[i] == 2 && time > Start_TasP_Campaign) cd4_cat2_person_days++;
    if (Status[i] == 1 && Treated[i] == 0 & CD4[i] == 3 && time > Start_TasP_Campaign) cd4_cat3_person_days++;
    if (Status[i] == 1 && Treated[i] == 0 & CD4[i] == 4 && time > Start_TasP_Campaign) cd4_cat4_person_days++;
    if (Status[i] == 1 && Treated[i] == 1 && time > Start_TasP_Campaign) tx_person_days++;
    if (Status[i] == -2 && time > Start_TasP_Campaign) cd4_cat5_person_days++;

  } // Status == 1
}


void UpdateViralLoads() 
{
  long i; 
  double V_Old;
  double Vpeak, Vtrans, max_decr_tx = -0.3;
  long t_2nd_phase = 33;
  double d_acute2 = 0.03;
  double Vmin = 0.0001; 
  int acute_up = 0, acute_1st = 0, acute_2nd = 0, chronic = 0, aids_up = 0, aids_max = 0;
  for (i=1;i<=N;i++){
    if (Status[i] == 1) {

      V_Old = V[i];
      Vpeak = pow(10.0, 4.639 + 0.495*LogSetPoint[i]); // Alternative w/out SPVL correlation would be Vpeak = V_peak;
      Vtrans = exp( (2.5*log(SetPoint[i]) + log(Vpeak)) / 3.5 );
      r0 = log(Vpeak/V0)/t_peak;
      d_acute[i] = log(Vpeak/Vtrans) / (t_2nd_phase - t_peak);
      d_acute2 = log(Vtrans/SetPoint[i]) / (t_acute - t_2nd_phase);

      acute_up = 0; acute_1st = 0; acute_2nd = 0; chronic = 0; aids_up = 0; aids_max = 0;
      if (i==10000000) printf("Agent %ld (Time %ld): V_Old = %lf, Vpeak = %lf, Vtrans = %lf, SPVL = %lf, r0 = %lf, treat = %ld, tx_days = %ld, s = %f\n", i,time,V_Old,Vpeak,Vtrans, SetPoint[i], r0, Treated[i], tx_days[i], s[i]); 
      if (i==10000000) printf("   d_acute = %lf, d_acute2 = %lf\n",d_acute[i], d_acute2);
      if (time <= Time_Inf[i] + t_peak ) acute_up = 1;
      if (time > Time_Inf[i] + t_peak && time <= Time_Inf[i] + t_2nd_phase) acute_1st = 1;
      if (time > Time_Inf[i] + t_2nd_phase && time <= Time_Inf[i] + t_acute) acute_2nd = 1;
      if (time > Time_Inf[i] + t_acute && ever_aids[i] == 0) chronic = 1;
      if (i==10000000) printf("time = %ld,  Time_Inf[i] - t_acute = %lf\n",time, Time_Inf[i] - t_acute);
      if (i==10000000) printf("\tTime_Inf = %lf, acute_up = %d, acute_1st = %d, acute_2nd = %d, chronic = %d, ever_aids = %d, t_acute = %lf\t",
                Time_Inf[i], acute_up, acute_1st, acute_2nd, chronic, ever_aids[i],t_acute);
    
      if (acute_up  == 1) V[i] = V0*exp(r0*(time-Time_Inf[i]));
      if (acute_1st == 1) V[i] = Vpeak*exp(-d_acute[i]*(time - t_peak - Time_Inf[i]));
      if (acute_2nd == 1) V[i] = Vtrans*exp(-d_acute2*(time - t_2nd_phase - Time_Inf[i]));
      if (chronic == 1)   V[i] = SetPoint[i]*exp(s[i]*(time - t_acute - tx_days[i] - Time_Inf[i])/365.0);
      if (CD4[i] == 4 || ever_aids[i] == 1) {
         if (ever_aids[i] == 0) {
           ever_aids[i] = 1;
           Time_Start_AIDS[i] = time;
           VL_Start_AIDS[i] = V[i];
         }
         V[i] = VL_Start_AIDS[i]*exp(VL_Incr_AIDS*(time - Time_Start_AIDS[i] - tx_days_aids[i])/365.0);
         if (V[i] > Max_VL_AIDS) V[i] = Max_VL_AIDS;
      }
      // Add a limits to the rate of VL increase (mainly for patients coming off of therapy)
      if (log(V[i]/V_Old) > r0) V[i] = V_Old*exp(r0);  
      
      if (Treated[i] == 1) {
         V[i] = VL_After_Treatment;
         if (chronic == 1) tx_days[i] = tx_days[i] + 1;
         if (ever_aids[i] == 1) tx_days_aids[i] = tx_days_aids[i] + 1;

         // Also add limits for rate of decrease after therapy
         if (log(V[i]/V_Old) < max_decr_tx )V[i] = V_Old*exp(max_decr_tx); 
      }
      if (i==10000000) {printf("Vnew = %lf\nType 'c' to continue",V[i]);  scanf("%s",junk); }
    }
  } //Status == 1
  if (V[i] < 0.0)  V[i] = Vmin;
}
/* current parameters
2 Progression_Model
1 CD4_Death
0 Gamma_Death
*/


// Alternative viral load update function under development
void UpdateViralLoadsAim3() 
{
  double K; // Carrying capacity
  double time_day;
  double step_size = 0.01; // 1/step_size must be an integer 
  double r, r_base;
  double deltaV, deltaI, deltaM, deltaL;
  double M_act = 0.01, L_act = 0.0005;
  double f = 0.98,  fM = 0.019999, fL = 1.0e-6;
  double h = 0.05, d = 0.6, dM = 0.04, dL = 0.001;
  double p = 100.0, pM = 10.0, pL = 1.0, c = 5.0;
  long i;
  r_base = (d + log(V_peak/V0)/t_peak)/(1-fM-fL);
  for (i=1;i<=N;i++){
    if (Status[i] == 1) { 
      // Get carrying capacity 
   
      if (time <= Time_Inf[i] + t_peak) {
        K = V_peak;
      } else {
        if (CD4[i] <= 3) {
          K = SetPoint[i]*exp(s[i]*(time - t_acute - Time_Inf[i])/365.0);
        } 
        if (CD4[i] == 4) {
          if (Time_Start_AIDS[i] > 1000000) {
            Time_Start_AIDS[i] = time;
            VL_Start_AIDS[i] = V[i];
          }
          K = VL_Start_AIDS[i] * exp(VL_Incr_AIDS*(time - Time_Start_AIDS[i]));
          if (K >  Max_VL_AIDS) K = Max_VL_AIDS; 
        }
        if (CD4[i] >= 5) K = 0.0;
      } 
      r = r_base*(K*K*K)/(K*K*K + V[i]*V[i]*V[i]);
      if (i == 1000000000) {
         printf("K[100] =  %lf, SPVL = %lf, V = %lf, r_base = %lf, r = %lf, Time_Start_AIDS = %ld,  ",K,SetPoint[i],V[i],r_base,r,Time_Start_AIDS[100]); 
         printf("    I = %lf, M = %lf, L = %lf, V = %lf\n",I[i],M[i],L[i],V[i]);
      }
      // Use carrying capacities to update viral loads
      for (time_day = 0; time_day <= 1; time_day = time_day + h) {
        r = r_base*(K*K*K)/(K*K*K + V[i]*V[i]*V[i]);
        deltaV = p*I[i] + pM*M[i] + pL*L[i] - c*V[i];
        deltaI = (1-fM - fL)*r*I[i] + M_act*M[i] + L_act*L[i] - d*I[i];
        deltaM = fM*r*I[i] - M_act*M[i] - dM*M[i];
        deltaL = fL*r*I[i] - L_act*L[i] - dL*L[i];
        if (i==10000000000) printf("    deltaI = %lf, deltaM = %lf, deltaL = %lf, deltaV = %lf,",deltaI,deltaM,deltaL,deltaV);
        V[i] = V[i] + h*deltaV; 
        I[i] = I[i] + h*deltaI; 
        M[i] = M[i] + h*deltaM; 
        L[i] = L[i] + h*deltaL; 
        if (V[i] < 0.0)  V[i] = 0.0;
        if (I[i] < 0.0)  I[i] = 0.0;
        if (M[i] < 0.0)  M[i] = 0.0;
        if (L[i] < 0.0)  L[i] = 0.0;
        if (i==1000000000000) printf("    I = %lf, M = %lf, L = %lf, V = %lf\n",I[i],M[i],L[i],V[i]);
  
      }
    }
  }
}

void CollectStatistics()
{
  /* Collect statistics on infecteds and susceptibles */
  long nInfected = 0, nInfectedTreated = 0;
  for (i = 1; i<= N; i++) {
    if (Status[i] == 0) {
      susc_days++; // Running total of susceptible-person days for incidence calculations 
      susc_days_last_5_years++;
      if (Age[i] < 25.0) susc_days_u25++;
      if (Age[i] >= 25.0 && Age[i] < 50.0) susc_days_Mid++;
      if (Age[i] >= 50.0) susc_days_o50++;
      if (Age[i] < 50.0) {
        susc_days_u50++; // Running total of susceptible-person days for incidence calculations 
        susc_days_last_5_years_u50++; // Running total of susceptible-person days for incidence calculations 
      }
    }
    if (Status[i] == 1) nInfected++;
    if (Status[i] == 1 && Treated[i] == 1) nInfectedTreated++; 
  }
  if (Time_All_Treated > 2*tfinal) {
    if (nInfectedTreated >= nInfected) {
      Time_All_Treated = time;
    }
  }
  if (nInfectedTreated < nInfected) {
    Last_Time_Not_All_Treated = time;
  }
  if (time == Start_TasP_Campaign+1) {
    if (nInfected >= 1) ActTxStartTasP = ((double) nInfectedTreated) / ( (double) nInfected);
                  else  ActTxStartTasP = 0.0;
  }
}

void TreatmentDropout() 
{ 
  long i;
  for (i=1; i<=N; i++) {
    if (Status[i] == 1 && Treated[i] == 1) {
      //if (time - Time_Treated[i] >= 30 && genrand_real2() < prob_dropout/365.0) {
      if (time - Time_Treated[i] >= 30 && genrand_real2() < ProbDropout[i]) {
        Treated[i] = 0;
        TargetedTreated[i] = 0;
      }
    } 
  }
}


void InitiateTreatmentGradual() 
{ 
  long i, jj, rand_i, newly_treated;
  double treatment_goal, num_new_treat, num_infected = 0.0, num_treated = 0.0, num_untreated_under_care = 0.0, num_infected_untreated=0.0;
  double rand_num;
  long num_new_treat_random, num_new_treat_targeted;
  double expected_num_new_treat_random, expected_num_new_treat_targeted;
  double frac_part;
  
  if (time <= Start_Spontaneous_Treatment) {
     treatment_goal_start = 0.0;
  }
  
  if (time >= Start_Spontaneous_Treatment) {
    /* Get number of untreated, infected agents (num_infected_untreated)  */
    num_treated = 0;
    for  (i = 1; i<=N; i++) {
      if (Status[i] == 1) num_infected++;
      if (Status[i] >= 0 && Treated[i] == 1) num_treated++; 
      if (Status[i] == 1 && Treated[i] == 0 && UnderCare[i] == 1 && Diagnosed[i] == 1 && (time - TimeDiag[i] > delay_diag_tx)) num_infected_untreated++; 
    }
  }
 
  if (time == Start_Spontaneous_Treatment) {
    treatment_goal_start = round(perc_tx_start * num_infected); // Use number of infected at start of the campaign to define the treatment goal
  }
 
  if (time >= Start_Spontaneous_Treatment) {
    /* Number of newly treated (num_new_treat) is a constant dependent on the "perc_tx_start" and "num_infected" at the start of the spontaneous treatment campaign */
    num_new_treat = treatment_goal_start/365.0;
    
    /* Divide the number of newly treated into random (num_new_treat_random) and target groups (num_new_treat_targeted) */ 
 
    if (time >= Start_TasP_Campaign) {
      expected_num_new_treat_random = PercentNonTargetedYear1 * num_new_treat;
      expected_num_new_treat_targeted =  (1.0 - PercentNonTargetedYear1) * num_new_treat;
    } else {
      expected_num_new_treat_random = num_new_treat;
      expected_num_new_treat_targeted =  0;
    }
     
    num_new_treat_random = trunc(expected_num_new_treat_random);
    frac_part = expected_num_new_treat_random - trunc(expected_num_new_treat_random);
    if (genrand_real2() < frac_part) num_new_treat_random = num_new_treat_random + 1;

    num_new_treat_targeted = trunc(expected_num_new_treat_targeted);
    frac_part = expected_num_new_treat_targeted - trunc(expected_num_new_treat_targeted);
    if (genrand_real2() < frac_part) num_new_treat_targeted = num_new_treat_targeted + 1;

    if (num_new_treat_random + num_new_treat_targeted < num_infected_untreated) {
      if (num_new_treat_random >= 1) TreatAgents(0,num_new_treat_random); 
      if (num_new_treat_targeted >= 1) TreatAgents(1,num_new_treat_targeted);
    }
     
    /* Error checking */
    for  (i = 1; i<=N; i++) {
      if ((Treated[i] == 1) && (Status[i] == 0 || UnderCare[i] == 0 || Diagnosed[i] == 0)) {
        printf("Warning at time %ld: treatment of agent %ld inconsistent with agent's status (%d)\n",time,i,Status[i]); 
        printf("Type 'c' to continue :");scanf("%s",junk);
      }
    }
 } // After start of treatment campaign
}

void InitiateTreatment() 
{ 
  long i, jj, rand_i, newly_treated;
  double treatment_goal, num_new_treat, num_infected = 0.0, num_treated = 0.0, num_untreated_under_care = 0.0, num_infected_untreated=0.0;
  double rand_num;
  long num_new_treat_random, num_new_treat_targeted;
  double expected_num_new_treat_random, expected_num_new_treat_targeted;
  
  if (time >= Start_Spontaneous_Treatment) {
    /* Get number of untreated, infected agents (num_infected_untreated)  */
    num_treated = 0;
    for  (i = 1; i<=N; i++) {
      if (Status[i] == 1) num_infected++;
      if (Status[i] >= 0 && Treated[i] == 1) num_treated++; 
      if (Status[i] == 1 && Treated[i] == 0 && UnderCare[i] == 1 && Diagnosed[i] == 1 && (time - TimeDiag[i] > delay_diag_tx)) num_infected_untreated++; 
    }
    treatment_goal = round(number_tx_spontaneous * (time - Start_Spontaneous_Treatment) / (Start_TasP_Campaign - Start_Spontaneous_Treatment));    
  }
 
  if (time == Start_TasP_Campaign) {
    treatment_goal_start = round(perc_tx_start * num_infected); // Use number of infected at start of the campaign to define the treatment goal
  }

  if (time >= Start_TasP_Campaign) {
    treatment_goal = treatment_goal_start*exp((time - Start_TasP_Campaign)  * yearly_incr_tx/365.0); // Update treatment goal based on yearly growth
  }
 
  if (time >= Start_Spontaneous_Treatment) {
    /* Number of newly treated (num_new_treat) will be the difference between the goal and the number of untreated */
    num_new_treat = treatment_goal - num_treated;
    if (num_infected_untreated < num_new_treat) num_new_treat = (double) num_infected_untreated; // Cannot exceed number untreated
    if (num_new_treat < 0 ) num_new_treat = 0.0; // Cannot be negative
    
    /* Divide the number of newly treated into random (num_new_treat_random) and target groups (num_new_treat_targeted) */ 
 
    if (time >= Start_TasP_Campaign) {
      expected_num_new_treat_random = PercentNonTargetedYear1 * num_new_treat;
      expected_num_new_treat_targeted =  (1.0 - PercentNonTargetedYear1) * num_new_treat;
    } else {
      expected_num_new_treat_random = num_new_treat;
      expected_num_new_treat_targeted =  0;
    }
    
    num_new_treat_random = round (expected_num_new_treat_random + genrand_real2() - 0.5); // Number will fluctuate around the expected number
    num_new_treat_targeted = round (expected_num_new_treat_targeted + genrand_real2() - 0.5);

    if (num_new_treat_random >= 1) TreatAgents(0,num_new_treat_random); 
    if (num_new_treat_targeted >= 1) TreatAgents(1,num_new_treat_targeted);
     
    /* Error checking */
    for  (i = 1; i<=N; i++) {
      if ((Treated[i] == 1) && (Status[i] == 0 || UnderCare[i] == 0 || Diagnosed[i] == 0)) {
        printf("Warning at time %ld: treatment of agent %ld inconsistent with agent's status (%d)\n",time,i,Status[i]); 
        printf("Type 'c' to continue :");scanf("%s",junk);
      }
    }
 } // After start of treatment campaign
}

void TreatAgents(int targeted_treatments, long num_newly_treated)
{
  double rand_num, baseline_treat_prob = 0.0;
  long i,jj;
  long newly_treated = 0, num_untreated_under_care = 0;
  
  newly_treated = 0; escape_counter = 0;
  do {
    /* Baseline treatment prioritization scores */
    for  (i = 1; i<=N; i++) {
      if (Status[i] == 1 && Treated[i] == 0 && UnderCare[i] == 1 && Diagnosed[i] == 1 && (time - TimeDiag[i] > delay_diag_tx)) {
        num_untreated_under_care++;
      }
    }
    baseline_treat_prob = 1.0 / num_untreated_under_care ; 
    
    /* Modified treatment probabilities based on TasP targeting priorities */
    if (time >= Start_TasP_Campaign && targeted_treatments == 1) {
       GetTreatmentProbsMultiplicative(baseline_treat_prob);  // returns treatment_prob[i] values for i = 1 to N
    } else {
       GetRandomTreatmentProbs(baseline_treat_prob); // returns treatment_prob[i] values for i = 1 to N
    }

    total_prob = 0.0;
    for (i=1; i<=N ; i++ ) total_prob = total_prob + treatment_prob[i]; 

    if (total_prob > 0.0) {
      /* Normalize prioritization scores for subsequent selection step h*/
      for (i=1; i<=N; i++) {
        rel_treatment_prob[i] = treatment_prob[i] / total_prob;
      }
        
      /* Generate cumulative distribution of treatment probs (needed for subsequent selection step) */
      cum_treatment_prob[0] = 0.0;
      for (i=1; i<=N; i++) {
        cum_treatment_prob[i] = cum_treatment_prob[i-1] + rel_treatment_prob[i];
      }
      if (cum_treatment_prob[N] > 1.00001) if (plt_warn==1) fprintf(WarningsFile,"Warning at time %ld: cum_treatment_prob[N] = %lf\n",time,cum_treatment_prob[N]);
   
      rand_num = genrand_real2();  // This will be the random patient.  Code below scans through to find the matching index.
      for (jj=1; jj<=N;jj++) {
        if ( (rand_num >= cum_treatment_prob[jj-1]) && (rand_num < cum_treatment_prob[jj])) {
          Treated[jj] = 1;
          Time_Treated[jj] = time;
          newly_treated++;
          pill_count++;
          TotalTreated++;
          TargetedTreated[jj] = 0;
          if (targeted_treatments == 1) {
            if ((cd4_low == 1 && CD4_nadir[jj] >= 2) || (under25 == 1 && Age[jj] <= 25.0) ||
                (under30 == 1 && Age[jj] < 30.0) || (under35 == 1 && Age[jj] <= 35.0) ||
                (vl_high == 1       && LogSetPoint[jj] >= 4.5) || 
                (cd4_super_low == 1  && CD4_nadir[jj] >= 4) || 
                (vl_very_high == 1  && LogSetPoint[jj] >= 5.0) ||
                (vl_super_high == 1 && LogSetPoint[jj] >= 5.5) ||
                (men27women23 == 1 && ((Sex[jj] == 2 && Age[jj] < 27.0) || (Sex[jj] == 1 & Age[jj] < 23.0))) ||
                (men == 1 && Sex[jj] == 2) || (women==1 && Sex[jj] == 1)) {
                TargetedTreated[jj] = 1;
            }
          }
          if (TargetedTreated[jj] == 0) NonTargetedTx++;
          if (time >= Start_TasP_Campaign && time < Start_TasP_Campaign + 365.0) {
             TotalTreatedYear1++;
             if (TargetedTreated[jj] == 0) NonTargetedTxYear1++;
          }
        }
      }
      escape_counter++;
    }       
  } while ((newly_treated < num_newly_treated) && (escape_counter < 4*N) && total_prob > 0.0);
  if (escape_counter > N) {
    printf("escape_counter = %ld at time %ld, targeted_treatments = %d, num_newly_treated_goal = %ld,total_prob = %lf,",escape_counter,time,targeted_treatments,num_newly_treated,total_prob);
    printf(" newly_treated = %ld, cum_treatment_prob[N] = %lf\n",newly_treated,cum_treatment_prob[N]);
  }
}

void GetTreatmentProbs(double baseline_treat_prob)
{
  long i;
  /* Modified prioritization scores under targeted TasP campaign */
  for  (i = 1; i<=N; i++) {
    treatment_prob[i] = 0.0;
    if (Status[i] == 1 && Treated[i] == 0 && UnderCare[i] == 1 && Diagnosed[i] == 1 && (time - TimeDiag[i] > delay_diag_tx)) {
      treatment_prob[i] = baseline_treat_prob; 
      if (under25 == 1 && time >= Start_TasP_Campaign) {
        if (Age[i] < 25.0) {
          treatment_prob[i] =  fold_incr_prob_with_targeting*baseline_treat_prob;       
        } 
      }
      if (under30 == 1 && time >= Start_TasP_Campaign) {
        if (Age[i] < 30.0) {
          treatment_prob[i] =  fold_incr_prob_with_targeting*baseline_treat_prob;       
        } 
      }
      if (under35 == 1 && time >= Start_TasP_Campaign) {
        if (Age[i] < 35.0) {
          treatment_prob[i] =  fold_incr_prob_with_targeting*baseline_treat_prob;       
        } 
      }
      if (cd4_low == 1 && time >= Start_TasP_Campaign) {
        if (CD4_nadir[i] >= 2) {
          treatment_prob[i] =  fold_incr_prob_with_targeting*baseline_treat_prob;       
        } 
      }
      if (cd4_super_low == 1 && time >= Start_TasP_Campaign) {
        if (CD4_nadir[i] >= 4) {
          treatment_prob[i] =  fold_incr_prob_with_targeting*baseline_treat_prob;       
        } 
      }
      if (vl_high == 1 && time >= Start_TasP_Campaign) {
        if (LogSetPoint[i] >= 4.5) {
          treatment_prob[i] =  fold_incr_prob_with_targeting*baseline_treat_prob;       
        } 
      }
      if (vl_very_high == 1 && time >= Start_TasP_Campaign) {
        if (LogSetPoint[i] >= 5.0) {
          treatment_prob[i] =  fold_incr_prob_with_targeting*baseline_treat_prob;       
        } 
      }
      if (vl_super_high == 1 && time >= Start_TasP_Campaign) {
        if (LogSetPoint[i] >= 5.5) {
          treatment_prob[i] =  fold_incr_prob_with_targeting*baseline_treat_prob;       
        } 
      }
      if (men27women23 == 1 && time >= Start_TasP_Campaign) {
        if ((Sex[i] == 2 & Age[i] < 27.0) || (Sex[i] == 1 && Age[i] < 23.0)){
          treatment_prob[i] =  fold_incr_prob_with_targeting*baseline_treat_prob;       
        } 
      }
      if (men == 1 && time >= Start_TasP_Campaign) {
        if (Sex[i] == 2) { 
          treatment_prob[i] =  fold_incr_prob_with_targeting*baseline_treat_prob;       
        } 
      }
      if (women == 1 && time >= Start_TasP_Campaign) {
        if (Sex[i] == 1) { 
          treatment_prob[i] =  fold_incr_prob_with_targeting*treatment_prob[i];       
        } 
      }
    } else {
      treatment_prob[i] = 0.0; // Note above that model does not consider PrEP
    }      
  } 
}

// Version in which persons have more than one favorable factor are especially likely to get treated
void GetTreatmentProbsMultiplicative(double baseline_treat_prob)
{
  long i;
  /* Modified prioritization scores under targeted TasP campaign */
  for  (i = 1; i<=N; i++) {
    treatment_prob[i] = 0.0;
    if (Status[i] == 1 && Treated[i] == 0 && UnderCare[i] == 1 && Diagnosed[i] == 1 && (time - TimeDiag[i] > delay_diag_tx)) {
      treatment_prob[i] = baseline_treat_prob; 
      if (under25 == 1 && time >= Start_TasP_Campaign) {
        if (Age[i] < 25.0) {
          treatment_prob[i] =  fold_incr_prob_with_targeting*treatment_prob[i];       
        } 
      }
      if (under30 == 1 && time >= Start_TasP_Campaign) {
        if (Age[i] < 30.0) {
          treatment_prob[i] =  fold_incr_prob_with_targeting*treatment_prob[i];       
        } 
      }
      if (under35 == 1 && time >= Start_TasP_Campaign) {
        if (Age[i] < 35.0) {
          treatment_prob[i] =  fold_incr_prob_with_targeting*treatment_prob[i];       
        } 
      }
      if (cd4_low == 1 && time >= Start_TasP_Campaign) {
        if (CD4_nadir[i] >= 2) {
          treatment_prob[i] =  fold_incr_prob_with_targeting*treatment_prob[i];       
        } 
      }
      if (cd4_super_low == 1 && time >= Start_TasP_Campaign) {
        if (CD4_nadir[i] >= 4) {
          treatment_prob[i] =  fold_incr_prob_with_targeting*treatment_prob[i];       
        } 
      }
      if (vl_high == 1 && time >= Start_TasP_Campaign) {
        if (LogSetPoint[i] >= 4.5) {
          treatment_prob[i] =  fold_incr_prob_with_targeting*treatment_prob[i];       
        } 
      }
      if (vl_very_high == 1 && time >= Start_TasP_Campaign) {
        if (LogSetPoint[i] >= 5.0) {
          treatment_prob[i] =  fold_incr_prob_with_targeting*treatment_prob[i];       
        } 
      }
      if (vl_super_high == 1 && time >= Start_TasP_Campaign) {
        if (LogSetPoint[i] >= 5.5) {
          treatment_prob[i] =  fold_incr_prob_with_targeting*treatment_prob[i];       
        } 
      }
      if (men27women23 == 1 && time >= Start_TasP_Campaign) {
        if ((Sex[i] == 2 & Age[i] < 27.0) || (Sex[i] == 1 && Age[i] < 23.0)){
          treatment_prob[i] =  fold_incr_prob_with_targeting*treatment_prob[i];       
        } 
      }
      if (men == 1 && time >= Start_TasP_Campaign) {
        if (Sex[i] == 2) { 
          treatment_prob[i] =  fold_incr_prob_with_targeting*treatment_prob[i];       
        } 
      }
      if (women == 1 && time >= Start_TasP_Campaign) {
        if (Sex[i] == 1) { 
          treatment_prob[i] =  fold_incr_prob_with_targeting*treatment_prob[i];       
        } 
      }
    } else {
      treatment_prob[i] = 0.0; // Note above that model does not consider PrEP
    }      
  } 
}

void GetRandomTreatmentProbs(double baseline_treat_prob)
{
  long i;
  /* Modified prioritization scores under targeted TasP campaign */
  for  (i = 1; i<=N; i++) {
    treatment_prob[i] = 0.0;
    if (Status[i] == 1 && Treated[i] == 0 && UnderCare[i] == 1 && Diagnosed[i] == 1 && (time - TimeDiag[i] > delay_diag_tx)) {
      treatment_prob[i] = baseline_treat_prob; 
    } 
  }
}

void GetCD4_BasedTreatmentProbs(double baseline_treat_prob)
{
  /* Modified prioritization scores under targeted TasP campaign */
  for  (i = 1; i<=N; i++) {
    treatment_prob[i] = 0.0;
  if (Status[i] == 1 && Treated[i] == 0 && UnderCare[i] == 1 && Diagnosed[i] == 1 && (time - TimeDiag[i] > delay_diag_tx)) {
      treatment_prob[i] = baseline_treat_prob; 
      if (time < Start_TasP_Campaign) {
        if (CD4[i] == 4) {
          treatment_prob[i] =  2.0*treatment_prob[i];       
        } 
        if (Sex[i] == 1) {
          treatment_prob[i] =  1.4*treatment_prob[i];       
        } 
        if (Age[i] >= 30) {
          treatment_prob[i] =  2.5*treatment_prob[i];       
        } 
      }
    } else {
      treatment_prob[i] = 0.0; // Note above that model does not consider PrEP
    }      
  }
}

void SimulateTransmission()
{
  long Infector, Recipient; // Potential infector and Recipient
  long i,j,temp_i,ii, jj, iii;
  long random_person;
  double average_age, condom_use_prob;
  double Adjusted_Prob_Sex;
  //printf("Starting SimulateTransmission\n");
  for (i=1;i<=N;i++ ) {
     HadSexAlready[i] = 0; // Refers to whether person had sex today
  }
  for (i=1;i<=N;i++) {
    RandomIndex[i] = i; 
  }
  // Now do a whole bunch of swaps to randomize the order
  for (ii=1; ii<=N; ii++) {
    jj = genrand_real2()*N + 1;
    temp_i = RandomIndex[ii];
    RandomIndex[ii] = RandomIndex[jj];
    RandomIndex[jj] = temp_i;
  }
  //printf("Finished swaps\n");

  for (random_person=1; random_person<=N; random_person++)
  {
     //i = RandomIndex[random_person]; // Select person i at random
     i = random_person;
     if (i > N || i < 1) printf("Trouble index out of range\n");
     if (NumLinks[i] >= 1) {
       random_partner = genrand_real2()*NumLinks[i] + 1;  // Person i selects one of his or her sex partners each day.  (Assumes coital dilution)
       j = Links[i][random_partner]; // That person is person j
       if ( ((Status[i] + Status[j] == 1) && (max(Status[i],Status[j]) == 1)) && (HadSexAlready[i] + HadSexAlready[j] == 0))
       {
         if (Status[i] == 1)
         {
           Infector = i;
           Recipient = j;
         } else {
           Infector = j;
           Recipient = i;
         }
         Adjusted_Prob_Sex = prob_sex;
         if (CD4[i] >=4) Adjusted_Prob_Sex = prob_sex_in_AIDS;
         if (age_dependent_sex == 1) Adjusted_Prob_Sex = Adjusted_Prob_Sex *  (MaxAgeSex - 0.5*(Age[Infector] + Age[Recipient])) / (75.0 - 15.0);
         if ((time - TimeLinked[i][random_partner]) <= HoneymoonDays) {
            Adjusted_Prob_Sex = IncreasedProbSexHoneymoonDays * Adjusted_Prob_Sex; // Higher probability of having sex in recently formed relationships 
            if (time >= 100000000 && time < 10100) printf("Time %ld: Adjusted_Prob_Sex increased to %lf because TimeLinked[%ld][%ld] = %ld (Dur1 = %lf, Dur2 = %lf)\n",
                                               time,Adjusted_Prob_Sex,i,random_partner,TimeLinked[i][random_partner],Duration[i],Duration[j]); 
         } else {
            if (time >= 100000000 && time < 10100) printf("Time %ld: Adjusted_Prob_Sex lef at %lf because TimeLinked[%ld][%ld] = %ld (Dur1 = %lf, Dur2 = %lf)\n",
                                               time,Adjusted_Prob_Sex,i,random_partner,TimeLinked[i][random_partner],Duration[i],Duration[j]); 
         } 
         if (Adjusted_Prob_Sex < 0.0) Adjusted_Prob_Sex = 0.0;
         
         if (time == 20000) {
           //printf("Time %ld: Adjusted_Prob_Sex = %lf (Dur[%ld] = %lf, Dur[%ld] = %lf, CD4[%ld] = %d, Age[%ld] = %lf, Age[%ld]= %lf, Status[%ld] = %d, Status[%ld] = %d\n",
           //      time, Adjusted_Prob_Sex,i, Duration[i], j, Duration[j],i, CD4[i], i, Age[i], j, Age[j], i, Status[i], j, Status[j]);
           //for (iii=N-400; iii<=N-360;iii++) {
           //  printf("Numlinks[%ld] = %ld\n",iii,NumLinks[iii]);
           //}
           //printf("Hit return to continue :"); scanf("%s",junk);     
         }
         if (time == 300000 && i == 30012) {
            printf("Time %ld, agent %ld,  Adj Prob Sex = %lf\n",time,i,Adjusted_Prob_Sex);
            printf("\tAge[%ld] = %lf, Age[%ld] = %lf, CD4[%ld] = %d,",i,Age[i],j,Age[j],i,CD4[i]);
            printf("TimeLinked[%ld][%ld] = %ld\n",i,random_partner,TimeLinked[i][random_partner]);
         }
         if (genrand_real2() < Adjusted_Prob_Sex) {
            HadSexAlready[Infector] = 1;
            HadSexAlready[Recipient] = 1;
             if (Inf_Rate_Function == 1) {
                if (V[Infector] > VL_After_Treatment)
                  Inf_Prob = InfRateBaseline*pow(log10(V[Infector]),InfRateExponent);
                else
                  Inf_Prob = 0.0; // power function gives negative values if V < 1.0 (Shouldn't be much of a problem in real-life)
             } else {
               if (V[Infector] > VL_After_Treatment) {
                 beta = HillCoeffInfRate;
                 Inf_Prob = MaxInfRate * pow(V[Infector],beta) / (pow(V[Infector],beta) + pow(VHalfMaxInfRate,beta));
               } else {
                 Inf_Prob = 0.0;
               }
            }
            // New code to have a reduced probability of transmission in vaccinated people (assuming donor's virus is sensitive)
            if ((Vaccinated[Recipient] == 1) && (ResistantToVaccine[Infector] == 0)) {
               Inf_Prob = RR_Vaccinated * Inf_Prob;
            }
            if (Sex[Recipient] == 1) {  
               Inf_Prob = RR_female_recipient * Inf_Prob; // Women thought to be more susceptible than men in general
            }
            if (Age[Recipient] < 45.0 ) {
              //if (time == 10000) printf("Time %ld: Inf_Prob before %lf",time,Inf_Prob); 
              Inf_Prob = Inf_Prob * pow(RR_youth, (45.0 - Age[Recipient])/10.0);  // Elevated risk for each decade below 45. Note that Evonet inverts RR: dat$param$trans_RR_youth^((age_vec_sus-dat$param$trans_base_age)/10)
              //if (time == 10000) printf(", after %lf\n",Inf_Prob); 
            }
            if (Sex[Recipient] == 2) {
              if (Circumcised[Recipient] == 1) Inf_Prob = RR_circum*Inf_Prob;
            }

            average_age = (Age[Infector] + Age[Recipient])/2.0;
            condom_use_prob = condom_prob16*pow(0.5,(average_age-15.0)/(age_condom_use_halves-15.0)); // Age 35: 0.5^[(35-15)/(35-15)] = 0.5, Age 50: 0.5^[(50-15)/(35-15)] = 0.3
            if (genrand_real2() < condom_use_prob) Inf_Prob = RR_condom*Inf_Prob;            
 
            randnum = genrand_real2();
            if (randnum < Inf_Prob) {
              //printf("%ld infects %ld at time %ld\n", Infector, Recipient, time);
              Time_Inf[Recipient] = time;
              V[Recipient] = V0;

              //ViralContribToLogSP0[Recipient] = ViralContribToLogSP0[Infector] + norm_rand(0,((time-Time_Inf[Infector])/365)*MutationVariance); 
              // Patient inherits previous patient's virus + mutational deviation scaled by time since infection
              ViralContribToLogSP0[Recipient] = ViralContribToLogSP0[Infector] + norm_rand(0, MutationVariance); 
              // Patient inherits previous patient's virus + mutational deviation
  
              EnvirContribToLogSP0[Recipient] = norm_rand(0, sqrt((1-H)*VarianceLogSP0)); 
              // Environmental component is independent
          
              LogSetPoint[Recipient] = ViralContribToLogSP0[Recipient] + EnvirContribToLogSP0[Recipient];
          
              SetPoint[Recipient] = pow(10.0,LogSetPoint[Recipient]);
  
              if (genrand_real2() < percent_under_care) {
                 UnderCare[Recipient] = 1;
              } else {
                 UnderCare[Recipient] = 0;
              }
  
              d_acute[Recipient] = 0.4; // Placeholder until exact value is calculated in UpdateViralLoads()
               
              ResistantToVaccine[Recipient] = ResistantToVaccine[Infector]; // Recipient inherits infectors vaccinne resistance status
               
              Generation[Recipient] = Generation[Infector] + 1;
              Status[Recipient] = 1;
              SpvlCat[Recipient]=SPVL_Category(LogSetPoint[Recipient] );
              CD4[Recipient]=initialCD4(SpvlCat[Recipient]);
              CD4_nadir[Recipient] = CD4[Recipient];
              CD4_initial_value[Recipient]=CD4[Recipient];
              CD4time[Recipient]=0;
              TimeInAIDS[Recipient]=0;
              //adding exponential distribution for time in each cd4 category
              //can tighten this up once it gets working....
              CD4_TimeToAIDS_exp[Recipient][1]=log( genrand_real2()) / (-1.0/CD4_lookup[SpvlCat[Recipient]][1]) ;
              CD4_TimeToAIDS_exp[Recipient][2]=log( genrand_real2()) / (-1.0/CD4_lookup[SpvlCat[Recipient]][2]) ;
              CD4_TimeToAIDS_exp[Recipient][3]=log( genrand_real2()) / (-1.0/CD4_lookup[SpvlCat[Recipient]][3]) ;
              CD4_TimeToAIDS_exp[Recipient][4]=log( genrand_real2()) / (-1.0/CD4_lookup[SpvlCat[Recipient]][4]) ;
 
              new_infections++; // Running total of number of people getting infected for incidence calculations 
              new_infections_last_5_years++; // Running total of number of people getting infected for incidence calculations 
              if (Age[Recipient] < 25.0)  new_infections_u25++;
              if (Age[Recipient] >= 25.0 && Age[Recipient] < 50.0) new_infections_Mid++;
              if (Age[Recipient] > 50.0)  new_infections_o50++;
              if (Age[Recipient] <= 50.0) {
                new_infections_u50++;  
                new_infections_last_5_years_u50++;  
              } 
              if(CD4_Exp_Flag==1) {
                double ll=0;
                int kk;
                for(kk=CD4[Recipient];kk<=3;kk++)
                {
                  ll+=CD4_TimeToAIDS_exp[Recipient][kk];
                }
                CD4_TimeToAIDS[Recipient]=ll;
              }
  
              Donors_ViralContribToLogSP0[Recipient] = ViralContribToLogSP0[Infector];
              Donors_EnvirContribToLogSP0[Recipient] = EnvirContribToLogSP0[Infector];
              Donors_d_acute[Recipient] = d_acute[Infector];
              Donors_Total_Time_Inf_At_Trans[Recipient] = time - Time_Inf[Infector];
              Donors_V[Recipient] = V[Infector];
              Donors_Generation[Recipient] = Generation[Infector];
              Donors_SetPoint[Recipient] = SetPoint[Infector];
              Donors_LogSetPoint[Recipient] = LogSetPoint[Infector];
              Donors_Index[Recipient] = Infector;
              NumRecipients[Infector] = NumRecipients[Infector] + 1;
            } // Infection event occurred
	 } // Couple had sex
       } // Status[i] + Status[j] == 1 (i.e., opporuntunity for new infection) & they hadn't already had sex
      } // Person i has at least one partner
   } // for i = 1 to N
}

void NaturalDeath()
{

   long i;
   for (i=1;i<=N;i++) {
     if (Status[i] >= 0) {
       randnum = genrand_real2();
       if (randnum < asmr[ ((long) Age[i]) ]) {
           Status[i] = -1;
           if (plt_R==1) fprintf(PatientRegistryOutputFile,"%ld\t%ld\t%ld\t%lf\t%9.6e\t%ld\t%ld\t%ld\t%ld\t%6.5e\t%6.5e\t%6.5e\t%6.5e\t%6.5e\t%6.5e\t%6.5e\t%6.5e\t%6.5e\t%6.5e\t%6.5e\t%6.5e\t%d\t%7.4f\t%7.4f\t%7.4f\t%7.4f\t%7.4f\t%d\t",
                         repl,time,i,Age[i],Duration[i],NumRecipients[i],Donors_Index[i],Generation[i],Donors_Generation[i],Time_Inf[i],V[i],Donors_V[i],LogSetPoint[i],Donors_LogSetPoint[i],ViralContribToLogSP0[i],Donors_ViralContribToLogSP0[i],EnvirContribToLogSP0[i],d_acute[i],Donors_d_acute[i],time-Time_Inf[i],Donors_Total_Time_Inf_At_Trans[i],CD4[i],CD4_TimeToAIDS[i],
                        CD4_TimeToAIDS_exp[i][1],CD4_TimeToAIDS_exp[i][2],CD4_TimeToAIDS_exp[i][3],CD4_TimeToAIDS_exp[i][4],CD4_initial_value[i]);
           if (plt_R==1) fprintf(PatientRegistryOutputFile,"DiedNatural\n");
           V[i] = 0.0;
       } // Person died of natural causes
     } // Status == 1 (i.e., person is alive)
   } // i = 1 to N
}


void DeathOfAIDSPatients()
{
  long i;
   NewAIDSDeaths = 0;

/* Current parameters
2 Progression_Model
1 CD4_Death
0 Gamma_Death
*/
    
   for (i=1;i<=N;i++) {
     if (CD4[i] == 5 && Status[i]==1 && Treated[i]!=1) {
     if (plt_R==1) fprintf(PatientRegistryOutputFile,"%ld\t%ld\t%ld\t%lf\t%9.6e\t%ld\t%ld\t%ld\t%ld\t%6.5e\t%6.5e\t%6.5e\t%6.5e\t%9.6e\t%9.5e\t%6.5e\t%6.5e\t%6.5e\t%6.5e\t%6.5e\t%6.5e\t%d\t%7.4f\t%7.4f\t%7.4f\t%7.4f\t%7.4f\t%d\t",
                         repl,time,i,Age[i],Duration[i],NumRecipients[i],Donors_Index[i],Generation[i],Donors_Generation[i],Time_Inf[i],V[i],Donors_V[i],LogSetPoint[i],Donors_LogSetPoint[i],ViralContribToLogSP0[i],Donors_ViralContribToLogSP0[i],EnvirContribToLogSP0[i],d_acute[i],Donors_d_acute[i],time-Time_Inf[i],Donors_Total_Time_Inf_At_Trans[i],CD4[i],CD4_TimeToAIDS[i],
                         CD4_TimeToAIDS_exp[i][1],CD4_TimeToAIDS_exp[i][2],CD4_TimeToAIDS_exp[i][3],CD4_TimeToAIDS_exp[i][4],CD4_initial_value[i]);
        V[i] = 0.0;
        Status[i] = -2;
        Time_AIDS_Death[i] = time;
        NewAIDSDeaths++;
        if (plt_R==1) fprintf(PatientRegistryOutputFile,"DiedAIDS\n");
     }
   }
   
   // Additional probabilities of death based on CD4 T-cell count whether treated or not (e.g., by cancer)
   double cd4_cat1_death_prob      = 0.0000112;    // prob. of death for cd4 cat1
   double cd4_cat2_death_prob      = 0.0000148;    // prob  of death for cd4 cat2
   double cd4_cat3_death_prob      = 0.0000333;    // prob of death for cd4 cat3
   double cd4_cat4_treated_death_prob = 0.0000760; // #prob death for cd4 cat4(aids) on tx
   for (i=1;i<=N;i++) {
     if (Status[i] == 1 && CD4[i] == 1) {
       if (genrand_real2() < cd4_cat1_death_prob) {
         V[i] = 0.0;
         Status[i] = -2;
         Time_AIDS_Death[i] = time;
         NewAIDSDeaths++;
       }
     }
     if (Status[i] == 1 && CD4[i] == 2) {
       if (genrand_real2() < cd4_cat2_death_prob) {
         V[i] = 0.0;
         Status[i] = -2;
         Time_AIDS_Death[i] = time;
         NewAIDSDeaths++;
       }
     }
     if (Status[i] == 1 && CD4[i] == 3) {
       if (genrand_real2() < cd4_cat3_death_prob) {
         V[i] = 0.0;
         Status[i] = -2;
         Time_AIDS_Death[i] = time;
         NewAIDSDeaths++;
       }
     }
     if (Status[i] == 1 && CD4[i] == 4 && Treated[i] == 1) {
       if (genrand_real2() < cd4_cat4_treated_death_prob) {
         V[i] = 0.0;
         Status[i] = -2;
         Time_AIDS_Death[i] = time;
         NewAIDSDeaths++;
       }
     }
  } // Additional sources of death 
}


void RemoveLinksToDeadPeople()
{
   long i, j, k, ii;
   /* Remove links to dead people */
   for (i = 1; i <= N; i++) {
      if ((Status[i] < 0) && (NumLinks[i] > 0)) {
        LinksToRemove = NumLinks[i];
        /*for (ii=1; ii<=NumLinks[ii]; ii++) {
          if (i < Links[i][ii]) {
             if (Status[i] == -1) fprintf(RelationshipRegistryFile,"%ld\t%ld\t%ld\t%ldDiedNat\n",i,Links[i][ii],time,i);
                             else fprintf(RelationshipRegistryFile,"%ld\t%ld\t%ld\t%ldDiedAIDS\n",i,Links[i][ii],time,i);
          } else {
             if (Status[i] == -1) fprintf(RelationshipRegistryFile,"%ld\t%ld\t%ld\t%ldDiedNat\n",Links[i][ii],i,time,i);
                             else fprintf(RelationshipRegistryFile,"%ld\t%ld\t%ld\t%ldDiedAIDS\n",Links[i][ii],i,time,i);
          }
        } */
        for (j = 1; j <= LinksToRemove; j++) {
           Widowed_partner = Links[i][j];
            /* Need to find which of the widowed partner's links points back to i */
           Link_back = 0;
           for (k = 1; k <= NumLinks[Widowed_partner]; k++) {
               if (Links[Widowed_partner][k] == i) {
                Link_back = k;
                //printf(" --- It Does!!! (so setting Link_back = %ld)",Link_back);
              }
              //printf("\n");
           }
           if (Link_back == 0) {
             printf("Trouble with unwinding links from the widow (%ld) of dead person %ld (NumLinks[%ld] = %ld: \n",Widowed_partner,i,i,NumLinks[i]);
             printf("   Couldn't find any link backs from the Widow's list!\n"); fflush(stdout); exit(1);
           }
           Links[Widowed_partner][Link_back] = Links[Widowed_partner][NumLinks[Widowed_partner]]; // Replace link to deceased with the last person on widowed persons contact list
           TimeLinked[Widowed_partner][Link_back] = TimeLinked[Widowed_partner][NumLinks[Widowed_partner]]; // Replace link to deceased with the last person on widowed persons contact list
           
           Links[Widowed_partner][NumLinks[Widowed_partner]] = 0; // Remove last person (as to not to count the last person twice)
           TimeLinked[Widowed_partner][NumLinks[Widowed_partner]] = 0; // Remove last person (as to not to count the last person twice)

           NumLinks[Widowed_partner] = NumLinks[Widowed_partner] - 1;

        } // j = 1 to NumLinks[i]
        NumLinks[i] = 0;
      }
   }
}


void SimulatePartnerShipDissolution()
{
  //printf("Inside SimulatePartnershipDissolution:\n"); fflush(stdout);
  long person1, person2, relationship, j;
  long counter;
  double Dur1, Dur2;
  breakups = 0;
  //printf("Inside SimulatePartnershipDissolution:\n"); fflush(stdout);
  for (person1 = 1; person1 <= N; person1++) {
    //printf("Testing if person %ld breaks up with anyone\n",person1); fflush(stdout);
    for (index_to_person2 = 1; index_to_person2 <= NumLinks[person1]; index_to_person2++) {
      person2 = Links[person1][index_to_person2];
      //printf("  Testing if person %ld breaks up with %ld\n",person1,person2); fflush(stdout);
  
      if (person2 >= person1) { // If person2 < person1, this relationship has already been tested.
         // Test if persons 1 and 2 break up
         Dur1 = Duration[person1];
         Dur2 = Duration[person2];
         if (HighRisk[person1] == 1) Dur1 = Dur1/2.0;
         if (HighRisk[person2] == 1) Dur2 = Dur2/2.0;
         PartnershipDuration = sqrt(Dur1*Dur2);  // Geometric Mean of each partner's duration
         ProbPartnersBreakUp = 1.0/PartnershipDuration; // Approximate method for determining probability of breakup
         if ((genrand_real2() < ProbPartnersBreakUp) || ((Status[person1] < 0) || (Status[person2] < 0)) ) { // Partnership broke up
           //if (time <= 3) printf("Targeting the relationship between %ld and %ld ", person1,person2);
           //if (time <= 3) printf("(index %ld from %ld's %ld sexual contacts)", index_to_person2, person1, NumLinks[person1]);
          
           TimeLastBreakUp[person1] = time;
           TimeLastBreakUp[person2] = time;
 
           /* Find the index in person2's list that points back to person 1 */
           found_link = 0;
           for (j=1;j<=NumLinks[person2];j++) {
              if (Links[person2][j] == person1) {
                found_link = 1;
                index_to_person1 = j;
              }
           }
         //if (time <= 3 && printOutput==1) printf(" (index %ld from %ld's %ld sexual contacts) for breakup\n", index_to_person1, person2, NumLinks[person2]);
            
           if (found_link == 0 && printOutput==1) {
              printf("Big trouble in breakups : couldn't find the corresponding link for person1 (%ld) in person2's (%ld) contact list\n",person1,person2);
              printf("   NumLinks[%ld] = %ld, Links[%ld][1] = %ld, NumLinks[%ld] = %ld, Links[%ld][1] = %ld, Links[%ld][2] = %ld, Status[%ld] = %d, Status[%ld] = %d\n",
                     person1,NumLinks[person1],person1,Links[person1][1],person2,NumLinks[person2],person2,Links[person2][1],person2,Links[person2][2],person1,Status[person1],person2,Status[person2]);
              fflush(stdout); exit(1);
           }
           breakups++;
  
          // Add this couple to the list of who broke up with whom.  But don't remove this link now.  (Doing so would change the upper limit of the for loop above)
          if (breakups < 3000000) {
             BreakUpList1[breakups] = person1;
             BreakUpList2[breakups] = index_to_person2;
             BreakUpList3[breakups] = person2;
             BreakUpList4[breakups] = index_to_person1;
           } else {
             printf("Error: BreakUpList1 and BreakupList2 arrays are not large enough to track of all of the breakups.  Type ctrl-C' to exit\n");  StopEarly = 1;
           }
        } // Persons 1 and 2 broke up
      } // Relationship not already tested
    } // partners of person1
  } // person1 from 1 to N

  // Flag links for deletion
  if (time <= -3 && printOutput==1) printf("Removing %ld partnerships from the list of Linked persons\n",breakups);   fflush(stdout);

  for (relationship = 1; relationship <= breakups; relationship++) {
	 person1 = BreakUpList1[relationship];
	 index_to_person2 = BreakUpList2[relationship];
	 person2 = BreakUpList3[relationship];
	 index_to_person1 = BreakUpList4[relationship];
         /*if (person1 < person2) fprintf(RelationshipRegistryFile,"%ld\t%ld\t%ld\tBreakup\n",person1,person2,time);
                           else fprintf(RelationshipRegistryFile,"%ld\t%ld\t%ld\tBreakup\n",person2,person1,time); */

	 /* Remove the link from person1 */
      Links[person1][index_to_person2] = -1; // Flag this link for deletion
      TimeLinked[person1][index_to_person2] = 0; // Zero out time that they got together

	  /* Remove the corresponding link index from person 2 */
      Links[person2][index_to_person1] = -1; // Flag this link for deletion
      TimeLinked[person2][index_to_person1] = 0; // Zero out time that they got together
  }

  // Skip over any links flagged for deletion.
  // For example, if person1 was originally linked to 1,2,3,4, but broke up with 1 and 3, the new list would be 2, 4, 3, 4 (with last two being orphans)
  //
   for (person1 = 1; person1 <= N; person1++) {
     RevisedNumLinks[person1] = 0;
	 counter = 0;
  	 for (index_to_person2 = 1; index_to_person2 <= NumLinks[person1]; index_to_person2++) {
		if (Links[person1][index_to_person2] > 0) {
		   counter++;
		   Links[person1][counter] = Links[person1][index_to_person2];
		   TimeLinked[person1][counter] = TimeLinked[person1][index_to_person2];
		}
		RevisedNumLinks[person1] = counter;
	 }
	 if (Status[person1] < 0) RevisedNumLinks[person1] = 0;
  }

  // Zero out any expired links (e.g., getting rid of 3 and 4 above), then update the number of links for each person
  for (person1 = 1; person1 <= N; person1++) {
    // Zero out any inactive links
    for (index_to_person2=RevisedNumLinks[person1]+1; index_to_person2<=NumLinks[person1]; index_to_person2++) {
      Links[person1][index_to_person2] = 0;
      TimeLinked[person1][index_to_person2] = 0;
    }
    NumLinks[person1] = RevisedNumLinks[person1];
  }
}


void SimulateNewPartnershipFormation()
{
  long i;
    /* Simulate process by which singles, widows, widowers, and jilted lovers look for new contacts */
   Active = 0.0;
   NLinks = 0.0;
   for (i=1;i<=N;i++) {
      if (Status[i] >= 0) {
        Active = Active + 1.0;
        NLinks = NLinks + (double) NumLinks[i];
      }
   }
   if (time_varying_mean_degree == 1) {
     LinksToAdd = GetPoisson(Active,d_rate);
   } else {
     ExpectedLinks = (long) ( Active* AverageDegree/2.0 + 0.5);
     dpLinksToAdd = ExpectedLinks - NLinks/2.0;
     LinksToAdd = (long) (dpLinksToAdd + 0.5);
  }
  if (LinksToAdd >= 1) AddNewLinks(LinksToAdd);
}

void CheckEdgeListForImpossibleValues()
{
   /* Do some elementary error checking before moving onto the next day*/
   for (i=1;i<=N;i++) {
     if ((NumLinks[i] < 0) || (NumLinks[i] > MaxLinks) ) {
       printf("Trouble: NumLinks[%ld] = %ld\n",i,NumLinks[i]);
       fflush(stdout); exit(1);
     }
     /*if ((NumLinks[i] >= 2) && (Concurrent[i] == 0)) {
       printf("Trouble: NumLinks[%ld] = %ld and Concurrent[%ld] = %ld\n",i,NumLinks[i],i,Concurrent[i]);
       fflush(stdout); exit(1);
     }*/
     for (j=1; j<=NumLinks[i];j++) {
        if ((Links[i][j] <= 0) || (Links[i][j] > N)) {
          printf("Trouble at end of day %ld: Links[%ld][%ld]= %ld (NumLinks[%ld] = %ld, N = %ld)\n",time,i,j,Links[i][j],i,NumLinks[i],N);
          fflush(stdout); exit(1);
        }
     }
     if (Status[i] < 0 &&  NumLinks[i] > 0) {
       printf("Trouble at end of day %ld: Status[%ld] = %d, but NumLinks[%ld] = %ld\n",time,i,Status[i],i,NumLinks[i]);
       fflush(stdout); exit(1);
     }
   }
}

void RecordStatusOfLivingHIVPatientsForRegistryFile()
{

 for (i=1;i<=N;i++) {
    if (Status[i] == 1) {
   if (plt_R==1) fprintf(PatientRegistryOutputFile,"%ld\t%ld\t%ld\t%lf\t%9.6e\t%ld\t%ld\t%ld\t%ld\t%6.5e\t%6.5e\t%6.5e\t%6.5e\t%9.6e\t%9.6e\t%6.5e\t%6.5e\t%6.5e\t%6.5e\t%6.5e\t%6.5e\t%d\t%7.4f\t%7.4f\t%7.4f\t%7.4f\t%7.4f\t%d\tAlive(HIV+)\n",
                         repl,time,i,Age[i],Duration[i],NumRecipients[i],Donors_Index[i],Generation[i],Donors_Generation[i],Time_Inf[i],V[i],Donors_V[i],LogSetPoint[i],Donors_LogSetPoint[i],ViralContribToLogSP0[i],Donors_ViralContribToLogSP0[i],EnvirContribToLogSP0[i],d_acute[i],Donors_d_acute[i],time-Time_Inf[i],Donors_Total_Time_Inf_At_Trans[i],CD4[i],CD4_TimeToAIDS[i],
                         CD4_TimeToAIDS_exp[i][1],CD4_TimeToAIDS_exp[i][2],CD4_TimeToAIDS_exp[i][3],CD4_TimeToAIDS_exp[i][4],CD4_initial_value[i]);

  }
 }

}


long max(long a, long b)
{
  if (a > b) return (a);
       else  return (b);
}


long GetPoisson(long num, double pois_prob)
{
   long numfound, i;
   numfound = 0;
   for (i = 1; i <= num; i++) {
     if (genrand_real2() < pois_prob) {
       numfound++;
     }
   }
   return(numfound);
}

void AddNewLinks(long LinksToAdd)
{ /* Note: While this is a clunky algorithm, it seems to work */
  long person1, person2; // Indices used to define pairings when building the sexual network
  long j, AcceptsNewPartner, SeeksNewPartner;
  double female_pref_older_men;
  double age_diff_50perc_less_likely = 0.15;
  double average_age, age_diff_with_offset, percent_diff_age, age_diff_index, rand_num;
  long male_age;
  int Unsuitable, failed_to_find_a_partner;
  long big_escape_counter = 0;
  long new_links_added = 0, attempts = 0;
  double Am, Af, prob_more_partners, age_specific_prob_conc;

  do {

    if (time >= 10000000) {printf("Seeking to add %ld new links at time %ld\n bescc = %ld\n Type c to continue : ",LinksToAdd,time,big_escape_counter); scanf("%s",junk); }
  
    big_escape_counter++;
    escape_counter = 0; failed_to_find_a_new_link = 0;
     
    /* Select the first person */ 
   
    do {
      person1 = genrand_real2() * N + 1;
      escape_counter++;

      SeeksNewPartner = 1;

      if (Status[person1] < 0) SeeksNewPartner = 0;  // Dead people don't seek partners
      
      if (NumLinks[person1] >= MaxLinks) { SeeksNewPartner = 0;  if (time >= 10000000) printf("  --> No longer seeking after test for NumLinks (=%ld) >= MaxLinks (%ld)\n",NumLinks[person1],MaxLinks); }

      if (HighRisk[person1] == 0) {
         if (genrand_real2() > reduced_linkage_low_risk) SeeksNewPartner = 0; // Low risk people have reduced chance of taking on a partner
      }
      if ( NumLinks[person1] >= 1) {
         if (Sex[person1] == 1) age_specific_prob_conc = female_concurrency*(45.0 - Age[person1])/(45.0 - 15.0);
                           else age_specific_prob_conc = male_concurrency*(45.0 - Age[person1])/(45.0 - 15.0);
        if (genrand_real2() > age_specific_prob_conc) SeeksNewPartner = 0;  // Already linked people are unliklely (assuming low concurrency values)  to seek a new partner
      }
      if (time == 10000000) printf("Inside AddNewLinks: person1 = %ld, HighRisk = %d, NumLinks = %ld, SeeksNewParnter = %ld\n",person1,HighRisk[person1],NumLinks[person1], SeeksNewPartner);

      if (escape_counter >= 50) {
        if (plt_warn==1) fprintf(WarningsFile,"Warning at time %ld: Unable to find anyone who could be linked to anyone (escape_counter = %ld)\nType Ctrl-C to quit : ",time,escape_counter);
        failed_to_find_a_new_link = 1;
      }

    } while ( SeeksNewPartner == 0 && failed_to_find_a_new_link == 0);
    // Logic is to go back and loop through with a new person1 if the current person1 is unsuitable (e.g., person1 is dead or already has the maximum number of partners)
    //   However, we exit this while loop if there were many, many failed attempts to find a suitable person (in which case we set failed_to_find_a_new_link = 1)
    
    escape_counter = 0;
    
    /* Select the second person */
   
    if (failed_to_find_a_new_link == 0) {
      failed_to_find_a_partner = 0;
      do {
        person2 = genrand_real2() * N + 1;
        Unsuitable = 0; // Assume that person2 is suitable until proven otherwise
        
        for (j=1; j<=NumLinks[person2]; j++) {
           if (Links[person2][j] == person1) Unsuitable = 1; // Don't link with someone you are already linked with
        }
       
        if (Status[person2] < 0) Unsuitable = 1; // Cannot link with a dead person
                                  
        if (person2 == person1) Unsuitable = 1; // Cannot link with one's self
   
        if (HighRisk[person2] == 0) {
          if (genrand_real2() > reduced_linkage_low_risk) Unsuitable = 1; // Low risk people have reduced chance of taking on a partner
        }
 
        if (Sex[person2] == Sex[person1]) Unsuitable = 1; // This version assumes heterosexual pairings
    
        if (NumLinks[person2] >= MaxLinks) Unsuitable = 1; 

        if (NumLinks[person2] >= 1) {
          if (Sex[person2] == 1)  age_specific_prob_conc = female_concurrency*(45.0 - Age[person2])/(45.0 - 15.0); // probability of having >1 partner decreases with age
                            else  age_specific_prob_conc = male_concurrency*(45.0 - Age[person2])/(45.0 - 15.0);
          if (genrand_real2() > age_specific_prob_conc) Unsuitable = 1;  // Already linked people are unlikely (assuming low concurrency values) to take on a new partner
        }
      
        if (age_dependent_mixing == 1 && Unsuitable == 0) {
           if (Sex[person1] == 2) {
              Am = Age[person1]; Af = Age[person2];
           } else {
              Am=Age[person2]; Af = Age[person1];
           }
           female_pref_older_men = 2.0;
           average_age = (Age[person1] + Age[person2])/2.0;
           age_diff_with_offset = Am - (Af + female_pref_older_men*(1.0 + (average_age-15.0)/15.0)); // Preference increases with age.  15 y.o. women pref men 17, women at 30 prefer men at 34, 

           age_diff_50perc_less_likely = 4.0; 
           age_diff_index =  exp(-fabs(age_diff_with_offset)/ age_diff_50perc_less_likely ); // =0 if same age (~ M-F age diff), = 0.5 if M - F - diff is 50% level, --> 0if very different ages
           // age_diff_index =  jmpow(fabs(age_diff_with_offset/age_diff_50perc_less_likely),-2.0); // =0 if same age (~ M-F age diff), = 0.5 if M - F - diff is 50% level, --> 0if very different ages
           rand_num = genrand_real2();
           if (rand_num > age_diff_index) Unsuitable = 1; // Probability of partnership decreases with age difference
        }
     
        escape_counter++;
        if (escape_counter >= 50) {
          failed_to_find_a_partner = 1; Unsuitable = 1;
         }
      } while ( Unsuitable == 1 && failed_to_find_a_partner == 0 );
      if (failed_to_find_a_partner == 0) {
        NumLinks[person1] = NumLinks[person1] + 1;
        NumLinks[person2] = NumLinks[person2] + 1;
        NumPartners[person1] = NumPartners[person1] + 1;
        NumPartners[person2] = NumPartners[person2] + 1;
        if (NumLinks[person1] > MaxLinks) printf("Big error in AddLinks for person1 at time %ld: NumLinks[%ld] = %ld\n",time,person1,NumLinks[person1]); 
        if (NumLinks[person2] > MaxLinks) printf("Big error in AddLinks for person2 at time %ld: NumLinks[%ld] = %ld\n",time,person2,NumLinks[person2]); 
        Links[person1][NumLinks[person1]] = person2;
        Links[person2][NumLinks[person2]] = person1;
        TimeLinked[person1][NumLinks[person1]] = time;
        TimeLinked[person2][NumLinks[person2]] = time;
        if (person1 == 3100002 || person2 == 3100002) {
          printf("Agent %ld linked to agent %ld at time %ld:\n",person1,person2,time);
          printf("New link established between persons %ld (S=%ld, A=%lf) and %ld (S=%ld, A=%lf), age_index= %lf,rand=%lf\n",person1,Sex[person1],Age[person1],person2,Sex[person2],Age[person2],age_diff_index,rand_num);
          printf("TimeLinked[%ld][%ld] = %ld\tTimeLinked[%ld][%ld] = %ld\n",person1,person2,TimeLinked[person1][NumLinks[person1]],person2,person1,TimeLinked[person2][NumLinks[person2]]);
        }
        /* if (person1 < person2) fprintf(RelationshipRegistryFile,"%ld\t%ld\t%ld\tNew\n",person1,person2,time);
                            else fprintf(RelationshipRegistryFile,"%ld\t%ld\t%ld\tNew\n",person2,person1,time); */
        new_links_added++;
      }
    }
    if (failed_to_find_a_new_link == 1 || failed_to_find_a_partner == 1) {
       //StopEarly = 1;
       num_cases_no_link_found++;
    }
    //printf("Type c to continue :"); scanf("%s",junk);
  } while ((new_links_added < LinksToAdd) && (big_escape_counter < 2*N));
  if (big_escape_counter >= 2*N) { if (plt_warn==1) fprintf(WarningsFile,"Warning at time %ld: Program gave up trying to find new links after %ld attempts\n",time,big_escape_counter); scanf("%s",junk); }
}

/* The following function calculates a^b (i.e., "a" raised to the power of "b"). */
double jmpow(double a, double b) {
  return(exp(b*log(a)));
}


/*************** Function GetGammaDelay ******************************/
/* This function returns a random number from the gamma distribution */
/* with a scale parameter of theta and a shape parameter of  k       */
/* Inputs:						  	     */
/*          theta -- the shape parameter			     */
/*          k     -- the scale parameter			     */
/* Ouput: A random number (delay time) 				     */
/*        with mean k*theta and variance k*theta^2   		     */
/* Dependencies: This function calls the MT random number generator  */
/*********************************************************************/

double GetGammaDelay(double k, double theta)
{
  double lnUsum, Final_Value;
  double V3m2, V3m1, V3m;
  double k_int, delta, epsilon, epsilon_m, nu_m, v0;
  long i,m,step, condition_met;
  step = 1;
  condition_met = 0;
  k_int = (int) k;
  delta = k - k_int;
  v0 = exp(1.0)/(exp(1.0) + delta);
  do {
    m = 1;
    V3m = genrand_real2();
    V3m1 = genrand_real2();
    V3m2 = genrand_real2();
    //printf("At step %ld: V3m-2 = %lf, V3m-1 = %lf, V3m = %lf\n",step,V3m2,V3m1,V3m);
    if (V3m2 <= v0) {
      epsilon_m = pow(V3m1,1/delta);
      nu_m = V3m*pow(epsilon_m,delta-1);
    } else {
      epsilon_m = 1.0 - log(V3m1);
      nu_m = V3m*exp(-epsilon_m);
    }
    if (nu_m > pow(epsilon_m,delta-1)*exp(-epsilon_m)) {
      m = m + 1;
    } else {
      condition_met = 1;
    }
    step++;
  } while (condition_met == 0);
  epsilon = epsilon_m;
  lnUsum = 0.0;
  for (i = 1; i<= k_int; i++) {
    lnUsum = lnUsum + log(genrand_real2());
  }
  Final_Value = theta*(epsilon - lnUsum);
  return(Final_Value);
}


void GetParameters(void)
{
 plt_out = 1; 
 double Nfloat;
 rewind(ParFile);

 fscanf(ParFile,"%d %s",&gradual_tx_incr, descript); fprintf(ParOut,"# %d %s\n",gradual_tx_incr, descript);
 fscanf(ParFile,"%d %s",&randomized_parameters, descript);    fprintf(ParOut,"# %d %s\n",randomized_parameters, descript);
 fscanf(ParFile,"%lf %s",&fold_incr_prob_with_targeting, descript);    fprintf(ParOut,"# %lf %s\n",fold_incr_prob_with_targeting, descript);
 fscanf(ParFile,"%lf %s",&lowest_tx, descript);    fprintf(ParOut,"# %lf %s\n",lowest_tx, descript);
 fscanf(ParFile,"%lf %s",&highest_tx, descript);    fprintf(ParOut,"# %lf %s\n",highest_tx, descript);
 fscanf(ParFile,"%lf %s",&tx_incr, descript);    fprintf(ParOut,"# %lf %s\n",tx_incr, descript);
 fscanf(ParFile,"%d %s",&first_tx_strategy, descript);    fprintf(ParOut,"# %d %s\n",first_tx_strategy, descript);
 fscanf(ParFile,"%d %s",&last_tx_strategy, descript);    fprintf(ParOut,"# %d %s\n",last_tx_strategy, descript);
 fscanf(ParFile,"%lf %s", &Nfloat, descript);
   N0 = (long) (Nfloat + 0.5);
      fprintf(ParOut,"# %ld %s\n",N0, descript);
   if (N0 > MaxN) {
     printf("Error: N0 (%ld) exceeds MaxN (%ld).  Hit ctrl-C to quit!\n",N0, MaxN);
     scanf("%s",junk);
   }
   OrigN = N0;
 fscanf(ParFile,"%ld %s",&Infected0, descript);                 fprintf(ParOut,"# %ld %s\n",Infected0, descript);
 fscanf(ParFile,"%d %s",&time_varying_mean_degree, descript);   fprintf(ParOut,"# %d %s\n",time_varying_mean_degree, descript);
 fscanf(ParFile,"%lf %s",&d_rate, descript);                    fprintf(ParOut,"# %lf %s\n",d_rate, descript);
 fscanf(ParFile,"%ld %s",&TrackedAgentsGroup1start, descript);  fprintf(ParOut,"# %ld %s\n",TrackedAgentsGroup1start, descript);
 fscanf(ParFile,"%ld %s",&TrackedAgentsGroup1end, descript);    fprintf(ParOut,"# %ld %s\n",TrackedAgentsGroup1end, descript);
 fscanf(ParFile,"%ld %s",&TrackedAgentsGroup2start, descript);  fprintf(ParOut,"# %ld %s\n",TrackedAgentsGroup2start, descript);
 fscanf(ParFile,"%ld %s",&TrackedAgentsGroup2end, descript);    fprintf(ParOut,"# %ld %s\n",TrackedAgentsGroup2end, descript);
 fscanf(ParFile,"%lf %s",&PercentVaccinated, descript);    fprintf(ParOut,"# %lf %s\n",PercentVaccinated, descript);
 fscanf(ParFile,"%lf %s",&RR_Vaccinated, descript);        fprintf(ParOut,"# %lf %s\n",RR_Vaccinated, descript);
 fscanf(ParFile,"%lf %s",&VaccineDuration, descript);      fprintf(ParOut,"# %lf %s\n",VaccineDuration, descript);
 fscanf(ParFile,"%lf %s",&PercentResistantToVaccine, descript);    fprintf(ParOut,"# %lf %s\n",PercentResistantToVaccine, descript);
 fscanf(ParFile,"%lf %s",&AverageDegree, descript);        fprintf(ParOut,"# %lf %s\n",AverageDegree, descript);
 fscanf(ParFile,"%lf %s",&MinDuration, descript);          fprintf(ParOut,"# %lf %s\n",MinDuration, descript);
 fscanf(ParFile,"%lf %s",&MaxDuration, descript);          fprintf(ParOut,"# %lf %s\n",MaxDuration, descript);
 fscanf(ParFile,"%lf %s",&prob_sex, descript);             fprintf(ParOut,"# %lf %s\n",prob_sex, descript);
 fscanf(ParFile,"%d %s", &age_dependent_sex, descript);    fprintf(ParOut,"# %d %s\n",age_dependent_sex, descript);
 fscanf(ParFile,"%d %s", &age_dependent_mixing, descript); fprintf(ParOut,"# %d %s\n",age_dependent_mixing, descript);
 fscanf(ParFile,"%ld %s" ,&MaxLinks, descript);            fprintf(ParOut,"# %ld %s\n",MaxLinks, descript);
 fscanf(ParFile,"%lf %s", &Start_Spontaneous_Treatment, descript);    fprintf(ParOut,"# %lf %s\n", Start_Spontaneous_Treatment, descript);
 fscanf(ParFile,"%lf %s", &Start_TasP_Campaign, descript);         fprintf(ParOut,"# %lf %s\n", Start_TasP_Campaign, descript);
 fscanf(ParFile,"%lf %s", &prob_dropout, descript);                fprintf(ParOut,"# %lf %s\n", prob_dropout, descript);
 fscanf(ParFile,"%lf %s", &PercentNonTargetedYear1, descript);     fprintf(ParOut,"# %lf %s\n", PercentNonTargetedYear1, descript);
 fscanf(ParFile,"%lf %s", &percent_under_care, descript);          fprintf(ParOut,"# %lf %s\n", percent_under_care, descript);
 fscanf(ParFile,"%lf %s", &VL_After_Treatment, descript);          fprintf(ParOut,"# %lf %s\n", VL_After_Treatment, descript);
 fscanf(ParFile,"%lf %s", &Start_SexReduction_Campaign, descript); fprintf(ParOut,"# %lf %s\n", Start_SexReduction_Campaign, descript);
 fscanf(ParFile,"%lf %s", &Reduction_Mean_Degree, descript);       fprintf(ParOut,"# %lf %s\n", Reduction_Mean_Degree, descript);
 fscanf(ParFile,"%ld %s", &Start_Condom_Campaign, descript);       fprintf(ParOut,"# %ld %s\n", Start_Condom_Campaign, descript);
 fscanf(ParFile,"%lf %s", &Percent_using_condums_after_condom_campaign, descript);    fprintf(ParOut,"# %lf %s\n", Percent_using_condums_after_condom_campaign, descript);
 fscanf(ParFile,"%lf %s", &Start_Faithfulness_Campaign, descript);    fprintf(ParOut,"# %lf %s\n", Start_Faithfulness_Campaign, descript);
 fscanf(ParFile,"%lf %s", &Increased_Duration, descript);  fprintf(ParOut,"# %lf %s\n", Increased_Duration, descript);
 fscanf(ParFile,"%ld %s", &Inf_Rate_Function, descript);   fprintf(ParOut,"# %ld %s\n", Inf_Rate_Function, descript);
 fscanf(ParFile,"%lf %s",&InfRateBaseline, descript);      fprintf(ParOut,"# %4.2e %s\n",InfRateBaseline, descript);
 fscanf(ParFile,"%lf %s",&InfRateExponent, descript);      fprintf(ParOut,"# %lf %s\n",InfRateExponent, descript);
 fscanf(ParFile,"%lf %s",&MaxInfRate, descript);           fprintf(ParOut,"# %lf %s\n",MaxInfRate, descript);
 fscanf(ParFile,"%lf %s",&VHalfMaxInfRate, descript);      fprintf(ParOut,"# %lf %s\n",VHalfMaxInfRate, descript);
 fscanf(ParFile,"%lf %s",&HillCoeffInfRate, descript);     fprintf(ParOut,"# %lf %s\n",HillCoeffInfRate, descript);
 fscanf(ParFile,"%lf %s",&BirthRate, descript);            fprintf(ParOut,"# %lf %s\n",BirthRate, descript);
 fscanf(ParFile,"%lf %s",&NaturalDeathRate, descript);     fprintf(ParOut,"# %lf %s\n",NaturalDeathRate, descript);
 fscanf(ParFile,"%lf %s",&shape_parameter, descript);      fprintf(ParOut,"# %lf shape_parameter(%s)\n",shape_parameter, descript);
 fscanf(ParFile,"%lf %s",&death_rate_constant, descript);  fprintf(ParOut,"# %lf death_rate_constant (%s)\n",death_rate_constant, descript);
 fscanf(ParFile,"%lf %s",&death_rate_exponent, descript);  fprintf(ParOut,"# %lf death_rate_exponent (%s)\n",D50, descript);
 fscanf(ParFile,"%lf %s",&Dmax, descript);    fprintf(ParOut,"# %lf Dmax (%s)\n",Dmax, descript);
 fscanf(ParFile,"%lf %s",&D50, descript);     fprintf(ParOut,"# %lf D50 (%s)\n",D50, descript);
 fscanf(ParFile,"%lf %s",&Dk, descript);      fprintf(ParOut,"# %lf Dk (%s)\n",Dk, descript);
 fscanf(ParFile,"%lf %s",&V0, descript);      fprintf(ParOut,"# %lf %s\n",V0, descript);
 fscanf(ParFile,"%lf %s",&V_peak, descript);  fprintf(ParOut,"# %lf %s\n",V_peak, descript);
 fscanf(ParFile,"%lf %s",&t_peak, descript);  fprintf(ParOut,"# %lf %s\n",t_peak, descript);
 fscanf(ParFile,"%lf %s",&t_acute, descript);  fprintf(ParOut,"# %lf %s\n",t_acute, descript);
 fscanf(ParFile,"%lf %s",&AverageLogSP0, descript);    fprintf(ParOut,"# %lf %s\n",AverageLogSP0, descript);
 fscanf(ParFile,"%lf %s",&VarianceLogSP0, descript);    fprintf(ParOut,"# %lf %s\n",VarianceLogSP0, descript);
 fscanf(ParFile,"%lf %s",&MutationVariance, descript);  fprintf(ParOut,"# %lf %s\n",MutationVariance, descript);
 fscanf(ParFile,"%ld %s",&tfinal, descript);          fprintf(ParOut,"# %ld %s\n",tfinal, descript);
 fscanf(ParFile,"%lf %s",&prog_rate, descript);       fprintf(ParOut,"# %lf %s\n",prog_rate, descript);
 fscanf(ParFile,"%lf %s",&Heritability, descript);    fprintf(ParOut,"# %lf %s\n",Heritability, descript);
 fscanf(ParFile,"%ld %s",&tprintHerit, descript);     fprintf(ParOut,"# %ld %s\n",tprintHerit, descript);
 fscanf(ParFile,"%ld %s",&VL_print_time, descript);   fprintf(ParOut,"# %ld %s\n",VL_print_time, descript);
 fscanf(ParFile,"%lf %s",&VL_print_lower, descript);  fprintf(ParOut,"# %lf %s\n",VL_print_lower, descript);
 fscanf(ParFile,"%lf %s",&VL_print_upper, descript);  fprintf(ParOut,"# %lf %s\n",VL_print_upper, descript);
 fscanf(ParFile,"%lf %s",&VL_print_interval, descript);   fprintf(ParOut,"# %lf %s\n",VL_print_interval, descript);
 fscanf(ParFile,"%ld %s",&random_number_seed, descript);  fprintf(ParOut,"# %ld %s\n",random_number_seed, descript);
 fscanf(ParFile,"%d %s",&printOutput, descript);      fprintf(ParOut,"# %d %s\n",printOutput, descript);
 fscanf(ParFile,"%lf %s",&perc_tx_start, descript);   fprintf(ParOut,"# %lf %s\n",perc_tx_start, descript);
 fscanf(ParFile,"%lf %s",&yearly_incr_tx, descript);  fprintf(ParOut,"# %lf %s\n",yearly_incr_tx, descript);
 fscanf(ParFile,"%ld %s",&replicates, descript);      fprintf(ParOut,"# %ld %s\n",replicates, descript);
 fscanf(ParFile,"%lf %s",&increased_duration_per_year, descript);    fprintf(ParOut,"# %lf %s\n",increased_duration_per_year, descript);
 fscanf(ParFile,"%lf %s",&rate_stops_being_concurrent, descript);    fprintf(ParOut,"# %lf %s\n",rate_stops_being_concurrent, descript);
 fscanf(ParFile,"%d %s",&plt_GS, descript);    fprintf(ParOut,"# %d %s\n",plt_GS, descript);
 fscanf(ParFile,"%d %s",&plt_NS, descript);    fprintf(ParOut,"# %d %s\n",plt_GS, descript);
 fscanf(ParFile,"%d %s",&plt_out, descript);   fprintf(ParOut,"# %d %s\n",plt_out, descript);
 fscanf(ParFile,"%d %s",&plt_H, descript);     fprintf(ParOut,"# %d %s\n",plt_H, descript);
 fscanf(ParFile,"%d %s",&plt_R, descript);     fprintf(ParOut,"# %d %s\n",plt_R, descript);
 fscanf(ParFile,"%d %s",&plt_AH, descript);    fprintf(ParOut,"# %d %s\n",plt_AH, descript);
 fscanf(ParFile,"%d %s",&plt_AM, descript);    fprintf(ParOut,"# %d %s\n",plt_AM, descript);
 fscanf(ParFile,"%d %s",&plt_AD, descript);      fprintf(ParOut,"# %d %s\n",plt_AD, descript);
 fscanf(ParFile,"%d %s",&plt_VL, descript);      fprintf(ParOut,"# %d %s\n",plt_VL, descript);
 fscanf(ParFile,"%d %s",&plt_spVL, descript);    fprintf(ParOut,"# %d %s\n",plt_spVL, descript);
 fscanf(ParFile,"%d %s",&plt_warn, descript);    fprintf(ParOut,"# %d %s\n",plt_warn, descript);
 fscanf(ParFile,"%lf %s",&MaxAgeSex, descript);  fprintf(ParOut,"# %lf %s\n",MaxAgeSex, descript);
 fscanf(ParFile,"%lf %s",&circum_prob, descript);       fprintf(ParOut,"# %lf %s\n",circum_prob, descript);
 fscanf(ParFile,"%lf %s",&RR_female_recipient, descript);    fprintf(ParOut,"# %lf %s\n",RR_female_recipient, descript);
 fscanf(ParFile,"%lf %s",&RR_circum, descript);        fprintf(ParOut,"# %lf %s\n",RR_circum, descript);
 fscanf(ParFile,"%lf %s",&RR_youth, descript);         fprintf(ParOut,"# %lf %s\n",RR_youth, descript);
 fscanf(ParFile,"%lf %s",&RR_condom, descript);        fprintf(ParOut,"# %lf %s\n",RR_condom, descript);
 fscanf(ParFile,"%lf %s",&condom_prob16, descript);    fprintf(ParOut,"# %lf %s\n",condom_prob16, descript);
 fscanf(ParFile,"%lf %s",&age_condom_use_halves, descript);    fprintf(ParOut,"# %lf %s\n",age_condom_use_halves, descript);
 fscanf(ParFile,"%lf %s",&number_tx_spontaneous, descript);    fprintf(ParOut,"# %lf %s\n",number_tx_spontaneous, descript);
 fscanf(ParFile,"%lf %s",&prob_sex_in_AIDS, descript);    fprintf(ParOut,"# %lf %s\n",prob_sex_in_AIDS, descript);
 fscanf(ParFile,"%lf %s",&male_concurrency, descript);    fprintf(ParOut,"# %lf %s\n",male_concurrency, descript);
 fscanf(ParFile,"%lf %s",&female_concurrency, descript);  fprintf(ParOut,"# %lf %s\n",female_concurrency, descript);
 fscanf(ParFile,"%lf %s",&LongevityFudgeFactor, descript);     fprintf(ParOut,"# %lf %s\n",LongevityFudgeFactor, descript);
 fscanf(ParFile,"%lf %s",&PopGrowth, descript);           fprintf(ParOut,"# %lf %s\n",PopGrowth, descript);
 fscanf(ParFile,"%ld %s",&HoneymoonDays, descript);       fprintf(ParOut,"# %ld %s\n",HoneymoonDays, descript);
 fscanf(ParFile,"%lf %s",&IncreasedProbSexHoneymoonDays, descript);    fprintf(ParOut,"# %lf %s\n",IncreasedProbSexHoneymoonDays, descript);
 fscanf(ParFile,"%lf %s",&Max_VL_AIDS, descript);         fprintf(ParOut,"# %lf %s\n",Max_VL_AIDS, descript);
 fscanf(ParFile,"%lf %s",&VL_Incr_AIDS, descript);        fprintf(ParOut,"# %lf %s\n",VL_Incr_AIDS, descript);
 fscanf(ParFile,"%lf %s",&prob_cd4_rebound, descript);    fprintf(ParOut,"# %lf %s\n",prob_cd4_rebound, descript);
 fscanf(ParFile,"%lf %s",&prob_cd4_rebound_further, descript); fprintf(ParOut,"# %lf %s\n",prob_cd4_rebound_further, descript);
 fscanf(ParFile,"%d %s",&enhanced_progression_age, descript);  fprintf(ParOut,"# %d %s\n",enhanced_progression_age, descript);
 fscanf(ParFile,"%lf %s",&daily_prob_diagnosis, descript);     fprintf(ParOut,"# %lf %s\n",daily_prob_diagnosis, descript);
 fscanf(ParFile,"%ld %s",&delay_diag_tx, descript);            fprintf(ParOut,"# %ld %s\n",delay_diag_tx, descript);
 fscanf(ParFile,"%ld %s",&delay_test_antibodies, descript);    fprintf(ParOut,"# %ld %s\n",delay_test_antibodies, descript);
 fscanf(ParFile,"%ld %s",&delay_test_antibodies, descript);    fprintf(ParOut,"# %ld %s\n",delay_test_antibodies, descript);
 fscanf(ParFile,"%lf %s",&life_expectancy, descript);     fprintf(ParOut,"# %lf %s\n",life_expectancy, descript);
 fscanf(ParFile,"%lf %s",&cost_hiv_treated, descript);    fprintf(ParOut,"# %lf %s\n",cost_hiv_treated, descript);
 fscanf(ParFile,"%lf %s",&cost_cd4_gt_350, descript);     fprintf(ParOut,"# %lf %s\n",cost_cd4_gt_350, descript);
 fscanf(ParFile,"%lf %s",&cost_cd4_200_350, descript);    fprintf(ParOut,"# %lf %s\n",cost_cd4_200_350, descript);
 fscanf(ParFile,"%lf %s",&cost_cd4_lt_200, descript);     fprintf(ParOut,"# %lf %s\n",cost_cd4_lt_200, descript);
 fscanf(ParFile,"%lf %s",&cost_died_AIDS, descript);      fprintf(ParOut,"# %lf %s\n",cost_died_AIDS, descript);
 
}


void PrintHeaders()
{
 if (plt_NS==1) fprintf(NetworkStatsOutputFile,"\tRepl\tTime\tSus\tPos\tAlv\tFem\tMale\tu25\tmidd\to50\tmd\tmd_u25\tmd_Mid\tmd50\tmd_f\tmd_m\tconc\tc_f\tc_m\tc_u25\tc2550\tc50\tTx\tTarTx\tdAIDS\tdNat\tsusc_d\tnew_inf\tsusc_d_u50\tnew_inf50");
 if (plt_NS==1) fprintf(NetworkStatsOutputFile,"\tf25\tm25\tfMid\tmMid\tf50\tm50,\tif\tim\ti25\tiMid\ti50\tif25\tim25\tifMid\timMid\tif50\tim50\t");
 if (plt_NS==1) fprintf(NetworkStatsOutputFile,"\tInf_Tx_u25\tInf_Tx_Mid\tInf_Tx50\tInf_Tx_f\tInf_Tx_m\tDur\tDur18\tDur25\tDurMid\tDur50\tTotTx\tNoTarg\t");
 if (plt_NS==1) fprintf(NetworkStatsOutputFile,"\tAvePartMen20\tAvePartMen25\tAvePartMen30\tAvePartMen40\tAvePartMen50\t");
 if (plt_NS==1) fprintf(NetworkStatsOutputFile,"\tAvePartWomen20\tAvePartWomen25\tAvePartWomen30\tAvePartWomen40\tAvePartWomen50\t");
 if (plt_NS==1) fprintf(NetworkStatsOutputFile,"\tconcF24\tconcM25\tconcFMid\tconcMMid\tconcF50\tconcM50\t");
 if (plt_NS==1) fprintf(NetworkStatsOutputFile,"\tmd_F24\tmd_M25\tmd_FMid\tmd_MMid\tmd_F50\tmd_M50\tPrev50\t");
 if (plt_NS==1) fprintf(NetworkStatsOutputFile,"\tMP_m25\tMP_mMid\tMP_m50\tMP_f25\tMP_fMid\tMP_f50\t");
 if (plt_NS==1) fprintf(NetworkStatsOutputFile,"\tsusc_d25\tsusc_dMid\tsusc_d_o50\tnewinfs25\tnewinfsMid\tnewinfs_o50\tAge\tAgeS\tAgeI\n");
 if (plt_NS==1) fprintf(RelLengthFile,"\ttime\tDur25\tDur25_30\tDur30\n");
 if (plt_GS==1) fprintf(GrandSummaryFile,"Strat\tu25\ttx\tActTx\tyearly_incr\trepl\tdAIDS\tHIV+\tAllTx\tNotAllTx\tUntreated\tNotTargeted\tPerDaysTx\tsusc_d\tnew_infs\tcd4_1days\tcd4_2days\tcd4_3days\tcd4_4days\tcd4_5days\ttx_days\t");
 if (plt_GS==1) fprintf(GrandSummaryFile,"D00\tD01\tD02\tD03\tD05\tD07\tD08\tD10\tD15\tD20\n");
 if (plt_out==1) fprintf(OutputFile,"\n\nRepl\ttime\tSusc\tInfected\tVaccinated\tVaccinatedInfected\tResistInfected\tInfectedTreated\tDead\tdAIDS\tdNat\tLogV\t\tVstd\t\tSPVL\t\tSet_std\t\td_ave\t\td_std\t\tLinks\t\tV1\t\tda1\t\tSts1\n");
 if (plt_H==1) fprintf(HeritOutputFile,"Repl\ttime\tIndex\tDuration\tDonorsIndex\tGen\tDonorsGen\tTimeInf\tVL\tTimeToAIDS\tDonorsTimeToAIDS\tDonorsVLatTrans\tLogSP0\tDonors_LogSP0\tVirContLogSP0\tDonors_VirContLogSP0\tEnvContLogSP0\tDonors_EnvContLogSP0\td_acute\tDonors_d_acute\tTotal_Time_Inf\tDonors_Total_Time_Inf_At_Trans\n");
 //printf("Initial Infected Population\nTime\tIndex\tGen\tVL\t\tLogSetPoint\tVContrSP0\tEnvContrSP0\td_acute\tTotal_Time_Inf\n");
 if (plt_R==1) fprintf(PatientRegistryOutputFile,"Repl\tTime\tIndex\tAge\tDuration\tNumRecipients\tDonorsIndex\tGen\tDonorsGen\tTimeInf\tVL\tDonorsVLatTrans\tLogSP0\tDonors_LogSP0\tVirContrLogSP0\tDonors_VirContLogSP0\tEnvContLogSP0\td_acute\tDonors_d_acute\tTotal_Time_Inf\tDonors_Total_Time_Inf_At_Trans\tCD4\tCD4_TimeToAIDS\tCD4_TimeToAIDS_exp_1\tCD4_TimeToAIDS_exp_2\tCD4_TimeToAIDS_exp_3\tCD4_TimeToAIDS_exp_4\tCD4_initial_value\tStatus\n");
 if ( (TrackedAgentsGroup1end - TrackedAgentsGroup1start > 0) || ( TrackedAgentsGroup2end - TrackedAgentsGroup2start > 0) )
     fprintf(AgentHistoryFile,"time\tAgent\tSPVL\t\tStatus\ttInf\t\ts\t\ttx\ttx_d\ttx_dA\tAge\tPdrop\t\tV\t\tCD4\n");


}

void PrintHeritabilityStats(void)
{
   long i;
   if ((labs(time - tprintHerit) < 0.01) || (StopEarly == 1) ) {  //  print data for heritabilities graphs
     for (i=1; i<= N; i++) {
       if ((Status[i] == 1) && (Generation[i] > 1)) {
         if (plt_H==1) fprintf(HeritOutputFile,"%ld\t%ld\t%6.5e\t%ld\t%ld\t%ld\t%6.5e\t%6.5e\t%6.5e\t%6.5e\t%6.5e\t%6.5e\t%6.5e\t%6.5e\t%6.5e\t%6.5e\t%6.5e\t%6.5e\t%6.5e\n",
           time,i,Duration[i],Donors_Index[i],Generation[i],Donors_Generation[i],Time_Inf[i],V[i],Donors_V[i],LogSetPoint[i],Donors_LogSetPoint[i],ViralContribToLogSP0[i],Donors_ViralContribToLogSP0[i],EnvirContribToLogSP0[i],Donors_EnvirContribToLogSP0[i],d_acute[i],Donors_d_acute[i],time-Time_Inf[i],Donors_Total_Time_Inf_At_Trans[i]);
       }
     }
   }
}

void PrintStats(void)
{
   long vaccinated, vaccinated_infected, resistant_to_vaccine;
   double md = 0.0, md_f = 0.0, md_m = 0.0;
   double md_u25 = 0.0, md_Mid = 0.0, md50 = 0.0;
   double md_f25 = 0.0, md_fMid = 0.0, md_f50 = 0.0;
   double md_m25 = 0.0, md_mMid = 0.0, md_m50 = 0.0;
   double concurrency = 0.0, conc_f = 0.0, conc_m = 0.0;
   double conc_u25 = 0.0, conc_Mid = 0.0, conc50 = 0.0;
   double conc_f25 = 0.0, conc_fMid = 0.0, conc_f50 = 0.0;
   double conc_m25 = 0.0, conc_mMid = 0.0, conc_m50 = 0.0;
   double SumLinks_alive = 0.0,  SumLinks_f = 0.0,  SumLinks_m = 0.0,  SumConc = 0.0,  SumConc_f = 0.0,  SumConc_m = 0.0;
   double  SumLinks_u25 = 0.0,  SumLinks_Mid = 0.0,  SumLinks50 = 0.0;
   double  SumLinks_f25 = 0.0,  SumLinks_fMid = 0.0,  SumLinks_f50 = 0.0;
   double  SumLinks_m25 = 0.0,  SumLinks_mMid = 0.0,  SumLinks_m50 = 0.0;
   double  SumConc_u25 = 0.0,  SumConc_Mid = 0.0, SumConc50 = 0.0;
   double  SumConc_f25 = 0.0,  SumConc_fMid = 0.0, SumConc_f50 = 0.0;
   double  SumConc_m25 = 0.0,  SumConc_mMid = 0.0, SumConc_m50 = 0.0;
   double num_alive = 0.0,  num_f = 0.0,  num_m = 0.0,  num25 = 0.0,  numMid = 0.0,  num50 = 0.0, NumUnlinked = 0.0, PercentUnlinked = 0.0;
   double f25 = 0.0, fMid = 0.0, f50 = 0.0, m25 = 0.0, mMid = 0.0, m50 = 0.0;
   double sa_f25 = 0.0, sa_fMid = 0.0, sa_f50 = 0.0, sa_m25 = 0.0, sa_mMid = 0.0, sa_m50 = 0.0;
   double inf_f = 0.0, inf_m = 0.0, inf25 = 0.0, infMid = 0.0, inf50 = 0.0, inf_f25 = 0.0, inf_m25 = 0.0, inf_fMid = 0.0, inf_mMid = 0.0, inf_f50 = 0.0, inf_m50 = 0.0;
   double MP_f25 = 0, MP_m25 = 0;
   double MP_fMid = 0, MP_mMid = 0;
   double MP_f50 = 0, MP_m50 = 0;
   long Inf_Tx = 0, Inf_Tx_m = 0, Inf_Tx_f = 0, Inf_Tx_u25 = 0, Inf_Tx_Mid = 0, Inf_Tx50 = 0; // HIV+ persons of various types that receive treatment
   long Inf_Tx_m25 = 0, Inf_Tx_mMid = 0, Inf_Tx_m50 = 0; // Number of HIV+ persons receiving treatment
   long Inf_Tx_f25 = 0, Inf_Tx_fMid = 0, Inf_Tx_f50 =0; // Number of HIV+ persons receiving treatment
   double DurSum = 0.0, DurSum18 = 0.0, DurSum25 = 0.0, DurSumMid = 0.0, DurSum50 = 0.0;
   double MeanDur = 0.0, MeanDur18 = 0.0, MeanDur25 = 0.0, MeanDurMid = 0.0, MeanDur50 = 0.0;
   long   DurCount = 0, DurCount18 = 0, DurCount25 = 0, DurCountMid = 0, DurCount50 = 0; 
   long NumUnder50 = 0, InfUnder50 = 00; 
   long MP_men25 = 0, MP_menMid = 0, MP_men50 = 0;
   long MP_women25 = 0, MP_womenMid = 0, MP_women50 = 0;
   long DiscordantCouples = 0;
   double PrevUnder50;
   Vsum = 0.0; d_sum = 0.0; G_sum = 0.0; set_sum = 0.0; Susceptible = 0; Dead = 0; DiedNat = 0; DiedAIDS = 0; TotalDiedAIDSAfterTasP = 0;
   Infected = 0; SumLinks = 0.0; AveLinks = 0.0;
   vaccinated = 0; vaccinated_infected = 0;resistant_to_vaccine = 0; InfectedTreated = 0; 
   long InfectedTreatedTargeted = 0;;

   long Num25partnerships = 0, Num25to30partnerships = 0, Num30partnerships = 0; 
   long SumLength25partnerships = 0, SumLength25to30partnerships = 0, SumLength30partnerships = 0; 
   double AveLength25partnerships = 0.0, AveLength25to30partnerships = 0.0, AveLength30partnerships = 0.0; 
   long num_diag = 0;
   double SumAge = 0.0, SumAgeSusc = 0.0, SumAgeInf= 0.0, Alive = 0.0;
   double AveAge = 0.0, AveAgeSusc = 0.0, AveAgeInf= 0.0;

   // Observed length of partnerships (without censoring)
   for (i =1 ; i<= N;i++) {
     if (Status[i] >= 0) {
       for (j=1;j<=NumLinks[i];j++) {
         if  (Age[i] < 25.0) {
           Num25partnerships++;     
           SumLength25partnerships = SumLength25partnerships  + (time - TimeLinked[i][j]);
         }
         if (Age[i] >= 25.0 && Age[i] < 30.0) {
           Num25to30partnerships++;     
           SumLength25to30partnerships = SumLength25to30partnerships  + (time - TimeLinked[i][j]);
         }
         if  (Age[i] >= 30.0) {
           Num30partnerships++;     
           SumLength30partnerships = SumLength30partnerships  + (time - TimeLinked[i][j]);
         }
       }
     }
   }
   AveLength25partnerships = ((double) SumLength25partnerships) / ((double) Num25partnerships); 
   AveLength25to30partnerships = ((double) SumLength25to30partnerships) / ((double) Num25to30partnerships); 
   AveLength30partnerships = ((double) SumLength30partnerships) / ((double) Num30partnerships); 
   fprintf(RelLengthFile,"%ld\t%lf\t%lf\t%lf\n",time,AveLength25partnerships,AveLength25to30partnerships,AveLength30partnerships);

   for (i =1 ; i<= N;i++) {
     if (Status[i] == 1 && Diagnosed[i] == 1) num_diag++;
     if (Status[i] == -2) {
       DiedAIDS++;
       if (Time_AIDS_Death[i] > Start_TasP_Campaign) TotalDiedAIDSAfterTasP++;
     }
     if (Status[i] == -1) DiedNat++;
     if ((Status[i] == -1)|| (Status[i] == -2)) Dead++;
     if (Status[i] == 0) {
        Susceptible++;
        SumAgeSusc = SumAgeSusc + Age[i];
        if (Vaccinated[i] == 1) vaccinated++;
     }
     if (Status[i] >= 0 && Age[i] < 50.0) NumUnder50++;
     if (Status[i] >= 0) {
       Alive = Alive + 1.0;
       SumAge = SumAge + Age[i];
       for (j=1;j<=NumLinks[i];j++) {
          if  (Status[i] != Status[Links[i][j]]) DiscordantCouples++;
       }
     }
     if (Status[i] == 1) {
       Infected++;
       SumAgeInf = SumAgeInf + Age[i];
       G_sum = G_sum + 1.0*Generation[i];
       d_sum = d_sum + d_acute[i];
       set_sum = set_sum + LogSetPoint[i];
       if (Treated[i] == 1) {
         InfectedTreated++;
       }
       if (TargetedTreated[i] == 1) InfectedTreatedTargeted++;
       if (Vaccinated[i] == 1) vaccinated_infected++; 
       if (ResistantToVaccine[i] == 1) resistant_to_vaccine++; 
       if (Age[i] < 50.0) InfUnder50++;
     }
     if (NumUnder50 >= 1) PrevUnder50 = ((double) InfUnder50) / ((double) NumUnder50);
     if ( Status[i] == 1) {
       Vsum = Vsum + log10(V[i]);
     }
     if (time == 1) printf("\n");
     SumLinks = SumLinks + (double) NumLinks[i];
     if (Status[i] >= 0) {
       if (NumLinks[i] == 0) NumUnlinked++;
       num_alive = num_alive + 1.0;
       //printf("First attempt to access NumLinks[%ld]\n",i);
       SumLinks_alive = SumLinks_alive + (double) NumLinks[i];
       if (NumLinks[i] >= 2) SumConc = SumConc + 1.0;
       if (Sex[i] == 1) {
         num_f = num_f + 1.0;
         SumLinks_f = SumLinks_f + (double) NumLinks[i];
         if (NumLinks[i] >= 2) SumConc_f = SumConc_f + 1.0;
         if (Status[i] == 1) {
           inf_f++;
           if (Treated[i] == 1) Inf_Tx_f++;
         }
       }
       if (Sex[i] == 2) {
         num_m = num_m + 1.0;
         SumLinks_m = SumLinks_m + (double) NumLinks[i];
         if (NumLinks[i] >= 2) SumConc_m = SumConc_m + 1.0;
         if (Status[i] == 1) {
           inf_m++;
           if (Treated[i] == 1) Inf_Tx_m++;
         }
       }
       if (Age[i] < 25.0) {
         num25 = num25 + 1.0;
         SumLinks_u25 = SumLinks_u25 + (double) NumLinks[i];
         if (NumLinks[i] >= 2) {
           SumConc_u25 = SumConc_u25 + 1.0;
           if (Sex[i] == 1) SumConc_f25 = SumConc_f25 + 1.0;
           if (Sex[i] == 2) SumConc_m25 = SumConc_m25 + 1.0;
         }
         if (Status[i] == 1) {
           inf25++; 
           if (Treated[i] == 1) Inf_Tx_u25++;
         }
         if (Sex[i] == 1) {
             f25++;
             SumLinks_f25 = SumLinks_f25 + (double) NumLinks[i];
             if ( NumLinks[i] + (NumPartners[i] - NumPartnersLastTime[i]) >= 1 ) { // Change to 1 for stats involving sexually actives only
               sa_f25++;  // Sexually active women under 25
               if ( NumLinks[i] + (NumPartners[i] - NumPartnersLastTime[i]) >= 2 ) {
                 MP_women25++;
               }
             }
         }
         if (Sex[i] == 2) {
             m25++;
             SumLinks_m25 = SumLinks_m25 + (double) NumLinks[i];
             if ( NumLinks[i] + (NumPartners[i] - NumPartnersLastTime[i]) >= 1 ) {
               sa_m25++;
               if ( NumLinks[i] + (NumPartners[i] - NumPartnersLastTime[i]) >= 2 ) {
                 MP_men25++;
               }
             }
         }
         if (Status[i] == 1 && Sex[i] == 1) {
           inf_f25++;
           if (Treated[i] == 1) Inf_Tx_f25++;
         }
         if (Status[i] == 1 && Sex[i] == 2) inf_m25++;
       }
       if (Age[i] >= 25.0 && Age[i] < 50.0) {
         numMid = numMid + 1.0;
         SumLinks_Mid = SumLinks_Mid + (double) NumLinks[i];
         if (NumLinks[i] >= 2) {
           SumConc_Mid = SumConc_Mid + 1.0;
           if (Sex[i] == 1) SumConc_fMid = SumConc_fMid + 1.0;
           if (Sex[i] == 2) SumConc_mMid = SumConc_mMid + 1.0;
         }
         if (Sex[i] == 1) {
            fMid++;
            if ( NumLinks[i] + (NumPartners[i] - NumPartnersLastTime[i]) >= 1 ){
              sa_fMid++;
              if ( NumLinks[i] + (NumPartners[i] - NumPartnersLastTime[i]) >= 2 ) {
                MP_womenMid++;
              }
            }
            SumLinks_fMid = SumLinks_fMid + (double) NumLinks[i];
         }
         if (Sex[i] == 2) {
            mMid++;
            if ( NumLinks[i] + (NumPartners[i] - NumPartnersLastTime[i]) >= 1 ) {
              sa_mMid++;
              if ( NumLinks[i] + (NumPartners[i] - NumPartnersLastTime[i]) >= 2 ) {
                MP_menMid++;
              }
            }
            SumLinks_mMid = SumLinks_mMid + (double) NumLinks[i];
         }
         if (Status[i] == 1) {
            infMid++;
            if (Treated[i] == 1) Inf_Tx_Mid++;
         }
         if (Status[i] == 1 && Sex[i] == 1) {
           inf_fMid++;
           if (Treated[i] == 1) Inf_Tx_fMid++;
         }
         if (Status[i] == 1 && Sex[i] == 2) {
           inf_mMid++;
           if (Treated[i] == 1) Inf_Tx_mMid++;
         }
       }
       if (Age[i] > 50.0) {
         num50 = num50 + 1.0;
         SumLinks50 = SumLinks50 + (double) NumLinks[i];
         if (NumLinks[i] >= 2) {
           SumConc50 = SumConc50 + 1.0;
           if (Sex[i] == 1) SumConc_f50 = SumConc_f50 + 1;
           if (Sex[i] == 2) SumConc_m50 = SumConc_m50 + 1;
         }
         if (Sex[i] == 1) {
           f50++;
           if ( NumLinks[i] + (NumPartners[i] - NumPartnersLastTime[i]) >= 1 ) {
             sa_f50++;
             if ( NumLinks[i] + (NumPartners[i] - NumPartnersLastTime[i]) >= 2 ) {
               MP_women50++;
             }
           }
           SumLinks_f50 = SumLinks_f50 + (double) NumLinks[i];
         }
         if (Sex[i] == 2) {
           m50++;
           if ( NumLinks[i] + (NumPartners[i] - NumPartnersLastTime[i]) >= 1 ) {
             sa_m50++;
             if ( NumLinks[i] + (NumPartners[i] - NumPartnersLastTime[i]) >= 2 ) {
               MP_men50++;
             }
           }
           SumLinks_m50 = SumLinks_m50 + (double) NumLinks[i];
         }
         if (Status[i] == 1) {
           inf50++;
           if (Treated[i] == 1) Inf_Tx50++;
         }
         if (Status[i] == 1 && Sex[i] == 1) {
           inf_f50++;
           if (Treated[i] == 1) Inf_Tx_f50++;
         }
         if (Status[i] == 1 && Sex[i] == 2) {
           inf_m50++;
           if (Treated[i] == 1) Inf_Tx_m50++;
         }
       }
       //printf("Done with first attempt to access NumLinks[%ld]\n",i);
     }
   }
   if (Alive > 0) AveAge = SumAge / Alive; 
             else AveAge = 0.0;
   if (Infected > 0) AveAgeInf = SumAgeInf / Infected;
                else AveAgeInf = 0.0;
   if (Susceptible > 0) AveAgeSusc = SumAgeSusc/ Susceptible;
                   else AveAgeSusc = 0.0;
   for (i=1;i<=N;i++) {
     if (Status[i] >= 0) NumPartnersLastTime[i] = NumPartners[i];
   }

   AveLinks = SumLinks / ((double) (Susceptible+Infected));
   PercentUnlinked = NumUnlinked / num_alive;
   if (fabs(SumLinks- SumLinks_alive) > 0.1) {
     printf("Error SumLinks (%lf) not equal to SumLinks_alive (%lf)\n",SumLinks,SumLinks_alive);
     printf("Type ctrl-c to stop"); scanf("%s",junk);
   }
   if (num_alive > 0) {md = SumLinks_alive / num_alive; } else {md = 0.0;}

   if (num_f > 0) {md_f = SumLinks_f / num_f; } else {md_f = 0.0;}
   if (num_m > 0) {md_m = SumLinks_m / num_m; } else {md_m = 0.0;}
 
   if (num25 > 0) {md_u25 = SumLinks_u25 / num25; } else {md_u25 = 0.0;}
   if (numMid > 0) {md_Mid = SumLinks_Mid / numMid; } else {md_Mid = 0.0;}
   if (num50 > 0) {md50 = SumLinks50 / num50; } else {md50 = 0.0;}
   
   if (f25 > 0) md_f25 = SumLinks_f25 / (double) f25;
           else md_f25 = 0.0;
   
   if (fMid > 0) md_fMid = SumLinks_fMid  /  (double)  fMid;
            else md_fMid = 0.0;
   
   if (f50 > 0) md_f50 = SumLinks_f50  /  (double)  f50;
           else md_f50 = 0.0;

   if (mMid > 0) md_mMid = SumLinks_mMid  / (double)  mMid; 
            else md_mMid = 0.0;

   if (m25 > 0) md_m25 = SumLinks_m25  / (double)  m25;
           else md_m25 = 0.0;
   
   if (m50 > 0) md_m50 = SumLinks_m50 / (double)  m50; 
           else md_m50 = 0.0;
 
   if (sa_f25 > 0) MP_f25 = ( (double) MP_women25 ) / ( (double) sa_f25);
              else MP_f25 = 0.0;
    
   if (sa_fMid > 0) MP_fMid = ( (double) MP_womenMid ) / ( (double) sa_fMid);
               else MP_fMid = 0.0;

   if (sa_f50 > 0) MP_f50 = ( (double) MP_women50 ) / ( (double) sa_f50);
              else MP_f50 = 0.0;
   
   if (sa_m25 > 0) MP_m25 = ( (double) MP_men25 ) / ( (double) sa_m25);
              else MP_m25 = 0.0;

   if (sa_mMid > 0) MP_mMid = ( (double) MP_menMid ) / ( (double) sa_mMid); 
               else MP_mMid = 0.0;
   
   if (sa_m50 > 0) MP_m50 = ( (double) MP_men50 ) / ( (double) sa_m50); 
              else MP_m50 = 0.0;
  
  
   if (num_alive > 0) {concurrency = SumConc / num_alive; } else {md = 0.0;}
   if (num_f > 0) {conc_f = SumConc_f / num_f; } else {conc_f = 0.0;}
   if (num_m > 0) {conc_m = SumConc_m / num_m; } else {conc_m = 0.0;}
 
   if (m25 > 0) {conc_u25 = SumConc_u25 / num25; } else {conc_u25 = 0.0;}
   if (mMid > 0) {conc_Mid = SumConc_Mid / numMid; } else {conc_Mid = 0.0;}
   if (50 > 0) {conc50 = SumConc50 / num50; } else {conc50 = 0.0;}
   
   if (f25 > 0) {conc_f25 = SumConc_f25 / f25; } else {conc_f25 = 0.0;}
   if (fMid > 0) {conc_fMid = SumConc_fMid / fMid; } else {conc_fMid = 0.0;}
   if (f50 > 0) {conc_f50 = SumConc_f50 / f50; } else {conc_f50 = 0.0;}
  
   if (m25 > 0) {conc_m25 = SumConc_m25 / m25; } else {conc_m25 = 0.0;}
   if (mMid > 0) {conc_mMid = SumConc_mMid / mMid; } else {conc_mMid = 0.0;}
   if (m50 > 0) {conc_m50 = SumConc_m50 / m50; } else {conc_m50 = 0.0;}
  
   
   double Sum20=0.0, Sum25=0.0, Sum30=0.0, Sum40 = 0.0, Sum50 = 0.0;
   double PartnerSum20=0.0, PartnerSum25=0.0, PartnerSum30=0.0, PartnerSum40 = 0.0, PartnerSum50 = 0.0;
   double AvePartnersMen20, AvePartnersMen25, AvePartnersMen30, AvePartnersMen40, AvePartnersMen50;
   double AvePartnersWomen20, AvePartnersWomen25, AvePartnersWomen30, AvePartnersWomen40, AvePartnersWomen50;
   
   Sum20=0.0; Sum25=0.0; Sum30=0.0; Sum40 = 0.0; Sum50 = 0.0;
   PartnerSum20=0.0; PartnerSum25=0.0; PartnerSum30=0.0; PartnerSum40 = 0.0; PartnerSum50 = 0.0;
   for (i=1; i<= N; i++) {
     if (Status[i] >= 0 && Sex[i] == 1 && NumPartners[i] >= 1) { // Last term means only collect statistics for those who ever had a partner
        if (Age[i] < 20) {
          Sum20 = Sum20 + 1.0;
          PartnerSum20 = PartnerSum20 + NumPartners[i];
        }
        if (Age[i] >= 20 && Age[i] < 25.0) {
          Sum25 = Sum25 + 1.0;
          PartnerSum25 = PartnerSum25 + NumPartners[i];
        }
        if (Age[i] >= 25 && Age[i] < 30.0) {
          Sum30 = Sum30 + 1.0;
          PartnerSum30 = PartnerSum30 + NumPartners[i];
        }
        if (Age[i] >= 30 && Age[i] < 40.0) {
          Sum40 = Sum40 + 1.0;
          PartnerSum40 = PartnerSum40 + NumPartners[i];
        }
        if (Age[i] >= 40 && Age[i] < 50.0) {
          Sum50 = Sum50 + 1.0;
          PartnerSum50 = PartnerSum50 + NumPartners[i];
        }
      }
    } 
    if (Sum20 > 0) AvePartnersWomen20 = ( (double) PartnerSum20 ) / Sum20; else AvePartnersWomen20 = 0.0;
    if (Sum25 > 0) AvePartnersWomen25 = ( (double) PartnerSum25 ) / Sum25; else AvePartnersWomen25 = 0.0;
    if (Sum30 > 0) AvePartnersWomen30 = ( (double) PartnerSum30 ) / Sum30; else AvePartnersWomen30 = 0.0;
    if (Sum40 > 0) AvePartnersWomen40 = ( (double) PartnerSum40 ) / Sum40; else AvePartnersWomen40 = 0.0;
    if (Sum50 > 0) AvePartnersWomen50 = ( (double) PartnerSum50 ) / Sum50; else AvePartnersWomen50 = 0.0;
   
    Sum20=0.0; Sum25=0.0; Sum30=0.0; Sum40 = 0.0; Sum50 = 0.0;
    PartnerSum20=0.0; PartnerSum25=0.0; PartnerSum30=0.0; PartnerSum40 = 0.0; PartnerSum50 = 0.0;
    for (i=1; i<= N; i++) {
     if (Status[i] >= 0 && Sex[i] == 2 && NumPartners[i] >= 1) { // Last term restricts stats to those who ever had a partner
        if (Age[i] < 20) {
          Sum20 = Sum20 + 1.0;
          PartnerSum20 = PartnerSum20 + NumPartners[i];
        }
        if (Age[i] >= 20 && Age[i] < 25.0) {
          Sum25 = Sum25 + 1.0;
          PartnerSum25 = PartnerSum25 + NumPartners[i];
        }
        if (Age[i] >= 25 && Age[i] < 30.0) {
          Sum30 = Sum30 + 1.0;
          PartnerSum30 = PartnerSum30 + NumPartners[i];
        }
        if (Age[i] >= 30 && Age[i] < 40.0) {
          Sum40 = Sum40 + 1.0;
          PartnerSum40 = PartnerSum40 + NumPartners[i];
        }
        if (Age[i] >= 40 && Age[i] < 50.0) {
          Sum50 = Sum50 + 1.0;
          PartnerSum50 = PartnerSum50 + NumPartners[i];
        }
      }
    } 
    if (Sum20 > 0) AvePartnersMen20 = ( (double) PartnerSum20 ) / Sum20; else AvePartnersMen20 = 0.0;
    if (Sum25 > 0) AvePartnersMen25 = ( (double) PartnerSum25 ) / Sum25; else AvePartnersMen25 = 0.0;
    if (Sum30 > 0) AvePartnersMen30 = ( (double) PartnerSum30 ) / Sum30; else AvePartnersMen30 = 0.0;
    if (Sum40 > 0) AvePartnersMen40 = ( (double) PartnerSum40 ) / Sum40; else AvePartnersMen40 = 0.0;
    if (Sum50 > 0) AvePartnersMen50 = ( (double) PartnerSum50 ) / Sum50; else AvePartnersMen50 = 0.0;

    for (i=1; i<=N; i++) {
     if (Status[i] >= 0) {
       DurSum = DurSum + Duration[i];
       DurCount++;
       if (Age[i] < 18.0) {
         DurSum18 = DurSum18 + Duration[i];
         DurCount18++;
       }
       if (Age[i] < 25.0) {
         DurSum25 = DurSum25 + Duration[i];
         DurCount25++;
       }
       if (Age[i] >= 25.0 && Age[i] < 50.0) {
         DurSumMid = DurSumMid + Duration[i];
         DurCountMid++;
       }
       if (Age[i] >= 50.0) {
         DurSum50 = DurSum50 + Duration[i];
         DurCount50++;
       }
     }
   } 
   if (DurCount > 0) MeanDur = DurSum/DurCount; 
                else MeanDur = 0.0;
   if (DurCount18 > 0) MeanDur18 = DurSum18/DurCount18; 
                  else MeanDur18 = 0.0;
   if (DurCount25 > 0) MeanDur25 = DurSum25/DurCount25; 
                  else MeanDur25 = 0.0;
   if (DurCountMid > 0) MeanDurMid = DurSumMid/DurCountMid; 
                  else MeanDurMid = 0.0;
   if (DurCount50 > 0) MeanDur50 = DurSum50/DurCount50; 
                  else MeanDur50 = 0.0;

   if (Infected >= 1) {
     Vave = Vsum/(1.0*Infected);
     d_ave = d_sum/(1.0*Infected);
     G_ave = G_sum/(1.0*Infected);
     set_ave = set_sum/(1.0*Infected);
   } else {
     Vave = 0.0;
     d_ave = 0.0;
     G_ave = 0.0;
     set_ave = 0.0;
   }
   Vstdsum = 0.0;  d_stdsum = 0.0; G_stdsum = 0.0; set_stdsum = 0.0; vcount = 0.0;
   for (i=1;i<=N;i++) {
     if (V[i] > epsilon) {
       vcount = vcount + 1.0;
        d_stdsum = d_stdsum + (d_acute[i] - d_ave)*(d_acute[i]-d_ave);
       set_stdsum = set_stdsum + (LogSetPoint[i] - set_ave)*(LogSetPoint[i]-set_ave);
       G_stdsum = G_stdsum + ((1.0*Generation[i]) - G_ave)* ((1.0*Generation[i]) - G_ave);
     }
     if ((Status[i] == 1) && ( time - Time_Inf[i] > 50.0)) {
       Vstdsum = Vstdsum + (log10(V[i]) - Vave)*(log10(V[i])-Vave);
     }

   }
   if (Infected >= 1) Vstd = sqrt(Vstdsum/Infected);
   else Vstd = 0.0;
   d_std = sqrt(d_stdsum/vcount); G_std = sqrt(G_stdsum/vcount); set_std = sqrt(set_stdsum/vcount);
   printcount++;
   if (printcount >= 30 && printOutput==1) {
    //printf("\nRepl\tTime\tN\tLinks\tUnInf\tInf\tTx\tNoTx\tDead\tdAIDS\tdNat\tLogV\tVstd\tSPVL\t\tSet_std\t\td_ave\t\td_std\tMD\n");
     printcount = 0;
   }
   if(printOutput==1)
   {
     printf("%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t",repl,time,N,(long)(SumLinks+0.5),Susceptible,Infected,num_diag,InfectedTreated,Infected-InfectedTreated,Dead,DiedAIDS,DiedNat);
     printf("%3.2lf\t%3.2lf\t%3.2lf\t%3.2lf\t%3.2lf\t%ld\t%ld\n",AveAge,Vave,set_ave,AveLinks,PercentUnlinked,new_infections,DiscordantCouples);
       //printf("Repl\tTime\tN\tLinks\tUnInf\tInf\tDiag\tTx\tNoTx\tDead\tdAIDS\tdNat\tLogV\tVstd\tSPVL\t\tSet_std\t\td_ave\t\td_std\t\tMD\tUnLink\tNewInfs\tDiscCouples\n");
   }
   // if (plt_out==1) fprintf(OutputFile,"\n\ntime\tN\tSusc\tInfected\tVaccinated\tVaccinatedInfected\tResistVaccine\tInfectedTreated\tDead\tdAIDS\tdNat\tLogV\tVstd\tSPVL\t\tSet_std\t\td_ave\t\td_std\t\tLinks\t\tV1\t\tda1\t\tSts1\n");
   //   long vaccinated, vaccinated_infected, resistant_to_vaccine;

   if (plt_out==1) fprintf(OutputFile,"%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t",repl,time,N,Susceptible,Infected,vaccinated, vaccinated_infected, resistant_to_vaccine,
                                                                 InfectedTreated,Dead,DiedAIDS,DiedNat);
   if (plt_out==1) fprintf(OutputFile,"%3.2lf\t%3.2lf\t%3.2e\t%3.2e\t%3.2e\t%3.2e\t",Vave,Vstd,set_ave,set_std,d_ave,d_std);
   if (plt_out==1) fprintf(OutputFile,"%3.2lf\n",AveLinks);
   //long Inf_Tx = 0, Inf_Tx_m = 0, Inf_Tx_f = 0, Inf_Tx_u25 = 0, Inf_Tx_Mid = 0, Inf_Tx50 = 0; // HIV+ persons of various types that receive treatment
   //long Inf_Tx_m25 = 0, Inf_Tx_mMid = 0, Inf_Tx_m50 = 0; // Number of HIV+ persons receiving treatment
   //long Inf_Tx_f25 = 0, Inf_Tx_fMid = 0, Inf_Tx_f50 =0; // Number of HIV+ persons receiving treatment
   if (time < 1) {
       NetStatsCount = 1;
   }
   if (plt_NS==1) fprintf(NetworkStatsOutputFile,"%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%3.2lf\t%3.2lf\t%3.2lf\t%3.2lf\t%3.2lf\t%3.2lf\t%3.2lf\t%3.2lf\t%3.2lf\t%3.2lf\t%3.2lf\t%3.2lf\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t",
                                   NetStatsCount,repl,time,Susceptible, Infected, (long) num_alive, (long) num_f, (long) num_m, (long) num25, (long) numMid, (long) num50, md, md_u25, md_Mid, md50, md_f, md_m,
                                   concurrency, conc_f, conc_m, conc_u25, conc_Mid, conc50,InfectedTreated,InfectedTreatedTargeted,DiedAIDS,DiedNat,susc_days,new_infections,susc_days_u50,new_infections_u50);
   if (plt_NS==1) fprintf(NetworkStatsOutputFile,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t",f25,m25,fMid,mMid,f50,m50,inf_f,inf_m,inf25,infMid,inf50,inf_f25,inf_m25,inf_fMid,inf_mMid,inf_f50,inf_m50);
 
   if (plt_NS==1) fprintf(NetworkStatsOutputFile,"%ld\t%ld\t%ld\t%ld\t%ld\t%lf\t%lf\t%lf\t%lf\t%lf\t%ld\t%ld\t",Inf_Tx_u25,Inf_Tx_Mid,Inf_Tx50,Inf_Tx_f,Inf_Tx_m,MeanDur,MeanDur18,MeanDur25,MeanDurMid,MeanDur50,TotalTreated,NonTargetedTx);
   if (plt_NS==1) fprintf(NetworkStatsOutputFile,"%lf\t%lf\t%lf\t%lf\t%lf\t",AvePartnersMen20,AvePartnersMen25,AvePartnersMen30,AvePartnersMen40,AvePartnersMen50);
   if (plt_NS==1) fprintf(NetworkStatsOutputFile,"%lf\t%lf\t%lf\t%lf\t%lf\t",AvePartnersWomen20,AvePartnersWomen25,AvePartnersWomen30,AvePartnersWomen40,AvePartnersWomen50);
   if (plt_NS==1) fprintf(NetworkStatsOutputFile,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t",conc_f25,conc_m25,conc_fMid,conc_mMid,conc_f50,conc_m50);
   if (plt_NS==1) fprintf(NetworkStatsOutputFile,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t",md_f25,md_m25,md_fMid,md_mMid,md_f50,md_m50,PrevUnder50);
   if (plt_NS==1) fprintf(NetworkStatsOutputFile,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t",MP_m25,MP_mMid,MP_m50,MP_f25,MP_fMid,MP_f50);
   if (plt_NS==1) fprintf(NetworkStatsOutputFile,"%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%lf\t%lf\t%lf\n",susc_days_u25,susc_days_Mid,susc_days_o50,new_infections_u25,new_infections_Mid,new_infections_o50,AveAge,AveAgeSusc,AveAgeInf);
   NetStatsCount++;
   fflush(stdout);
}

double nrand(void)
{
  return(sqrt(-2.0*log(genrand_real2()))*cos(2.0*3.14159265358979*genrand_real2()));
}

/*
   A C-program for MT19937, with initialization improved 2002/1/26.
   Coded by Takuji Nishimura and Makoto Matsumoto.

   Before using, initialize the state by using init_genrand(seed)
   or init_by_array(init_key, key_length).

   Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
   All rights reserved.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

     1. Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.

     2. Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.

     3. The names of its contributors may not be used to endorse or promote
        products derived from this software without specific prior written
        permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


   Any feedback is very welcome.
   http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
   email: m-mat @ math.sci.hiroshima-u.ac.jp (remove space)
*/


/* Period parameters */
#define N 624
#define M397 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */

static unsigned long mt[N]; /* the array for the state vector  */
static int mti=N+1; /* mti==N+1 means mt[N] is not initialized */

/* initializes mt[N] with a seed */
void init_genrand(unsigned long SSS)
{
    mt[0]= SSS & 0xffffffffUL;
    for (mti=1; mti<N; mti++) {
        mt[mti] =
	    (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti);
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array mt[].                        */
        /* 2002/01/09 modified by Makoto Matsumoto             */
        mt[mti] &= 0xffffffffUL;
        /* for >32 bit machines */
    }
}

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
/* slight change for C++, 2004/2/26 */
void init_by_array(unsigned long init_key[], int key_length)
{
    int i, j, k;
    init_genrand(19650218UL);
    i=1; j=0;
    k = (N>key_length ? N : key_length);
    for (; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1664525UL))
          + init_key[j] + j; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++; j++;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
        if (j>=key_length) j=0;
    }
    for (k=N-1; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1566083941UL))
          - i; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
    }

    mt[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */
}

/* generates a random number on [0,0xffffffff]-interval */
unsigned long genrand_int32(void)
{
    unsigned long y;
    static unsigned long mag01[2]={0x0UL, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if (mti >= N) { /* generate N words at one time */
        int kk;

        if (mti == N+1)   /* if init_genrand() has not been called, */
            init_genrand(5489UL); /* a default initial seed is used */

        for (kk=0;kk<N-M397;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+M397] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (;kk<N-1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(M397-N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[N-1] = mt[M397-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

        mti = 0;
    }

    y = mt[mti++];

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
}

/* generates a random number on [0,0x7fffffff]-interval */
long genrand_int31(void)
{
    return (long)(genrand_int32()>>1);
}

/* generates a random number on [0,1]-real-interval */
double genrand_real1(void)
{
    return genrand_int32()*(1.0/4294967295.0);
    /* divided by 2^32-1 */
}

/* generates a random number on [0,1)-real-interval */
double genrand_real2(void)
{
    return genrand_int32()*(1.0/4294967296.0);
    /* divided by 2^32 */
}

/* generates a random number on (0,1)-real-interval */
double genrand_real3(void)
{
    return (((double)genrand_int32()) + 0.5)*(1.0/4294967296.0);
    /* divided by 2^32 */
}

/* generates a random number on [0,1) with 53-bit resolution*/
double genrand_res53(void)
{
    unsigned long a=genrand_int32()>>5, b=genrand_int32()>>6;
    return(a*67108864.0+b)*(1.0/9007199254740992.0);
}
/* These real versions are due to Isaku Wada, 2002/01/09 added */
