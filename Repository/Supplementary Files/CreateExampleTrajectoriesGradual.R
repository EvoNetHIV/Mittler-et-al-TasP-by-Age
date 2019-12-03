#note, graphing routine starts line 163, need to change pathways below

re_read_data <- FALSE
par(mfrow=c(3,3))
legend_text_size <- 0.65

if (re_read_data == TRUE) {
  run1_title <- "'Random'"
  load("H:/TasP2/Exp20_Gradual_95_care_N10000/txt_opto_Exp20_paramsGradual_95_N10000_TasP60perc__random_prop0.2.RData")
  run1_1=evomodel
  remove(evomodel)
  load("H:/TasP2/Exp20_Gradual_95_care_N10000/txt_opto_Exp20_paramsGradual_95_N10000_TasP60perc__random_prop0.3.RData")
  run1_2=evomodel
  remove(evomodel)
  load("H:/TasP2/Exp20_Gradual_95_care_N10000/txt_opto_Exp20_paramsGradual_95_N10000_TasP60perc__random_prop0.4.RData")
  run1_3=evomodel
  remove(evomodel)
  run2_title <- "'Under Age 30'"
  load("H:/TasP2/Exp20_Gradual_95_care_N10000/txt_opto_Exp20_paramsGradual_95_N10000_TasP60perc__under25_under30_random_prop0.2.RData")
  run2_1=evomodel
  remove(evomodel)
  load("H:/TasP2/Exp20_Gradual_95_care_N10000/txt_opto_Exp20_paramsGradual_95_N10000_TasP60perc__under25_under30_random_prop0.3.RData")
  run2_2=evomodel
  remove(evomodel)
  load("H:/TasP2/Exp20_Gradual_95_care_N10000/txt_opto_Exp20_paramsGradual_95_N10000_TasP60perc__under25_under30_random_prop0.4.RData")
  run2_3=evomodel
  remove(evomodel)
  run1_1_ps <- run1_1$popsumm
  run1_2_ps <- run1_2$popsumm
  run1_3_ps <- run1_3$popsumm
  run2_1_ps <- run2_1$popsumm
  run2_2_ps <- run2_2$popsumm
  run2_3_ps <- run2_3$popsumm
}
nsteps <- run1_1$param[[1]]$n_steps
nyears <- (nsteps/365)


max_infected <- 1400
max_incidence <- 2.0
max_death_rate <- 0.18

##############################

extract_evonet_data <-function(model_obj,overlay){
  #plots incidence rate
  #input netsim object
  #output: plot of incidence rate by year
  if (overlay==TRUE) graph_col = "red"
               else  graph_col = "black"
  
  model <- model_obj
  #model <- model_obj 
  total_sims <- model$param[[1]]$nsims
  result_list <- vector('list',length=total_sims)
  
  for(nsim in 1:total_sims)
  {
    popsumm_freq <- model$param[[nsim]]$popsumm_frequency
    nsteps <- model$param[[nsim]]$n_steps
    
    if(nsteps>=365)
      steps_per_year <- floor(365/model$param[[nsim]]$popsumm_frequency)
    else
      steps_per_year <- floor(nsteps/model$param[[nsim]]$popsumm_frequency)
    
    sus <- model$popsumm[[nsim]]$susceptibles
    inf <- model$popsumm[[nsim]]$new_infections
    nsteps <- model$param[[nsim]]$n_steps
    nyears <- (nsteps/365)
    
    if(nyears<1)
      tli <- floor(seq(1,length(sus),length=2))
    if(nyears>=1)
      tli <- floor(seq(1,length(sus),length=floor(nyears)+1))   
    
    number_sus <- rep(NA_real_,length(tli)-1)
    total_new_inf <- rep(NA_real_,length(tli)-1)
    
    #browser()
    for(ii in 1:(length(tli)-1))
    {
      total_new_inf[ii] <- sum(inf[(tli[ii]+1):(tli[ii+1])])
      
      if(popsumm_freq>1){
        number_sus[ii] <- mean(sus[(tli[ii]+1):(tli[ii+1])])
      }
      #note: this gives exact same results as above
      if(popsumm_freq==1){
        number_sus[ii] <- sum(sus[(tli[ii]+1):(tli[ii+1])])/365
      }
    }
    
    if(sum(number_sus)==0){next}
    
    # scalar accounts for round-off error
    scalar<-365/(diff(tli)*popsumm_freq)
    inc <- scalar*100*total_new_inf/number_sus
    inc[number_sus==0]<-0
    result_list[[nsim]] <- inc
    
  }
  result_mat <- do.call(rbind,result_list)
  ymax=max(unlist(result_mat),na.rm=T)
  for(jj in 1:total_sims) {
     if(jj==1 && overlay==FALSE){
       if(nsteps>=365){
         plot((1:nyears)-20,result_list[[jj]],type='l',xlim=c(-10,25),ylim=c(0,1.6),xlab="Years before/after TasP campaingn",ylab="Incidence [rate per 100 person years]",
             #xlim=c(10,nyears),
             lty=2,col=graph_col,axes=F)
         labseq=seq(0,nyears,by=10)-20
         #labseq=labseq[-1]
         axis(1,at=labseq,labels=labseq)
       }
      
       if(nsteps<365){
          plot(nyears,result_list[[jj]],type='l',xlim=c(-10,25),ylim=c(0,ymax),xlab="Years before/after TasP campaingn",ylab="Incidence [rate per 100 person years]",
             lty=2,col=graph_col,axes=F)
         axis(1,at=1:nyears,labels=1:nyears)
       }
       axis(2);box()
      } else {
       lines((1:nyears)-20,result_list[[jj]],type='l',col=graph_col,ylim=c(0,max(inc)),lty=2)
      }
  }
  lines((1:nyears)-20,colMeans(result_mat),type='l',lwd=2,col=graph_col)
  
 
 
  # Now calculate mean and standard deviations from the last 10 years
  
  mean_incid <- 0 ; mean_prev <- 0; mean_pills_taken <- 0; mean_not_prioritized_yr10 <- 0; mean_died_AIDS <- 0
  for(jj in 1:total_sims)
  {
    last_recorded_step <- length(model$popsumm[[jj]]$prevalence)
    after_rampup <- 20
    
    mean_incid <- mean_incid + mean(result_list[[jj]][40:45])
    mean_prev <- mean_prev + model$popsumm[[jj]]$prevalence[last_recorded_step]
    mean_pills_taken <- mean_pills_taken + sum(model$popsumm[[jj]]$total_pills_taken[last_recorded_step])
    mean_not_prioritized_yr10 <- mean_not_prioritized_yr10 + model$param[[jj]]$num_randomly_chosen_start_campaign
    mean_died_AIDS <- mean_died_AIDS + sum(model$popsumm[[jj]]$aids_deaths[after_rampup:last_recorded_step])
  }
  mean_incid <- mean_incid/total_sims
  mean_prev  <- mean_prev/total_sims
  mean_pills_taken <- mean_pills_taken/total_sims
  mean_not_prioritized_yr10 <- mean_not_prioritized_yr10/total_sims
  mean_died_AIDS <- mean_died_AIDS/total_sims
  
  ss_incid <- 0 ; ss_prev <- 0; ss_pills_taken <- 0; ss_not_prioritized_yr10 <- 0; ss_died_AIDS <- 0
  for(jj in 1:total_sims)
  {
    ss_incid <- ss_incid + (mean(result_list[[jj]][40:45]) - mean_incid)^2
    ss_prev <- ss_prev + (model$popsumm[[jj]]$prevalence[last_recorded_step] - mean_prev)^2
    ss_pills_taken <- ss_pills_taken + (model$popsumm[[jj]]$total_pills_taken[last_recorded_step] - mean_pills_taken)^2
    ss_not_prioritized_yr10 <- ss_not_prioritized_yr10 + (model$param[[jj]]$num_randomly_chosen_start_campaign - mean_not_prioritized_yr10)^2
    ss_died_AIDS <- ss_died_AIDS + (sum(model$popsumm[[jj]]$aids_deaths[after_rampup:last_recorded_step]) - mean_died_AIDS)^2
  }
  sd_incid <- sqrt(ss_incid/(total_sims-1))
  sd_prev <- sqrt(ss_prev/(total_sims-1))
  sd_pills_taken <- sqrt(ss_pills_taken/(total_sims-1))
  sd_not_prioritized_yr10 <- sqrt(ss_not_prioritized_yr10/(total_sims-1))
  sd_died_AIDS <- sqrt(ss_died_AIDS/(total_sims-1))
  cat("Incidence last 10 years:",mean_incid,"(",sd_incid,")\n")
  cat("Prevalence final time point:",mean_prev,"(",sd_prev,")\n")
  cat("Total pills taken:",mean_pills_taken,"(",sd_pills_taken,")\n")
  cat("Not prioritized year 10:",mean_not_prioritized_yr10,"(",sd_not_prioritized_yr10,")\n")
  cat("Died of AIDS after ramp-up: ",mean_died_AIDS,"(",sd_died_AIDS,")\n")
  yearly_inc <- 1 + model$param[[1]]$yearly_incr_tx 
  evonet_results <- list("proportion_treated_yr10" = model$param[[1]]$proportion_treated/( (yearly_inc^(1/365))^(365*30) ),
                         "proportion_treated_end" = model$param[[1]]$proportion_treated,
                         "under_care" = model$param[[1]]$prob_care,
                         "incid" = mean_incid, "incid_low"  = mean_incid  - sd_incid,
                         "incid_high" = mean_incid  + sd_incid,
                         "prev" = mean_prev,  "prev_low"  = mean_prev - sd_prev,
                         "prev_high" = mean_prev + sd_prev,
                         "pills_taken" = mean_pills_taken, "pills_taken_low"  = mean_pills_taken - sd_pills_taken,
                         "pills_taken_high" = mean_pills_taken + sd_pills_taken,
                         "not_prioritized_yr10" = mean_not_prioritized_yr10, 
                         "not_prioritized_yr10_low"  = mean_not_prioritized_yr10 - sd_not_prioritized_yr10,
                         "not_prioritized_yr10_high" = mean_not_prioritized_yr10 + sd_not_prioritized_yr10,
                         "died_AIDS" = mean_died_AIDS, "died_AIDS_low"  = mean_died_AIDS - sd_died_AIDS,
                         "died_AIDS_high" = mean_died_AIDS + sd_died_AIDS
  )
  
  
  
  return(evonet_results)
}
#####################################

#par(mfrow=c(3,3))
nreps = 16

###

###
#Plot number infected and number treated for strategy # 1_1
labseq=seq(0,nyears,by=10)-20
plot(-20+run1_1_ps[[1]]$timestep/365,run1_1_ps[[1]]$total_infections_alive,xlim=c(-10,25),ylim=c(0,1700),
     #xlim=c(10,45),
     type="l",lwd=0.5,lty=2,xlab="Years before/after TasP campaign",ylab="Number Infected and Treated",col="black",axes=F)
axis(1,at=labseq,labels=labseq)
lines(-20+run1_1_ps[[1]]$timestep/365,run1_1_ps[[1]]$no_treated  ,type="l",lwd=0.5,lty=2,col="blue")
mean_no_treated_run1_1_ps <- run1_1_ps[[1]]$no_treated/nreps
mean_infected_run1_1_ps <- run1_1_ps[[1]]$total_infections_alive/nreps
for (jj in 2:nreps) {
  lines(-20+run1_1_ps[[jj]]$timestep/365,run1_1_ps[[jj]]$total_infections_alive,type="l",lwd=0.5,lty=2,col="black")
  lines(-20+run1_1_ps[[jj]]$timestep/365,run1_1_ps[[jj]]$no_treated  ,type="l",lwd=0.5,lty=2,col="blue")
  #lines(run1_1_ps[[jj]]$timestep/365,run1_1_ps[[jj]]$no_prioritized_treated  ,type="l",lwd=0.5,lty=2,col="orange")
  mean_no_treated_run1_1_ps <- mean_no_treated_run1_1_ps + run1_1_ps[[jj]]$no_treated/nreps
  mean_infected_run1_1_ps <- mean_infected_run1_1_ps + run1_1_ps[[jj]]$total_infections_alive/nreps
}
lines(-20+run1_1_ps[[jj]]$timestep/365,mean_infected_run1_1_ps,type="l",lwd=2,col="black")
lines(-20+run1_1_ps[[jj]]$timestep/365,mean_no_treated_run1_1_ps,type="l",lwd=2,col="blue")
axis(2); box()
#title("Slow, linear")

#Plot number infected and number treated for strategy # 2_1

lines(-20+run2_1_ps[[1]]$timestep/365,run2_1_ps[[1]]$total_infections_alive,type="l",lwd=0.5,col="red")
lines(-20+run2_1_ps[[1]]$timestep/365,run2_1_ps[[1]]$no_treated  ,type="l",lwd=0.5,lty=2,col="green")
mean_no_treated_run2_1_ps <- run2_1_ps[[1]]$no_treated/nreps
mean_infected_run2_1_ps <- run2_1_ps[[1]]$total_infections_alive/nreps
mean_prioritizied_run2_1_ps <- run2_1_ps[[1]]$prioritized_tx/nreps
for (jj in 2:nreps) {
  lines(-20+run2_1_ps[[jj]]$timestep/365,run2_1_ps[[jj]]$total_infections_alive,type="l",lwd=0.5,lty=2,col="red")
  lines(-20+run2_1_ps[[jj]]$timestep/365,run2_1_ps[[jj]]$no_treated  ,type="l",lwd=0.5,lty=2,col="green")
  lines(-20+run2_1_ps[[jj]]$timestep/365,run2_1_ps[[jj]]$prioritized_tx  ,type="l",lwd=0.5,lty=2,col="purple")
  #lines(run2_1_ps[[jj]]$timestep/365,run2_1_ps[[jj]]$no_prioritized_treated  ,type="l",lwd=0.5,lty=2,col="orange")
  mean_no_treated_run2_1_ps <- mean_no_treated_run2_1_ps + run2_1_ps[[jj]]$no_treated/nreps
  mean_infected_run2_1_ps <- mean_infected_run2_1_ps + run2_1_ps[[jj]]$total_infections_alive/nreps
  mean_prioritizied_run2_1_ps <- mean_prioritizied_run2_1_ps + run2_1_ps[[jj]]$prioritized_tx/nreps
}
lines(-20+run2_1_ps[[jj]]$timestep/365,mean_infected_run2_1_ps,type="l",lwd=2,col="red")
lines(-20+run2_1_ps[[jj]]$timestep/365,mean_no_treated_run2_1_ps,type="l",lwd=2,col="green")
lines(-20+run2_1_ps[[jj]]$timestep/365,mean_prioritizied_run2_1_ps,type="l",lwd=2,col="purple")
#title(run2_1_title)
if (1==1)legend(-10,1770,c(paste("Infected ",run1_title,sep=""),
                                        paste("Infected ",run2_title,sep=""),
                                        paste("Treated ",run1_title,sep=""),
                                        paste("Treated ",run2_title,sep=""),
                                        paste("Targeted ",run2_title,sep=""))
                ,lwd=c(2,2,2,2,2),col=c("black","red","blue","green","purple"),bty='n',cex=legend_text_size,
                y.intersp = 1)

#Plot number infected and number treated for strategy # 1_2
plot(-20+run1_2_ps[[1]]$timestep/365,run1_2_ps[[1]]$total_infections_alive,xlim=c(-10,25),ylim=c(0,1700),
     #xlim=c(10,45),
     type="l",lwd=0.5,lty=2,xlab="Years before/after TasP campaign",ylab="Number Infected and Treated",col="black",axes="F")
axis(1,at=labseq,labels=labseq)
lines(-20+run1_2_ps[[1]]$timestep/365,run1_2_ps[[1]]$no_treated  ,type="l",lwd=0.5,lty=2,col="blue")
mean_no_treated_run1_2_ps <- run1_2_ps[[1]]$no_treated/nreps
mean_infected_run1_2_ps <- run1_2_ps[[1]]$total_infections_alive/nreps
for (jj in 2:nreps) {
  lines(-20+run1_2_ps[[jj]]$timestep/365,run1_2_ps[[jj]]$total_infections_alive,type="l",lwd=0.5,lty=2,col="black")
  lines(-20+run1_2_ps[[jj]]$timestep/365,run1_2_ps[[jj]]$no_treated  ,type="l",lwd=0.5,lty=2,col="blue")
  #lines(run1_2_ps[[jj]]$timestep/365,run1_2_ps[[jj]]$no_prioritized_treated  ,type="l",lwd=0.5,lty=2,col="orange")
  mean_no_treated_run1_2_ps <- mean_no_treated_run1_2_ps + run1_2_ps[[jj]]$no_treated/nreps
  mean_infected_run1_2_ps <- mean_infected_run1_2_ps + run1_2_ps[[jj]]$total_infections_alive/nreps
}
lines(-20+run1_2_ps[[jj]]$timestep/365,mean_infected_run1_2_ps,type="l",lwd=2,col="black")
lines(-20+run1_2_ps[[jj]]$timestep/365,mean_no_treated_run1_2_ps,type="l",lwd=2,col="blue")
#title("Slow, linear")
axis(2); box()
#Plot number infected and number treated for strategy # 2_2

lines(-20+run2_2_ps[[1]]$timestep/365,run2_2_ps[[1]]$total_infections_alive,type="l",lwd=0.5,col="red")
lines(-20+run2_2_ps[[1]]$timestep/365,run2_2_ps[[1]]$no_treated  ,type="l",lwd=0.5,lty=2,col="green")
mean_no_treated_run2_2_ps <- run2_2_ps[[1]]$no_treated/nreps
mean_infected_run2_2_ps <- run2_2_ps[[1]]$total_infections_alive/nreps
mean_prioritizied_run2_2_ps <- run2_2_ps[[1]]$prioritized_tx/nreps
for (jj in 2:nreps) {
  lines(-20+run2_2_ps[[jj]]$timestep/365,run2_2_ps[[jj]]$total_infections_alive,type="l",lwd=0.5,lty=2,col="red")
  lines(-20+run2_2_ps[[jj]]$timestep/365,run2_2_ps[[jj]]$no_treated  ,type="l",lwd=0.5,lty=2,col="green")
  lines(-20+run2_2_ps[[jj]]$timestep/365,run2_2_ps[[jj]]$prioritized_tx  ,type="l",lwd=0.5,lty=2,col="purple")
  #lines(run2_2_ps[[jj]]$timestep/365,run2_2_ps[[jj]]$no_prioritized_treated  ,type="l",lwd=0.5,lty=2,col="orange")
  mean_no_treated_run2_2_ps <- mean_no_treated_run2_2_ps + run2_2_ps[[jj]]$no_treated/nreps
  mean_infected_run2_2_ps <- mean_infected_run2_2_ps + run2_2_ps[[jj]]$total_infections_alive/nreps
  mean_prioritizied_run2_2_ps <- mean_prioritizied_run2_2_ps + run2_2_ps[[jj]]$prioritized_tx/nreps
}
lines(-20+run2_2_ps[[jj]]$timestep/365,mean_infected_run2_2_ps,type="l",lwd=2,col="red")
lines(-20+run2_2_ps[[jj]]$timestep/365,mean_no_treated_run2_2_ps,type="l",lwd=2,col="green")
lines(-20+run2_2_ps[[jj]]$timestep/365,mean_prioritizied_run2_2_ps,type="l",lwd=2,col="purple")
#title(run2_2_title)
if (1==2)legend(-10,1.03*max_infected,c(paste("Infected ",run1_title,sep=""),
                                        paste("Infected ",run2_title,sep=""),
                                        paste("Treated ",run1_title,sep=""),
                                        paste("Treated ",run2_title,sep=""),
                                        paste("Targeted ",run2_title,sep=""))
                ,lwd=c(2,2,2,2,2),col=c("black","red","blue","green","purple"),bty='n',cex=legend_text_size,
                y.intersp = 1)


#Plot number infected and number treated for strategy # 1_3
plot(-20+run1_3_ps[[1]]$timestep/365,run1_3_ps[[1]]$total_infections_alive,xlim=c(-10,25),ylim=c(0,1700),
     #xlim=c(10,45),axes="F",
     type="l",lwd=0.5,lty=2,xlab="Years before/after TasP campaign",ylab="Number Infected and Treated",col="black")
axis(1,at=labseq,labels=labseq)
lines(-20+run1_3_ps[[1]]$timestep/365,run1_3_ps[[1]]$no_treated  ,type="l",lwd=0.5,lty=2,col="blue")
mean_no_treated_run1_3_ps <- run1_3_ps[[1]]$no_treated/nreps
mean_infected_run1_3_ps <- run1_3_ps[[1]]$total_infections_alive/nreps
for (jj in 2:nreps) {
  lines(-20+run1_3_ps[[jj]]$timestep/365,run1_3_ps[[jj]]$total_infections_alive,type="l",lwd=0.5,lty=2,col="black")
  lines(-20+run1_3_ps[[jj]]$timestep/365,run1_3_ps[[jj]]$no_treated  ,type="l",lwd=0.5,lty=2,col="blue")
  #lines(run1_3_ps[[jj]]$timestep/365,run1_3_ps[[jj]]$no_prioritized_treated  ,type="l",lwd=0.5,lty=2,col="orange")
  mean_no_treated_run1_3_ps <- mean_no_treated_run1_3_ps + run1_3_ps[[jj]]$no_treated/nreps
  mean_infected_run1_3_ps <- mean_infected_run1_3_ps + run1_3_ps[[jj]]$total_infections_alive/nreps
}
lines(-20+run1_3_ps[[jj]]$timestep/365,mean_infected_run1_3_ps,type="l",lwd=2,col="black")
lines(-20+run1_3_ps[[jj]]$timestep/365,mean_no_treated_run1_3_ps,type="l",lwd=2,col="blue")
axis(2); box()#title("Slow, linear")

#Plot number infected and number treated for strategy # 2_3

lines(-20+run2_3_ps[[1]]$timestep/365,run2_3_ps[[1]]$total_infections_alive,type="l",lwd=0.5,col="red")
lines(-20+run2_3_ps[[1]]$timestep/365,run2_3_ps[[1]]$no_treated  ,type="l",lwd=0.5,lty=2,col="green")
mean_no_treated_run2_3_ps <- run2_3_ps[[1]]$no_treated/nreps
mean_infected_run2_3_ps <- run2_3_ps[[1]]$total_infections_alive/nreps
mean_prioritizied_run2_3_ps <- run2_3_ps[[1]]$prioritized_tx/nreps
for (jj in 2:nreps) {
  lines(-20+run2_3_ps[[jj]]$timestep/365,run2_3_ps[[jj]]$total_infections_alive,type="l",lwd=0.5,lty=2,col="red")
  lines(-20+run2_3_ps[[jj]]$timestep/365,run2_3_ps[[jj]]$no_treated  ,type="l",lwd=0.5,lty=2,col="green")
  lines(-20+run2_3_ps[[jj]]$timestep/365,run2_3_ps[[jj]]$prioritized_tx  ,type="l",lwd=0.5,lty=2,col="purple")
  #lines(run2_3_ps[[jj]]$timestep/365,run2_3_ps[[jj]]$no_prioritized_treated  ,type="l",lwd=0.5,lty=2,col="orange")
  mean_no_treated_run2_3_ps <- mean_no_treated_run2_3_ps + run2_3_ps[[jj]]$no_treated/nreps
  mean_infected_run2_3_ps <- mean_infected_run2_3_ps + run2_3_ps[[jj]]$total_infections_alive/nreps
  mean_prioritizied_run2_3_ps <- mean_prioritizied_run2_3_ps + run2_3_ps[[jj]]$prioritized_tx/nreps
}
lines(-20+run2_3_ps[[jj]]$timestep/365,mean_infected_run2_3_ps,type="l",lwd=2,col="red")
lines(-20+run2_3_ps[[jj]]$timestep/365,mean_no_treated_run2_3_ps,type="l",lwd=2,col="green")
lines(-20+run2_3_ps[[jj]]$timestep/365,mean_prioritized_run2_3_ps,type="l",lwd=2,col="purple")
#title(run2_3_title)
if (1==2)legend(-10,1.03*max_infected,c(paste("Infected ",run1_title,sep=""),
                               paste("Infected ",run2_title,sep=""),
                               paste("Treated ",run1_title,sep=""),
                               paste("Treated ",run2_title,sep=""),
                               paste("Targeted ",run2_title,sep=""))
       ,lwd=c(2,2,2,2,2),col=c("black","red","blue","green","purple"),bty='n',cex=legend_text_size,
       y.intersp = 1)


# Plot incidence for under strategy 1_1
extract_evonet_data(run1_1,overlay=FALSE)
#title(run1_title)

# Plot incidence for strategy 2_1
extract_evonet_data(run2_1,overlay=TRUE)
legend(-10,1.68,c(run1_title, run2_title), cex=legend_text_size,
      lwd=c(2,2),col=c("black","red"),bty='n',y.intersp = 1)

# Plot incidence for under strategy 1_2
extract_evonet_data(run1_2,overlay=FALSE)
#title(run1_title)

# Plot incidence for strategy 2_2
extract_evonet_data(run2_2,overlay=TRUE)
#legend(10,1.5,c(run1_title, run2_title), cex=legend_text_size,
#       lwd=c(2,2),col=c("black","red"),bty='n',y.intersp = 1)

# Plot incidence for under strategy 1_3
extract_evonet_data(run1_3,overlay=FALSE)
#title(run1_title)

# Plot incidence for strategy 2_3
extract_evonet_data(run2_3,overlay=TRUE)
#legend(10,1.5,c(run1_title, run2_title), cex=legend_text_size,
#       lwd=c(2,2),col=c("black","red"),bty='n',y.intersp = 1)

#Plot AIDS Death Rates for  Strategy 1_1
plot(-20+run1_1_ps[[1]]$timestep/365,run1_1_ps[[1]]$aids_deaths/run1_1_ps[[1]]$total_infections_alive,
     xlim=c(-10,25),
     ylim=c(0,0.16),type="l",lwd=0.5,lty=2,xlab="Years before/after TasP campaign",ylab="AIDS Death Rate",col="red",axes="F")
axis(1,at=labseq,labels=labseq)
mean_aids_deaths_run1_1_ps <- (run1_1_ps[[1]]$aids_deaths/run1_1_ps[[jj]]$total_infections_alive)/nreps
for (jj in 2:nreps) {
  lines(-20+run1_1_ps[[jj]]$timestep/365,run1_1_ps[[jj]]$aids_deaths/run1_1_ps[[jj]]$total_infections_alive,
        type="l",lwd=0.5,lty=2,col="black")
  mean_aids_deaths_run1_1_ps <- mean_aids_deaths_run1_1_ps + (run1_1_ps[[jj]]$aids_deaths/run1_1_ps[[jj]]$total_infections_alive)/nreps
}
lines(-20+run1_1_ps[[jj]]$timestep/365,mean_aids_deaths_run1_1_ps,type="l",lwd=2,col="black")
axis(2); box()

#Plot AIDS Death Rates for  Strategy 2_1
lines(-20+run2_1_ps[[1]]$timestep/365,run2_1_ps[[1]]$aids_deaths/run2_1_ps[[1]]$total_infections_alive,
      lwd=0.5,lty=2,col="red")
mean_aids_deaths_run2_1_ps <- (run2_1_ps[[1]]$aids_deaths/run2_1_ps[[1]]$total_infections_alive)/nreps
for (jj in 2:nreps) {
  lines(-20+run2_1_ps[[jj]]$timestep/365,run2_1_ps[[jj]]$aids_deaths/run2_1_ps[[jj]]$total_infections_alive,type="l",lwd=0.5,lty=2,col="red")
  mean_aids_deaths_run2_1_ps <- mean_aids_deaths_run2_1_ps + (run2_1_ps[[jj]]$aids_deaths/run2_1_ps[[jj]]$total_infections_alive)/nreps
}
lines(-20+run2_1_ps[[jj]]$timestep/365,mean_aids_deaths_run2_1_ps,type="l",lwd=2,col="red")
#lines(-20+run1_ps[[jj]]$timestep/365,mean_aids_deaths_run1_ps,type="l",lwd=0.5,col="black")
#title(run2_1_title)
legend(-10,0.168,c(run1_title, run2_title),cex=legend_text_size,
       lwd=c(2,2),col=c("black","red"),bty='n',y.intersp = 1)

#Plot AIDS Death Rates for  Strategy 1_2
plot(-20+run1_2_ps[[1]]$timestep/365,run1_2_ps[[1]]$aids_deaths/run1_2_ps[[1]]$total_infections_alive,
     xlim=c(-10,25),axes="F",
     ylim=c(0,0.16),type="l",lwd=0.5,lty=2,xlab="Years before/after TasP campaign",ylab="AIDS Death Rate",col="red")
axis(1,at=labseq,labels=labseq)

mean_aids_deaths_run1_2_ps <- (run1_2_ps[[1]]$aids_deaths/run1_2_ps[[jj]]$total_infections_alive)/nreps
for (jj in 2:nreps) {
  lines(-20+run1_2_ps[[jj]]$timestep/365,run1_2_ps[[jj]]$aids_deaths/run1_2_ps[[jj]]$total_infections_alive,
        type="l",lwd=0.5,lty=2,col="black")
  mean_aids_deaths_run1_2_ps <- mean_aids_deaths_run1_2_ps + (run1_2_ps[[jj]]$aids_deaths/run1_2_ps[[jj]]$total_infections_alive)/nreps
}
lines(-20+run1_2_ps[[jj]]$timestep/365,mean_aids_deaths_run1_2_ps,type="l",lwd=2,col="black")
axis(2); box()

#Plot AIDS Death Rates for  Strategy 2_2
lines(-20+run2_2_ps[[1]]$timestep/365,run2_2_ps[[1]]$aids_deaths/run2_2_ps[[1]]$total_infections_alive,
      lwd=0.5,lty=2,col="red")
mean_aids_deaths_run2_2_ps <- (run2_2_ps[[1]]$aids_deaths/run2_2_ps[[1]]$total_infections_alive)/nreps
for (jj in 2:nreps) {
  lines(-20+run2_2_ps[[jj]]$timestep/365,run2_2_ps[[jj]]$aids_deaths/run2_2_ps[[jj]]$total_infections_alive,type="l",lwd=0.5,lty=2,col="red")
  mean_aids_deaths_run2_2_ps <- mean_aids_deaths_run2_2_ps + (run2_2_ps[[jj]]$aids_deaths/run2_2_ps[[jj]]$total_infections_alive)/nreps
}
lines(-20+run2_2_ps[[jj]]$timestep/365,mean_aids_deaths_run2_2_ps,type="l",lwd=2,col="red")
#lines(-20+run1_ps[[jj]]$timestep/365,mean_aids_deaths_run1_ps,type="l",lwd=0.5,col="black")
#title(run2_2_title)
#legend(10,0.15,c(run1_title, run2_title),cex=legend_text_size,
#       lwd=c(2,2),col=c("black","red"),bty='n',y.intersp = 1)

#Plot AIDS Death Rates for  Strategy 1_3
plot(-20+run1_3_ps[[1]]$timestep/365,run1_3_ps[[1]]$aids_deaths/run1_3_ps[[1]]$total_infections_alive,
     xlim=c(-10,25),axes="F",
     ylim=c(0,0.16),type="l",lwd=0.5,lty=2,xlab="Years before/after TasP campaign",ylab="AIDS Death Rate",col="red")
axis(1,at=labseq,labels=labseq)
mean_aids_deaths_run1_3_ps <- (run1_3_ps[[1]]$aids_deaths/run1_3_ps[[jj]]$total_infections_alive)/nreps
for (jj in 2:nreps) {
  lines(-20+run1_3_ps[[jj]]$timestep/365,run1_3_ps[[jj]]$aids_deaths/run1_3_ps[[jj]]$total_infections_alive,
        type="l",lwd=0.5,lty=2,col="black")
  mean_aids_deaths_run1_3_ps <- mean_aids_deaths_run1_3_ps + (run1_3_ps[[jj]]$aids_deaths/run1_3_ps[[jj]]$total_infections_alive)/nreps
}
lines(-20+run1_3_ps[[jj]]$timestep/365,mean_aids_deaths_run1_3_ps,type="l",lwd=2,col="black")
axis(2); box()

#Plot AIDS Death Rates for  Strategy 2_3
lines(-20+run2_3_ps[[1]]$timestep/365,run2_3_ps[[1]]$aids_deaths/run2_3_ps[[1]]$total_infections_alive,
      lwd=0.5,lty=2,col="red")
mean_aids_deaths_run2_3_ps <- (run2_3_ps[[1]]$aids_deaths/run2_3_ps[[1]]$total_infections_alive)/nreps
for (jj in 2:nreps) {
  lines(-20+run2_3_ps[[jj]]$timestep/365,run2_3_ps[[jj]]$aids_deaths/run2_3_ps[[jj]]$total_infections_alive,type="l",lwd=0.5,lty=2,col="red")
  mean_aids_deaths_run2_3_ps <- mean_aids_deaths_run2_3_ps + (run2_3_ps[[jj]]$aids_deaths/run2_3_ps[[jj]]$total_infections_alive)/nreps
}
lines(-20+run2_3_ps[[jj]]$timestep/365,mean_aids_deaths_run2_3_ps,type="l",lwd=2,col="red")
#lines(-20+run1_ps[[jj]]$timestep/365,mean_aids_deaths_run1_ps,type="l",lwd=0.5,col="black")
#title(run2_3_title)
#legend(10,0.15,c(run1_title, run2_title), cex=legend_text_size,
#       lwd=c(2,2),col=c("black","red"),bty='n',y.intersp = 1)



