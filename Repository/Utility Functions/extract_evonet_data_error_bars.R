#example: incidence_rate_plot_v2(evomodel)

extract_evonet_data <-function(model_obj){
  #plots incidence rate
  #input netsim object
  #output: plot of incidence rate by year
  
  
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
  for(jj in 1:total_sims)
  {
    #browser()
    if(jj==1){
      if(nsteps>=365){
        plot(1:nyears,result_list[[jj]],type='o',ylim=c(0,ymax),ylab=NA,xlab="years",
             pch=1,lty=2,col="darkblue",axes=F)
        axis(1,at=1:nyears,labels=1:nyears)
      }
      
      if(nsteps<365){
        plot(nyears,result_list[[jj]],type='o',ylim=c(0,ymax),ylab=NA,xlab="years",
             pch=16,lty=2,col="darkblue",axes=F)
        axis(1,at=nyears,labels=round(nyears,1))
      }
      axis(2);box()
    }else{
      lines(1:nyears,result_list[[jj]],type='o',col="darkblue",ylim=c(0,max(inc)),pch=1,lty=2)
    }
  }
  lines(1:nyears,colMeans(result_mat),type='o',pch=16,lwd=2,col="darkblue")
  mtext("incidence rate per 100 person years",side=3,line=2.7,col="blue")

  # Now calculate mean and standard deviations from the last 10 years
  last_recorded_step <- length(model$popsumm[[jj]]$prevalence)
  after_rampup <- 20
  
  mean_incid <- 0 ; mean_prev <- 0; mean_pills_taken <- 0; mean_not_prioritized_yr10 <- 0; mean_died_AIDS <- 0
  mean_died_AIDS_d00 <- 0; mean_died_AIDS_d03 <- 0; mean_died_AIDS_d05 <- 0; mean_died_AIDS_d07 <- 0
  d00_vec <- 1.00 ^ (0:(last_recorded_step - after_rampup +1))
  d03_vec <- 0.97 ^ (0:(last_recorded_step - after_rampup +1))
  d05_vec <- 0.95 ^ (0:(last_recorded_step - after_rampup +1))
  d07_vec <- 0.93 ^ (0:(last_recorded_step - after_rampup +1))
  
  for(jj in 1:total_sims)
  {
    last_recorded_step <- length(model$popsumm[[jj]]$prevalence)
    after_rampup <- 20
    
    mean_incid <- mean_incid + mean(result_list[[jj]][40:45])
    mean_prev <- mean_prev + model$popsumm[[jj]]$prevalence[last_recorded_step]
    mean_pills_taken <- mean_pills_taken + sum(model$popsumm[[jj]]$total_pills_taken[last_recorded_step])
    mean_not_prioritized_yr10 <- mean_not_prioritized_yr10 + model$param[[jj]]$num_randomly_chosen_start_campaign
    mean_died_AIDS <- mean_died_AIDS + sum(model$popsumm[[jj]]$aids_deaths[after_rampup:last_recorded_step])
    mean_died_AIDS_d00 <- mean_died_AIDS_d00 + sum(model$popsumm[[jj]]$aids_deaths[after_rampup:last_recorded_step]*d00_vec)
    mean_died_AIDS_d03 <- mean_died_AIDS_d03 + sum(model$popsumm[[jj]]$aids_deaths[after_rampup:last_recorded_step]*d03_vec)
    mean_died_AIDS_d05 <- mean_died_AIDS_d05 + sum(model$popsumm[[jj]]$aids_deaths[after_rampup:last_recorded_step]*d05_vec)
    mean_died_AIDS_d07 <- mean_died_AIDS_d07 + sum(model$popsumm[[jj]]$aids_deaths[after_rampup:last_recorded_step]*d07_vec)
  }
  mean_incid <- mean_incid/total_sims
  mean_prev  <- mean_prev/total_sims
  mean_pills_taken <- mean_pills_taken/total_sims
  mean_not_prioritized_yr10 <- mean_not_prioritized_yr10/total_sims
  mean_died_AIDS <- mean_died_AIDS/total_sims
  mean_died_AIDS_d00 <- mean_died_AIDS_d00/total_sims
  mean_died_AIDS_d03 <- mean_died_AIDS_d03/total_sims
  mean_died_AIDS_d05 <- mean_died_AIDS_d05/total_sims
  mean_died_AIDS_d07 <- mean_died_AIDS_d07/total_sims
 
  ss_incid <- 0 ; ss_prev <- 0; ss_pills_taken <- 0; ss_not_prioritized_yr10 <- 0; ss_died_AIDS <- 0
  ss_died_AIDS_d00 <- 0; ss_died_AIDS_d03 <- 0; ss_died_AIDS_d05 <- 0; ss_died_AIDS_d07 <- 0
  for(jj in 1:total_sims)
  {
    ss_incid <- ss_incid + (mean(result_list[[jj]][40:45]) - mean_incid)^2
    ss_prev <- ss_prev + (model$popsumm[[jj]]$prevalence[last_recorded_step] - mean_prev)^2
    ss_pills_taken <- ss_pills_taken + (model$popsumm[[jj]]$total_pills_taken[last_recorded_step] - mean_pills_taken)^2
    ss_not_prioritized_yr10 <- ss_not_prioritized_yr10 + (model$param[[jj]]$num_randomly_chosen_start_campaign - mean_not_prioritized_yr10)^2
    ss_died_AIDS <- ss_died_AIDS + (sum(model$popsumm[[jj]]$aids_deaths[after_rampup:last_recorded_step]) - mean_died_AIDS)^2
    ss_died_AIDS_d00 <- ss_died_AIDS_d00 + (sum(model$popsumm[[jj]]$aids_deaths[after_rampup:last_recorded_step]*d00_vec) - mean_died_AIDS_d00)^2
    ss_died_AIDS_d03 <- ss_died_AIDS_d03 + (sum(model$popsumm[[jj]]$aids_deaths[after_rampup:last_recorded_step]*d03_vec) - mean_died_AIDS_d03)^2
    ss_died_AIDS_d05 <- ss_died_AIDS_d05 + (sum(model$popsumm[[jj]]$aids_deaths[after_rampup:last_recorded_step]*d05_vec) - mean_died_AIDS_d05)^2
    ss_died_AIDS_d07 <- ss_died_AIDS_d07 + (sum(model$popsumm[[jj]]$aids_deaths[after_rampup:last_recorded_step]*d07_vec) - mean_died_AIDS_d07)^2
  }
  sd_incid <- sqrt(ss_incid/(total_sims-1))
  sd_prev <- sqrt(ss_prev/(total_sims-1))
  sd_pills_taken <- sqrt(ss_pills_taken/(total_sims-1))
  sd_not_prioritized_yr10 <- sqrt(ss_not_prioritized_yr10/(total_sims-1))
  sd_died_AIDS <- sqrt(ss_died_AIDS/(total_sims-1))
  sd_died_AIDS_d00 <- sqrt(ss_died_AIDS_d00/(total_sims-1))
  sd_died_AIDS_d03 <- sqrt(ss_died_AIDS_d03/(total_sims-1))
  sd_died_AIDS_d05 <- sqrt(ss_died_AIDS_d05/(total_sims-1))
  sd_died_AIDS_d07 <- sqrt(ss_died_AIDS_d07/(total_sims-1))
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
                                                       "died_AIDS_high" = mean_died_AIDS + sd_died_AIDS,
                         "died_AIDS_d00"      = mean_died_AIDS_d00,
                         "died_AIDS_d00_low"  = mean_died_AIDS_d00 - sd_died_AIDS_d00,
                         "died_AIDS_d00_high" = mean_died_AIDS_d00 + sd_died_AIDS_d00,
                         "died_AIDS_d03"      = mean_died_AIDS_d03,
                         "died_AIDS_d03_low"  = mean_died_AIDS_d03 - sd_died_AIDS_d03,
                         "died_AIDS_d03_high" = mean_died_AIDS_d03 + sd_died_AIDS_d03,
                         "died_AIDS_d05"      = mean_died_AIDS_d05,
                         "died_AIDS_d05_low"  = mean_died_AIDS_d05 - sd_died_AIDS_d05,
                         "died_AIDS_d05_high" = mean_died_AIDS_d05 + sd_died_AIDS_d05,
                         "died_AIDS_d07"      = mean_died_AIDS_d07,
                         "died_AIDS_d07_low"  = mean_died_AIDS_d07 - sd_died_AIDS_d07,
                         "died_AIDS_d07_high" = mean_died_AIDS_d07 + sd_died_AIDS_d07
  )
  
  
  
  return(evonet_results)
}
