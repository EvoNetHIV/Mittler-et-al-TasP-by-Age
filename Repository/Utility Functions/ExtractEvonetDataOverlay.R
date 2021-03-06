
extract_evonet_data <-function(model_obj,max_incidence,overlay,graphcol){
  #plots incidence rate
  #input netsim object
  #output: plot of incidence rate by year
  
  
  model <- model_obj
  #model <- model_obj 
  total_sims <- model$param[[1]]$nsims
  result_list <- vector('list',length=total_sims)
  yylab = "Incidence"
  if (graphcol >=2) yylab = ""
  
  
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
  if (overlay==0) {
    for(jj in 1:total_sims)
    {
      #browser()
      if(jj==1){
        if(nsteps>=365){
            plot((1:nyears)-20,result_list[[jj]],type='l',xlim=c(-10,25),ylim=c(0,max_incidence),xlab="",ylab="",cex.lab=1.1,
               #xlim=c(10,nyears),
               lty=2,col="red",axes=F)
          title(ylab=yylab, mgp=c(1.6,1,0),cex.lab=1.1)
          labseq=seq(0,nyears,by=5)-20
          #labseq=labseq[-1]
          axis(1,at=labseq,labels=labseq)
        }
        
        if(nsteps<365){
          plot(nyears,result_list[[jj]],type='l',xlim=c(-10,25),ylim=c(0,max_incidence),xlab="",ylab="",cex.lab=1.1,
               lty=2,col="red",axes=T)
          #axis(1,at=1:nyears,labels=1:nyears)
          title(ylab=yylab, mgp=c(1.6,1,0),cex.lab=1.1)
        }
        axis(2);box()
      } else {
        lines((1:nyears)-20,result_list[[jj]],type='l',col="red",ylim=c(0,max(inc)),lty=2)
      }
    }
    lines((1:nyears)-20,colMeans(result_mat),type='l',lwd=2,col="red")
  } else {
    if (overlay==1) lines((1:nyears)-20,colMeans(result_mat),type='l',lwd=1,col="black")
    else lines((1:nyears)-20,colMeans(result_mat),type='l',lwd=0.5,col="green")
  }
  
  #mtext("incidence",side=3,line=2.7,col="blue")
  
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