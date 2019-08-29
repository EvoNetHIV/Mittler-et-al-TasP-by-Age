reload <- FALSE
if (reload == TRUE) {
  load("H:/TasP2/Exp20_Gradual_95_care_N10000/txt_opto_Exp20_paramsGradual_95_N10000_TasP60perc__random_prop0.2.RData")
  u35 <- evomodel$popsumm
  u35_full <- evomodel
  load("H:/TasP2/Exp20_Gradual_95_care_N10000/txt_opto_Exp20_paramsGradual_95_N10000_TasP60perc__under25_under30_random_prop0.2.RData")
  cd4 <- evomodel$popsumm
  cd4_full <- evomodel
  par(mfrow=c(3,2))
  nreps = 16
}

if (1==2) {
  #Plot number infected and number treated for CD4 strategy
plot(cd4[[1]]$timestep/365,cd4[[1]]$total_infections_alive,ylim=c(0,1600),type="l",lwd=0.5,lty=2,xlab="Years",ylab="Number Infected and Treated",col="red")
lines(cd4[[1]]$timestep/365,cd4[[1]]$no_treated  ,type="l",lwd=0.5,lty=2,col="blue")
mean_no_treated_cd4 <- cd4[[1]]$no_treated/nreps
mean_infected_cd4 <- cd4[[1]]$total_infections_alive/nreps
for (jj in 2:nreps) {
  lines(cd4[[jj]]$timestep/365,cd4[[jj]]$total_infections_alive,type="l",lwd=0.5,lty=2,col="red")
  lines(cd4[[jj]]$timestep/365,cd4[[jj]]$no_treated  ,type="l",lwd=0.5,lty=2,col="blue")
  #lines(cd4[[jj]]$timestep/365,cd4[[jj]]$no_prioritized_treated  ,type="l",lwd=0.5,lty=2,col="orange")
  mean_no_treated_cd4 <- mean_no_treated_cd4 + cd4[[jj]]$no_treated/nreps
  mean_infected_cd4 <- mean_infected_cd4 + cd4[[jj]]$total_infections_alive/nreps
}
lines(cd4[[jj]]$timestep/365,mean_infected_cd4,type="l",lwd=2,col="red")
lines(cd4[[jj]]$timestep/365,mean_no_treated_cd4,type="l",lwd=2,col="blue")

#Plot number infected and number treated for under 25 strategy
plot(u35[[1]]$timestep/365,u35[[1]]$total_infections_alive,ylim=c(0,1600),type="l",lwd=0.5,lty=2,xlab="Years",ylab="Number Infected and Treated",col="red")
lines(u35[[1]]$timestep/365,u35[[1]]$no_treated  ,type="l",lwd=0.5,lty=2,col="blue")
mean_no_treated_u35 <- u35[[1]]$no_treated/nreps
mean_infected_u35 <- u35[[1]]$total_infections_alive/nreps
for (jj in 2:nreps) {
  lines(u35[[jj]]$timestep/365,u35[[jj]]$total_infections_alive,type="l",lwd=0.5,lty=2,col="red")
  lines(u35[[jj]]$timestep/365,u35[[jj]]$no_treated  ,type="l",lwd=0.5,lty=2,col="blue")
  #lines(u35[[jj]]$timestep/365,u35[[jj]]$no_prioritized_treated  ,type="l",lwd=0.5,lty=2,col="orange")
  mean_no_treated_u35 <- mean_no_treated_u35 + u35[[jj]]$no_treated/nreps
  mean_infected_u35 <- mean_infected_u35 + u35[[jj]]$total_infections_alive/nreps
}
lines(u35[[jj]]$timestep/365,mean_infected_u35,type="l",lwd=2,col="red")
lines(u35[[jj]]$timestep/365,mean_no_treated_u35,type="l",lwd=2,col="blue")
}
# Plot incidence for CD4 strategy
extract_evonet_data(cd4_full)

# Plot incidence for under 25 strategy
extract_evonet_data_lines(u35_full)

# Cumulative AIDS Deaths
if (1==2) {
  plot(cd4[[1]]$timestep/365,cumsum(cd4[[1]]$aids_deaths),ylim=c(0,1200),type="l",lwd=0.5,lty=2,xlab="Years",ylab="AIDS Deaths",col="black")
  lines(u35[[1]]$timestep/365,cumsum(u35[[1]]$aids_deaths),type="l",lwd=0.5,lty=2,col="blue")
  mean_aids_deaths_cd4 <- cumsum(cd4[[1]]$aids_deaths)/nreps
  mean_aids_deaths_u35 <- cumsum(u35[[1]]$aids_deaths)/nreps
  for (jj in 2:nreps) {
    lines(cd4[[jj]]$timestep/365,cumsum(cd4[[jj]]$aids_deaths),type="l",lwd=0.5,lty=2,col="black")
    lines(u35[[jj]]$timestep/365,cumsum(u35[[jj]]$aids_deaths),type="l",lwd=0.5,lty=2,col="blue")
    mean_aids_deaths_cd4 <- mean_aids_deaths_cd4 + cumsum(cd4[[jj]]$aids_deaths)/nreps
    mean_aids_deaths_u35 <- mean_aids_deaths_u35 + cumsum(u35[[jj]]$aids_deaths)/nreps
    
  }
  lines(cd4[[jj]]$timestep/365,mean_aids_deaths_cd4,type="l",lwd=2,col="black")
  lines(u35[[jj]]$timestep/365,mean_aids_deaths_u35,type="l",lwd=2,col="blue")
  
  cum_aids_deaths_cd4 <- cd4[[1]]$aids_deaths
}

#Plot AIDS Death Rates for CD4 Strategy
aids_death_rate <- cd4[[1]]$aids_deaths 


plot(cd4[[1]]$timestep/365,cd4[[1]]$aids_deaths,ylim=c(0,70),type="l",lwd=0.5,lty=2,xlab="Years",ylab="AIDS Deaths",col="red")
mean_aids_deaths_cd4 <- cd4[[1]]$aids_deaths/nreps
for (jj in 2:nreps) {
  lines(cd4[[jj]]$timestep/365,cd4[[jj]]$aids_deaths,type="l",lwd=0.5,lty=2,col="red")
  mean_aids_deaths_cd4 <- mean_aids_deaths_cd4 + cd4[[jj]]$aids_deaths/nreps
 }
lines(cd4[[jj]]$timestep/365,mean_aids_deaths_cd4,type="l",lwd=2,col="red")

#Plot AIDS Deaths for under 25 Strategy (with CD4 overlay for comparison)
plot(u35[[1]]$timestep/365,u35[[1]]$aids_deaths,ylim=c(0,70),type="l",lwd=0.5,lty=2,xlab="Years",ylab="AIDS Deaths",col="red")
mean_aids_deaths_u35 <- u35[[1]]$aids_deaths/nreps
for (jj in 2:nreps) {
  lines(u35[[jj]]$timestep/365,u35[[jj]]$aids_deaths,type="l",lwd=0.5,lty=2,col="red")
   mean_aids_deaths_u35 <- mean_aids_deaths_u35 + u35[[jj]]$aids_deaths/nreps
  
}
lines(u35[[jj]]$timestep/365,mean_aids_deaths_u35,type="l",lwd=2,col="red")
lines(cd4[[jj]]$timestep/365,mean_aids_deaths_cd4,type="l",lwd=1,col="black")



