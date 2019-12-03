#load("H:/TasP2/Exp1_Base10000_second_try/txt_opto_Exp1_paramsBase10000__random_prop0.RData")

par(mfrow=c(4,2))
inf_start <- 20*365
inf_end <- 25*365
age_groups_lower <- c(15,20,25,30,35,40,45,50)
age_groups_upper <- c(20,25,30,35,40,45,50,55)

#tiff("NewSuppFigInfectorsAgesFixedPart2.tif", res=600, compression = "lzw", height=8, width=7.322917, units="in")
#for (ii in 1: 4) {
for (ii in 5: 8) {
#for (ii in 1: length(age_groups)) {
  
  if (ii >= 5) break_dist <- c(15,20,25,30,35,40,45,50,55,60,65,70,75,80)
          else break_dist <- c(15,20,25,30,35,40,45,50,55,60,65)
  
  for (jj in 1:16) {
    ppp <- evomodel$pop[[jj]]
    base_crit <-  ppp$Time_Inf < 7300 & ppp$Time_Inf < 9125 &
      (ppp$Status == 1 | ppp$Status == -2) & is.na(ppp$Donors_Index) == FALSE
    
    fem_inf  <- which(base_crit & ppp$sex == "f" &
                     ppp$age_infection > age_groups_lower[ii] &
                     ppp$age_infection < age_groups_upper[ii])
    male_inf <- which(base_crit & ppp$sex == "m" &
                     ppp$age_infection > age_groups_lower[ii] &
                     ppp$age_infection < age_groups_upper[ii])
    cat("For rep",jj,"ages:",round(age_groups_lower[ii]),"-",round(age_groups_upper[ii]),":",
          length(male_inf),"M",length(fem_inf),"F\n")
    fh <- hist(ppp$Donors_age[fem_inf],breaks=break_dist,plot=FALSE)
    mh <- hist(ppp$Donors_age[male_inf],breaks=break_dist,plot=FALSE)
    if(ii==1) fh_young_women <- fh
    if (jj==1) {
      fh_tot <- fh
      mh_tot <- mh
    } else{
      fh_tot$counts <- fh_tot$counts + fh$counts
      mh_tot$counts <- mh_tot$counts + mh$counts
    }
  }
  fh_tot$counts <- fh_tot$counts/16
  mh_tot$counts <- mh_tot$counts/16
  
  ccat = cut(fh_tot$breaks, c(0,6+5*ii,11+5*ii,90))
  main_title <- paste("Females:",round(age_groups_lower[ii]),"-",round(age_groups_upper[ii]-1))
  if (ii == 18) main_title <- "Females > 50"
  plot(fh_tot, xlab = "Age of male infector",
       main=main_title,col=c("white","red","white")[ccat])
  
  ccat = cut(mh_tot$breaks, c(0,6+5*ii,11+5*ii,90))
  main_title <- paste("Males:",round(age_groups_lower[ii]),"-",round(age_groups_upper[ii]-1))
  if (ii == 18) main_title <- "Males > 50"
  plot(mh_tot, xlab = "Age of female infector",
       main=main_title,col=c("white","red","white")[ccat])
 
}
#dev.off()

