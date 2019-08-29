partner_list=curr_model$partner_list[[jrep]]
ppp <- curr_model$pop[[jrep]]

no_partners=unlist(lapply(partner_list,length))
no_partners[which(is.na(partner_list))]=0

arrival_time=evomodel$pop[[jrep]]$arrival_time
arrival_time[is.na(arrival_time)]=0
time_death=evomodel$pop[[jrep]]$Time_Death
time_death[which(is.na(time_death))]=curr_model$param[[jrep]]$n_steps
life_days=time_death-arrival_time
life_years=life_days/365
life_decades = life_days/3650
partners_per_year=no_partners/life_years
partners_per_decade=no_partners/life_decades

short_life <- which(life_decades < 1)
partner_data <- partners_per_decade[-short_life]
time_window <- which(time_death < curr_model$param[[jrep]]$n_steps)

pd_subset <- which (ppp$age > lower_age & ppp$age < upper_age & ppp$Status >=0)
partner_data <- no_partners[pd_subset]

cat("For",main_title,"mean=",mean(partner_data),"Var = ",var(partner_data),
     "V/m=",var(partner_data)/mean(partner_data),"\n")

partner_data_pois = rpois(length(partner_data),mean(partner_data))
hist_data <- hist(partner_data,plot=FALSE,type="h")
hist_data_pois <-  hist(partner_data_pois,plot=FALSE,type="h")
if (graph_hists == TRUE) {
  hist(partner_data,main=main_title,type="h",xlab ="Number of Partners")
  hist(partner_data_pois,breaks=nbreaks,main=paste("Poisson ",main_title,sep=""),xlab ="Number of Partners",type="h")
}#plot(log(hist_data$mids),log(hist_data$counts),main=main_title,xlab = "Log Partners", ylab= "Log Density")
counts_no_zeros <- hist_data$counts
zeros <- which(counts_no_zeros == 0)
counts_no_zeros[zeros] <- 1
counts_no_zeros_pois <- hist_data_pois$counts
zeros <- which(counts_no_zeros_pois == 0)
counts_no_zeros_pois[zeros] <- 1
max_val <- max(log(counts_no_zeros))
max_val2 <- max(log(rev(cumsum(rev(counts_no_zeros)))))
min_val2 <- min(log(rev(cumsum(rev(counts_no_zeros))))-max_val2)
if (cum_only == FALSE) {
  plot(log(hist_data$mids),log(counts_no_zeros),main=main_title,ylim=c(0,max_val),
     xlab = "Log Partners", ylab= "Log Density",type="p")

  lines(log(hist_data_pois$mids),log(counts_no_zeros_pois),type="p",col="green")
}

plot(log(hist_data$mids),log(rev(cumsum(rev(counts_no_zeros))))-max_val2,
     main=main_title,ylim=c(min_val2,0),
     xlab = "Log Partners", ylab= "Log Density",col="blue",type="p",pch=0)

lines(c(log(min(hist_data$mids)),log(hist_data_pois$mids)),
        c(0,log(rev(cumsum(rev(counts_no_zeros_pois))))-max_val2),
        type="l",col="blue")

 