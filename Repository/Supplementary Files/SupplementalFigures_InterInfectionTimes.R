#load("H:/TasP2/Exp1_Base10000_second_try/txt_opto_Exp1_paramsBase10000__random_prop0.RData")
ppp <- evomodel$pop[[1]]
par(mfrow=c(1,1))
inf_start <- 20*365
inf_end <- 45*365
infected <-  which(ppp$Time_Inf > inf_start & ppp$Time_Inf < inf_end &
              (ppp$Status == 1 | ppp$Status == -2) & is.na(ppp$Donors_Index) == FALSE)

infected_with_nas <-  which(ppp$Time_Inf > inf_start & ppp$Time_Inf < inf_end &
                     (ppp$Status == 1 | ppp$Status == -2))   # For cross-checking.

main_tile = "Histogram of inter-infection times"
histdata <- hist(ppp$Donors_Total_Time_Inf_At_Trans[infected]/365, plot=FALSE, right=FALSE)
hist(ppp$Donors_Total_Time_Inf_At_Trans[infected]/365,main="",
     xlab = "Inter-infection time (years)")

v1 <- round(mean(ppp$Donors_Total_Time_Inf_At_Trans[infected]/365),1)
v2 <- round(median(ppp$Donors_Total_Time_Inf_At_Trans[infected]/365),1)
box_text <- paste("Mean:",v1,"\nMedian:",v2)
                    
#text(20,1000,box_text,cex=0.9)
legend(7,max(histdata$counts),box.lty=0,
       leg=paste(c("Mean ","Median "),c(v1,v2),sep="= "))
cat(length(infected)/length(infected_with_nas))


