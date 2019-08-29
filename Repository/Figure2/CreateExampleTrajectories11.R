#note, graphing routine starts line 163, need to change pathways below

re_read_data <- FALSE  # Set to false if tinkering with data already loaded into memory

if (re_read_data == TRUE) {
  load("H:/TasP2/Exp39_N100000/txt_opto_Exp39_paramsBase20000_nodropFig1__CD4_nadir_under500_random_prop0.5.RData")
  run1_title <- "'CD4 Under 500'"
  run1=evomodel
  remove(evomodel)
  load("H:/TasP2/Exp39_N100000/txt_opto_Exp39_paramsBase20000_nodropFig1__random_prop0.5.RData")
  run2_title <- "'Random (Untargeted)'"
  run2=evomodel
  remove(evomodel)
  load("H:/TasP2/Exp39_N100000//txt_opto_Exp39_paramsBase20000_nodropFig1__under25_under30_random_prop0.5.RData")
  run3_title <- "'Under Age 30'"
  run3=evomodel
  remove(evomodel)
  run1_ps <- run1$popsumm
  run2_ps <- run2$popsumm
  run3_ps <- run3$popsumm
  
}

max_infected <- 1500
max_incidence <- 1.24
max_death_rate <- 0.145
max_prev <- 0.123

par(mfrow=c(4,3))
nreps = 16


#####################################
# First row
yylab1="Prevalence"
#yylab2="Percent treated"
yylab2=""

par(mfrow=c(3,3),mgp=c(1.2,.2,0),tcl = -0.2,mar=c(2,2.5,1,1),oma=c(2,3.2,2,2))
for(ii in 1:3){
  if(ii==1)mod=run1
  else if(ii==2)mod=run2
  else mod=run3
  out=summary_popsumm_mat(mod)
  ps=out$popsumm_mats
  timev=(ps$timestep$timestep[[1]]/365)-20
  
  mean_tx_perc <- rep(0,length=length(timev))
  mean_targeted <- rep(0,length=length(timev))
  mean_prev <- rep(0,length=length(timev))
  
  yvec=seq(0,.2,length=length(timev))
  if (ii>=2) yylab1 = ""
  plot(timev,yvec,type='n',ylab="",xlab="",cex.lab=1.1,
       xlim=c(-10,25),ylim=c(0,max_prev)) #cex.lab=0.8,cex.axis=0.8)
  title(ylab=yylab1, mgp=c(1.6,1,0),cex.lab=1.1)
  for(jj in 1:16) {
    lines(timev,ps$prevalence$values[jj,],col="red",lty=2)
    mean_prev <- mean_prev + ps$prevalence$values[jj,]
  }
  lines(timev,mean_prev/16,lwd=2, col="red")
  
  yvec=seq(0,1.0,length=length(timev))
  par(new=T)
  #plot(timev,yvec,type='n',ylab="",xlab="",xlim=c(-10,25),
  #     axes=F,ylim=c(0,100))
 # axis(4)
  for(jj in 1:16){
    lines(timev,
          ps$prevalence$values[jj,]*ps$no_treated$values[jj,]/ps$total_infections_alive$values[jj,],lty=2,
          col="blue")
    lines(timev,
          ps$prevalence$values[jj,]*(ps$prioritized_tx$values[jj,]/ps$no_treated$values[jj,])*(ps$no_treated$values[jj,]/ps$total_infections_alive$values[jj,]),lty=2,
          col="purple")
    mean_tx_perc <- mean_tx_perc + ps$prevalence$values[jj,]*ps$no_treated$values[jj,]/ps$total_infections_alive$values[jj,]/16
    mean_targeted <- mean_targeted + ps$prevalence$values[jj,]*(ps$prioritized_tx$values[jj,]/ps$no_treated$values[jj,])*(ps$no_treated$values[jj,]/ps$total_infections_alive$values[jj,])/16
  }
  
  lines(timev,mean_targeted,lwd=2, col="purple")
  lines(timev,mean_tx_perc,lwd=2, col="blue")

  if(ii==1) { 
    tt=run1_title
    fig_label_with_offsets("(a)",x_off=0.77,y_off= 0.94)
  }
  if(ii==2) {
    tt=run2_title
    fig_label_with_offsets("(b)",x_off=0.87,y_off= 0.94)
  }
  if (ii==3) {
    tt=run3_title
    fig_label_with_offsets("(c)",x_off=1.05,y_off= 0.94)
  }
  #mtext(tt,side=3,line=1)
  #mtext(yylab1,side=2,line=2,cex=.7)
  #mtext("Percent Treated",side=4,line=1.3,cex=.7)
  title(tt,cex.main=0.90)
  legend(-12,1.05*max_prev,c("HIV+","Treated HIV+","Targeted Treated HIV+"),bty="n",col=c("red","blue","purple"),
         cex=.75,lwd=c(2,2,2),lty=1)
  
  
}


# Plot incidence for under strategy 1
extract_evonet_data(run1,max_incidence,overlay=0,graphcol=1)
title(run1_title,cex.main=0.90)
fig_label_with_offsets("(d)",x_off=0.77)
#legend(-12,1.06*max_incidence,run1_title,cex=0.8,
 #      lwd=c(2),col=c("red"),bty='n',y.intersp = 1)

# Plot incidence for strategy 2
extract_evonet_data(run2,max_incidence,overlay=0, graphcol=2)
title(run2_title,cex.main=0.90)
fig_label_with_offsets("(e)",x_off=0.87)
#extract_evonet_data(run1,max_incidence,overlay=1)
#legend(-12,1.06*max_incidence,c(run2_title, paste("Overlay: Mean ",run1_title,sep="")),cex=0.8,
#       lwd=c(2,1),col=c("red","black"),bty='n',y.intersp = 1)

# Plot incidence for strategy 3
#extract_evonet_data(run1,max_incidence,overlay=1
extract_evonet_data(run3,max_incidence,overlay=0, graphcol=3)
title(run3_title,cex.main=0.90)
fig_label_with_offsets("(f)",x_off=1.05)
#legend(-12,1.06*max_incidence,c(run3_title, paste("Overlay: Mean ",run1_title,sep="")),cex=0.8,
#       lwd=c(2,1),col=c("red","black"),bty='n',y.intersp = 1)


#extract_evonet_data(run1,overlay=1)
#extract_evonet_data(run3,overlay=2)


#Plot AIDS Death Rates
seq1 <- seq(1,221,by=5)
seq2 <- seq(2,226,by=5)
seq3 <- seq(3,226,by=5)
seq4 <- seq(4,226,by=5)
seq5 <- seq(5,226,by=5)
xvals <- -20+(run1_ps[[1]]$timestep[seq3] + run1_ps[[1]]$timestep[seq4])/(2*365)
pfreq <- run1$param[[1]]$popsumm_frequency/365

for (ii in 1:3) {
  if (ii==1) { run_ps <- run1_ps; yylab = "AIDS Death Rate"}
  if (ii==2) { run_ps <- run2_ps; yylab = ""}
  if (ii==3) { run_ps <- run3_ps; yylab = ""}
  mean_aids_deaths <- rep(0,length(seq1))
  for (jj in 1:nreps) {
    yvals <- run_ps[[jj]]$aids_deaths[seq1]/run_ps[[jj]]$total_infections_alive[seq1] +
             run_ps[[jj]]$aids_deaths[seq2]/run_ps[[jj]]$total_infections_alive[seq2] +
             run_ps[[jj]]$aids_deaths[seq3]/run_ps[[jj]]$total_infections_alive[seq3] +
             run_ps[[jj]]$aids_deaths[seq4]/run_ps[[jj]]$total_infections_alive[seq4] +
             run_ps[[jj]]$aids_deaths[seq5]/run_ps[[jj]]$total_infections_alive[seq5] 
    yvals <- yvals/(5*pfreq)
    if (jj == 1) {
       plot(xvals,yvals, xlim=c(-10,25),
         ylim=c(0,max_death_rate),type="l",lwd=0.5,lty=2,xlab="",ylab="",cex.lab = 1.1,col="red")
       title(ylab=yylab, mgp=c(1.6,1,0),cex.lab=1.1)
    } else {
      lines(xvals,yvals, type="l",lwd=0.5,lty=2,col="red")
    }
    mean_aids_deaths <- mean_aids_deaths + yvals/nreps
  }
  lines(xvals,mean_aids_deaths,type="l",lwd=2,col="red")
  
  if (ii == 1) {
    title(run1_title,cex.main= 0.9)
    legend(-12,1.06*max_death_rate,run1_title,cex=0.75, lwd=c(2),col=c("red"),bty='n',y.intersp = 1)
    fig_label_with_offsets("(g)",x_off=0.77,y_off=1.07)
    mean_aids_deaths_run1_ps <- mean_aids_deaths
  }
  if (ii == 2) {
    lines(xvals,mean_aids_deaths_run1_ps,type="l",lwd=0.5,col="black")
    title(run2_title,cex.main= 0.9)
    fig_label_with_offsets("(h)",x_off=0.87,y_off=1.07)
    legend(-12,1.06*max_death_rate,c(run2_title,paste("Overlay: Mean ",run1_title,sep="")),cex=0.75, lwd=c(2,0.5),col=c("red","black"),bty='n',y.intersp = 1)
  }
  if (ii == 3) {
    lines(xvals,mean_aids_deaths_run1_ps,type="l",lwd=0.5,col="black")
    title(run3_title,cex.main= 0.9)
    legend(-12,1.06*max_death_rate,c(run3_title,paste("Overlay: Mean ",run1_title,sep="")),cex=0.75, lwd=c(2,0.5),col=c("red","black"),bty='n',y.intersp = 1)
    fig_label_with_offsets("(i)",x_off=1.05,y_off=1.07)
  }
}

mtext("Years before/after TasP campaign",outer=T,side=1,cex=0.8)
