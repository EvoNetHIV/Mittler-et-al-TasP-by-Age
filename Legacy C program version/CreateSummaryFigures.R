SimDir = "Exp22_100_Targeting_After_Random"  # 
gs <- read.table(paste(SimDir,"/GrandSummary.txt",sep=""),header=TRUE)
u25 <- gs[,"u25"]
perc_tx <- gs[,"tx"] 
yearly_incr_tx <- gs[,"yearly_incr"]
DaysTx <- gs[,"PerDaysTx"]
NonTargetedTx <- gs[,"NotTargeted"]
AllTreated <- gs[,"AllTx"]
strategy <- gs[,"Strat"]
repl <- gs[,"repl"]
dAIDS <- gs[,"dAIDS"]
D00 <- gs[,"D00"]
D03 <- gs[,"D03"]
D05 <- gs[,"D05"]
D07 <- gs[,"D07"]

HIV_pos <- gs[,"HIV."]
Untreated <- gs[,"Untreated"]
total_strategies <- unique(strategy)
strategies <- c(1:5)  # Must include 1 as first element
replicates <- unique(repl)
treatments <- unique(perc_tx)
final_susc_d <- gs[,"susc_d"]
final_new_infs <- gs[,"new_infs"]
Final_Incid <- 100*final_new_infs / (final_susc_d / 365)
line_width = 1

par(mfrow=c(2,2))
under25_index <- which(gs[,"u25"]==1 & gs[,"repl"] ==1)
rand_index <- which(gs[,"u25"]==1)
jit <- 0.05
#perc_treated <- c(0,10,20,30,40,50,60,70,80,90,100)
perc_treated <- c(7.0,7.2,7.4,7.6,7.8,8.0,8.2,8.4,8.6,8.8,9)/100
#perc_treated <- c(40)
legend_vector <- c("Random","Under age 25","Under age 30","Under age 35","CD4<500","Under25, CD4<500","SPVL > 5.5","SPVL > 5.0", "SPVL > 4.5","M27, W23","Men","Women")
color_vector <- c("black","darkred","red","orange","green","violet","darkblue","blue","lightblue","purple","brown","grey")
pch_vector <- c(0,1,2,3,4,5,6)

legend_vector <- legend_vector[strategies]
color_legend_vector <- color_vector[strategies]
pch_legend_vector <- pch_vector[strategies]

#### Plot final incidence ###
gs_mean <- matrix(0,nrow=length(total_strategies),ncol=length(treatments))
gs_sum <- matrix(0,nrow=length(total_strategies),ncol=length(treatments))
gs_SS  <- matrix(0,nrow=length(total_strategies),ncol=length(treatments))
gs_lower  <- matrix(0,nrow=length(total_strategies),ncol=length(treatments))
gs_upper  <- matrix(0,nrow=length(total_strategies),ncol=length(treatments))
for (ii in (strategies-min(strategies)+ 1)) {
  for (jj in 1:length(treatments)) {
    for (kk in 1: length(replicates)) {
      index <- which(gs[,"Strat"]==ii & gs[,"tx"] == treatments[jj] & gs[,"repl"]==replicates[kk])
      gs_sum[ii,jj] = gs_sum[ii,jj] + Final_Incid[index]  # This line references quantity to be plotted
    }
    gs_mean[ii,jj] <- gs_sum[ii,jj]/length(replicates)
    for (kk in 1: length(replicates)) {
      index <- which(gs[,"Strat"]==ii & gs[,"tx"] == treatments[jj] & gs[,"repl"]==replicates[kk])
      gs_SS[ii,jj] <-  gs_SS[ii,jj] + (Final_Incid[index] - gs_mean[ii,jj])^2   # This line references quantity to be plotted
    }
    #gs_lower[ii,jj] <- gs_mean[ii,jj] - 2.131 * sqrt(gs_SS[ii,jj])/(length(replicates)-1)
    #gs_upper[ii,jj] <- gs_mean[ii,jj] + 2.131 * sqrt(gs_SS[ii,jj])/(length(replicates)-1)
    gs_lower[ii,jj] <- gs_mean[ii,jj] -  sqrt((gs_SS[ii,jj])/(length(replicates)-1))
    gs_upper[ii,jj] <- gs_mean[ii,jj] + sqrt((gs_SS[ii,jj])/(length(replicates)-1))
  }
}
plot(perc_treated-4*jit,gs_mean[1,],col=color_vector[1],type="b",xlim=c(0,100),lwd=line_width,pch=pch_vector[1],
     ylim=c(0,max(gs_upper)),
     xlab="TasP Target [Percent Treated]",
     ylab="Incidence Last 10 Years")
arrows(perc_treated-4*jit, gs_lower[1,], perc_treated-4*jit, gs_upper[1,], code=3, length=0.02, angle = 90, col=color_vector[1])
for (ii in strategies[-1]) {
  lines(perc_treated-(5-ii)*jit,gs_mean[ii,],col=color_vector[ii],type = "b", lwd=line_width,pch = pch_vector[ii])
  arrows(perc_treated-(5-ii)*jit, gs_lower[ii,], perc_treated-(5-ii)*jit, gs_upper[ii,], code=3, length=0.02, angle = 90, col=color_vector[ii])
}


legend( x="topright", 
        legend=legend_vector,
        #        legend=c("CD4<500","Random","Under age 25"),
        col=  color_legend_vector, 
        text.col= color_legend_vector, 
        lwd=line_width,bty="n",cex = 0.8,pt.cex=1.0,
        pch=pch_legend_vector) # 0 1 6 5

#    ylab="Incidence Last 10 Years")
#    ylab="Person-years of therapy")
#    ylab="Percent Treated Non-targeted")

#Final_Incid[index]
#DaysTx[index] 
#NonTargetedTx[index]

### Plot AIDS Deaths ###
gs_mean <- matrix(0,nrow=length(total_strategies),ncol=length(treatments))
gs_sum <- matrix(0,nrow=length(total_strategies),ncol=length(treatments))
gs_SS  <- matrix(0,nrow=length(total_strategies),ncol=length(treatments))
gs_lower  <- matrix(0,nrow=length(total_strategies),ncol=length(treatments))
gs_upper  <- matrix(0,nrow=length(total_strategies),ncol=length(treatments))
for (ii in strategies) {
  for (jj in 1:length(treatments)) {
    for (kk in 1: length(replicates)) {
      index <- which(gs[,"Strat"]==ii & gs[,"tx"] == treatments[jj] & gs[,"repl"]==replicates[kk])
      gs_sum[ii,jj] = gs_sum[ii,jj] + dAIDS[index]  # This line references quantity to be plotted
    }
    gs_mean[ii,jj] <- gs_sum[ii,jj]/length(replicates)
    for (kk in 1: length(replicates)) {
      index <- which(gs[,"Strat"]==ii & gs[,"tx"] == treatments[jj] & gs[,"repl"]==replicates[kk])
      gs_SS[ii,jj] <-  gs_SS[ii,jj] + (dAIDS[index] - gs_mean[ii,jj])^2   # This line references quantity to be plotted
    }
    #gs_lower[ii,jj] <- gs_mean[ii,jj] - 2.131 * sqrt(gs_SS[ii,jj])/(length(replicates)-1)
    #gs_upper[ii,jj] <- gs_mean[ii,jj] + 2.131 * sqrt(gs_SS[ii,jj])/(length(replicates)-1)
    gs_lower[ii,jj] <- gs_mean[ii,jj] - sqrt((gs_SS[ii,jj])/(length(replicates)-1))
    gs_upper[ii,jj] <- gs_mean[ii,jj] + sqrt((gs_SS[ii,jj])/(length(replicates)-1))
  }
}
plot(perc_treated-4*jit,gs_mean[1,],col=color_vector[1],type="b",xlim=c(0,100),lwd=line_width,pch=pch_vector[1],
     ylim=c(0,max(gs_upper)),
     xlab="TasP Target [Percent Treated]",
     ylab="Total AIDS Deaths")
arrows(perc_treated-4*jit, gs_lower[1,], perc_treated-4*jit, gs_upper[1,], code=3, length=0.02, angle = 90, col=color_vector[1])
for (ii in strategies[-1]) {
  lines(perc_treated-(5-ii)*jit,gs_mean[ii,],col=color_vector[ii],type = "b", lwd=line_width,pch = pch_vector[ii])
  arrows(perc_treated-(5-ii)*jit, gs_lower[ii,], perc_treated-(5-ii)*jit, gs_upper[ii,], code=3, length=0.02, angle = 90, col=color_vector[ii])
}

if (1==2) {legend( x="topright", 
                   legend=legend_vector,
                   #        legend=c("CD4<500","Random","Under age 25"),
                   col=  color_legend_vector, 
                   text.col= color_legend_vector, 
                   lwd=line_width,bty="n",cex = 1,pt.cex=1.0,
                   pch=pch_legend_vector) # 0 1 6 5
}
#### Plot person-years of therapy ####
gs_mean <- matrix(0,nrow=length(total_strategies),ncol=length(treatments))
gs_sum <- matrix(0,nrow=length(total_strategies),ncol=length(treatments))
gs_SS  <- matrix(0,nrow=length(total_strategies),ncol=length(treatments))
gs_lower  <- matrix(0,nrow=length(total_strategies),ncol=length(treatments))
gs_upper  <- matrix(0,nrow=length(total_strategies),ncol=length(treatments))
for (ii in strategies) {
  for (jj in 1:length(treatments)) {
    for (kk in 1: length(replicates)) {
      index <- which(gs[,"Strat"]==ii & gs[,"tx"] == treatments[jj] & gs[,"repl"]==replicates[kk])
      gs_sum[ii,jj] = gs_sum[ii,jj] + DaysTx[index]  # This line references quantity to be plotted
    }
    gs_mean[ii,jj] <- gs_sum[ii,jj]/length(replicates)
    for (kk in 1: length(replicates)) {
      index <- which(gs[,"Strat"]==ii & gs[,"tx"] == treatments[jj] & gs[,"repl"]==replicates[kk])
      gs_SS[ii,jj] <-  gs_SS[ii,jj] + (DaysTx[index] - gs_mean[ii,jj])^2   # This line references quantity to be plotted
    }
    #gs_lower[ii,jj] <- gs_mean[ii,jj] - 2.131 * sqrt(gs_SS[ii,jj])/(length(replicates)-1)
    #gs_upper[ii,jj] <- gs_mean[ii,jj] + 2.131 * sqrt(gs_SS[ii,jj])/(length(replicates)-1)
    gs_lower[ii,jj] <- gs_mean[ii,jj] - sqrt((gs_SS[ii,jj])/(length(replicates)-1))
    gs_upper[ii,jj] <- gs_mean[ii,jj] + sqrt((gs_SS[ii,jj])/(length(replicates)-1))
  }
}
plot(perc_treated-4*jit,gs_mean[1,],col=color_vector[1],type="b",xlim=c(0,100),lwd=line_width,pch=pch_vector[1],
     ylim=c(0,max(gs_upper)),
     xlab="TasP Target [Percent Treated]",
     ylab="Person-Years of Therapy")
arrows(perc_treated-4*jit, gs_lower[1,], perc_treated-4*jit, gs_upper[1,], code=3, length=0.02, angle = 90, col=color_vector[1])
for (ii in strategies[-1]) {
  lines(perc_treated-(5-ii)*jit,gs_mean[ii,],col=color_vector[ii],type = "b", lwd=line_width,pch = pch_vector[ii])
  arrows(perc_treated-(5-ii)*jit, gs_lower[ii,], perc_treated-(5-ii)*jit, gs_upper[ii,], code=3, length=0.02, angle = 90, col=color_vector[ii])
}
if (1==2) {legend( x="bottomright", 
                   legend=legend_vector,
                   #        legend=c("CD4<500","Random","Under age 25"),
                   col=  color_legend_vector, 
                   text.col= color_legend_vector, 
                   lwd=line_width,bty="n",cex = 1,pt.cex=1.2,
                   pch=pch_legend_vector) # 0 1 6 5
}

#### Plot percent non-targeted ###
perc_treated <- perc_treated
gs_mean <- matrix(0,nrow=length(total_strategies),ncol=length(treatments))
gs_sum <- matrix(0,nrow=length(total_strategies),ncol=length(treatments))
gs_SS  <- matrix(0,nrow=length(total_strategies),ncol=length(treatments))
gs_lower  <- matrix(0,nrow=length(total_strategies),ncol=length(treatments))
gs_upper  <- matrix(0,nrow=length(total_strategies),ncol=length(treatments))
for (ii in strategies) {
  for (jj in 1:length(treatments)) {
    for (kk in 1: length(replicates)) {
      index <- which(gs[,"Strat"]==ii & gs[,"tx"] == treatments[jj] & gs[,"repl"]==replicates[kk])
      gs_sum[ii,jj] = gs_sum[ii,jj] + 100*NonTargetedTx[index]  # This line references quantity to be plotted
    }
    gs_mean[ii,jj] <- gs_sum[ii,jj]/length(replicates)
    for (kk in 1: length(replicates)) {
      index <- which(gs[,"Strat"]==ii & gs[,"tx"] == treatments[jj] & gs[,"repl"]==replicates[kk])
      gs_SS[ii,jj] <-  gs_SS[ii,jj] + (100*NonTargetedTx[index] - gs_mean[ii,jj])^2   # This line references quantity to be plotted
    }
    #gs_lower[ii,jj] <- gs_mean[ii,jj] - 2.131 * sqrt(gs_SS[ii,jj])/(length(replicates)-1)
    #gs_upper[ii,jj] <- gs_mean[ii,jj] + 2.131 * sqrt(gs_SS[ii,jj])/(length(replicates)-1)
    gs_lower[ii,jj] <- gs_mean[ii,jj] - sqrt((gs_SS[ii,jj])/(length(replicates)-1))
    gs_upper[ii,jj] <- gs_mean[ii,jj] + sqrt((gs_SS[ii,jj])/(length(replicates)-1))
  }
}
plot(perc_treated-4*jit,gs_mean[1,],col=color_vector[1],type="b",xlim=c(0,100),lwd=line_width,pch=pch_vector[1],
     ylim=c(0,max(gs_upper)),
     xlab="TasP Target [Percent Treated]",
     ylab="Percent Treated Non-targeted")
arrows(perc_treated-4*jit, gs_lower[1,], perc_treated-4*jit, gs_upper[1,], code=3, length=0.02, angle = 90, col=color_vector[1])
for (ii in strategies[-1]) {
  lines(perc_treated-(5-ii)*jit,gs_mean[ii,],col=color_vector[ii],type = "b", lwd=line_width,pch = pch_vector[ii])
  arrows(perc_treated-(5-ii)*jit, gs_lower[ii,], perc_treated-(5-ii)*jit, gs_upper[ii,], code=3, length=0.02, angle = 90, col=color_vector[ii])
}

if (1==2) {legend( x="bottomright", 
                   legend=legend_vector,
                   #        legend=c("CD4<500","Random","Under age 25"),
                   col=  color_legend_vector, 
                   text.col= color_legend_vector, 
                   lwd=line_width,bty="n",cex = 1,pt.cex=1.2,
                   pch=pch_legend_vector) # 0 1 6 5
}


### Plot DALYs discount 0 ###
gs_mean <- matrix(0,nrow=length(total_strategies),ncol=length(treatments))
gs_sum <- matrix(0,nrow=length(total_strategies),ncol=length(treatments))
gs_SS  <- matrix(0,nrow=length(total_strategies),ncol=length(treatments))
gs_lower  <- matrix(0,nrow=length(total_strategies),ncol=length(treatments))
gs_upper  <- matrix(0,nrow=length(total_strategies),ncol=length(treatments))
for (ii in strategies) {
  for (jj in 1:length(treatments)) {
    for (kk in 1: length(replicates)) {
      index <- which(gs[,"Strat"]==ii & gs[,"tx"] == treatments[jj] & gs[,"repl"]==replicates[kk])
      gs_sum[ii,jj] = gs_sum[ii,jj] + D00[index]  # This line references quantity to be plotted
    }
    gs_mean[ii,jj] <- gs_sum[ii,jj]/length(replicates)
    for (kk in 1: length(replicates)) {
      index <- which(gs[,"Strat"]==ii & gs[,"tx"] == treatments[jj] & gs[,"repl"]==replicates[kk])
      gs_SS[ii,jj] <-  gs_SS[ii,jj] + (D00[index] - gs_mean[ii,jj])^2   # This line references quantity to be plotted
    }
    #gs_lower[ii,jj] <- gs_mean[ii,jj] - 2.131 * sqrt(gs_SS[ii,jj])/(length(replicates)-1)
    #gs_upper[ii,jj] <- gs_mean[ii,jj] + 2.131 * sqrt(gs_SS[ii,jj])/(length(replicates)-1)
    gs_lower[ii,jj] <- gs_mean[ii,jj] - sqrt((gs_SS[ii,jj])/(length(replicates)-1))
    gs_upper[ii,jj] <- gs_mean[ii,jj] + sqrt((gs_SS[ii,jj])/(length(replicates)-1))
  }
}
plot(perc_treated-4*jit,gs_mean[1,],col=color_vector[1],type="b",xlim=c(min(perc_treated),max(perc_treated)),lwd=line_width,pch=pch_vector[1],
     ylim=c(0,max(gs_upper)),
     xlab="TasP Target [Percent Treated]",
     ylab="DALYS (0% discount rate)")
arrows(perc_treated-4*jit, gs_lower[1,], perc_treated-4*jit, gs_upper[1,], code=3, length=0.02, angle = 90, col=color_vector[1])
for (ii in strategies[-1]) {
  lines(perc_treated-(5-ii)*jit,gs_mean[ii,],col=color_vector[ii],type = "b", lwd=line_width,pch = pch_vector[ii])
  arrows(perc_treated-(5-ii)*jit, gs_lower[ii,], perc_treated-(5-ii)*jit, gs_upper[ii,], code=3, length=0.02, angle = 90, col=color_vector[ii])
}

if (1==1) {legend( x="topright", 
                   legend=legend_vector,
                   #        legend=c("CD4<500","Random","Under age 25"),
                   col=  color_legend_vector, 
                   text.col= color_legend_vector, 
                   lwd=line_width,bty="n",cex = 1,pt.cex=1.0,
                   pch=pch_legend_vector) # 0 1 6 5
}

### Plot DALYS, Discount 3% ###
gs_mean <- matrix(0,nrow=length(total_strategies),ncol=length(treatments))
gs_sum <- matrix(0,nrow=length(total_strategies),ncol=length(treatments))
gs_SS  <- matrix(0,nrow=length(total_strategies),ncol=length(treatments))
gs_lower  <- matrix(0,nrow=length(total_strategies),ncol=length(treatments))
gs_upper  <- matrix(0,nrow=length(total_strategies),ncol=length(treatments))
for (ii in strategies) {
  for (jj in 1:length(treatments)) {
    for (kk in 1: length(replicates)) {
      index <- which(gs[,"Strat"]==ii & gs[,"tx"] == treatments[jj] & gs[,"repl"]==replicates[kk])
      gs_sum[ii,jj] = gs_sum[ii,jj] + D03[index]  # This line references quantity to be plotted
    }
    gs_mean[ii,jj] <- gs_sum[ii,jj]/length(replicates)
    for (kk in 1: length(replicates)) {
      index <- which(gs[,"Strat"]==ii & gs[,"tx"] == treatments[jj] & gs[,"repl"]==replicates[kk])
      gs_SS[ii,jj] <-  gs_SS[ii,jj] + (D03[index] - gs_mean[ii,jj])^2   # This line references quantity to be plotted
    }
    #gs_lower[ii,jj] <- gs_mean[ii,jj] - 2.131 * sqrt(gs_SS[ii,jj])/(length(replicates)-1)
    #gs_upper[ii,jj] <- gs_mean[ii,jj] + 2.131 * sqrt(gs_SS[ii,jj])/(length(replicates)-1)
    gs_lower[ii,jj] <- gs_mean[ii,jj] - sqrt((gs_SS[ii,jj])/(length(replicates)-1))
    gs_upper[ii,jj] <- gs_mean[ii,jj] + sqrt((gs_SS[ii,jj])/(length(replicates)-1))
  }
}
plot(perc_treated-4*jit,gs_mean[1,],col=color_vector[1],type="b",xlim=c(min(perc_treated),max(perc_treated)),lwd=line_width,pch=pch_vector[1],
     ylim=c(0,max(gs_upper)),
     xlab="TasP Target [Percent Treated]",
     ylab="DALYS (3% discount rate))")
arrows(perc_treated-4*jit, gs_lower[1,], perc_treated-4*jit, gs_upper[1,], code=3, length=0.02, angle = 90, col=color_vector[1])
for (ii in strategies[-1]) {
  lines(perc_treated-(5-ii)*jit,gs_mean[ii,],col=color_vector[ii],type = "b", lwd=line_width,pch = pch_vector[ii])
  arrows(perc_treated-(5-ii)*jit, gs_lower[ii,], perc_treated-(5-ii)*jit, gs_upper[ii,], code=3, length=0.02, angle = 90, col=color_vector[ii])
}

### Plot DALYS Discount 5% ###
gs_mean <- matrix(0,nrow=length(total_strategies),ncol=length(treatments))
gs_sum <- matrix(0,nrow=length(total_strategies),ncol=length(treatments))
gs_SS  <- matrix(0,nrow=length(total_strategies),ncol=length(treatments))
gs_lower  <- matrix(0,nrow=length(total_strategies),ncol=length(treatments))
gs_upper  <- matrix(0,nrow=length(total_strategies),ncol=length(treatments))
for (ii in strategies) {
  for (jj in 1:length(treatments)) {
    for (kk in 1: length(replicates)) {
      index <- which(gs[,"Strat"]==ii & gs[,"tx"] == treatments[jj] & gs[,"repl"]==replicates[kk])
      gs_sum[ii,jj] = gs_sum[ii,jj] + D05[index]  # This line references quantity to be plotted
    }
    gs_mean[ii,jj] <- gs_sum[ii,jj]/length(replicates)
    for (kk in 1: length(replicates)) {
      index <- which(gs[,"Strat"]==ii & gs[,"tx"] == treatments[jj] & gs[,"repl"]==replicates[kk])
      gs_SS[ii,jj] <-  gs_SS[ii,jj] + (D05[index] - gs_mean[ii,jj])^2   # This line references quantity to be plotted
    }
    #gs_lower[ii,jj] <- gs_mean[ii,jj] - 2.131 * sqrt(gs_SS[ii,jj])/(length(replicates)-1)
    #gs_upper[ii,jj] <- gs_mean[ii,jj] + 2.131 * sqrt(gs_SS[ii,jj])/(length(replicates)-1)
    gs_lower[ii,jj] <- gs_mean[ii,jj] - sqrt((gs_SS[ii,jj])/(length(replicates)-1))
    gs_upper[ii,jj] <- gs_mean[ii,jj] + sqrt((gs_SS[ii,jj])/(length(replicates)-1))
  }
}
plot(perc_treated-4*jit,gs_mean[1,],col=color_vector[1],type="b",xlim=c(min(perc_treated),max(perc_treated)),lwd=line_width,pch=pch_vector[1],
     ylim=c(0,max(gs_upper)),
     xlab="TasP Target [Percent Treated]",
     ylab="DALYS (5% discount rate))")
arrows(perc_treated-4*jit, gs_lower[1,], perc_treated-4*jit, gs_upper[1,], code=3, length=0.02, angle = 90, col=color_vector[1])
for (ii in strategies[-1]) {
  lines(perc_treated-(5-ii)*jit,gs_mean[ii,],col=color_vector[ii],type = "b", lwd=line_width,pch = pch_vector[ii])
  arrows(perc_treated-(5-ii)*jit, gs_lower[ii,], perc_treated-(5-ii)*jit, gs_upper[ii,], code=3, length=0.02, angle = 90, col=color_vector[ii])
}


### Plot DALYs Discount 7% ###
gs_mean <- matrix(0,nrow=length(total_strategies),ncol=length(treatments))
gs_sum <- matrix(0,nrow=length(total_strategies),ncol=length(treatments))
gs_SS  <- matrix(0,nrow=length(total_strategies),ncol=length(treatments))
gs_lower  <- matrix(0,nrow=length(total_strategies),ncol=length(treatments))
gs_upper  <- matrix(0,nrow=length(total_strategies),ncol=length(treatments))
for (ii in strategies) {
  for (jj in 1:length(treatments)) {
    for (kk in 1: length(replicates)) {
      index <- which(gs[,"Strat"]==ii & gs[,"tx"] == treatments[jj] & gs[,"repl"]==replicates[kk])
      gs_sum[ii,jj] = gs_sum[ii,jj] + D07[index]  # This line references quantity to be plotted
    }
    gs_mean[ii,jj] <- gs_sum[ii,jj]/length(replicates)
    for (kk in 1: length(replicates)) {
      index <- which(gs[,"Strat"]==ii & gs[,"tx"] == treatments[jj] & gs[,"repl"]==replicates[kk])
      gs_SS[ii,jj] <-  gs_SS[ii,jj] + (D07[index] - gs_mean[ii,jj])^2   # This line references quantity to be plotted
    }
    #gs_lower[ii,jj] <- gs_mean[ii,jj] - 2.131 * sqrt(gs_SS[ii,jj])/(length(replicates)-1)
    #gs_upper[ii,jj] <- gs_mean[ii,jj] + 2.131 * sqrt(gs_SS[ii,jj])/(length(replicates)-1)
    gs_lower[ii,jj] <- gs_mean[ii,jj] - sqrt((gs_SS[ii,jj])/(length(replicates)-1))
    gs_upper[ii,jj] <- gs_mean[ii,jj] + sqrt((gs_SS[ii,jj])/(length(replicates)-1))
  }
}
plot(perc_treated-4*jit,gs_mean[1,],col=color_vector[1],type="b",xlim=c(min(perc_treated),max(perc_treated)),lwd=line_width,pch=pch_vector[1],
     ylim=c(0,max(gs_upper)),
     xlab="TasP Target [Percent Treated]",
     ylab="DALYS (7% discount rate)")
arrows(perc_treated-4*jit, gs_lower[1,], perc_treated-4*jit, gs_upper[1,], code=3, length=0.02, angle = 90, col=color_vector[1])
for (ii in strategies[-1]) {
  lines(perc_treated-(5-ii)*jit,gs_mean[ii,],col=color_vector[ii],type = "b", lwd=line_width,pch = pch_vector[ii])
  arrows(perc_treated-(5-ii)*jit, gs_lower[ii,], perc_treated-(5-ii)*jit, gs_upper[ii,], code=3, length=0.02, angle = 90, col=color_vector[ii])
}

