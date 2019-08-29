par(pty = "s")
par(mfrow=c(2,2))

for (repl in 1:4) {

ave_age_still_have_sex <- 100

#age of partners
elr <- evomodel$el[[repl]][[1]]  # Compacted notation for edgelist (not sure why we need double [[1]][[1]]'s)
agevec <- evomodel$pop[[repl]]$age  # Ages according to the dat$pop structure
sexvec <- evomodel$pop[[repl]]$sex
id_trans <- evomodel$attr[[repl]]$id  # Translates btw dat$pop IDs and dat$el IDs
age.el <- cbind(agevec[id_trans[elr[,1]]],agevec[id_trans[elr[,2]]])  # Need to insert translation for fast edgelist
#plot(age.el[,1],age.el[,2],xaxt = "n", xlab="Age partner 1",ylab="Age partner 2") 


par(pty = "s")


#parter ages by sex
sx1=sexvec[id_trans[elr[,1]]]
sx2=sexvec[id_trans[elr[,2]]]
ix1=which(sx1=="f")
ix2=which(sx2=="f")
femvec=c(age.el[ix1,1],age.el[ix2,2])
malevec=c(age.el[ix1,2],age.el[ix2,1])
ave_age_vec <- (femvec + malevec)/2
ave_age_still_sex <- which(ave_age_vec < ave_age_still_have_sex)
one_one_points <- c(16,85)
#plot(one_one_points,one_one_points,type="l",col="blue",lwd=2,xaxt = "n",xlab="Female Age",ylab="Male Age")
#plot(femvec,malevec,xaxt = "n",xlab="Female Age",ylab="Male Age",xlim=c(10,100),ylim=c(10,100))
#axis(1,at=c(20,30,40,50,60,70,80,90),labels = c("20","30","40","50","60","70","80","90"))
#axis(2,at=c(20,30,40,50,60,70,80,90),labels = c("20","30","40","50","60","70","80","90"))
#lines(femvec,malevec,type ="p")
#lines(one_one_points,one_one_points,type="l",col="blue",lwd=2)

plot(femvec[ave_age_still_sex],malevec[ave_age_still_sex],xaxt = "n",xlab="Female Age",ylab="Male Age")
axis(1,at=c(20,30,40,50,60,70,80,90),labels = c("20","30","40","50","60","70","80","90"))
axis(2,at=c(20,30,40,50,60,70,80,90),labels = c("20","30","40","50","60","70","80","90"))
lines(lowess(femvec,malevec),col="red",lwd=2)
lines(one_one_points,one_one_points,type="l",col="blue",lwd=2)

#plot(malevec[ave_age_still_sex],femvec[ave_age_still_sex],xaxt = "n",xlab="Male Age",ylab="Female Age")
#axis(1,at=c(20,30,40,50,60,70,80,90),labels = c("20","30","40","50","60","70","80","90"))
#axis(2,at=c(20,30,40,50,60,70,80,90),labels = c("20","30","40","50","60","70","80","90"))
#lines(lowess(malevec,femvec),col="red",lwd=2)
#lines(one_one_points,one_one_points,type="l",col="blue",lwd=2)

}
