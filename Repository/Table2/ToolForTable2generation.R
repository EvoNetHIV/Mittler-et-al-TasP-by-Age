GetTable2Data <- function(tx_file_loc,tx_exp_num) {
  
tx_percent <- c("0","0.4","0.5","0.6","0.7","0.8")
par(mfrow=c(3,2))
tx_type = "random"
for (jj in 1:length(tx_percent)) {
  file_name = paste(tx_file_loc,"txt_opto_",tx_exp_num,"_",tx_type,"_prop",tx_percent[jj],".RData",sep="")
  cat(".")
  #cat(paste("About to load '",file_name,"'\n",sep=""))
  load(file_name)
  if (jj==1) rand <- extract_evonet_data(evomodel)
      else   rand <- cbind(rand,extract_evonet_data(evomodel))
  rm(evomodel)
}

tx_type = "under25_random"
for (jj in 1:length(tx_percent)) {
  file_name = paste(tx_file_loc,"txt_opto_",tx_exp_num,"_",tx_type,"_prop",tx_percent[jj],".RData",sep="")
  cat(".")
  #cat(paste("About to load '",file_name,"'\n",sep=""))
  load(file_name)
  if (jj==1) u25_rand <- extract_evonet_data(evomodel)
  else   u25_rand <- cbind(u25_rand,extract_evonet_data(evomodel))
  rm(evomodel)
}

tx_type = "under25_under30_random"
for (jj in 1:length(tx_percent)) {
  file_name = paste(tx_file_loc,"txt_opto_",tx_exp_num,"_",tx_type,"_prop",tx_percent[jj],".RData",sep="")
  cat(".")
  #cat(paste("About to load '",file_name,"'\n",sep=""))
  load(file_name)
  if (jj==1) u30_rand <- extract_evonet_data(evomodel)
  else   u30_rand <- cbind(u30_rand,extract_evonet_data(evomodel))
  rm(evomodel)
}

tx_type = "under25_under30_under35_random"
for (jj in 1:length(tx_percent)) {
  file_name = paste(tx_file_loc,"txt_opto_",tx_exp_num,"_",tx_type,"_prop",tx_percent[jj],".RData",sep="")
  cat(".")
  #cat(paste("About to load '",file_name,"'\n",sep=""))
  load(file_name)
  if (jj==1) u35_rand <- extract_evonet_data(evomodel)
  else   u35_rand <- cbind(u35_rand,extract_evonet_data(evomodel))
  rm(evomodel)
}

tx_type = "S6.0_S5.5_S5.0_S4.5_S4.0_random"
for (jj in 1:length(tx_percent)) {
  file_name = paste(tx_file_loc,"txt_opto_",tx_exp_num,"_",tx_type,"_prop",tx_percent[jj],".RData",sep="")
  cat(".")
  #cat(paste("About to load '",file_name,"'\n",sep=""))
  load(file_name)
  if (jj==1) spvl_rand <- extract_evonet_data(evomodel)
  else   spvl_rand <- cbind(spvl_rand,extract_evonet_data(evomodel))
  rm(evomodel)
}

tx_type = "CD4_nadir_under500_random"
for (jj in 1:length(tx_percent)) {
  file_name = paste(tx_file_loc,"txt_opto_",tx_exp_num,"_",tx_type,"_prop",tx_percent[jj],".RData",sep="")
  cat(".")
  #cat(paste("About to load '",file_name,"'\n",sep=""))
  load(file_name)
  if (jj==1) cd4_rand <- extract_evonet_data(evomodel)
  else   cd4_rand <- cbind(cd4_rand,extract_evonet_data(evomodel))
  rm(evomodel)
}
ave_incid_0tx <- (unlist(u35_rand[4,1]) + unlist(u30_rand[4,1]) + unlist(u25_rand[4,1]) + unlist(spvl_rand[4,1]) + unlist(cd4_rand[4,1]) + unlist(rand[4,1]))/6
cat("\nFor ",tx_file_loc,"(ave incidence no tx",ave_incid_0tx,") : \n")
cat(" u35:", min(30 + 10*which(unlist( u35_rand[4,2:6])/unlist( u35_rand[4,1]) < 0.05)))
cat(" u30:", min(30 + 10*which(unlist( u30_rand[4,2:6])/unlist( u30_rand[4,1]) < 0.05)))
cat(" u25:", min(30 + 10*which(unlist( u25_rand[4,2:6])/unlist( u25_rand[4,1]) < 0.05)))
cat(" spvl:", min(30 + 10*which(unlist(spvl_rand[4,2:6])/unlist(spvl_rand[4,1]) < 0.05)))
cat(" cd4:", min(30 + 10*which(unlist( cd4_rand[4,2:6])/unlist( cd4_rand[4,1]) < 0.05)))
cat(" rand:", min(30 + 10*which(unlist(     rand[4,2:6])/unlist(     rand[4,1]) < 0.05)))
cat("\n")

}




