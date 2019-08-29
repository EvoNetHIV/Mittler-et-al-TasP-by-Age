par(mfrow=c(3,2))
tx_percent <- c("0","0.3","0.4","0.45","0.5","0.55","0.6","0.65","0.7","0.75","0.8","0.85","0.9","1")
tx_file_loc <- "H:/TasP2/Exp5_Equal_Rel_Durs_Lower_Mean_Duration/"
tx_exp_num <- "Exp5_2paramsEqualRelDursN10000_"
tx_type = "random"
for (jj in 1:length(tx_percent)) {
  file_name = paste(tx_file_loc,"txt_opto_",tx_exp_num,"_",tx_type,"_prop",tx_percent[jj],".RData",sep="")
  cat(paste("About to load '",file_name,"'\n",sep=""))
  load(file_name)
  if (jj==1) rand <- extract_evonet_data(evomodel)
      else   rand <- cbind(rand,extract_evonet_data(evomodel))
  rm(evomodel)
}

tx_type = "under25_random"
for (jj in 1:length(tx_percent)) {
  file_name = paste(tx_file_loc,"txt_opto_",tx_exp_num,"_",tx_type,"_prop",tx_percent[jj],".RData",sep="")
  cat(paste("About to load '",file_name,"'\n",sep=""))
  load(file_name)
  if (jj==1) u25_rand <- extract_evonet_data(evomodel)
  else   u25_rand <- cbind(u25_rand,extract_evonet_data(evomodel))
  rm(evomodel)
}

tx_type = "under25_under30_random"
for (jj in 1:length(tx_percent)) {
  file_name = paste(tx_file_loc,"txt_opto_",tx_exp_num,"_",tx_type,"_prop",tx_percent[jj],".RData",sep="")
  cat(paste("About to load '",file_name,"'\n",sep=""))
  load(file_name)
  if (jj==1) u30_rand <- extract_evonet_data(evomodel)
  else   u30_rand <- cbind(u30_rand,extract_evonet_data(evomodel))
  rm(evomodel)
}

tx_type = "under25_under30_under35_random"
for (jj in 1:length(tx_percent)) {
  file_name = paste(tx_file_loc,"txt_opto_",tx_exp_num,"_",tx_type,"_prop",tx_percent[jj],".RData",sep="")
  cat(paste("About to load '",file_name,"'\n",sep=""))
  load(file_name)
  if (jj==1) u35_rand <- extract_evonet_data(evomodel)
  else   u35_rand <- cbind(u35_rand,extract_evonet_data(evomodel))
  rm(evomodel)
}

tx_type = "S6.0_S5.5_S5.0_S4.5_S4.0_random"
for (jj in 1:length(tx_percent)) {
  file_name = paste(tx_file_loc,"txt_opto_",tx_exp_num,"_",tx_type,"_prop",tx_percent[jj],".RData",sep="")
  cat(paste("About to load '",file_name,"'\n",sep=""))
  load(file_name)
  if (jj==1) spvl_rand <- extract_evonet_data(evomodel)
  else   spvl_rand <- cbind(spvl_rand,extract_evonet_data(evomodel))
  rm(evomodel)
}

tx_type = "CD4_nadir_under500_random"
for (jj in 1:length(tx_percent)) {
  file_name = paste(tx_file_loc,"txt_opto_",tx_exp_num,"_",tx_type,"_prop",tx_percent[jj],".RData",sep="")
  cat(paste("About to load '",file_name,"'\n",sep=""))
  load(file_name)
  if (jj==1) cd4_rand <- extract_evonet_data(evomodel)
  else   cd4_rand <- cbind(cd4_rand,extract_evonet_data(evomodel))
  rm(evomodel)
}





