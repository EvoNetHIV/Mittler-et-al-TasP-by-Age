param_set="Base10000"
param_set_name=paste("params",param_set,sep="")
param_file=paste(param_set_name,".R",sep="")

setwd("H://TasP2/Exp1_Base10000_second_try")
source(param_file) # Default parameters for Tas-By-Age study


# Biological parameters specific to this simulation
exp_number = paste("Exp1_",param_set_name,"_",sep="")

percent_treated = list(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)
nw <- nw_setup(evoparams) # Sets up the initial network
nw_saved_john <- nw
cat("Finished setting up the model\n")

for (ii in 1:6){
  evoparams$tx_type <- TasP_criteria[ii] # See "TasP_Params_March2018.R" for list of criteria

   for(jj in 1:length(percent_treated)) {
    model_name = paste("txt_opto_",exp_number,"_",paste(TasP_criteria[ii][[1]],collapse="_"),"_prop",as.character(percent_treated[jj]),sep="")
    evoparams$proportion_treated  = percent_treated[jj][[1]]
    
    cat("About to run evorun: (evoparams$mean_test_interval_under25 = ",evoparams$mean_test_interval_under25,"\n")
    
    nw <- nw_saved_john # Don't assume that evorun leaves nw alone
     evomodel <- evorun(modules,evoparams,nw) # Runs the simulation (Line above b/c cannot be sure that evorun doesn't change nw!)
    closeAllConnections()
    
    assign(model_name,evomodel)
    file_name <- paste(model_name,".RData",sep="")
    save(evomodel,file = file.path(getwd(),file_name) )
    evoplot(model=evomodel, name = paste(model_name,".pdf",sep=""), path= file.path(getwd())) # Saves a pdf with standard plots
    remove(evomodel); remove(list=model_name)
    
    #remove(evomodel); remove(list=model_name)
  } # Percent treated
} # TasP criterion


