# Functions

# Prepare data
get_data_for_stan_model_mab <- function(threshold.seropositive,n.age.classes,m0){
  
  df_2006 <- readRDS(file=paste0("data/seroprev_all_ages_year_2006_thres_",threshold.seropositive,".rds"))
  df_2011 <- readRDS(file=paste0("data/seroprev_all_ages_year_2011_thres_",threshold.seropositive,".rds"))
  df_2016 <- readRDS(file=paste0("data/seroprev_all_ages_year_2016_thres_",threshold.seropositive,".rds"))
  df_2017 <- readRDS(file=paste0("data/seroprev_all_ages_year_2017_thres_",threshold.seropositive,".rds"))
  
  # Correct NAs
  df_2006$Seropositive[is.na(df_2006$Total)] = 0
  df_2006$Total[is.na(df_2006$Total)] = 0
  
  df_2011$Seropositive[is.na(df_2011$Total)] = 0
  df_2011$Total[is.na(df_2011$Total)] = 0
  
  df_2016$Seropositive[is.na(df_2016$Total)] = 0
  df_2016$Total[is.na(df_2016$Total)] = 0
  
  df_2017$Seropositive[is.na(df_2017$Total)] = 0
  df_2017$Total[is.na(df_2017$Total)] = 0
  
  age_classes_over1_v <- paste0(seq(1,n.age.classes)," yo")
  age_classes_under1_v <- paste0(seq(0,11)," mo")
  
  ent68_dat <- list(year1 = 2006,
                    year2 = 2011,
                    year3 = 2017,
                    year3_1 = 2016,
                    N = n.age.classes, 
                    N1 = 12,
                    m0 = m0,
                    Decimal_age = df_2006$Decimal_age[df_2006$Age_class %in% age_classes_over1_v],
                    Decimal_age_under1 = df_2006$Decimal_age[df_2006$Age_class %in% age_classes_under1_v],
                    n_2006 = df_2006$Total[df_2006$Age_class %in% age_classes_over1_v],
                    z_2006 = df_2006$Seropositive[df_2006$Age_class %in% age_classes_over1_v],
                    n_2011 = df_2011$Total[df_2011$Age_class %in% age_classes_over1_v],
                    z_2011 = df_2011$Seropositive[df_2011$Age_class %in% age_classes_over1_v],
                    n_2017 = df_2017$Total[df_2017$Age_class %in% age_classes_over1_v],
                    z_2017 = df_2017$Seropositive[df_2017$Age_class %in% age_classes_over1_v],
                    n1_2006 = df_2006$Total[df_2006$Age_class %in% age_classes_under1_v],
                    z1_2006 = df_2006$Seropositive[df_2006$Age_class %in% age_classes_under1_v],
                    n1_2011 = df_2011$Total[df_2011$Age_class %in% age_classes_under1_v],
                    z1_2011 = df_2011$Seropositive[df_2011$Age_class %in% age_classes_under1_v],
                    n1_2017 = df_2016$Total[df_2016$Age_class %in% age_classes_under1_v],
                    z1_2017 = df_2016$Seropositive[df_2016$Age_class %in% age_classes_under1_v]
                    )
  
  return(ent68_dat)
}


# Return name of stan program
get_model_name_str <- function(model_abb){
  if(!model_abb %in% c("mab_const","mab_rw")){
    stop("The model does not exist")
  }else{
    if(model_abb == "mab_const") return("model_ent68_mab_const.stan")
    if(model_abb == "mab_rw") return("model_ent68_mab_rw.stan")
  }
}


# Set function to generate initial values
set_init_f <- function(model_abb){

  if(model_abb %in% c("mab_const","mab_rw")){
    init_f1 <- "random"
  }
  
  return(init_f1)
}


# Given a fit object, extract vector of parameters
get_model_parameters <- function(fit){
  
  idx_omega = which(fit@model_pars == "omega")
  model_pars = fit@model_pars[1:idx_omega]
  
  # get dimension of parameters
  pars_dims = fit@par_dims
  # create vector of parameters to be returned
  pars_v <- c()
  for(i in 1:length(model_pars)){
    ind_par_dim = pars_dims[[i]]
    if(length(ind_par_dim) == 0){
      pars_v = c(pars_v,model_pars[i])
    }else{
      if(ind_par_dim > 1){
        new_pars = paste0(model_pars[i],"[", 1:ind_par_dim, "]")
        pars_v = c(pars_v,new_pars)
      }
    }
  }
    
  return(pars_v)
}



