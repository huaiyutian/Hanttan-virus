library(zoo)


data.clean <- function(dataset){
  
  # We load the data and create our main dataframe PR, which stands for "Puerto Rico"
  
  ################# Add Gaussian noise #####  perAAvirus
  for (i in 1:length(monthdata$perAAvirus)){
    if (is.na(monthdata$perAAvirus[i])){
      monthdata$perAAvirus[i] <- 0
    }
  }
  monthdata$perAAvirus[monthdata$perAAvirus == 0] <- 0.01
  # set.seed(123)  
  # mean_noise <- 0
  # sd_noise <- 0.5
  # monthdata$perAAvirus <- monthdata$perAAvirus+ abs(rnorm(length(monthdata), mean = mean_noise, sd = sd_noise))
  # monthdata$perAAvirus <- abs(monthdata$perAAvirus)
  
  ############ Add Gaussian noise #####   perAA
  for (i in 1:length(monthdata$perAA)){
    if (is.na(monthdata$perAA[i])){
      monthdata$perAA[i] <- 0
    }
  }
  monthdata$perAA[monthdata$perAA == 0] <- 0.01

  ##############AAdensity
  monthdata$AAdensity[monthdata$AAdensity< 0] <- 0
  monthdata$AAdensity[monthdata$AAdensity == 0] <- 0.01

  ############## RFdensity
  monthdata$RFdensity[monthdata$RFdensity< 0] <- 0
  # 
  ############## RNdensity
  monthdata$RNdensity[monthdata$RNdensity< 0] <- 0

  PR <- monthdata 
  return(PR)
}


#######################################################################

###################Generate control variables for 12 months #############

data.control <- function(PR){
  for (i in 1:nrow(PR)){
    if (PR$month[[i]] == 1){PR$jan[[i]] <- 1
    }else{PR$jan[[i]] <- 0}
  }
  PR$jan <- do.call(rbind, PR$jan)


  for (i in 1:nrow(PR)){
    if (PR$month[[i]] == 2){PR$feb[[i]] <- 1
    }else{PR$feb[[i]] <- 0}
  }
  PR$feb <- do.call(rbind, PR$feb)

  for (i in 1:nrow(PR)){
    if (PR$month[[i]] == 3){PR$mar[[i]] <- 1
    }else{PR$mar[[i]] <- 0}
  }
  PR$mar <- do.call(rbind, PR$mar)

  for (i in 1:nrow(PR)){
    if (PR$month[[i]] == 4){PR$apr[[i]] <- 1
    }else{PR$apr[[i]] <- 0}
  }
  PR$apr <- do.call(rbind, PR$apr)

  for (i in 1:nrow(PR)){
    if (PR$month[[i]] == 5){PR$may[[i]] <- 1
    }else{PR$may[[i]] <- 0}
  }
  PR$may <- do.call(rbind, PR$may)

  for (i in 1:nrow(PR)){
    if (PR$month[[i]] == 6){PR$jun[[i]] <- 1
    }else{PR$jun[[i]] <- 0}
  }
  PR$jun <- do.call(rbind, PR$jun)

  for (i in 1:nrow(PR)){
    if (PR$month[[i]] == 7){PR$jul[[i]] <- 1
    }else{PR$jul[[i]] <- 0}
  }
  PR$jul <- do.call(rbind, PR$jul)

  for (i in 1:nrow(PR)){
    if (PR$month[[i]] == 8){PR$aug[[i]] <- 1
    }else{PR$aug[[i]] <- 0}
  }
  PR$aug <- do.call(rbind, PR$aug)

  for (i in 1:nrow(PR)){
    if (PR$month[[i]] == 9){PR$sep[[i]] <- 1
    }else{PR$sep[[i]] <- 0}
  }
  PR$sep <- do.call(rbind, PR$sep)

  for (i in 1:nrow(PR)){
    if (PR$month[[i]] == 10){PR$oct[[i]] <- 1
    }else{PR$oct[[i]] <- 0}
  }
  PR$oct <- do.call(rbind, PR$oct)

  for (i in 1:nrow(PR)){
    if (PR$month[[i]] == 11){PR$nov[[i]] <- 1
    }else{PR$nov[[i]] <- 0}
  }
  PR$nov <- do.call(rbind, PR$nov)

  for (i in 1:nrow(PR)){
    if (PR$month[[i]] == 12){PR$dec[[i]] <- 1
    }else{PR$dec[[i]] <- 0}
  }
  PR$dec <- do.call(rbind, PR$dec)

  return(PR)
}

# Add lagged and averaged lagged variables
data.process <- function(PR){
  for (i in 1:3){
    PR[[paste("perAA.lag_", i, sep="")]] <- data.table::shift(PR$perAA, i)
  }

  for (i in 1:60){
    PR[[paste("avg_patch_size.lag_", i, sep="")]] <- data.table::shift(PR$avg_patch_size, i)
  }

  for (i in 1:24){
    PR[[paste("avg_temp.lag_", i, sep="")]] <- data.table::shift(PR$avg_temp, i)
  }

  for (i in 1:24){
    PR[[paste("rainfall.lag_", i, sep="")]] <- data.table::shift(PR$rainfall, i)
  }

  for (i in 1:3){
    PR[[paste("perAAvirus.lag_", i, sep="")]] <- data.table::shift(PR$perAAvirus, i)
  }

  for (i in 1:60){
    PR[[paste("AAdensity.lag_", i, sep="")]] <- data.table::shift(PR$AAdensity, i)
  }

  for (i in 1:12){
    PR[[paste("CTdensity.lag_", i, sep="")]] <- data.table::shift(PR$CTdensity, i)
  }

  for (i in 1:16){
    PR[[paste("RNdensity.lag_", i, sep="")]] <- data.table::shift(PR$RNdensity, i)
  }

  for (i in 1:12){
    PR[[paste("avg_relative_humidity.lag_", i, sep="")]] <- data.table::shift(PR$AAdensity, i)
  }
  for (i in 24:48){
    PR[[paste("growth_rate.lag_", i, sep="")]] <- data.table::shift(PR$avg_patch_size, i)
  }
  for (i in 24:48){
    PR[[paste("num_patch_cul.lag_", i, sep="")]] <- data.table::shift(PR$num_patch_cul, i)
  }

  return(PR)
}
