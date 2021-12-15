# Housekeeping
rm(list=ls()) #Clear all
graphics.off() #Close all

# Loading required package(s)
require(tseries); require(ggplot2); require(rugarch)

#Loading daily prices of S&P 500
SP500 <- get.hist.quote("^GSPC", provider="yahoo"); head(SP500)
SP500_logreturns <- diff(log(SP500$Close), lag=1)

#Disregarding data before 03/01/2000
SP500_logreturns <- SP500_logreturns[2275:length(SP500_logreturns)]

#### Stylized facts ####
## Plot of raw data
autoplot(SP500$Close) + labs(x="Time", y="Close Price")

## (i) Volatility clustering
autoplot(SP500_logreturns) + labs(x="Time", y="Log Returns")

## (ii) Autocorrelation
ggacfhomemade <- function(series, alpha=0.05, na.action){
  acfdata <- acf(series, na.action=na.action, plot=FALSE)
  df_acfdata <- with(acfdata, data.frame(lag, acf))
  
  lim <- qnorm((1+(1-alpha))/2) / sqrt(acfdata$n.used)
  
  ggplot(data = df_acfdata, mapping = aes(x = lag, y = acf)) +
    geom_hline(aes(yintercept = 0)) +
    geom_segment(mapping = aes(xend = lag, yend = 0))   +
    geom_hline(aes(yintercept = lim), linetype = 2, color = 'blue') +
    geom_hline(aes(yintercept = -lim), linetype = 2, color = 'blue') +
    labs(x = "Lag", y = "ACF")
}
ggacfhomemade(SP500_logreturns, na.action=na.exclude)
ggacfhomemade(abs(SP500_logreturns), na.action=na.exclude)
ggacfhomemade(SP500_logreturns^2, na.action=na.exclude)

## (iii) Leptokurticity
ggqqhomemade <- function(series){
  df_data <- data.frame("series"=series)
  ggplot(df_data, aes(sample=series)) + stat_qq() + stat_qq_line() +
    labs(x="Theoretical Quantiles", y="Sample Quantiles")
}
gghistogramhomemade <- function(series){
  df_data <- data.frame("series"=series)
  ggplot(df_data, aes(series)) + geom_histogram(aes(y=..density..), bins=nclass.FD(series)) +
    stat_function(fun=dnorm, color="red", args=list(mean = mean(series), sd = sd(series))) +    #theoretical Gaussian pdf
    geom_density(color="blue")
}
ggqqhomemade(SP500_logreturns)
gghistogramhomemade(SP500_logreturns)

## (iv) Leverage effects
# Under construction

## Testing GARCH effects using LM test

#### Estimation of GARCH models using data up and including 30/04/2021 ####
## Input:
#'tsdata' is a time series object of the data to be modelled in the GARCH framework
#'x' is a vector setting the maximum considered order of the GARCH model
#'y' is a list of the types of GARCH models to be considered
#'z' is a list of the distributions to be considered
#'solver' specifies the type of solver used in ugarchfit
#'out.sample' specifies the amount of data to be left out in the fitting procedure (used later for forecasting)
gmorder <- function(x){
  output.int <- expand.grid(0:x[1],0:x[2])
  n <- nrow(output.int)
  output <- NULL
  for(i in 1:n){
    output[[i]] <- c(output.int[i,1],output.int[i,2])
  }
  return(output[-1])
}
gmisfit <- function(tsdata,x,y,z,solver,out.sample){
  output <- list(); h <- 1; print(h)
  res.gmorder <- gmorder(x)
  for(i in 1:length(y)){
    print(i)
    if(y[[i]]!="TGARCH" & y[[i]]!="AVGARCH" & y[[i]]!="fiGARCH"){    #Not submodel or fiGARCH
      for(j in 1:length(z)){
        for(k in 1:length(res.gmorder)){
          spec <- ugarchspec(variance.model=list(model=y[[i]], garchOrder=res.gmorder[[k]]), 
                             distribution.model=z[[j]], mean.model=list(armaOrder=c(0,0), include.mean = FALSE))
          fit <- ugarchfit(spec=spec, data=tsdata, out.sample=out.sample, solver=solver)
          if(fit@fit$convergence==0){       #Solver converges
            output[[h]] <- fit; h <- h+1; print(h)
          } else {                          #Solver failed to converge
            output[[h]] <- list("FAILTED TO CONVERGE", res.gmorder[[k]], y[[i]], z[[j]]); h <- h+1; print(h)
          }
        }
      }
    } else if(y[[i]]=="fiGARCH"){                  #fiGARCH
      print("fiGARCH")
      remove.ord.figarch <- rep(0,length(res.gmorder))
      for(l in 1:length(res.gmorder)){
        if(res.gmorder[[l]][1]==0 || res.gmorder[[l]][2]==0){
          remove.ord.figarch[l] <- l
        }
      }
      res.gmorder.figarch <- res.gmorder[-remove.ord.figarch]
      for(j in 1:length(z)){
        for(k in 1:length(res.gmorder.figarch)){
          spec <- ugarchspec(variance.model=list(model=y[[i]], garchOrder=res.gmorder.figarch[[k]]), 
                             distribution.model=z[[j]], mean.model=list(armaOrder=c(0,0)))
          fit <- ugarchfit(spec=spec, data=tsdata, out.sample=out.sample, solver=solver)
          if(fit@fit$convergence==0){       #Solver converges
            output[[h]] <- fit; h <- h+1; print(h)
          } else {                          #Solver failed to converge
            output[[h]] <- list("FAILTED TO CONVERGE", res.gmorder.figarch[[k]], y[[i]], z[[j]]); h <- h+1; print(h)
          }
        }
      }
    } else {                                      #Submodel
      print("Submodel")
      for(j in 1:length(z)){
        for(k in 1:length(res.gmorder)){
          spec <- ugarchspec(variance.model=list(model="fGARCH", submodel=y[[i]], garchOrder=res.gmorder[[k]]), 
                             distribution.model=z[[j]], mean.model=list(armaOrder=c(0,0)))
          fit <- ugarchfit(spec=spec, data=tsdata, out.sample=out.sample, solver=solver)
          if(fit@fit$convergence==0){       #Solver converges
            output[[h]] <- fit; h <- h+1; print(h)
          } else {                          #Solver failed to converge
            output[[h]] <- list("FAILTED TO CONVERGE", res.gmorder[[k]], y[[i]], z[[j]]); h <- h+1; print(h)
          }
        }
      }
    }
  }
  return(output)
}

# Setting input
x1 <- c(3,3)
y1 <- list("sGARCH", "eGARCH", "TGARCH", "gjrGARCH", "apARCH", "AVGARCH", "fiGARCH"); y2 <- list("sGARCH")
z1 <- list("norm", "std", "ged"); z2 <- list("norm")

#Note that the 'out.sample = 156' corresponds leaving out all data after 30/04/2021
allgms <- gmisfit(SP500_logreturns, x1, y1, z1, solver="hybrid", out.sample=156)

#### Model selection according to AIC, BIC, SIC, and HQIC ####
#Important: This procedure only works if z is a list of three distributions.
modelselect.part1 <- function(models,y,z){
  output <- list(); s <- 1
  for(i in 1:length(y)){
    onetypemodel <- list(); h <- 1
    if(y[[i]]!="TGARCH" & y[[i]]!="AVGARCH"){                         #Not submodel
      for(k in 1:length(models)){
        if(models[[k]]@model$modeldesc$vmodel == y[[i]]){
          onetypemodel[[h]] <- models[[k]]; h <- h+1
        }
      }
      if(length(onetypemodel)==0){print(i); print("idiot?")}
      onetypemodel.norm <- list(); h1 <- 1
      onetypemodel.std <- list(); h2 <- 1
      onetypemodel.ged <- list(); h3 <- 1
      for(j in 1:length(onetypemodel)){
        if(onetypemodel[[j]]@model$modeldesc$distribution == z[[1]]){
          onetypemodel.norm[[h1]] <- onetypemodel[[j]]; h1 <- h1+1 
        } else if(onetypemodel[[j]]@model$modeldesc$distribution == z[[2]]){
          onetypemodel.std[[h2]] <- onetypemodel[[j]]; h2 <- h2+1 
        } else if(onetypemodel[[j]]@model$modeldesc$distribution == z[[3]]){
          onetypemodel.ged[[h3]] <- onetypemodel[[j]]; h3 <- h3+1
        }
      }
      output[[s]] <- list(onetypemodel.norm,onetypemodel.std,onetypemodel.ged); s <- s+1
    } else {                                                        #Submodel
      for(k in 1:length(models)){
        if(isTRUE(models[[k]]@model$modeldesc$vsubmodel == y[[i]])){
          onetypemodel[[h]] <- models[[k]]; h <- h+1
        }
      }
      if(length(onetypemodel)==0){print(i); print("idiot?")}
      onetypemodel.norm <- list(); h1 <- 1
      onetypemodel.std <- list(); h2 <- 1
      onetypemodel.ged <- list(); h3 <- 1
      for(j in 1:length(onetypemodel)){
        if(onetypemodel[[j]]@model$modeldesc$distribution == z[[1]]){
          onetypemodel.norm[[h1]] <- onetypemodel[[j]]; h1 <- h1+1 
        } else if(onetypemodel[[j]]@model$modeldesc$distribution == z[[2]]){
          onetypemodel.std[[h2]] <- onetypemodel[[j]]; h2 <- h2+1 
        } else if(onetypemodel[[j]]@model$modeldesc$distribution == z[[3]]){
          onetypemodel.ged[[h3]] <- onetypemodel[[j]]; h3 <- h3+1
        }
      }
      output[[s]] <- list(onetypemodel.norm,onetypemodel.std,onetypemodel.ged); s <- s+1
    }
  }
  return(output)
}
modelselect.part2 <- function(output.modelselect.part1){
  output <- list()
  for(i in 1:length(output.modelselect.part1)){
    output[[i]] <- list()
    onetypealldist <- output.modelselect.part1[[i]]
    for(j in 1:length(onetypealldist)){
      onetypeonedist <- onetypealldist[[j]]
      minaic <- Inf; minbic <- Inf; minhqic <- Inf; minsic <- Inf
      for(k in 1:length(onetypeonedist)){
        onetypeonedistmodel <- onetypeonedist[[k]]
        if(infocriteria(onetypeonedistmodel)[1]<minaic){
          minaic <- infocriteria(onetypeonedistmodel)[1]
          aicmodel <- onetypeonedistmodel
        }
        if(infocriteria(onetypeonedistmodel)[2]<minbic){
          minbic <- infocriteria(onetypeonedistmodel)[2]
          bicmodel <- onetypeonedistmodel
        }
        if(infocriteria(onetypeonedistmodel)[3]<minsic){
          minsic <- infocriteria(onetypeonedistmodel)[3]
          sicmodel <- onetypeonedistmodel
        }
        if(infocriteria(onetypeonedistmodel)[4]<minhqic){
          minhqic <- infocriteria(onetypeonedistmodel)[4]
          hqicmodel <- onetypeonedistmodel
        }
      }
      output[[i]][[j]] <- list(aicmodel,bicmodel,sicmodel,hqicmodel)
    }
  }
  return(output)
}
modelselect <- function(models,y,z){
  output.part1 <- modelselect.part1(models,y,z)
  output <- modelselect.part2(output.part1)
  return(output)
}

ms <- modelselect(allgms, y1, z1)

#Testing procedure above for one type of model with same distribution
#Input: 'gm' is a list of one type of models with same distribution
insample1 <- function(gm){
  output <- list()
  l <- length(gm)
  order <- rep(0,l); pq <- rep(0,l); distribution <- rep(0,l); akaika <- rep(0,l); 
  bayes <- rep(0,l); shibata <- rep(0,l); hannanquinn <- rep(0,l)
  for(i in 1:length(gm)){
    order <- unname(c(gm[[i]]@model$modelinc["alpha"],gm[[i]]@model$modelinc["beta"]))
    pq[i] <- paste0("(",order[1],",",order[2],")")
    distribution[i] <- gm[[i]]@model$modeldesc$distribution
    akaika[i] <- infocriteria(gm[[i]])[1]
    bayes[i] <- infocriteria(gm[[i]])[2]
    shibata[i] <- infocriteria(gm[[i]])[3]
    hannanquinn[i] <- infocriteria(gm[[i]])[4]
  }
  df <- data.frame(pq, distribution, akaika, bayes, shibata, hannanquinn)
  output[[1]] <- df
  output[[2]] <- order(df$akaika)[1:3]
  output[[3]] <- order(df$bayes)[1:3]
  output[[4]] <- order(df$shibata)[1:3]
  output[[5]] <- order(df$hannanquinn)[1:3]
  names(output) <- c("Table(0.1)", "Akaika top three", "Bayes top three", "Shibata top three", "Hannan-Quinn top three")
  return(output)
}
gm.test <- gmisfit(SP500_logreturns, x1, y2, z2, solver="hybrid", out.sample=156)
insample1.test <- insample1(gm.test)[[1]]

#### Forecasting ####
#Input:
#'fitten' is one fitted model from 'ugarchfit'
#'n.ahead' determines the forecast horizon
#'sigma.out' if true only forecasted sigmas are returned
hm.recursive.gmfit <- function(fitten, data, n.ahead, solver){
  if(fitten@model$modeldesc$vmodel!="fGARCH"){                #Not submodel
    modelS <- fitten@model$modeldesc$vmodel
    garchOrderS <- unname(c(fitten@model$modelinc["alpha"],fitten@model$modelinc["beta"]))
    distribution.modelS <- fitten@model$modeldesc$distribution
    out.sampleS <- fitten@model$n.start
    
    output <- list(); h <- 1; #print(h)
    output[[1]] <- fitten; h <- h+1; #print(h)
    
    for(i in 1:(out.sampleS/n.ahead-1)){
      spec <- ugarchspec(variance.model=list(model=modelS, garchOrder=garchOrderS), distribution.model=distribution.modelS, 
                         mean.model=list(armaOrder=c(0,0), include.mean=FALSE))
      output[[h]] <- ugarchfit(spec=spec, data=data, solver=solver, out.sample=(out.sampleS-i))
      h <- h+1; #print(h)
    }
  } else {
    modelS <- fitten@model$modeldesc$vsubmodel
    garchOrderS <- unname(c(fitten@model$modelinc["alpha"],fitten@model$modelinc["beta"]))
    distribution.modelS <- fitten@model$modeldesc$distribution
    out.sampleS <- fitten@model$n.start
    
    output <- list(); h <- 1; #print(h)
    output[[1]] <- fitten; h <- h+1; #print(h)
    
    for(i in 1:(out.sampleS/n.ahead-1)){
      spec <- ugarchspec(variance.model=list(model="fGARCH", submodel=modelS, garchOrder=garchOrderS), distribution.model=distribution.modelS, 
                         mean.model=list(armaOrder=c(0,0), include.mean=FALSE))
      output[[h]] <- ugarchfit(spec=spec, data=data, solver=solver, out.sample=(out.sampleS-i))
      h <- h+1; #print(h)
    }
  }  
  return(output)
}
homemadeugarchforecast.recursive <- function(fitten, data, n.ahead, solver, sigma.out=FALSE){
  gm.rec <- hm.recursive.gmfit(fitten=fitten, data, n.ahead=n.ahead, solver=solver)
  print("Finished recursive fitting")
  output <- list()
  for(i in 1:length(gm.rec)){
    output[[i]] <- ugarchforecast(gm.rec[[i]], n.ahead=n.ahead, solver=solver)
    #print(i)
  }
  if(sigma.out == FALSE){     #Not only return sigma forecasts
    return(output)
  } else {                    #Only return sigma forecasts
    output.sigma <- rep(0, length(output))
    for(i in 1:length(output)){
      output.sigma[i] <- unname(c(output[[i]]@forecast$sigmaFor[n.ahead]))
    }
    return(output.sigma)
  }
}
final.homemadeugarchforecast.recursive <- function(output.modelselect, data, n.ahead, solver){
  all <- unlist(output.modelselect); all <- all[!duplicated(all)]
  output <- list()
  for(i in 1:length(all)){
    output[[i]] <- homemadeugarchforecast.recursive(all[[i]], data, n.ahead, solver)
    print(i)
  }
  return(output)
}

forecastres.1 <- final.homemadeugarchforecast.recursive(ms, data=SP500_logreturns, n.ahead=1, solver="hybrid")

#Testing procedure above using 'ugarchroll' from 'rugarch' package
forecastreshomemade <- homemadeugarchforecast.recursive(ms[[3]][[1]][[1]], n.ahead=1, data=SP500_logreturns, 
                                                        solver="hybrid", sigma.out=TRUE)
getspec(ms[[3]][[1]][[1]])
testspec <- ugarchspec(variance.model=list(model="sGARCH", garchOrder=c(2,2)),
                       mean.model=list(armaOrder=c(0,0), include.mean=FALSE), distribution.model="norm")
forecastrestest <- ugarchroll(testspec, data=SP500_logreturns, n.ahead=1, forecast.length=156,
                              refit.every=1, refit.window="recursive", solver="hybrid")

length(forecastrestest@forecast$density$Sigma)
length(forecastreshomemade)
forecastrestest@forecast$density$Sigma == forecastreshomemade

#Extracting the forecasted sigma
forecastedsigmahomemade.one <- function(oneforecastedmodel, n.ahead){
  output <- list()
  for(i in 1:length(allforecastedmodels)){
    output.sigma <- rep(0, length(allforecastedmodels[[i]]))
    for(j in 1:length(allforecastedmodels[[i]])){
      output.sigma[j] <- unname(c(allforecastedmodels[[i]]@forecast$sigmaFor[n.ahead]))
    }
    output[[i]] <- output.sigma
  }
  return(output)
}
forecastedsigmahomemade.all <- function(allforecastedmodels, n.ahead){
  output <- list()
  for(i in 1:length(allforecastedmodels)){
    oneforecastedmodels <- allforecastedmodels[[i]]
    output.sigma <- rep(0,length(oneforecastedmodels))
    for(j in 1:length(oneforecastedmodels)){
      output.sigma[j] <- unname(c(oneforecastedmodels[[j]]@forecast$sigmaFor[n.ahead]))
    }
    output[[i]] <- output.sigma
  }
  return(output)
}

#### Diebold-Mariano test ####

