rm(list=ls())

dataFileName <- 'analysisData.txt'
source('prepareData.R')
source('MLEfunctions.R')

paramInit <- c(pCom=0.01, pHouse=0.1, dHouse=1, phiV=0.3, phiA=0.8, piV=0.995, piA=0.995)

getModelResults <- function(fixedParams){
	optParams <- setdiff(names(paramInit),names(fixedParams))
	
	optFn <- function(x) -getLL(x['pCom'], Inf, x['pHouse'], x['dHouse'], x['phiV'], x['phiA'], x['piV'], x['piA'])

	if('piV' %in% names(fixedParams)){
		sol <- getSol(paramInit[optParams], function(x) optFn(c(x, fixedParams)))
	}else{
		getModelSol <- function(piV) getSol(paramInit[setdiff(optParams,'piV')],
								function(x) optFn(c(x, fixedParams, piV=piV)))
		piVsol <- optimize(function(x) getModelSol(x)$value, c(0.998,1), tol=1e-10)$minimum
		sol <- getModelSol(piVsol)
		solPiV1 <- getModelSol(1)
		if(solPiV1$value < sol$value){
			piVsol <- 1
			sol <- solPiV1
		}

		sol$par <- c(sol$par,piV=piVsol)
	}
	sol$par <- c(sol$par,fixedParams)[names(paramInit)]
	sol
}

fixP <- vector('list',18)
modelNames <- rep('',18)

modelNames[1] <- 'none'; 					fixP[[1]] <- {}
modelNames[2] <- 'dH=0';					fixP[[2]] <- c(dHouse=0)
modelNames[3] <- 'dH=Inf';					fixP[[3]] <- c(dHouse=Inf)
modelNames[4] <- 'piV=1';					fixP[[4]] <- c(piV=1)
modelNames[5] <- 'phiA=0.893';				fixP[[5]] <- c(phiA=0.893)
modelNames[6] <- 'piA=0.996';					fixP[[6]] <- c(piA=0.996)
modelNames[7] <- 'piA=1';					fixP[[7]] <- c(piA=1)
modelNames[8] <- 'dH=0;piV=1';				fixP[[8]] <- c(dHouse=0,piV=1)
modelNames[9] <- 'dH=0;phiA=0.893';				fixP[[9]] <- c(dHouse=0,phiA=0.893)
modelNames[10] <- 'dH=0;piA=0.996';				fixP[[10]] <- c(dHouse=0,piA=0.996)
modelNames[11] <- 'piV=1;phiA=0.893';			fixP[[11]] <- c(piV=1, phiA=0.893)
modelNames[12] <- 'piV=1;piA=0.996';			fixP[[12]] <- c(piV=1, piA=0.996)
modelNames[13] <- 'phiA=0.893;piA=0.996';			fixP[[13]] <- c(phiA=0.893, piA=0.996)
modelNames[14] <- 'dH=0;piV=1;phiA=0.893';		fixP[[14]] <- c(dHouse=0, piV=1, phiA=0.893)
modelNames[15] <- 'dH=0;piV=1;piA=0.996';			fixP[[15]] <- c(dHouse=0, piV=1, piA=0.996)
modelNames[16] <- 'dH=0;phiA=0.893;piA=0.996';		fixP[[16]] <- c(dHouse=0, phiA=0.893, piA=0.996)
modelNames[17] <- 'phiA=0.893;piV=1;piA=0.996';		fixP[[17]] <- c(phiA=0.893, piV=1, piA=0.996)
modelNames[18] <- 'dH=0;phiA=0.893;piV=1;piA=0.996';	fixP[[18]] <- c(dHouse=0, phiA=0.893, piV=1, piA=0.996)


getPval <- function(LL, optLL, df) pchisq(2*(optLL-LL), df=df, lower.tail=FALSE)
getAIC <- function(LL, numParam) 2*(numParam - LL)

freeParams <- rep(NA,length(modelNames))
table2opt <- matrix(NA,length(modelNames),7)
logLik <- rep(NA,length(modelNames))
pVal <- rep(NA,length(modelNames))
aic <- rep(NA,length(modelNames))

for(i in seq_along(modelNames)){
	freeParams[i] <- 7 - length(fixP[[i]])
	sol <- getModelResults(fixP[[i]])
	table2opt[i,] <- sol$par
	logLik[i] <- -sol$value
	if(i==1){
		optLL <- logLik[1]
	}else{
		pVal[i] <- getPval(logLik[i], optLL, length(fixP[[i]]))
	}
	aic[i] <- getAIC(logLik[i], freeParams[i])
}

colnames(table2opt) <- names(paramInit)

deltaAIC <- aic-aic[1]

table2 <- cbind(table2opt,
		    logLik = logLik,
		    Pvalue = pVal,
		    freeParams = freeParams,
		    deltaAIC = deltaAIC)

rownames(table2) <- modelNames

print(table2)
