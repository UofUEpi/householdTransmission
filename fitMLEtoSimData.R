rm(list=ls())

numSims <- 500

source('MLEfunctions.R')

getSolDcInf <- function(piVfix) getSol(c(pCom=0.01,pHouse=0.1,dHouse=1,phiV=0.3,phiA=0.8,piA=0.995),
							      function(x) -getLL(x[1],Inf,x[2],x[3],x[4],x[5],piVfix,x[6]))

getSolDcInfDh0 <- function(piV) getSol(c(pCom=0.01,pHouse=0.1,phiV=0.3,phiA=0.8,piA=0.995),
				    function(x) -getLL(x[1],Inf,x[2],0,x[3],x[4],piV,x[5]))

getSolDcInfDhInf <- function(piV) getSol(c(pCom=0.01,pHouse=0.1,phiV=0.3,phiA=0.8,piA=0.995),
				    function(x) -getLL(x[1],Inf,x[2],Inf,x[3],x[4],piV,x[5]))

solBoots <- matrix(0,numSims,7)

for(b in 1:numSims){
	dataFileName <- paste0('simData/simData',b,'.txt')
	source('prepareData.R')
	piVsolBoot <- optimize(function(x) getSolDcInf(x)$value,c(0.99,1))$minimum
	solBoot <- getSolDcInf(piVsolBoot)
	solBootPar <- c(solBoot$par[1:5],piVsolBoot,solBoot$par[6])

	#Check if optimum occurs at a boundary of the dHouse parameter
	if(solBoot$par['dHouse'] < 0.1){
		piVsolBootDh0 <- optimize(function(x) getSolDcInfDh0(x)$value,c(0.99,1))
		if(piVsolBootDh0$objective < solBoot$value){
			solBoot <- getSolDcInfDh0(piVsolBootDh0$minimum)
			solBootPar <- c(solBoot$par[1:2],0,solBoot$par[3:4],piVsolBootDh0$minimum,solBoot$par[5])
		}
	}else if(solBoot$par['dHouse'] > 10){
		piVsolBootDhInf <- optimize(function(x) getSolDcInfDhInf(x)$value,c(0.99,1))
		if(piVsolBootDhInf$objective < solBoot$value){
			solBoot <- getSolDcInfDhInf(piVsolBootDhInf$minimum)
			solBootPar <- c(solBoot$par[1:2],Inf,solBoot$par[3:4],piVsolBootDhInf$minimum,solBoot$par[5])
		}
	}
	solBoots[b,] <- solBootPar
}

colnames(solBoots) <- c('pCom','pHouse','dHouse','phiV','phiA','piV','piA')
 
write.table(solBoots,'bootstrapResults.txt',quote=FALSE,row.names=FALSE)
