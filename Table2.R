rm(list=ls())

dataFileName <- 'analysisData.txt'
source('prepareData.R')
source('MLEfunctions.R')

getSolDcInf <- function(piVfix) getSol(c(pCom=0.01,pHouse=0.1,dHouse=1,phiV=0.3,phiA=0.8,piA=0.995),
							      function(x) -getLL(x[1],Inf,x[2],x[3],x[4],x[5],piVfix,x[6]))

piVsolDcInf <- optimize(function(x) getSolDcInf(x)$value,c(0.999,1),tol=1e-10)$minimum
solDcInf <- getSolDcInf(piVsolDcInf)
solDcInf$par <- c(solDcInf$par[1:5],piV=piVsolDcInf,solDcInf$par[6])

getSolDcInfDh0 <- function(piV) getSol(c(pCom=0.01,pHouse=0.1,phiV=0.3,phiA=0.8,piA=0.995),
				    function(x) -getLL(x[1],Inf,x[2],0,x[3],x[4],piV,x[5]))

piVsolDcInfDh0 <- optimize(function(x) getSolDcInfDh0(x)$value,c(0.999,1),tol=1e-10)
solDcInfDh0 <- getSolDcInfDh0(piVsolDcInfDh0$minimum)
solDcInfDh0$par <- c(solDcInfDh0$par[1:4],piV=piVsolDcInfDh0$minimum,solDcInfDh0$par[5])


solDcInfPiV1 <- getSol(c(pCom=0.01,pHouse=0.1,dHouse=1,phiV=0.3,phiA=0.8,piA=0.995), 
			     function(x) -getLL(x[1],Inf,x[2],x[3],x[4],x[5],1,x[6]))

solDcInfDhInf <- getSol(c(pCom=0.01,pHouse=0.1,phiV=0.3,phiA=0.8,piV=0.995,piA=0.995), 
				function(x) -getLL(x[1],Inf,x[2],Inf,x[3],x[4],x[5],x[6]))

solDcInfDh0PiV1 <- getSol(c(pCom=0.01,pHouse=0.1,phiV=0.3,phiA=0.8,piA=0.995),
				  function(x) -getLL(x[1],Inf,x[2],0,x[3],x[4],1,x[5]))

solDcInfPiA1PiV1 <- getSol(c(pCom=0.01,pHouse=0.1,dHouse=1,phiV=0.3,phiA=0.8), 
			         function(x) -getLL(x[1],Inf,x[2],x[3],x[4],x[5],1,1))

getSolDcInfPiA1 <- function(piV) getSol(c(pCom=0.01,pHouse=0.1,dHouse=1,phiV=0.3,phiA=0.8),
			    function(x) -getLL(x[1],Inf,x[2],x[3],x[4],x[5],piV,1))

#Applying getSolDcInfPiA1 for any piV<1 produces lower likelihood than piV=1, so: 
solDcInfPiA1 <- solDcInfPiA1PiV1
solDcInfPiA1$par <- c(solDcInfPiA1$par,piV=1)

getSolDcInfPiA996 <- function(piV) getSol(c(pCom=0.01,pHouse=0.1,dHouse=1,phiV=0.3,phiA=0.8),
			    function(x) -getLL(x[1],Inf,x[2],x[3],x[4],x[5],piV,0.996))

piVsolDcInfPiA996 <- optimize(function(x) getSolDcInfPiA996(x)$value,c(0.999,1))
solDcInfPiA996 <- getSolDcInfPiA996(piVsolDcInfPiA996$minimum)
solDcInfPiA996$par <- c(solDcInfPiA996$par,piV=piVsolDcInfPiA996$minimum)

PvalDcInfDh0 <- (1-pchisq(2*(solDcInfDh0$value-solDcInf$value), df=1))
PvalDcInfDhInf <- (1-pchisq(2*(solDcInfDhInf$value-solDcInf$value), df=1))
PvalDcInfPiV1 <- (1-pchisq(2*(solDcInfPiV1$value-solDcInf$value),df=1))
PvalDcInfPiA1 <- (1-pchisq(2*(solDcInfPiA1$value-solDcInf$value),df=1))
PvalDcInfPiA996 <- (1-pchisq(2*(solDcInfPiA996$value-solDcInf$value),df=1))
PvalDcInfPiA1PiV1 <- (1-pchisq(2*(solDcInfPiA1PiV1$value-solDcInf$value), df=2))

table2opt <- rbind(solDcInf$par,
			 c(solDcInfDh0$par[1:2],dHouse=0,solDcInfDh0$par[3:6]),
			 c(solDcInfDhInf$par[1:2],dHouse=Inf,solDcInfDhInf$par[3:6]),
			 c(solDcInfPiV1$par[1:5],piV=1,solDcInfPiV1$par[6]),
			 c(solDcInfPiA996$par,piA=0.996),
			 c(solDcInfPiA1$par,piA=1))

table2 <- cbind(table2opt,
		    logLik = -c(solDcInf$value,solDcInfDh0$value,solDcInfDhInf$value,
				    solDcInfPiV1$value,solDcInfPiA996$value,solDcInfPiA1$value),
		    pValue = c(NA,PvalDcInfDh0,PvalDcInfDhInf,PvalDcInfPiV1,
				   PvalDcInfPiA996,PvalDcInfPiA1))

rownames(table2) <- c('none','dH=0','dH=Inf','piV=1','piA=0.996','piA=1')

print(table2)