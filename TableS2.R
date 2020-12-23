rm(list=ls())

dataFileName <- 'analysisData.txt'
source('prepareData.R')
source('MLEfunctions.R')

optAll <- c(pCom = 0.005198326, dCom = 2.353515625, pHouse = 0.277844661,
		dHouse = 0.300781250, phiV = 0.758066320, phiA = 0.864529682,
		piV = 0.999216578, piA = 0.993014304) 

optDcInf <- c(pCom = 0.003851334, dCom = Inf, pHouse = 0.354171115,
		  dHouse = 0.535584255, phiV = 0.762914827, phiA = 0.852187169,
		  piV = 0.999311905, piA = 0.992963410) 

solAllValue <- -do.call(getLL,as.list(optAll))
solDcInfValue <- -do.call(getLL,as.list(optDcInf))

getSolDh0 <- function(pHfix) getSol(c(pCom=0.01,dCom=1,phiV=0.3,phiA=0.8,piV=0.999,piA=0.995), 
						function(x) -getLL(x[1],x[2],pHfix,0,x[3],x[4],x[5],x[6]))
pHsolDh0 <- optimize(function(x) getSolDh0(x)$value,c(0,0.5), tol=1e-10)
solDh0 <- getSolDh0(pHsolDh0$minimum)
solDh0$par <- c(solDh0$par[1:2],pHouse = pHsolDh0$minimum,solDh0$par[3:6])

getSolDhInf <- function(pHfix) getSol(c(pCom=0.01,dCom=1,phiV=0.3,phiA=0.8,piV=0.999,piA=0.995), 
			             	  function(x) -getLL(x[1],x[2],pHfix,Inf,x[3],x[4],x[5],x[6]))
pHsolDhInf <- optimize(function(x) getSolDhInf(x)$value,c(0,0.1),tol=1e-10)
solDhInf <- getSolDhInf(pHsolDhInf$minimum)
solDhInf$par <- c(solDhInf$par[1:2],pHouse = pHsolDhInf$minimum,solDhInf$par[3:6])

solPh0 <- getSolDhInf(0)

if(solPh0$value < solDhInf$value){
	solDhInf <- solPh0
	solDhInf$par <- c(solDhInf$par[1:2],pHouse = 0,solDhInf$par[3:6])
}

PvalDcInf <- (1-pchisq(2*(solDcInfValue-solAllValue), df=1))	
PvalDh0 <- (1-pchisq(2*(solDh0$value-solAllValue), df=1))
PvalDhInf <- (1-pchisq(2*(solDhInf$value-solAllValue), df=1))
PvalPh0 <- (1-pchisq(2*(solPh0$value-solAllValue), df=2))

TableS2opt <- rbind(optAll,
			  optDcInf,
			  c(solDh0$par[1:3],dHouse=0,solDh0$par[4:7]),
			  c(solDhInf$par[1:3],dHouse=Inf,solDhInf$par[4:7]),
			  c(solPh0$par[1:2],pHouse=0,dHouse=NA,solPh0$par[3:6]))

TableS2 <- cbind(TableS2opt,
	     	     logLik = -c(solAllValue,solDcInfValue,solDh0$value,solDhInf$value,solPh0$value),
		     pValue = c(NA,PvalDcInf,PvalDh0,PvalDhInf,PvalPh0))

rownames(TableS2) <- c('none','dC=Inf','dH=0','dH=Inf','pH=0')

print(TableS2)
