rm(list=ls())

dataFileName <- 'analysisData.txt'
source('prepareData.R')
source('MLEfunctions.R')

getSolDcInf <- function(piVfix) getSol(c(pCom=0.01,pHouse=0.1,dHouse=1,phiV=0.3,phiA=0.8,piA=0.995),
							      function(x) -getLL(x[1],Inf,x[2],x[3],x[4],x[5],piVfix,x[6]))
piVsolDcInf <- optimize(function(x) getSolDcInf(x)$value,c(0.999,1), tol=1e-10)$minimum

solDcInf <- getSolDcInf(piVsolDcInf)
solDcInf$par <- c(solDcInf$par[1:5],piV=piVsolDcInf,solDcInf$par['piA'])

optAll <- solDcInf$par

llRatioStat <- function(pCom=optAll['pCom'], dCom=Inf, pHouse=optAll['pHouse'],
				dHouse=optAll['dHouse'], phiV=optAll['phiV'], phiA=optAll['phiA'], 
				piV=optAll['piV'], piA=optAll['piA'])
	-2*(solDcInf$value + getLL(pCom,dCom,pHouse,dHouse,phiV,phiA,piV,piA))

pCI <- 0.95
	
ciFun_pCom <- function(pC) llRatioStat(pCom=pC) - qchisq(pCI,1)
ciFun_pHouse <- function(pH) llRatioStat(pHouse=pH) - qchisq(pCI,1)
ciFun_phiV <- function(pV) llRatioStat(phiV=pV) - qchisq(pCI,1)
ciFun_phiA <- function(pA) llRatioStat(phiA=pA) - qchisq(pCI,1)
ciFun_dHouse <- function(dH) llRatioStat(dHouse=dH) - qchisq(pCI,1)
ciFun_piV <- function(piV) llRatioStat(piV=piV) - qchisq(pCI,1)
ciFun_piA <- function(piA) llRatioStat(piA=piA) - qchisq(pCI,1)

getCI <- function(ciFun, parName, interval)
	c(uniroot(ciFun, c(interval[1],optAll[parName]),tol=1e-10)$root,
	  uniroot(ciFun, c(optAll[parName],interval[2]),tol=1e-10)$root)

pHouseCI <- getCI(ciFun_pHouse,'pHouse',c(0.001,0.999))
pComCI <- getCI(ciFun_pCom,'pCom',c(0.001,0.999))
phiVCI <- getCI(ciFun_phiV,'phiV',c(0.001,0.999))
phiACI <- getCI(ciFun_phiA,'phiA',c(0.001,0.999))
dHouseCI <- getCI(ciFun_dHouse,'dHouse',c(0,1000))
piACI <- getCI(ciFun_piA,'piA',c(0.5,1))
piVCI <- getCI(ciFun_piV,'piV',c(0.5,1))

cis <- rbind(pCom = pComCI, pHouse = pHouseCI, dHouse = dHouseCI, 
	 	 phiV = phiVCI, phiA = phiACI, piV = piVCI, piA = piACI)

table1MLE <- data.frame(mle = optAll, ciLow = cis[,1], ciHigh = cis[,2])

print(table1MLE)
