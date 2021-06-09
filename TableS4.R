rm(list=ls())

dataFileName <- 'analysisData.txt'
source('prepareData.R')
source('MLEfunctions.R')

getSolAll <- function(dCom,dHouse) getSol(c(pCom=0.01,pHouse=0.35,phiV=0.8,phiA=0.8,piV=0.995,piA=0.995),
			function(x) -getLL(x[1],dCom,x[2],dHouse,x[3],x[4],x[5],x[6]))

dCdHsol <- optim(fn = function(x) getSolAll(x[1],x[2])$value, par=c(2.5,0.5))

solAll <- getSolAll(dCdHsol$par[1],dCdHsol$par[2]) 
solAll$par <- c(solAll$par[1],dCom=dCdHsol$par[1],solAll$par[2],dHouse=dCdHsol$par[2],solAll$par[3:6])

optAll <- solAll$par

llRatioStat <- function(pCom=optAll['pCom'], dCom=optAll['dCom'], pHouse=optAll['pHouse'],
				dHouse=optAll['dHouse'], phiV=optAll['phiV'], phiA=optAll['phiA'], 
				piV=optAll['piV'], piA=optAll['piA'])
	-2*(solAll$value + getLL(pCom,dCom,pHouse,dHouse,phiV,phiA,piV,piA))

pCI <- 0.95
	
ciFun_pCom <- function(pC) llRatioStat(pCom=pC) - qchisq(pCI,1)
ciFun_pHouse <- function(pH) llRatioStat(pHouse=pH) - qchisq(pCI,1)
ciFun_phiV <- function(pV) llRatioStat(phiV=pV) - qchisq(pCI,1)
ciFun_phiA <- function(pA) llRatioStat(phiA=pA) - qchisq(pCI,1)
ciFun_dHouse <- function(dH) llRatioStat(dHouse=dH) - qchisq(pCI,1)
ciFun_dCom <- function(dC) llRatioStat(dCom=dC) - qchisq(pCI,1)
ciFun_piV <- function(piV) llRatioStat(piV=piV) - qchisq(pCI,1)
ciFun_piA <- function(piA) llRatioStat(piA=piA) - qchisq(pCI,1)

getCI <- function(ciFun, parName, interval)
	c(uniroot(ciFun, c(interval[1],optAll[parName]),tol=1e-10)$root,
	  uniroot(ciFun, c(optAll[parName],interval[2]),tol=1e-10)$root)

getCIoneSide <- function(ciFun,interval) uniroot(ciFun, interval, tol=1e-10)$root

pHouseCI <- getCI(ciFun_pHouse,'pHouse',c(0.001,0.999))
pComCI <- getCI(ciFun_pCom,'pCom',c(0.001,0.999))
phiVCI <- getCI(ciFun_phiV,'phiV',c(0.001,0.999))
phiACI <- getCI(ciFun_phiA,'phiA',c(0.001,0.999))
dHouseCI <- c(0,getCIoneSide(ciFun_dHouse,c(0.001,1000)))
dComCI <- getCI(ciFun_dCom,'dCom',c(0,1000))
piACI <- getCI(ciFun_piA,'piA',c(0.5,1))
piVCI <- getCI(ciFun_piV,'piV',c(0.5,1))

cis <- rbind(pCom = pComCI, dCom = dComCI, pHouse = pHouseCI, dHouse = dHouseCI, 
	 	 phiV = phiVCI, phiA = phiACI, piV = piVCI,piA = piACI)

TableS4 <- cbind(mle = optAll, low = cis[,1], high = cis[,2])

print(TableS4)
