rm(list=ls())

#Results from Table 1:
pH_mle <- 0.363258707
dH_mle <- 0.433838990

dataFileName <- 'analysisData.txt'
source('prepareData.R')
source('MLEfunctions.R')

#Solve for optimum & confidence interval of pH_binom (via likelihood ratio test):
solBinom <- getSol(c(pCom=0.01,pHouse=0.1,phiV=0.3,phiA=0.8,piV=0.995,piA=0.995), 
					function(x) -getLL(x[1],Inf,x[2],Inf,x[3],x[4],x[5],x[6]))
optBinom <- solBinom$par
pH_binom <- as.numeric(optBinom['pHouse'])

llRatioStatBinom <- function(pH)
	-2*(solBinom$value + getLL(optBinom['pCom'],Inf,pH,Inf,
					   optBinom['phiV'],optBinom['phiA'],optBinom['piV'],optBinom['piA']))

ciFun <- function(pH) llRatioStatBinom(pH) - qchisq(0.95,1)

pH_binomLow <- uniroot(ciFun, c(0.001,optBinom['pHouse']),tol=1e-10)$root
pH_binomHigh <- uniroot(ciFun, c(optBinom['pHouse'],0.999),tol=1e-10)$root

boots <- read.table('bootstrapResults.txt',header=TRUE)

nt <- 1:9

dbetabinom <- function (y, size, p, d){
    a <- d * p
    b <- d * (1 - p)
    exp(lbeta(y + a, size - y + b) - lbeta(a, b) + lchoose(size, y))
}

noTci <- matrix(0,length(nt),2)
allTci <- matrix(0,length(nt),2)

getBndry = function(boot,opt,p){
	f <- splinefun(1:length(boot),sort(boot),method='hyman')
	z0 <- qnorm(sum(boot>opt)/length(boot))
	zetalow <- pnorm(2*z0+qnorm((1-p)/2))
	zetahigh <- pnorm(2*z0+qnorm((1+p)/2))
	c(f(zetalow*length(boot)),f(zetahigh*length(boot)))
}

dbb <- function(k,n,pH,dH){
	out <- 0
	if(dH == 0){
		if(k == 0){
			out <- 1-pH
		}else if(k==n){
			out <- pH
		} 
	}else if(dH == Inf){
		out <- dbinom(k,n,pH)
	}else{
		out <- dbetabinom(k,n,pH,dH)
	}
	out
}
dbb <- Vectorize(dbb,c('pH','dH'))

noTransmProbMLE <- dbb(0,nt,pH_mle,dH_mle)
allTransmProbMLE <- dbb(nt,nt,pH_mle,dH_mle)

noTransmProbBinom <- dbb(0,nt,pH_binom,Inf)
allTransmProbBinom <- dbb(nt,nt,pH_binom,Inf)

noTransmProbBinomLow <- dbb(0,nt,pH_binomHigh,Inf)
noTransmProbBinomHigh <- dbb(0,nt,pH_binomLow,Inf)
allTransmProbBinomLow <- dbb(nt,nt,pH_binomLow,Inf)
allTransmProbBinomHigh <- dbb(nt,nt,pH_binomHigh,Inf)

for(nti in nt){
	noTransmBoot <- dbb(0,nti,boots[,'pHouse'],boots[,'dHouse'])
	allTransmBoot <- dbb(nti,nti,boots[,'pHouse'],boots[,'dHouse'])
	noTci[nti,] <- getBndry(noTransmBoot,noTransmProbMLE[nti],0.95)
	allTci[nti,] <- getBndry(allTransmBoot,allTransmProbMLE[nti],0.95)
}

noTransmBootTable <- data.frame(mle = noTransmProbMLE,
			         mleLow = noTci[,1],
			         mleHigh = noTci[,2],
			         binom = noTransmProbBinom,
				   binomLow = noTransmProbBinomLow,
				   binomHigh = noTransmProbBinomHigh)

allTransmBootTable <- data.frame(mle = allTransmProbMLE,
			          mleLow = allTci[,1],
			          mleHigh = allTci[,2],
			          binom = allTransmProbBinom,
				    binomLow = allTransmProbBinomLow,
				    binomHigh = allTransmProbBinomHigh)

rownames(noTransmBootTable) <- 2:10
rownames(allTransmBootTable) <- 2:10

Table3 <- list(transmitToNone = 100*noTransmBootTable,
		   transmitToAll = 100*allTransmBootTable)

print(Table3)
