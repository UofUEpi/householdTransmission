rm(list=ls())

d <- read.table('analysisData.txt',header=TRUE)
dnas <- aggregate(list(freq=d$freq),list(n=d$n,a=d$a,s=d$s),sum)

getFinalSize <- function(n,init,a,b){
	uninf <- n - init
	inf <- init
	while(inf > 0){
		p <- rbeta(inf,a,b)
		newInf <- 0
		for(i in 1:inf){
			transm <- rbinom(1,uninf,p[i])
			uninf <- uninf - transm
			newInf <- newInf + transm
		}
		inf <- newInf
	}
	n - uninf
}

getFinalSize <- Vectorize(getFinalSize,'init')

getSimNAS <- function(pC,pH,dH){
	a <- dH*pH
	b <- dH*(1-pH)

	out <- {}
	for(row in 1:nrow(dnas)){
		inits <- rbinom(dnas$freq[row],dnas$n[row],pC)	
		fs <- getFinalSize(dnas$n[row],inits,a,b)
		nOut <- cbind(dnas$n[row],dnas$a[row],dnas$s[row],0:max(fs),tabulate(fs+1))
		out <- rbind(out,nOut)	
	}

	out <- out[out[,5]>0,]
	data.frame(n = out[,1], a = out[,2], s = out[,3], k = out[,4], freq = out[,5])
}

rmultinomial <- function (n, size, prob) {
    if (n == 1) return(t(rmultinom(n, size, prob)))
    prob <- matrix(prob, nrow = 1)
    ss <- rep(size, length.out = n)
    prob <- matrix(t(prob), ncol = ncol(prob), nrow = n, byrow = TRUE)
    res <- sapply(1:n, function(x) rmultinom(1, ss[x], prob[x,]))
    res <- matrix(res, ncol = ncol(prob), byrow = TRUE)
    return(res)
}

getSimTestNAS <- function(pC,pH,dH,phiV,phiA,piV,piA){

	pI <- c(phiV*phiA, phiV*(1-phiA), (1-phiV)*phiA, (1-phiV)*(1-phiA))
	pU <- c((1-piV)*(1-piA), (1-piV)*piA, piV*(1-piA), piV*piA)

	sNAS <- getSimNAS(pC,pH,dH)

	out <- {}
	for(row in 1:nrow(sNAS)){
		n <- sNAS$n[row]; a <- sNAS$a[row]; s <- sNAS$s[row]; k <- sNAS$k[row]; fr <- sNAS$freq[row]
		aI <- rhyper(fr,k,n-k,a)
		aU <- a-aI
		sI <- rhyper(fr,k-aI,n-k-aU,s)
		sU <- s-sI

		aTst <- rmultinomial(fr, aI, pI) + rmultinomial(fr, aU, pU)
		sP <- rbinom(fr,sI,phiV) + rbinom(fr,sU,1-piV)
		outRows <- cbind(n,a,s,aTst,sP,1)
		out <- rbind(out,outRows)			
	}
	df <- data.frame(n = out[,1], a = out[,2], s = out[,3], 
	                 aPP = out[,4], aPN = out[,5], aNP = out[,6], sP = out[,8], freq = out[,9])
	
	aggregate(list(freq=df$freq),
		    list(n=df$n,a=df$a,s=df$s,aPP=df$aPP,aPN=df$aPN,aNP=df$aNP,sP=df$sP),sum)
}

#Before running, create folder named "simData" in active folder
for(sim in 1:500){
	simData <- getSimTestNAS(0.004078679, 0.363258707, 0.433838990, 0.724447576, 0.855594551, 0.999356963, 0.993141399)	
	fileName <- paste0('simData/simData',sim,'.txt')
	write.table(simData,fileName,quote=FALSE,row.names=FALSE)
}
