rm(list=ls())

pHouse <- 0.35342938
dHouse <- 0.54052428

dbetabinom <- function (y, size, p, d){
    a <- d * p
    b <- d * (1 - p)
    exp(lbeta(y + a, size - y + b) - lbeta(a, b) + lchoose(size, y))
}


num <- c(1376,3223,1385,1228,887,567,258,87,50,17,7,5,1)

N <- 13
i <- 1:N
h <- num/sum(num)
mu <- sum(i*h)
ci <- i*h/mu

sig2 <- sum(i^2*h)-mu^2
meanC <- mu + sig2/mu - 1

Rh <- pHouse * meanC

getRhStar <- function(pHouse,dHouse){
	M <- 12
	M2 <- M^2
	if(dHouse == 0){
		RhStar <- pHouse * meanC
	}else{
		ph <- array(0,c(M,M+1,M))

		if(dHouse == Inf){ for(j in 0:(M-1)) ph[1,j+1,(j+1):M] <- dbinom(j,(j+1):M,pHouse)
		}else			 for(j in 0:(M-1)) ph[1,j+1,(j+1):M] <- dbetabinom(j,(j+1):M,pHouse,dHouse)

		for(i in 2:M)
			for(j in 0:(M-i)) 
				for(k in (j+1):(M-i+1)) 
						ph[i,j+1,k] <- sum(ph[i-1,1:(j+1),k] * ph[1 + M*(j:0) + M2*((k-1):(k-j-1))])
		pt <- ph

		for(j in 1:(M-1))
			for(i in 1:(M-j))
				for(k in (j+1):(M-i+1))
					pt[i,j+1,k] <- sum(ph[i,(j+1):2,k] * pt[j:1 + (0:(j-1))*M + ((k-j-1):(k-2))*M2])
	
		for(k in 1:M) pt[1,k+1,k] <- 1-sum(pt[1,1:k,k])

		RhTot <- rep(0,N)

		for(i in 2:N)
			RhTot[i] <- sum((1:(i-1))*pt[1,2:i,i-1])
		
		RhStar <- sum(ci*RhTot)
	}
	RhStar
}


RhStar <- getRhStar(pHouse,dHouse)

boots <- read.table('bootstrapResults.txt',header=TRUE)

RhBoot <- boots$pHouse * meanC
RhStarBoot <- Vectorize(getRhStar)(boots$pHouse,boots$dHouse)

getBndry = function(boot,opt,p){
	f <- splinefun(1:length(boot),sort(boot),method='hyman')
	z0 <- qnorm(sum(boot>opt)/length(boot))
	zetalow <- pnorm(2*z0+qnorm((1-p)/2))
	zetahigh <- pnorm(2*z0+qnorm((1+p)/2))
	c(f(zetalow*length(boot)),f(zetahigh*length(boot)))
}

RhCI <- getBndry(RhBoot,Rh,0.95)
RhStarCI <- getBndry(RhStarBoot,RhStar,0.95)

RcThresh <- 1/(1+RhStar)
RcThreshCI <- rev(1/(1+RhStarCI))

RhResults <- rbind(c(Rh,RhCI),c(RhStar,RhStarCI),c(RcThresh,RcThreshCI))
rownames(RhResults) <- c('Rh','RhStar','RcThreshold')
colnames(RhResults) <- c('mle','low','high')

print(RhResults)
