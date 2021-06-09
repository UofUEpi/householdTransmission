rm(list=ls())

pH <- 0.363258707
dH <- 0.433838990

getk <- function(pH,dH){
	ifelse(dH==0, 0, ifelse(dH==Inf,Inf,
	uniroot(function(k) ((1-pH)*dH+1)^(1/k)*(2-(1-pH)^(1/k))-(dH+1)^(1/k) ,c(dH/10,dH))$root))
}

k <- getk(pH,dH)

boot <- read.table('bootstrapResults.txt',header=TRUE)

kboot <- Vectorize(getk)(boot$pHouse,boot$dHouse)

getBndry = function(boot,opt,p){
	f <- splinefun(1:length(boot),sort(boot),method='hyman')
	z0 <- qnorm(sum(boot>opt)/length(boot))
	zetalow <- pnorm(2*z0+qnorm((1-p)/2))
	zetahigh <- pnorm(2*z0+qnorm((1+p)/2))
	c(f(zetalow*length(boot)),f(zetahigh*length(boot)))
}

kbootAdj <- 1-exp(-kboot)
kAdj <- 1-exp(-k)

kci <- -log(1-getBndry(kbootAdj,kAdj,0.95))

gammaDispersionResults <- c(mle = k, low = kci[1], high = kci[2])

print(gammaDispersionResults)




