rm(list=ls())

optPar <- c(pCom = 0.003851334, pHouse = 0.354171115, dHouse = 0.535584255,
	      phiV = 0.762914827, phiA = 0.852187169, piV = 0.999311905, piA = 0.992963410)

boots <- read.table('bootstrapResults.txt',header=TRUE)

bootsAdj <- boots
bootsAdj[,'dHouse'] <- 1-exp(-boots[,'dHouse'])

optParAdj <- optPar
optParAdj['dHouse'] <- 1-exp(-optPar['dHouse'])

getBndry = function(boot,opt,p){
	f <- splinefun(1:length(boot),sort(boot),method='hyman')
	z0 <- qnorm(sum(boot>opt)/length(boot))
	zetalow <- pnorm(2*z0+qnorm((1-p)/2))
	zetahigh <- pnorm(2*z0+qnorm((1+p)/2))
	c(f(zetalow*length(boot)),f(zetahigh*length(boot)))
}

cis <- matrix(0,length(optPar),2)
rownames(cis) <- names(optPar)
colnames(cis) <- c('low','high')

for(col in names(optPar)) cis[col,] <- getBndry(bootsAdj[,col],optParAdj[col],0.95)

cis['dHouse',] <- -log(1-cis['dHouse',])

Table1Bootstrap <- cbind(median = apply(boots,2,median),cis)
 
print(Table1Bootstrap)