rm(list=ls())

optPar <- c(pCom = 0.004078679, pHouse = 0.363258707, dHouse = 0.433838990,
		phiV = 0.724447576, phiA = 0.855594551, piV = 0.999356963, piA = 0.993141399)

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
