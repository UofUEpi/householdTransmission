rm(list=ls())

source('table1MLE.R')

llrs <- Vectorize(llRatioStat)

plotVals <- rbind(pCom = seq(0.0026,0.0056,len=30),
			pHouse = seq(0.2,0.55,len=30),
			dHouse = exp(seq(0,log(5),len=30))-1,
			phiV = seq(0.57,0.87,len=30),
			phiA = seq(0.71,0.96,len=30),
			piA = seq(0.9909,0.995,len=30),
			piV = seq(0.9984,1,len=30))

plotLim <- rbind(pCom = c(0.0026,0.0056),
		     pHouse = c(0.2,0.55),
		     dHouse = c(0,4),
		     phiV = c(0.57,0.87),
		     phiA = c(0.71,0.96),
		     piA = c(0.9909,0.995),
		     piV = c(0.9985,1))

lab <- c(pCom = expression(italic(p[c])),
	   pHouse = expression(italic(p[h])),
	   dHouse = expression(italic(d[h])),
	   phiV = expression(italic(phi[V])),
	   phiA = expression(italic(phi[A])),
	   piA = expression(italic(pi[A])),
	   piV = expression(italic(pi[V])))

getZ <- function(xName,yName,f) outer(plotVals[xName,],plotVals[yName,], f)

z_pCpH <- getZ('pCom', 'pHouse', function(pC,pH) llrs(pCom=pC,pHouse=pH))
z_pCdH <- getZ('pCom', 'dHouse', function(pC,dH) llrs(pCom=pC,dHouse=dH))
z_pCphiV <- getZ('pCom', 'phiV', function(pC,phiV) llrs(pCom=pC,phiV=phiV))
z_pCphiA <- getZ('pCom', 'phiA', function(pC,phiA) llrs(pCom=pC,phiA=phiA))
z_pHdH <- getZ('pHouse', 'dHouse', function(pH,dH) llrs(pHouse=pH,dHouse=dH))
z_pHphiV <- getZ('pHouse', 'phiV', function(pH,phiV) llrs(pHouse=pH,phiV=phiV))
z_pHphiA <- getZ('pHouse', 'phiA', function(pH,phiA) llrs(pHouse=pH,phiA=phiA))
z_dHphiV <- getZ('dHouse', 'phiV', function(dH,phiV) llrs(dHouse=dH,phiV=phiV))
z_dHphiA <- getZ('dHouse', 'phiA', function(dH,phiA) llrs(dHouse=dH,phiA=phiA))
z_phiVphiA <- getZ('phiV', 'phiA', function(phiV,phiA) llrs(phiV=phiV,phiA=phiA))
z_pCpiA <- getZ('pCom', 'piA', function(pC,piA) llrs(pCom=pC, piA=piA))
z_pHpiA <- getZ('pHouse', 'piA', function(pH,piA) llrs(pHouse=pH, piA=piA))
z_dHpiA <- getZ('dHouse', 'piA', function(dH,piA) llrs(dHouse=dH, piA=piA))
z_phiVpiA <- getZ('phiV', 'piA', function(phiV,piA) llrs(phiV=phiV, piA=piA))
z_phiApiA <- getZ('phiA', 'piA', function(phiA,piA) llrs(phiA=phiA, piA=piA))
z_pCpiV <- getZ('pCom', 'piV', function(pC,piV) llrs(pCom=pC, piV=piV))
z_pHpiV <- getZ('pHouse', 'piV', function(pH,piV) llrs(pHouse=pH, piV=piV))
z_dHpiV <- getZ('dHouse', 'piV', function(dH,piV) llrs(dHouse=dH, piV=piV))
z_phiVpiV <- getZ('phiV', 'piV', function(phiV,piV) llrs(phiV=phiV, piV=piV))
z_phiApiV <- getZ('phiA', 'piV', function(phiA,piV) llrs(phiA=phiA, piV=piV))
z_piApiV <- getZ('piA', 'piV', function(piA,piV) llrs(piA=piA, piV=piV))

boots <- read.table('bootstrapResults.txt',header=TRUE)
bootsAdj <- boots
bootsAdj[boots[,'dHouse']==Inf,'dHouse'] <- 1000

plotContour <- function(xName,yName,z,qlevel){
	contour(plotVals[xName,],plotVals[yName,],z,levels=qlevel,drawlabels=FALSE,
		  xlab=lab[xName],ylab=lab[yName],xlim=plotLim[xName,],ylim=plotLim[yName,])
	points(optAll[xName],optAll[yName],pch=19)
	lines(rep(optAll[xName],2), cis[yName,], lty=2)
	lines(cis[xName,], rep(optAll[yName],2), lty=2)
	points(bootsAdj[,xName],bootsAdj[,yName],pch='.')
}

q2 <- qchisq(pCI,2)

dev.new(width=7,height=5)
par(mfrow=c(2,3),mar=c(5,4,2,2)+0.1)
plotContour('pCom','pHouse',z_pCpH,q2)
plotContour('pCom','dHouse',z_pCdH,q2)
plotContour('pCom','phiV',z_pCphiV,q2)
plotContour('pCom','phiA',z_pCphiA,q2)
plotContour('pCom','piA',z_pCpiA,q2)
plotContour('pCom','piV',z_pCpiV,q2)

dev.new(width=7,height=5)
par(mfrow=c(2,3),mar=c(5,4,2,2)+0.1)
plotContour('pHouse','pCom',t(z_pCpH),q2)
plotContour('pHouse','dHouse',z_pHdH,q2)
plotContour('pHouse','phiV',z_pHphiV,q2)
plotContour('pHouse','phiA',z_pHphiA,q2)
plotContour('pHouse','piA',z_pHpiA,q2)
plotContour('pHouse','piV',z_pHpiV,q2)

dev.new(width=7,height=5)
par(mfrow=c(2,3),mar=c(5,4,2,2)+0.1)
plotContour('dHouse','pCom',t(z_pCdH),q2)
plotContour('dHouse','pHouse',t(z_pHdH),q2)
plotContour('dHouse','phiV',z_dHphiV,q2)
plotContour('dHouse','phiA',z_dHphiA,q2)
plotContour('dHouse','piA',z_dHpiA,q2)
plotContour('dHouse','piV',z_dHpiV,q2)

dev.new(width=7,height=5)
par(mfrow=c(2,3),mar=c(5,4,2,2)+0.1)
plotContour('phiV','pCom',t(z_pCphiV),q2)
plotContour('phiV','pHouse',t(z_pHphiV),q2)
plotContour('phiV','dHouse',t(z_dHphiV),q2)
plotContour('phiV','phiA',z_phiVphiA,q2)
plotContour('phiV','piA',z_phiVpiA,q2)
plotContour('phiV','piV',z_phiVpiV,q2)

dev.new(width=7,height=5)
par(mfrow=c(2,3),mar=c(5,4,2,2)+0.1)
plotContour('phiA','pCom',t(z_pCphiA),q2)
plotContour('phiA','pHouse',t(z_pHphiA),q2)
plotContour('phiA','dHouse',t(z_dHphiA),q2)
plotContour('phiA','phiV',t(z_phiVphiA),q2)
plotContour('phiA','piA',z_phiApiA,q2)
plotContour('phiA','piV',z_phiApiV,q2)

dev.new(width=7,height=5)
par(mfrow=c(2,3),mar=c(5,4,2,2)+0.1)
plotContour('piA','pCom',t(z_pCpiA),q2)
plotContour('piA','pHouse',t(z_pHpiA),q2)
plotContour('piA','dHouse',t(z_dHpiA),q2)
plotContour('piA','phiV',t(z_phiVpiA),q2)
plotContour('piA','phiA',t(z_phiApiA),q2)
plotContour('piA','piV',z_piApiV,q2)

dev.new(width=7,height=5)
par(mfrow=c(2,3),mar=c(5,4,2,2)+0.1)
plotContour('piV','pCom',t(z_pCpiV),q2)
plotContour('piV','pHouse',t(z_pHpiV),q2)
plotContour('piV','dHouse',t(z_dHpiV),q2)
plotContour('piV','phiV',t(z_phiVpiV),q2)
plotContour('piV','phiA',t(z_phiApiV),q2)
plotContour('piV','piA',t(z_piApiV),q2)





