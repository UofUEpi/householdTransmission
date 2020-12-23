cf <- read.table(dataFileName, header=TRUE)

vseq <- Vectorize(seq.default, c('from', 'to'))
vrep <- Vectorize(rep.int)

maxSize <- max(cf$n)

cf_aNN <- cf$a-cf$aPP-cf$aPN-cf$aNP
cf_sN <- cf$s-cf$sP
cf_q <- cf$n-cf$a-cf$s

ul <- cf$aPP + 1
vl <- cf$aPN + 1
wl <- cf$aNP + 1
xl <- cf_aNN + 1
yl <- cf$sP + 1
zl <- cf_sN + 1
kl <- cf_q + 1 
numTerms <- ul*vl*wl*xl*yl*zl*kl

u <- rep(unlist(vseq(0,cf$aPP)),rep(vl*wl*xl*yl*zl*kl,ul))
v <- rep(unlist(vrep(vseq(0,cf$aPN),ul)),rep(wl*xl*yl*zl*kl,ul*vl))
w <- rep(unlist(vrep(vseq(0,cf$aNP),ul*vl)),rep(xl*yl*zl*kl,ul*vl*wl))
x <- rep(unlist(vrep(vseq(0,cf_aNN),ul*vl*wl)),rep(yl*zl*kl,ul*vl*wl*xl))
y <- rep(unlist(vrep(vseq(0,cf$sP),ul*vl*wl*xl)),rep(zl*kl,ul*vl*wl*xl*yl))
z <- rep(unlist(vrep(vseq(0,cf_sN),ul*vl*wl*xl*yl)),rep(kl,ul*vl*wl*xl*yl*zl))
Ia <- u+v+w+x
Is <- y+z
k <- Ia+Is + unlist(vrep(vseq(0,cf_q),ul*vl*wl*xl*yl*zl))

stretch <- function(x) rep(x,numTerms)
n <- stretch(cf$n)
a <- stretch(cf$a)
s <- stretch(cf$s)
sP <- stretch(cf$sP)
aPP <- stretch(cf$aPP)
aPN <- stretch(cf$aPN)
aNP <- stretch(cf$aNP)
aNN <- stretch(cf_aNN)
sN <- stretch(cf_sN)
Mi <- k+1+(n-1)*(maxSize+1)

#true positive (tp) / false negative (fn) / true negative (tn) / false positive (fp) exponents
tpV <- u+v+y
fnV <- w+x+z
tpA <- u+w
fnA <- v+x
tnV <- aNP+aNN+sN-fnV
fpV <- aPP+aPN+sP-tpV
tnA <- aPN+aNN-fnA
fpA <- aPP+aNP-tpA 

AinfCoef <- factorial(Ia)/factorial(u)/factorial(v)/factorial(w)/factorial(x)
AuninfCoef <- factorial(a-Ia)/factorial(aPP-u)/factorial(aPN-v)/factorial(aNP-w)/factorial(aNN-x)
Bcoef <- choose(Is,y) * choose(s-Is,sN-z)
Hcoef <- dhyper(Ia,k,n-k,a) * dhyper(Is,k-Ia,n-k-a+Ia,s) 

coef <- AinfCoef*AuninfCoef*Bcoef*Hcoef

