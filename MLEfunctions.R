dbetabinom <- function (y, size, p, d){
    a <- d * p
    b <- d * (1 - p)
    exp(lbeta(y + a, size - y + b) - lbeta(a, b) + lchoose(size, y))
}

getProbsMaster <- function(pCom, dCom, pHouse, dHouse, maxSize){

	if(dCom == 0){
		m <- matrix(0, maxSize+1, maxSize, dimnames = list(0:maxSize, 1:maxSize))
		m[1,1:maxSize] <- 1-pCom
		m[row(m)-col(m) == 1] <- pCom
	}else if(dCom == Inf & dHouse == Inf){
		m <- matrix(0, maxSize+1, maxSize+1, dimnames = list(0:maxSize, 0:maxSize))
		m[1,1] <- 1
		for(k in 1:maxSize){
			m[1:k,k+1] <- choose(k,0:(k-1)) * diag(m)[1:k] * (1-pCom)^(k:1) * (1-pHouse)^((0:(k-1))*(k:1))
			m[k+1,k+1] <- 1 - sum(m[,k+1])
		}
		m <- m[,-1]
	}else if(dHouse == 0){
		m <- matrix(0, maxSize+1, maxSize, dimnames = list(0:maxSize, 1:maxSize))
		noTransm <- (1-pHouse)^(0:(maxSize-1))
		k <- sequence(1:maxSize)
		n <- rep(1:maxSize, 1:maxSize)
		if(dCom == Inf) db <- dbinom(k-1,n,pCom) else db <- dbetabinom(k-1,n,pCom,dCom)
		m[upper.tri(m,diag=TRUE)] <- noTransm[k] * db	
		m[row(m)-col(m) == 1] <- 1 - colSums(m) 
	}else{
		m <- matrix(0, maxSize+1, maxSize, dimnames = list(0:maxSize, 1:maxSize))
		
		if(dCom == Inf){
			m[1,1:maxSize] <- (1-pCom)^(1:maxSize)
			dbCom <- dbinom(sequence(1:maxSize),rep(1:maxSize,1:maxSize),pCom)
		}else{
			m[1,1:maxSize] <- dbetabinom(0,1:maxSize,pCom,dCom)
			dbCom <- dbetabinom(sequence(1:maxSize),rep(1:maxSize,1:maxSize),pCom,dCom)
		}
		m[2,1] <- 1 - m[1,1]

		pc <- matrix(0,maxSize,maxSize)
		pc[upper.tri(pc,diag=TRUE)] <- dbCom

		if(pHouse == 0){
			m[2:(maxSize+1),1:maxSize] <- pc
		}else{
			M <- maxSize-1
			M2 <- M^2
			ph <- array(0,c(M,M,M))

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
	
			for(n in 2:maxSize){
				for(k in 1:(n-1)) m[k+1,n] <- sum(pc[k:1,n] * pt[k:1 + (0:(k-1))*M + ((n-k-1):(n-2))*M2])
				m[n+1,n] <- 1 - sum(m[,n])
			}
		}
	}
	m
}

getP <- function(pCom,dCom,pHouse,dHouse,phiV,phiA,piV,piA){
	p <- 0
	if(phiA > 0 & phiA <= 1 & phiV > 0 & phiV <= 1 & pCom > 0 & pCom < 1 & pHouse >= 0 & pHouse < 1 & dHouse >= 0 & dCom >= 0 & 
	   piV <= 1 & piA <= 1 & piV > 0 & piA > 0){
	
		M <- getProbsMaster(pCom,dCom,pHouse,dHouse,maxSize)

		x <- M[Mi] * coef * phiV^tpV * (1-phiV)^fnV * phiA^tpA * (1-phiA)^fnA *
					  piV^tnV * (1-piV)^fpV * piA^tnA * (1-piA)^fpA

		p <- tapply(x,rep(seq_along(numTerms),numTerms),sum)
	}
	p
}

getLL <- function(pCom,dCom,pHouse,dHouse,phiV,phiA,piV,piA) 
	sum(cf$freq * log(getP(pCom,dCom,pHouse,dHouse,phiV,phiA,piV,piA)))

getSol <- function(initpar,fn)
	optim(initpar,fn,control=list(abstol = 1e-10, maxit=10000))


