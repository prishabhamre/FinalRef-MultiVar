class10code.pck <-
c("class10code.pck", "Fisher.disc.space", "Perm.cov.test", "gen.inv1.2"
)
Fisher.disc.space <-
function(dat0,idcol,num.vec,do.plot=T){
mat0<-convert.data(dat0,idcol)
id1<-mat0[,idcol]
ud1<-unique(id1)
#mean calculation and SS calculation
mean.mat<-NULL
nvec<-NULL
SSB<-0
SSW<-0
residmat<-mat0
mumat<-mat0
detvec<-NULL
for(i in 1:length(ud1)){
	I1<-(id1==ud1[i])
	m0<-mat0[I1,-idcol]
#print(m0)
	mu0<-apply(m0,2,mean)
	musubt<-function(vec){
	vec-mu0
	}
	m1<-apply(m0,1,musubt)
	SSW<-SSW+m1%*%t(m1)
	detvec<-c(detvec,mdet(m1%*%t(m1)))
	mean.mat<-rbind(mean.mat,c(mu0))
	nvec<-c(nvec,length(m1[1,]))
	residmat[I1,-idcol]<-t(m1)
	for(i in 1:nvec[i]){
		mumat[I1,-idcol][i,]<-mu0
	}
}
mutot<-(nvec%*%mean.mat)/sum(nvec)
musubt<-function(vec){vec-mutot}
b0<-apply(mean.mat,1,musubt)
#print(dim(b0))
#print(length(nvec))

for(i in 1:length(nvec)){
	b0[,i]<-b0[,i]*sqrt(nvec[i])
}
SSB<-b0%*%t(b0)
mat.toexp<-gen.inv1.2(SSW)%*%SSB%*%gen.inv1.2(SSW)
z1<-eigen(mat.toexp)
#print(z1)
n1<-length(z1$val)
matdisp<-mat0[,-idcol]%*%z1$vec[,1:num.vec]
if(do.plot){
par(mfrow=c(num.vec,num.vec))
for(i in 1:num.vec){
for(j in 1:num.vec){
plot(matdisp[,i],matdisp[,j],type="n",main=paste("vec",i,"vs vec",j))
for(k in 1:length(ud1)){
	I1<-(id1==ud1[k])
	text(matdisp[I1,i],matdisp[I1,j],k,col=k)
}
}
}
}
lambda1<-mdet(SSW)/detvec
llambda<-(sum(log(lambda1)))


list(mat.out=cbind(matdisp,id1),llambda=llambda)
}
Perm.cov.test <-
function(dat0,idcol,numvec){
mat0<-convert.data(dat0,idcol)
n1<-length(mat0[,1])
nperm<-10000
llambda0<-Fisher.disc.space(mat0,idcol,numvec,T)$llambda
llambdavec<-NULL
for(i in 1:nperm){
if((i/500)==floor(i/500)){print(i)}
idvec<-mat0[,idcol]
nv1<-sample(n1,n1)
mat1<-mat0[nv1,]
mat1[,idcol]<-idvec
llambdavec<-c(llambdavec,Fisher.disc.space(mat1,idcol,numvec,F)$llambda)
}
plambda<-sum(llambdavec>llambda0)/nperm
plambda
	
}
gen.inv1.2 <-
function(mat, thresh = 1e-10)
{
	v1 <- sum(is.na(mat))
	v2 <- sum(is.inf(mat))
	if((v1 + v2) > 0.5) {
		print(mat)
	}
	e1 <- eigen(mat, symmetric = T)
	val <- Re(e1$val)
	vec <- Re(e1$vec)
	val1 <- val/max(val)
	#
	#	print("normalized eigen values")
	#	print(val1)	#
	#	n1 <- length(val1)
	#	plot(c(1:n1), abs(val1), log = "y", xlab = "eigen rank", ylab
	#		 = "log10 of value")
	I1 <- val1 > thresh
	I3 <- is.na(I1)
	if(sum(I3) < 0.5) {
		val2 <- val[I1]
		I2 <- I1
		if(sum(I2) > 1.5) {
			ret <- vec[, I1] %*% diag(1/sqrt(val2)) %*% t(vec[, I1])
		}
		else {
			v1 <- as.matrix(vec[, I1], length(c(vec[, I1])), 1)
			ret <- (1/val2) * v1 %*% t(v1)
		}
	}
	else {
		ret <- diag(length(I1)) * 0
	}
	ret
}
