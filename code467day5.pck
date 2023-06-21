code467day5.pck <-
c("code467day5.pck", "Hotellings.twosample.conf", "Bootstrap.twosample.simconf", 
"gui.bootstrap.twosample.simconf", "Manova.calculation.multiway.test0", 
"Manova.calculation.multiway.test.portmanteau", "Perm.calc.multiwayManova1", 
"gui.multiway.manova.test.portmanteau", "T620", "P620", "int.list.620", 
"mdet")
Hotellings.twosample.conf <-
function(dat0,idcol,conmat,alpha=.05,mu00=0){
mat0<-convert.data(dat0,idcol)
id1<-mat0[,idcol]
ud1<-unique(id1)
#mean calculation and SS calculation
mean.mat<-NULL
nvec<-NULL
SSW<-list()

for(i in 1:length(ud1)){
	I1<-(id1==ud1[i])
	m0<-mat0[I1,-idcol]
	mu0<-apply(m0,2,mean)
	musubt<-function(vec){
	vec-mu0
	}
	m1<-apply(m0,1,musubt)
	SSW[[i]]<-m1%*%t(m1)
	mean.mat<-rbind(mean.mat,c(mu0))
	nvec<-c(nvec,length(m1[1,]))
}
mudiff<-mean.mat[1,]-mean.mat[2,]
S<-SSW[[1]]/(nvec[1]^2)+SSW[[2]]/(nvec[2]^2)
H<-t(mudiff-mu00)%*%gen.inv1(S)%*%(mudiff-mu00)
vv<-diag(conmat%*%S%*%conmat)
p<-length(mudiff)
denom<-0
for(i in 1:2){
	denom<-denom+(1/nvec[i])*tr(((SSW[[i]]/(nvec[i]^2))%*%gen.inv1(S))^2)
}
nu<-(p*(p+1))/denom
H1<-H*(nu-p+1)/(nu*p)
qval<-(nu*p)/(nu-p+1)*qf(1-alpha,p,nu-p+1)
mmudiff<-conmat%*%mudiff
confmat<-cbind(mmudiff-sqrt(qval*vv),mmudiff,mmudiff+sqrt(qval*vv))
list(H=H,confmat=confmat,conf=100*(1-alpha))
}
Bootstrap.twosample.simconf <-
function(zdata,idcol,conmat,mu00=0,nboot=10000,alpha=.05)
{
Hconf<-Hotellings.twosample.conf(zdata,idcol,conmat,alpha,mu00)$confmat
mmu0<-Hconf[,2]
mmumat<-NULL
v1<-zdata[,idcol]
vid<-unique(v1)
nvec<-NULL
for(i in 1:2){
	nvec<-c(nvec,sum(v1==vid[i]))
}
I1<-v1==vid[1]
I2<-v1==vid[2]
zbdata<-zdata
for(i in 1:nboot){
if(floor(i/500)==(i/500)){print(i)}
bootsamp1<-sample(nvec[1],replace=T)
bootsamp2<-sample(nvec[2],replace=T)
zbdata[I1,]<-(zdata[I1,])[bootsamp1,]
zbdata[I2,]<-(zdata[I2,])[bootsamp2,]
mmub<-Hotellings.twosample.conf(zbdata,idcol,conmat,alpha,mu00)$confmat[,2]
mmumat<-rbind(mmumat,c(mmub-mmu0))
}
cv1<-cov(mmumat)
z1<-NULL
cv2<-gen.inv1(cv1)
print(dim(mmumat))
print(dim(cv2))
for(i in 1:nboot){
z1<-c(z1,t(mmumat[i,])%*%cv2%*%(mmumat[i,]))
}
o1<-order(z1)
q1<-quantile(z1,1-alpha)
I1<-(z1<=q1)
mubound<-mmumat[I1,]
print(dim(mubound))
muconf<-NULL
for(j in 1:length(mmu0)){
muconf<-rbind(muconf,c(mmu0[j]-max(mubound[,j]),mmu0[j],mmu0[j]-min(mubound[,j])))
}

list(bootconf=muconf, Hotelling=Hconf)
}
gui.bootstrap.twosample.simconf <-
function(){
library(tcltk)
#function(zdata,idcol,conmat,mu00=0,nboot=10000,alpha=.05)

inputs <- function(){

   x <- tclVar("data")
   y <- tclVar("idcol")
   z <- tclVar("conmat")
  w<-tclVar("mu00")
   zx<-tclVar("10000")
   zy<-tclVar(".05")

   tt <- tktoplevel()
   tkwm.title(tt,"Choose parameters for new function                   ")
   x.entry <- tkentry(tt, textvariable=x)
   y.entry <- tkentry(tt, textvariable=y)
   z.entry <- tkentry(tt, textvariable=z)
 w.entry <- tkentry(tt, textvariable=w)
   zx.entry<-tkentry(tt, textvariable=zx)
   zy.entry<-tkentry(tt, textvariable=zy) 
   reset <- function()
    {
     tclvalue(x)<-""
     tclvalue(y)<-""
     tclvalue(z)<-""
tclvalue(w)<-""
      tclvalue(zx)<-""
	tclvalue(zy)<-""

 }

   reset.but <- tkbutton(tt, text="Reset", command=reset)

   submit <- function() {
     x <- tclvalue(x)
     y <- tclvalue(y)
     z <- tclvalue(z)
     w<-tclvalue(w)
     zx<-tclvalue(zx)
     zy<-tclvalue(zy)
     e <- parent.env(environment())
     e$x <- x
     e$y <- y
     e$z <- z
e$w<-w
     e$zx<-zx
     e$zy<-zy
        tkdestroy(tt)
   }

   submit.but <- tkbutton(tt, text="start", command=submit)
   tkgrid(tklabel(tt,text="Input data matrix"),columnspan=2)
   tkgrid(tklabel(tt,text="data"), x.entry, pady=10, padx =30)


   tkgrid(tklabel(tt,text="Input idcolumn"),columnspan=2)
   tkgrid(tklabel(tt,text="col #"), y.entry, pady=10, padx =30)


   tkgrid(tklabel(tt,text="contrasts"),columnspan=2)
   tkgrid(tklabel(tt,text="contrast matrix"), z.entry, pady=10, padx =30)

	tkgrid(tklabel(tt,text="Hypothesized mean difference"),columnspan=2)
   tkgrid(tklabel(tt,text="mu00"), w.entry, pady=10, padx =30)
   
   tkgrid(tklabel(tt,text="Nboot"),columnspan=2)
   tkgrid(tklabel(tt,text="N bootstrap samples"), zx.entry, pady=10, padx =30)

   tkgrid(tklabel(tt,text="alpha"),columnspan=2)
   tkgrid(tklabel(tt,text="alpha for 100(1-alpha) interval"), zy.entry, pady=10, padx =30)

  
   tkgrid(submit.but, reset.but)

  tkwait.window(tt)
  return(c(x,y,z,w,zx,zy))
}
#Now run the function like:
predictor_para <- inputs()
print(predictor_para)
xmat<-eval(parse(text=predictor_para[1]))
xidcol<-eval(parse(text=predictor_para[2]))
xconmat<-eval(parse(text=predictor_para[3]))
xmu0<-eval(parse(text=predictor_para[4]))
xnboot<-eval(parse(text=predictor_para[5]))
xalpha<-eval(parse(text=predictor_para[6]))

Bootstrap.twosample.simconf(xmat,xidcol,xconmat,xmu0,xnboot,xalpha)

}
Manova.calculation.multiway.test0 <-
function(dat0,idcol,rid=0){
#This program calculates SS for the full interaction. We assume that there are no more than 9 levels for any
#factor. It is assumed this is higher than 2 way and every cell has more than 1 observation
if(abs(rid[1])>.2){
mat0<-convert.data(dat0,c(idcol,rid))
}
else{
mat0<-convert.data(dat0,c(idcol))
}
nway<-length(idcol)
idvec0<-0
for(i in 1:nway){
idvec0<-idvec0+mat0[,idcol[i]]*10^(i-1)
}
mat1<-cbind(idvec0,mat0[,-c(idcol,rid)])
dum1<-Manova.calculation.1way.test(mat1,1)
SSE<-dum1$SSW
SSB<-dum1$SSB
SST<-dum1$SST
mean.mat<-dum1$mumat
residmat<-dum1$resid
list(lambda=mdet(SSE)/(mdet(SSE+SSB)),SSE=dum1$SSW,SSB=dum1$SSB,SST=dum1$SST,mean.mat=mean.mat,mumatA=dum1$mumatA,residmat=residmat)
}
Manova.calculation.multiway.test.portmanteau <-
function(dat0,idcol,int.str){
#list(lambda=mdet(SSW)/mdet(SSW+SSB),SSW=SSW,SSB=SSB,SST=SST,mumat=mean.mat)
#This program calculates SS for each factor and combination of factors for all terms listed in the int.str list
#
n1<-length(int.str)
out.str<-rep(list(),(n1+1))
SSElist<-rep(list(),(n1+1))
SSBlist<-rep(list(),(n1+1))
SSTlist<-rep(list(),(n1+1))
mulist<-rep(list(),(n1+1))
muAlist<-rep(list(),(n1+1))
residlist<-rep(list(),(n1+1))
lambdavec<-NULL
clambdavec<-NULL
n1<-length(int.str)
#SSE and full interaction
dum1<-Manova.calculation.multiway.test0(dat0,idcol)
out.str[[1]]<-dum1
lambdavec<-NULL
SSElist[[1]]<-dum1$SSE
SSBlist[[1]]<-dum1$SSB
SSTlist[[1]]<-dum1$SST
mulist[[1]]<-dum1$mean.mat
residlist[[1]]<-dum1$residmat
muAlist[[1]]<-dum1$mumatA
SSBtot<-0
dat00<-convert.data(dat0,idcol)
for(j in 1:n1){
	idcoltemp<-idcol[int.str[[j]]]

throw.away<-(idcol[-int.str[[j]]])
	 dum2<-Manova.calculation.multiway.test0(dat00,idcoltemp,rid=throw.away)
	lambdavec<-c(lambdavec,mdet(dum1$SSE)/mdet(dum1$SSE+dum2$SSB))
	clambdavec<-c(clambdavec,det(dum1$SSE)/det(dum1$SSE+dum2$SSB))
	out.str[[(j+1)]]<-list(int=int.str[[j]],raw=dum2)
	SSElist[[(j+1)]]<-(dum2$SSE)
	SSBlist[[(j+1)]]<-(dum2$SSB)
	SSTlist[[(j+1)]]<-dum2$SST
	SSBtot<-SSBtot+dum2$SSB
	mulist[[(j+1)]]<-dum2$mean.mat
	residlist[[(j+1)]]<-dum2$residmat
	muAlist[[(j+1)]]<-dum2$mumatA

}
lambda.int<-mdet(dum1$SSE)/mdet(dum1$SST-SSBtot)
clambda.int<-det(dum1$SSE)/det(dum1$SST-SSBtot)
out<-list(Lambda=c(lambdavec,lambda.int),cLambda=c(clambdavec,clambda.int),outinf=out.str,SSE=SSElist,SSB=SSBlist,SST=SSTlist,mu=mulist,muAlist=muAlist,resid=residlist)
out
}
Perm.calc.multiwayManova1 <-
function(mat,idcol,nperm=10000){
#function(dat0,idcol,int.str)
#NOTE IDCOL SHOULD BE IN ORDER THEY APPEAR IN MATRIX
n00<-length(idcol)
int.str<-rep(list(),n00)
for(i in 1:n00){
int.str[[i]]<-i
}
dum<-Manova.calculation.multiway.test.portmanteau(mat,idcol,int.str)
lambda0<-dum$Lambda
clambda0<-dum$cLambda
m1<-mat[,-idcol]
idvec<-mat[,idcol]
n0<-length(idvec[1,])
n1<-length(m1[,1])
lambdamat<-NULL
clambdamat<-NULL

print(lambda0)
print(clambda0)
for(i in 1:nperm){
	nsamp<-sample(n1)
	mp<-m1[nsamp,]
	matp<-cbind(idvec,mp)
	dumP<-Manova.calculation.multiway.test.portmanteau(matp,c(1:n0),int.str)
	lambdamat<-rbind(lambdamat,dumP$Lambda)
	clambdamat<-rbind(clambdamat,dumP$cLambda)
if((i/500)==floor(i/500)){print(c(i))}
}
lambcomp<-function(vec){1*(vec<lambda0)}
dumcomp<-t(apply(lambdamat,1,lambcomp))
pvec<-apply(dumcomp,2,sum)/nperm
clambcomp<-function(vec){1*(vec<clambda0)}
cdumcomp<-t(apply(clambdamat,1,clambcomp))
cpvec<-apply(cdumcomp,2,sum)/nperm

dum$p<-pvec
dum$cp<-cpvec
dum	
}
gui.multiway.manova.test.portmanteau <-
function(){
library(tcltk)
inputs <- function(){

   x <- tclVar("T620")
   y <- tclVar("c(1:3)")
     w<-tclVar("10000")
   tt <- tktoplevel()
   tkwm.title(tt,"Choose parameters for new function                   ")
   x.entry <- tkentry(tt, textvariable=x)
   y.entry <- tkentry(tt, textvariable=y)
     w.entry<-tkentry(tt,textvariable=w)
     reset <- function()
    {
     tclvalue(x)<-""
     tclvalue(y)<-""
          tclvalue(w)<-""
       }

   reset.but <- tkbutton(tt, text="Reset", command=reset)

   submit <- function() {
     x <- tclvalue(x)
     y <- tclvalue(y)
         w<-tclvalue(w)
           e <- parent.env(environment())
     e$x <- x
     e$y <- y
    
     e$w<-w
        tkdestroy(tt)
   }

   submit.but <- tkbutton(tt, text="start", command=submit)
   tkgrid(tklabel(tt,text="Input data matrix"),columnspan=2)
   tkgrid(tklabel(tt,text="data"), x.entry, pady=10, padx =30)


   tkgrid(tklabel(tt,text="Input id column matrix"),columnspan=2)
   tkgrid(tklabel(tt,text="idcol"), y.entry, pady=10, padx =30)

  


   tkgrid(tklabel(tt,text="Input number of permutation tries"),columnspan=2)
   tkgrid(tklabel(tt,text="nperm"), w.entry, pady=10, padx =30)

  
   tkgrid(submit.but, reset.but)

  tkwait.window(tt)
  return(c(x,y,w))
}
#Now run the function like:
predictor_para <- inputs()
print(predictor_para)
mat<-eval(parse(text=predictor_para[1]))
idcol<-eval(parse(text=predictor_para[2]))
#int.str<-eval(parse(text=predictor_para[3]))
nperm<-eval(parse(text=predictor_para[3]))

Perm.calc.multiwayManova1(mat,idcol,nperm)

}
T620 <-
structure(list(V1 = c("Low", "Low", "Low", "Low", "Low", "Low", 
"Low", "Low", "High", "High", "High", "High", "High", "High", 
"High", "High"), V2 = c("Simple", "Simple", "Simple", "Simple", 
"Complex", "Complex", "Complex", "Complex", "Simple", "Simple", 
"Simple", "Simple", "Complex", "Complex", "Complex", "Complex"
), V3 = c("Novice", "Novice", "Guru", "Guru", "Novice", "Novice", 
"Guru", "Guru", "Novice", "Novice", "Guru", "Guru", "Novice", 
"Novice", "Guru", "Guru"), V4 = c(3, 2.2999999999999998, 1.7, 
1.2, 6.7000000000000002, 7.0999999999999996, 5.5999999999999996, 
4.5, 4.5, 4.7000000000000002, 3.1000000000000001, 3, 7.9000000000000004, 
6.9000000000000004, 5, 5.2999999999999998), V5 = c(6.2999999999999998, 
5.2999999999999998, 2.1000000000000001, 1.6000000000000001, 12.6, 
12.800000000000001, 8.8000000000000007, 9.1999999999999993, 9.5, 
10.699999999999999, 6.2999999999999998, 5.5999999999999996, 15.6, 
14.9, 10.4, 10.4), V6 = c(9.3000000000000007, 7.5999999999999996, 
3.7999999999999998, 2.7999999999999998, 19.300000000000001, 19.899999999999999, 
14.4, 13.699999999999999, 14, 15.4, 9.4000000000000004, 8.5999999999999996, 
23.5, 21.800000000000001, 15.4, 15.699999999999999)), class = "data.frame", row.names = c(NA, 
-16L))
P620 <-
structure(list(V1 = c("Low", "Low", "Low", "Low", "Low", "Low", 
"Low", "Low", "High", "High", "High", "High", "High", "High", 
"High", "High"), V2 = c("Simple", "Simple", "Simple", "Simple", 
"Complex", "Complex", "Complex", "Complex", "Simple", "Simple", 
"Simple", "Simple", "Complex", "Complex", "Complex", "Complex"
), V3 = c("Novice", "Novice", "Guru", "Guru", "Novice", "Novice", 
"Guru", "Guru", "Novice", "Novice", "Guru", "Guru", "Novice", 
"Novice", "Guru", "Guru"), V4 = c(6.7000000000000002, 4.5, 3, 
4.5, 7.0999999999999996, 3.1000000000000001, 1.2, 5.5999999999999996, 
3, 6.9000000000000004, 2.2999999999999998, 7.9000000000000004, 
5, 4.7000000000000002, 1.7, 5.2999999999999998), V5 = c(12.6, 
9.1999999999999993, 5.5999999999999996, 9.5, 12.800000000000001, 
6.2999999999999998, 1.6000000000000001, 8.8000000000000007, 6.2999999999999998, 
14.9, 5.2999999999999998, 15.6, 10.4, 10.699999999999999, 2.1000000000000001, 
10.4), V6 = c(19.300000000000001, 13.699999999999999, 8.5999999999999996, 
14, 19.899999999999999, 9.4000000000000004, 2.7999999999999998, 
14.4, 9.3000000000000007, 21.800000000000001, 7.5999999999999996, 
23.5, 15.4, 15.4, 3.7999999999999998, 15.699999999999999)), row.names = c(NA, 
-16L), class = "data.frame")
int.list.620 <-
list(1, 2, 3, c(1, 2), c(1, 3), c(2, 3))
mdet <-
function(mat){
z1<-eigen(mat)$val
I1<-z1>1e-12
exp(sum(log(z1[I1])))
}
