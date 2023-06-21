code467day4.pck <-
c("code467day4.pck", "orthogonalize", "explore", "gui.explore", 
"gui.orth", "Manova.calculation.1way.test", "Perm.calc.Manova1way", 
"convert.data", "mdet", "gui.1waymanovatest", "Hotellings.twosample", 
"Hotellings.twosample.perm.test", "gui.twosample.perm", "tr")
orthogonalize <-
function(mat){
n1<-length(mat[,1])
v1<-mat[1,]
u1<-v1/sqrt(sum(v1^2))
m0<-t(as.matrix(u1,length(u1),1))
print(m0)
for(i in 2:n1){
	v0<-mat[i,]
	vv0<-v0-t(m0)%*%m0%*%v0
	u0<-vv0/sqrt(sum(vv0^2))
	m0<-rbind(m0,c(u0))
}
m0
}
explore <-
function(mat,dat,do.pairs=T,do.3d=F,do.animate=F,vec3d=c(1,2,3)){
library(scatterplot3d)
library(tourr)
dat0<-dat%*%t(mat)
n1<-length(dat0[1,])
dimnames(dat0)[[2]]<-c(1:n1)
if(do.pairs){
pairs(dat0)
}else{
	if(do.3d){
par(mfrow=c(1,1))
	scatterplot3d(dat0[,vec3d])
	}else{
par(mfrow=c(1,1))
	animate(dat0)
}
}
dat0
}
gui.explore <-
function(){
library(tcltk)
inputs <- function(){
#function(mat,dat,do.pairs=T,do.3d=F,do.animate=F,vec3d=c(1,2,3))

   x <- tclVar("matdir")
   y <- tclVar("dat")
   z <- tclVar("T")
   w<-tclVar("F")
zx<-tclVar("F")
zy<-tclVar("c(1,2,3)")

   tt <- tktoplevel()
   tkwm.title(tt,"Choose parameters for new function                   ")
   x.entry <- tkentry(tt, textvariable=x)
   y.entry <- tkentry(tt, textvariable=y)
   z.entry <- tkentry(tt, textvariable=z)
   w.entry <- tkentry(tt, textvariable=w)
   zx.entry <- tkentry(tt, textvariable=zx)
   zy.entry <- tkentry(tt, textvariable=zy)
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
     w <- tclvalue(w)
     zx <- tclvalue(zx)
     zy <- tclvalue(zy)
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
   tkgrid(tklabel(tt,text="Input direction matrix"),columnspan=2)
   tkgrid(tklabel(tt,text="mat"), x.entry, pady=10, padx =30)


   tkgrid(tklabel(tt,text="Input data matrix"),columnspan=2)
   tkgrid(tklabel(tt,text="data"), y.entry, pady=10, padx =30)

   tkgrid(tklabel(tt,text="Pairs plot?"),columnspan=2)
   tkgrid(tklabel(tt,text="Pairs?"), z.entry, pady=10, padx =30)

   tkgrid(tklabel(tt,text="Tour"),columnspan=2)
   tkgrid(tklabel(tt,text="Animate?"), w.entry, pady=10, padx =30)
   
   tkgrid(tklabel(tt,text="3d plot?"),columnspan=2)
   tkgrid(tklabel(tt,text="3d?"), zx.entry, pady=10, padx =30)

   tkgrid(tklabel(tt,text="Rows dmat for 3d"),columnspan=2)
   tkgrid(tklabel(tt,text="rows?"), zy.entry, pady=10, padx =30)

  
   tkgrid(submit.but, reset.but)

  tkwait.window(tt)
  return(c(x,y,z,w,zx,zy))
}
#Now run the function like:
predictor_para <- inputs()
print(predictor_para)
mat<-eval(parse(text=predictor_para[1]))
dat<-eval(parse(text=predictor_para[2]))
do.pairs<-eval(parse(text=predictor_para[3]))
do.animate<-eval(parse(text=predictor_para[4]))
do.3d<-eval(parse(text=predictor_para[5]))
zvec<-eval(parse(text=predictor_para[6]))
explore(mat,dat,do.pairs,do.3d,do.animate,zvec)

}
gui.orth <-
function(){
library(tcltk)
inputs <- function(){
#function(m1,ntest=99,Q=.01){

   x <- tclVar("mat")
    
   tt <- tktoplevel()
   tkwm.title(tt,"Choose parameters for new function                   ")
   x.entry <- tkentry(tt, textvariable=x)
        reset <- function()
    {
     tclvalue(x)<-""
            }

   reset.but <- tkbutton(tt, text="Reset", command=reset)

   submit <- function() {
     x <- tclvalue(x)
     e <- parent.env(environment())
     e$x <- x
         tkdestroy(tt)
   }

   submit.but <- tkbutton(tt, text="start", command=submit)
   tkgrid(tklabel(tt,text="Input row matrix of directions"),columnspan=2)
   tkgrid(tklabel(tt,text="data"), x.entry, pady=10, padx =30)


  
   tkgrid(submit.but, reset.but)

  tkwait.window(tt)
  return(c(x))
}
#Now run the function like:
predictor_para <- inputs()
print(predictor_para)
mat<-eval(parse(text=predictor_para[1]))
orthogonalize(mat)

}
Manova.calculation.1way.test <-
function(dat0,idcol){
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
T0<-apply(mat0[,-idcol],1,musubt)
SST<-(T0)%*%t(T0)
list(lambda=mdet(SSW)/mdet(SSW+SSB),SSW=SSW,SSB=SSB,SST=SST,mumat=mean.mat,mumatA=mumat,resid=residmat)
}
Perm.calc.Manova1way <-
function(mat,idcol,nperm=10000){
dum<-Manova.calculation.1way.test(mat,idcol)
lambda0<-dum$lambda
m1<-mat[,-idcol]
idvec<-mat[,idcol]
n1<-length(m1[,1])
lambdavec<-NULL
for(i in 1:nperm){
	nsamp<-sample(n1)
	mp<-m1[nsamp,]
	matp<-cbind(idvec,mp)
	lambdavec<-c(lambdavec,Manova.calculation.1way.test(matp,1)$lambda)
if((i/500)==floor(i/500)){print(c(i))}

}
pval<-sum(lambdavec<lambda0)/nperm
dum$p<-pval
#dum$lambdavec<-lambdavec
dum	
}
convert.data <-
function(dat,idcol){
n1<-length(idcol)
#print(idcol)
dat00<-as.matrix(dat[,-idcol])
n0<-length(dat[,1])
n2<-length(dat[1,])
dat0<-matrix(rep(0,n0*n2),n0,n2)
zmatnumid<-c(1:n2)[-idcol]
dat0[,zmatnumid]<-dat00
for(i in 1:n1){
	idv<-dat[,idcol[i]]
	un1<-unique(idv)
	n2<-length(un1)
	vec0<-rep(0,length(idv))
	for(j in 1:n2){
		Id0<-(idv==un1[j])
		vec0[Id0]<-j
	}
	dat0[,idcol[i]]<-vec0
	}	
dat0
}
mdet <-
function(mat){
z1<-eigen(mat)$val
I1<-z1>1e-12
exp(sum(log(z1[I1])))
}
gui.1waymanovatest <-
function(){
library(tcltk)
inputs <- function(){
#function(m1,ntest=99,Q=.01){

   x <- tclVar("mat")
   y <- tclVar("idcol")
   z <- tclVar("10000")
  
   tt <- tktoplevel()
   tkwm.title(tt,"Choose parameters for new function                   ")
   x.entry <- tkentry(tt, textvariable=x)
   y.entry <- tkentry(tt, textvariable=y)
   z.entry <- tkentry(tt, textvariable=z)
     reset <- function()
    {
     tclvalue(x)<-""
     tclvalue(y)<-""
     tclvalue(z)<-""
       }

   reset.but <- tkbutton(tt, text="Reset", command=reset)

   submit <- function() {
     x <- tclvalue(x)
     y <- tclvalue(y)
     z <- tclvalue(z)
           e <- parent.env(environment())
     e$x <- x
     e$y <- y
     e$z <- z
        tkdestroy(tt)
   }

   submit.but <- tkbutton(tt, text="start", command=submit)
   tkgrid(tklabel(tt,text="Input data matrix"),columnspan=2)
   tkgrid(tklabel(tt,text="data"), x.entry, pady=10, padx =30)


   tkgrid(tklabel(tt,text="Input id column"),columnspan=2)
   tkgrid(tklabel(tt,text="idcol"), y.entry, pady=10, padx =30)


   tkgrid(tklabel(tt,text="Input number of permutation tries"),columnspan=2)
   tkgrid(tklabel(tt,text="nperm"), z.entry, pady=10, padx =30)

  
   tkgrid(submit.but, reset.but)

  tkwait.window(tt)
  return(c(x,y,z))
}
#Now run the function like:
predictor_para <- inputs()
print(predictor_para)
mat<-eval(parse(text=predictor_para[1]))
idcol<-eval(parse(text=predictor_para[2]))
nperm<-eval(parse(text=predictor_para[3]))
Perm.calc.Manova1way(mat,idcol,nperm)

}
Hotellings.twosample <-
function(dat0,idcol,mu00){
mat0<-convert.data(dat0,idcol)
id1<-mat0[,idcol]
ud1<-unique(id1)
#mean calculation and SS calculation
mean.mat<-NULL
nvec<-NULL
SSW<-list()
SS2<-0
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
p<-length(mudiff)
denom<-0
for(i in 1:2){
	denom<-denom+(1/nvec[i])*tr(((SSW[[i]]/(nvec[i]^2))%*%gen.inv1(S))^2)
}
nu<-(p*(p+1))/denom
H1<-H*(nu-p+1)/(nu*p)
pval.standard<-1-pf(H1,p,nu-p+1)
list(H=H,H1=H1,P=pval.standard,mu=mudiff-mu00,df=c(p,nu-p+1))
}
Hotellings.twosample.perm.test <-
function(dat0,idcol,mu0=0,nperm=10000){
dum<-Hotellings.twosample(dat0,idcol,mu0)
H1<-dum$H1
mat0<-convert.data(dat0,idcol)
idvec<-mat0[,idcol]
n0<-length(idvec)
H1vec<-NULL
for(i in 1:nperm){
	if((i/500)==floor(i/500)){print(i)}
	nvec<-sample(n0)
	matp<-mat0[nvec,-idcol]
	datp<-cbind(idvec,matp)
	H1p<-Hotellings.twosample(datp,1,mu0)$H1
	H1vec<-c(H1vec,H1p)
}
pval<-sum(H1vec>c(H1))/nperm
dum$permP<-pval
dum
}
gui.twosample.perm <-
function(){
library(tcltk)
inputs <- function(){
#function(m1,ntest=99,Q=.01){

   x <- tclVar("mat")
   y <- tclVar("1")
   z <- tclVar("0")
   w<-tclVar("10000")
   tt <- tktoplevel()
   tkwm.title(tt,"Choose parameters for new function                   ")
   x.entry <- tkentry(tt, textvariable=x)
   y.entry <- tkentry(tt, textvariable=y)
   z.entry <- tkentry(tt, textvariable=z)
   w.entry <- tkentry(tt, textvariable=w)

     reset <- function()
    {
     tclvalue(x)<-""
     tclvalue(y)<-""
     tclvalue(z)<-""
 tclvalue(w)<-""

       }

   reset.but <- tkbutton(tt, text="Reset", command=reset)

   submit <- function() {
     x <- tclvalue(x)
     y <- tclvalue(y)
     z <- tclvalue(z)
	w<-tclvalue(w)
           e <- parent.env(environment())
     e$x <- x
     e$y <- y
     e$z <- z
	e$w<-w
        tkdestroy(tt)
   }

   submit.but <- tkbutton(tt, text="start", command=submit)
   tkgrid(tklabel(tt,text="Input data matrix"),columnspan=2)
   tkgrid(tklabel(tt,text="data"), x.entry, pady=10, padx =30)


   tkgrid(tklabel(tt,text="Input id column"),columnspan=2)
   tkgrid(tklabel(tt,text="idcol"), y.entry, pady=10, padx =30)

   tkgrid(tklabel(tt,text="Input mu0"),columnspan=2)
   tkgrid(tklabel(tt,text="Hypothesized difference"), z.entry, pady=10, padx =30)


   tkgrid(tklabel(tt,text="Input number of permutation tries"),columnspan=2)
   tkgrid(tklabel(tt,text="nperm"), w.entry, pady=10, padx =30)

  
   tkgrid(submit.but, reset.but)

  tkwait.window(tt)
  return(c(x,y,z,w))
}
#Now run the function like:
predictor_para <- inputs()
print(predictor_para)
mat<-eval(parse(text=predictor_para[1]))
idcol<-eval(parse(text=predictor_para[2]))
mu0<-eval(parse(text=predictor_para[3]))

nperm<-eval(parse(text=predictor_para[4]))
Hotellings.twosample.perm.test(mat,idcol,mu0,nperm)

}
tr <-
function(mat){sum(diag(mat))}
