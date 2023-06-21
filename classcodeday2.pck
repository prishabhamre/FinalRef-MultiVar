classcodeday2.pck <-
c("classcodeday2.pck", "Bootstrap.1sample", "Bootstrap.simconf", 
"HotellingTcalc", "HotellingTconf", "gen.inv1", "gui.bootstrap1sample", 
"gui.bootstrapsimconf")
Bootstrap.1sample <-
function(zdata,mean0=0,nboot=10000)
{
meandat<-apply(zdata,2,mean)
if(length(mean0)!=length(meandat)){
mean0<-rep(mean0,length(meandat))
}
T0calc<-HotellingTcalc(zdata,mean0)
Tvec<-NULL
nsamp<-length(zdata[,1])
for(i in 1:nboot){
if(floor(i/500)==(i/500)){print(i)}
bootsampind<-sample(nsamp,replace=T)
zbdata<-zdata[bootsampind,]
Tb<-HotellingTcalc(zbdata,meandat)$T
Tvec<-c(Tvec,Tb)
}
#print(length(Tvec))
#print(T0calc)
Pval<-sum(1*(c(Tvec)>c(T0calc$T)))/nboot
npval<-T0calc$pval
list(bootpval=Pval,normtheoryP=npval,T2=T0calc$T)
}
Bootstrap.simconf <-
function(zdata,conmat,nboot=10000,alpha=.05)
{
meandat<-apply(zdata,2,mean)
Hconf<-HotellingTconf(zdata,conmat)$confmat
mmu0<-Hconf[,2]
mmumat<-NULL
nsamp<-length(zdata[,1])
for(i in 1:nboot){
if(floor(i/500)==(i/500)){print(i)}
bootsampind<-sample(nsamp,replace=T)
zbdata<-zdata[bootsampind,]
mmub<-HotellingTconf(zbdata,conmat)$confmat[,2]
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
HotellingTcalc <-
function(data,mean0){
meanx<-apply(data,2,mean)
n<-length(data[,1])
cov<-var(data)/n
p<-length(data[1,])
Hstat<-t(meanx-mean0)%*%solve(cov)%*%(meanx-mean0)
HT<-((n-p)/((n-1)*p))*Hstat
pval<-1-pf(HT,p,n-p)
list(T=HT,pval=pval)
	
}
HotellingTconf <-
function(data,conmat,alpha=.05){
meanx<-apply(data,2,mean)
n<-length(data[,1])
cov<-(1/n)*var(data)
p<-length(data[1,])
n1<-length(conmat[,1])
qval<-((n-1)*p/(n-p))*qf(1-alpha,p,n-p)
conf.int<-NULL
vv<-diag(conmat%*%cov%*%t(conmat))
mmu<-conmat%*%(meanx)
confmat<-cbind(mmu-sqrt(qval*vv),mmu,mmu+sqrt(qval*vv))
list(alpha=alpha,confmat=confmat)
}
gen.inv1 <-
function(mat, thresh = 0)
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
			ret <- vec[, I1] %*% diag(1/val2) %*% t(vec[, I1])
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
gui.bootstrap1sample <-
function(){
library(tcltk)
inputs <- function(){

   x <- tclVar("data")
   y <- tclVar("0")
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


   tkgrid(tklabel(tt,text="Input H0 mean"),columnspan=2)
   tkgrid(tklabel(tt,text="mean vector name"), y.entry, pady=10, padx =30)


   tkgrid(tklabel(tt,text="Nboot"),columnspan=2)
   tkgrid(tklabel(tt,text="nboot"), z.entry, pady=10, padx =30)

  
   tkgrid(submit.but, reset.but)

  tkwait.window(tt)
  return(c(x,y,z))
}
#Now run the function like:
predictor_para <- inputs()
print(predictor_para)
xmat<-eval(parse(text=predictor_para[1]))
xmean<-eval(parse(text=predictor_para[2]))
xnboot<-eval(parse(text=predictor_para[3]))
Bootstrap.1sample(xmat,xmean,xnboot)

}
gui.bootstrapsimconf <-
function(){
library(tcltk)
inputs <- function(){

   x <- tclVar("data")
   y <- tclVar("con.mat")
   z <- tclVar("10000")
  w<-tclVar("alpha")

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


   tkgrid(tklabel(tt,text="Input contrast matrix"),columnspan=2)
   tkgrid(tklabel(tt,text="matrix name"), y.entry, pady=10, padx =30)


   tkgrid(tklabel(tt,text="Nboot"),columnspan=2)
   tkgrid(tklabel(tt,text="nboot"), z.entry, pady=10, padx =30)

	tkgrid(tklabel(tt,text="alpha for 100(1-alpha) conf"),columnspan=2)
   tkgrid(tklabel(tt,text="alpha"), w.entry, pady=10, padx =30)

  
   tkgrid(submit.but, reset.but)

  tkwait.window(tt)
  return(c(x,y,z,w))
}
#Now run the function like:
predictor_para <- inputs()
print(predictor_para)
xmat<-eval(parse(text=predictor_para[1]))
xconmat<-eval(parse(text=predictor_para[2]))
xnboot<-eval(parse(text=predictor_para[3]))
xalpha<-eval(parse(text=predictor_para[4]))
Bootstrap.simconf(xmat,xconmat,xnboot,xalpha)

}
