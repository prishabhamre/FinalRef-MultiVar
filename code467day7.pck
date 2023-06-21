code467day7.pck <-
c("code467day7.pck", "matrix.2ndorder.make", "Residplot2", "multismooth", 
"multireg", "multitest", "gui.selectres", "FuninT.pck", "Fun.in.Tstats"
)
matrix.2ndorder.make <-
function(x, only.quad=T){
x0<-x
dimn<-dimnames(x)[[2]] #extract the names of the variables
num.col<-length(x[1,]) # how many columns
for(i in 1:num.col){
# if we are doing all 2nd order
if(!only.quad){
for(j in i:num.col){
x0<-cbind(x0,x[,i]*x[,j])
dimn<-c(dimn,paste(dimn[i],dimn[j],sep=""))
#create interaction dimnames

}
}
else{
#in here only if doing only squared terms
x0<-cbind(x0,x[,i]*x[,i])
dimn<-c(dimn,paste(dimn[i],"2",sep="")) # squared dimmension names
}
}
dimnames(x0)[[2]]<-dimn
x0
}
Residplot2 <-
function(str,mat,idcol,i,varsel=0)
{
if(varsel[1]==0){
	n1<-length(str$resid[[i]][1,])
	varsel<-c(1:n1)
}

	m3<-str$resid[[i]]
	m1<-m3[,c(1,varsel)]
	m2<-convert.data(mat,idcol)
	m0<-cbind(m2[,idcol],m1[,-1])
	pairs(m0,pch=as.character(m1[,1]),col=m1[,1])
	assign("m99a3", m0, envir = .GlobalEnv)
	gui.selectres()
}
multismooth <-
function(mat,ycol,xcol,idcol,do.plot=T){
#note xcol and ycol are residual columns from a MANOVA so already centered, 
#here we are exploring the residuals of ymat regressed on xmat either as a
#multivariate linear model or as a multivariate smooth
xid<-mat[,idcol]
u1<-unique(xid)
SSE<-0
SST<-0
xmat0<-mat
for(i in 1:length(u1)){
I1<-xid==u1[i]
nx<-length(xcol)
ny<-length(ycol)
xmat<-mat[I1,xcol]
ymat<-mat[I1,ycol]
SST<-SST+t(ymat)%*%ymat
for(k in 1:ny){
y0<-ymat[,k]
smst.str<-list()
for(j in 1:nx){
smst<-smooth.spline(xmat[,j],y0)
smst.str[[j]]<-smst
resid<-c(y0)-c(predict(smst,xmat[,j])$y)
y0<-resid
}
xmat0[I1,ycol[k]]<-resid
}
}
residmat<-xmat0[,ycol]
SSE<-t(residmat)%*%residmat
if(do.plot){
pairs(xmat0,pch=as.character(xmat0[,idcol]),col=xmat0[,idcol])
}
list(lambda=mdet(SSE)/mdet(SST),SSE=SSE,SSB=SST-SSE,SST=SST,outmat=xmat0,smstr=smst.str)
}
multireg <-
function(mat,ycol,xcol,idcol,do.plot=F){
#note xcol and ycol are residual columns from a MANOVA so already centered, 
#here we are exploring the residuals of ymat regressed on xmat either as a
#multivariate linear model or as a multivariate smooth
v1<-mat[,idcol]
u1<-unique(v1)
SST<-0
SSE<-0
xmat0<-mat
for(i in 1:length(u1)){
I1<-v1==u1[i]
nx<-length(xcol)
ny<-length(ycol)
xmat<-mat[I1,xcol]
ymat<-mat[I1,ycol]
SST<-SST+t(ymat)%*%ymat
ls.str<-lsfit(xmat,ymat)
res<-ls.str$res
SSE<-SSE+t(res)%*%res
xmat0[I1,ycol]<-res
}
if(do.plot){
pairs(xmat0,pch=as.character(xmat0[,idcol]),col=xmat0[,idcol])
}
SSB<-SST-SSE
lambda<-mdet(SSE)/mdet(SST)
list(lambda=lambda,SSE=SSE,SSB=SSB,SST=SST,outmat=xmat0,lsout=ls.str)
}
multitest <-
function(mat,ycol,xcol,idcol){
reg.out<-multireg(mat,ycol,xcol,idcol,T)
print("type in Y and hit return twice when ready to move on")
duh<-scan(,what="character")
smooth.out<-multismooth(reg.out$outmat,ycol,xcol,idcol,T)
lambdar0<-reg.out$lambda
lambdas0<-smooth.out$lambda
lambdarvec<-NULL
lambdasvec<-NULL
n1<-length(mat[,1])
maty<-mat[,ycol]
mat0<-mat
xid<-mat[,idcol]
u1<-unique(xid)
#END OF BORROWED
for(i1 in 1:10000){
if((i1/500)==floor(i1/500)){print(i1)}
v1<-c(1:n1)
v2<-v1
for(j1 in 1:length(u1)){
I1<-xid==u1[j1]
v2a<-sample(v2[I1])
v1[I1]<-v2a
}
mat0[,ycol]<-maty[v1,]
preg.out<-multireg(mat0,ycol,xcol,idcol,F)
lambdarvec<-c(lambdarvec,preg.out$lambda)
lambdasvec<-multismooth(preg.out$outmat,ycol,xcol,idcol,F)$lambda
}
preg<-sum(lambdarvec<lambdar0)/10000
ps<-sum(lambdasvec<lambdas0)/10000
list(preg=preg,psmooth=ps,reg=reg.out,smooth=smooth.out)
}
gui.selectres <-
function(){
library(tcltk)
inputs <- function(){
#function(mat,dat,do.pairs=T,do.3d=F,do.animate=F,vec3d=c(1,2,3))

   x <- tclVar("m99a3")
   y <- tclVar("ycol")
   z <- tclVar("xcol")
   w<-tclVar("idcol")

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
   tkgrid(tklabel(tt,text="Matrix from residplot"),columnspan=2)
   tkgrid(tklabel(tt,text="dont touch"), x.entry, pady=10, padx =30)


   tkgrid(tklabel(tt,text="Input Y variables from plot"),columnspan=2)
   tkgrid(tklabel(tt,text="columns"), y.entry, pady=10, padx =30)

   tkgrid(tklabel(tt,text="Input X variables from plot"),columnspan=2)
   tkgrid(tklabel(tt,text="columns"), z.entry, pady=10, padx =30)

    tkgrid(tklabel(tt,text="Input idcolumnd from plot"),columnspan=2)
   tkgrid(tklabel(tt,text="1 columne"), w.entry, pady=10, padx =30)

   

  
   tkgrid(submit.but, reset.but)

  tkwait.window(tt)
  return(c(x,y,z,w))
}
#Now run the function like:
predictor_para <- inputs()
print(predictor_para)
mat<-eval(parse(text=predictor_para[1]))
ycol<-eval(parse(text=predictor_para[2]))
xcol<-eval(parse(text=predictor_para[3]))
idcol<-eval(parse(text=predictor_para[4]))
multitest(mat,ycol,xcol,idcol)
}
FuninT.pck <-
c("FuninT.pck", "Fun.in.Tstats")
Fun.in.Tstats <-
function(Bx=.1,By=1,n=50){
x1<-rnorm(n)
x2<-(x1*Bx+rnorm(n))
x2<-x2/(sd(x2))
y<-x1+x2+rnorm(n)+1
par(mfrow=c(2,2))
plot(x1,x2,main=paste("cor=",cor(x1,x2)))
lstr<-lsfit(cbind(x1,x2),y)
print(lstr$coef)
yhat<-lstr$coef[1]+cbind(x1,x2)%*%c(lstr$coef[-1])
mse<-sum(lstr$resid^2)/(length(y)-3)
plot(yhat,y)
plot(yhat,lstr$resid)
Press1<-lspluspress(cbind(x1,x2),y)$PRESS
Press2<-lspluspress(x2,y)$PRESS
Press3<-lspluspress(x1,y)$PRESS
ls.print((lstr))
ls.print(lsfit(x2,y))
print(paste("PRESSALL=",Press1,"PRESS1=",Press3,"PRESS2=",Press2))
}
