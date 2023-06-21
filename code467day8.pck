code467day8.pck <-
c("code467day8.pck", "gui.princomp", "eigenvalue.inference")
gui.princomp <-
function(){
library(tcltk)
inputs <- function(){
#function(dat,alpha=.05,scaled=F)
   x <- tclVar("T620")
   y <- tclVar(".05")
   z <- tclVar("F")
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


   tkgrid(tklabel(tt,text="alpha"),columnspan=2)
   tkgrid(tklabel(tt,text="alpha"), y.entry, pady=10, padx =30)

   tkgrid(tklabel(tt,text="Scale cov?"),columnspan=2)
   tkgrid(tklabel(tt,text="F"), z.entry, pady=10, padx =30)

tkgrid(tklabel(tt,text="Nboot"),columnspan=2)
   tkgrid(tklabel(tt,text="10000"), w.entry, pady=10, padx =30)


 
   tkgrid(submit.but, reset.but)

  tkwait.window(tt)
  return(c(x,y,z,w))
}
#Now run the function like:
predictor_para <- inputs()
print(predictor_para)
mat<-eval(parse(text=predictor_para[1]))
alpha<-eval(parse(text=predictor_para[2]))
xscale<-eval(parse(text=predictor_para[3]))
nboot<-eval(parse(text=predictor_para[4]))

eigenvalue.inference(mat,alpha,xscale,nboot)

}
eigenvalue.inference <-
function(dat,alpha=.05,scaled=F,nboot=10000){
v1<-var(dat)
n1<-length(dat[,1])
if(scaled){
v1<-cor(dat)
}
e1<-eigen(v1)
p1<-length(e1$val)
val0<-e1$val
vec0<-e1$vec
alpha0<-alpha/2
alphabon<-alpha0/p1
z1<-qnorm(alpha0)
zbon<-qnorm(alphabon)
denL1<-1-z1*sqrt(2/n1)
denLbon1<-1-zbon*sqrt(2/n1)
denL2<-1+z1*sqrt(2/n1)
denLbon2<-1+zbon*sqrt(2/n1)
bootmat<-NULL
bootmatd<-NULL
for(i in 1:nboot){
	if((i/500)==floor(i/500)){print(i)}
	vbn<-sample(n1,replace=T)
	bdat<-dat[vbn,]
	vb1<-var(bdat)
	if(scaled){vb1<-cor(bdat)}
	eb1<-eigen(vb1)
	valb<-eb1$val
	bootmat<-cbind(bootmat,val0/valb)
	bootmatd<-cbind(bootmatd,diff(valb))
}
my.quantile<-function(x){quantile(x,c(alphabon,alpha0,1-alpha0,1-alphabon))}
vz<-apply(bootmat,1,my.quantile)
vsd<-sqrt((2/n1)*(val0[-1]^2+val0[-p1]^2))
vt<-diff(val0)/vsd
Bootstrap<-t(vz)*val0
Normal<-cbind(val0/denLbon1,val0/denL1,val0,val0/denL2,val0/denLbon2)
vec.mat.list<-list()
for(i in 1:p1){
	lambda0<-val0[i]
	E0<-0
	for(j in 1:p1){
		if(j!=i){
		E0<-E0+(lambda0*val0[j]/(sqrt(n1)*((lambda0-val0[j])^2)))*vec0%*%t(vec0)
		}
	}
	vec.mat.list[[i]]<-E0
}
	
list(eigen=e1,test.notunique=pt(vt,n1-1),bootstrap.val=Bootstrap,normal.val=Normal,proportion.var=rbind(val0,cumsum(val0)/sum(val0)),vec.cov=vec.mat.list)
}
