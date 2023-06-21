code467day6.pck <-
c("code467day6.pck", "Manova.calculation.multiway.test", "Perm.calc.multiwayManova", 
"gui.multiway.manova.test.int", "resid.meanplot", "Residplot1", 
"I620", "mquiz")
Manova.calculation.multiway.test <-
function(dat0,idcol,int.str){
#list(lambda=mdet(SSW)/mdet(SSW+SSB),SSW=SSW,SSB=SSB,SST=SST,mumat=mean.mat)
#This program calculates SS for each factor and combination of factors for all terms listed in the int.str list
# int.str used to identify interactions to study, including 1st order.
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
SSBtot1<-0
SSBtot2<-0
SSBtot3<-0
dat00<-convert.data(dat0,idcol)
# the assumption is that the interactions in int.str go from lowest to highest So all previous SSB must be subtracted from SST 
for(j in 1:n1){
	idcoltemp<-idcol[int.str[[j]]]
throw.away<-(idcol[-int.str[[j]]])
	 dum2<-Manova.calculation.multiway.test0(dat00,idcoltemp,rid=throw.away)
if(length(int.str[[j]])==1){
	lambdavec<-c(lambdavec,mdet(dum1$SSE)/mdet(dum1$SSE+dum2$SSB))
      clambdavec<-c(clambdavec,det(dum1$SSE)/det(dum1$SSE+dum2$SSB))
	out.str[[(j+1)]]<-list(int=int.str[[j]],raw=dum2)
	SSElist[[(j+1)]]<-(dum2$SSE)
	SSBlist[[(j+1)]]<-(dum2$SSB)
	SSTlist[[(j+1)]]<-dum2$SST
	SSBtot1<-SSBtot1+dum2$SSB
	mulist[[(j+1)]]<-dum2$mean.mat
      residlist[[(j+1)]]<-dum2$residmat
	muAlist[[(j+1)]]<-dum2$mumatA

	
}else{
	if( length(int.str[[j]])==2){
		vec0<-int.str[[j]]
		jvec<-vec0+1
		SSB<-dum2$SST-SSBlist[[jvec[1]]]-SSBlist[[jvec[2]]]-dum2$SSE
		lambdavec<-c(lambdavec,mdet(dum1$SSE)/mdet(dum1$SSE+SSB))
		clambdavec<-c(clambdavec,det(dum1$SSE)/det(dum1$SSE+SSB))
		out.str[[(j+1)]]<-list(int=int.str[[j]],raw=dum2)
		SSElist[[(j+1)]]<-(dum2$SSE)
		SSBlist[[(j+1)]]<-SSB
		SSTlist[[(j+1)]]<-dum2$SST
      	SSBtot2<-SSBtot2+SSB
		mulist[[(j+1)]]<-dum2$mean.mat
            residlist[[(j+1)]]<-dum2$residmat
		muAlist[[(j+1)]]<-dum2$mumatA

	}else{
		if( length(int.str[[j]])==3){
			jvec<-int.str[[j]]+1+lenght(idcol)*(length(idcol)-1)/2
			SSB<-dum2$SST-dum2$SSE
			for(k in 1:3){
			SSB<-SSB-SSBlist[[jvec[k]]]
			}
			lambdavec<-c(lambdavec,mdet(dum1$SSE)/mdet(dum1$SSE+SSB))
			clambdavec<-c(clambdavec,det(dum1$SSE)/det(dum1$SSE+SSB))
			out.str[[(j+1)]]<-list(int=int.str[[j]],raw=dum2)
			SSElist[[(j+1)]]<-(dum2$SSE)
			SSBlist[[(j+1)]]<-SSB
			SSTlist[[(j+1)]]<-dum2$SST
      		SSBtot3<-SSBtot3+SSB
			mulist[[(j+1)]]<-dum2$mean.mat
        		residlist[[(j+1)]]<-dum2$residmat
			muAlist[[(j+1)]]<-dum2$mumatA

		}
	}
}
}
	SSBtottot<-SSBtot1+SSBtot2+SSBtot3
lambda.int<-mdet(dum1$SSE)/mdet(dum1$SST-SSBtottot)
clambda.int<-det(dum1$SSE)/det(dum1$SST-SSBtottot)
out<-list(Lambda=c(lambdavec,lambda.int),cLambda=c(clambdavec,clambda.int),outinf=out.str,SSE=SSElist,SSB=SSBlist,SST=SSTlist,mulist=mulist,muAlist=muAlist,residlist=residlist)
out
}
Perm.calc.multiwayManova <-
function(mat,idcol,int.str,nperm=10000){
#function(dat0,idcol,int.str)
#NOTE IDCOL SHOULD BE IN ORDER THEY APPEAR IN MATRIX
dum<-Manova.calculation.multiway.test(mat,idcol,int.str)
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
	dumP<-Manova.calculation.multiway.test(matp,c(1:n0),int.str)
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
gui.multiway.manova.test.int <-
function(){
library(tcltk)
inputs <- function(){

   x <- tclVar("T620")
   y <- tclVar("c(1:3)")
   z<-tclVar("I620")
     w<-tclVar("10000")
   tt <- tktoplevel()
   tkwm.title(tt,"Choose parameters for new function                   ")
   x.entry <- tkentry(tt, textvariable=x)
   y.entry <- tkentry(tt, textvariable=y)
   z.entry<-tkentry(tt,textvariable=z)
     w.entry<-tkentry(tt,textvariable=w)
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
	z<-tclvalue(z)
         w<-tclvalue(w)
           e <- parent.env(environment())
     e$x <- x
     e$y <- y
     e$z<-z
     e$w<-w
        tkdestroy(tt)
   }

   submit.but <- tkbutton(tt, text="start", command=submit)
   tkgrid(tklabel(tt,text="Input data matrix"),columnspan=2)
   tkgrid(tklabel(tt,text="data"), x.entry, pady=10, padx =30)


   tkgrid(tklabel(tt,text="Input id column matrix"),columnspan=2)
   tkgrid(tklabel(tt,text="idcol"), y.entry, pady=10, padx =30)

  
   tkgrid(tklabel(tt,text="Input effects list"),columnspan=2)
   tkgrid(tklabel(tt,text="intlist"), z.entry, pady=10, padx =30)


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
int.str<-eval(parse(text=predictor_para[3]))
nperm<-eval(parse(text=predictor_para[4]))

Perm.calc.multiwayManova(mat,idcol,int.str,nperm)

}
resid.meanplot <-
function(str,mat,idcol,i,varsel=0)
{
if(varsel[1]==0){
	n1<-length(str$resid[[i]][1,])
	varsel<-c(1:n1)
}
	m1<-str$resid[[i]][,c(1,varsel)]
	m3<-str$muAlist[[i]][,c(1,varsel)]
	ma<-cbind(m1,m3[,-1]+m1[,-1])
	m2<-convert.data(mat,idcol)
	m0<-cbind(m2[,idcol],ma[,-1])
	pairs(m0,pch=as.character(ma[,1]),col=ma[,1])
}
Residplot1 <-
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
}
I620 <-
list(1, 2, 3, c(1, 2), c(1, 3), c(2, 3))
mquiz <-
structure(list(`c(rep(1, 48), rep(0, 47), rep(1, 23), rep(0, 26))` = c(1, 
1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), origin = c(1, 
1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3), mpg = c(17, 
15, 21, 10, 19, 14, 13, 18, 21, 14, 13, 22, 14, 13, 18, 12, 19, 
26, 16, 14, 15, 14, 13, 23, 18, 14.5, 24, 18, 26.5, 13, 18.5, 
16, 25.5, 36.100000000000001, 20.199999999999999, 25.100000000000001, 
18.100000000000001, 22.300000000000001, 17.600000000000001, 15.5, 
26.399999999999999, 34.399999999999999, 29.899999999999999, 20.199999999999999, 
24, 22, 27, 28, 18, 14, 15, 10, 28, 17, 14, 22, 20, 13, 13, 13, 
13, 16, 11, 21, 15, 15, 25, 16, 18, 15, 16, 20, 17.5, 22, 29, 
24.5, 13, 17.5, 17.5, 15.5, 30.5, 20.5, 19.199999999999999, 30, 
17, 19.199999999999999, 28.800000000000001, 32.100000000000001, 
28, 23.5, 28, 27, 34, 27, 31, 27, 27, 35, 28, 22, 31, 32, 24, 
32, 33.5, 22, 39.399999999999999, 27.199999999999999, 23.899999999999999, 
31.800000000000001, 37.200000000000003, 37, 40.799999999999997, 
32.700000000000003, 37, 32.899999999999999, 24.199999999999999, 
38, 24, 25, 31, 24, 23, 27, 20, 20, 32, 31, 29, 24, 28, 19, 26, 
27.5, 21.100000000000001, 38.100000000000001, 29.800000000000001, 
32.200000000000003, 39.100000000000001, 37.700000000000003, 32.399999999999999, 
25.399999999999999, 34, 32), cylinders = c(8, 8, 6, 8, 6, 8, 
8, 6, 4, 8, 8, 4, 8, 8, 6, 8, 4, 4, 8, 8, 6, 8, 8, 4, 6, 8, 6, 
6, 4, 8, 6, 8, 4, 4, 6, 4, 8, 4, 8, 8, 4, 4, 4, 6, 4, 6, 4, 4, 
8, 8, 8, 8, 4, 6, 8, 4, 4, 8, 8, 8, 8, 6, 8, 4, 8, 6, 4, 6, 6, 
8, 6, 8, 8, 6, 4, 4, 8, 8, 6, 8, 4, 6, 8, 4, 8, 8, 6, 4, 4, 6, 
4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4, 4, 4, 
4, 4, 6, 4, 4, 6, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 4, 4, 4, 6, 
4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4), displacement = c(302, 429, 
200, 360, 250, 351, 400, 250, 122, 351, 302, 122, 302, 351, 250, 
400, 122, 122, 302, 302, 250, 351, 302, 140, 171, 351, 200, 250, 
140, 302, 250, 351, 140, 98, 200, 140, 302, 140, 302, 351, 140, 
98, 98, 200, 140, 232, 140, 120, 307, 454, 400, 307, 140, 250, 
350, 140, 140, 350, 307, 350, 400, 250, 400, 140, 350, 250, 140, 
250, 250, 350, 250, 262, 305, 250, 85, 98, 350, 305, 250, 350, 
98, 200, 305, 98, 305, 267, 173, 98, 151, 173, 112, 112, 112, 
151, 119, 97, 97, 72, 97, 108, 79, 83, 119, 85, 85, 146, 85, 
119, 119, 85, 86, 119, 85, 168, 85, 119, 146, 91, 113, 113, 71, 
113, 120, 97, 97, 156, 71, 76, 97, 134, 97, 156, 97, 134, 134, 
89, 134, 108, 79, 89, 108, 168, 108, 144), horsepower = c(140, 
198, 85, 215, 88, 153, 170, 88, 86, 153, 140, 86, 137, 158, 88, 
167, 85, 80, 140, 140, 72, 148, 129, 83, 97, 152, 81, 78, 72, 
130, 98, 149, 89, 66, 85, 88, 139, 88, 129, 142, 88, 65, 65, 
88, 92, 112, 86, 79, 130, 220, 150, 200, 90, 100, 165, 72, 90, 
165, 130, 145, 150, 100, 150, 72, 145, 100, 75, 100, 105, 145, 
105, 110, 140, 105, 52, 60, 145, 145, 110, 170, 63, 95, 145, 
68, 130, 125, 115, 70, 90, 110, 88, 88, 88, 90, 82, 88, 88, 69, 
92, 94, 67, 61, 97, 70, 70, 97, 70, 97, 97, 65, 65, 92, 65, 132, 
65, 100, 120, 67, 95, 95, 65, 95, 97, 88, 88, 122, 65, 52, 75, 
96, 75, 108, 75, 95, 95, 60, 90, 75, 58, 62, 75, 116, 70, 96), 
    weight = c(3449, 4341, 2587, 4615, 3302, 4154, 4746, 3139, 
    2226, 4129, 4294, 2395, 4042, 4363, 3021, 4906, 2310, 2451, 
    4141, 4638, 3158, 4657, 3169, 2639, 2984, 4215, 3012, 3574, 
    2565, 3870, 3525, 4335, 2755, 1800, 2965, 2720, 3205, 2890, 
    3725, 4054, 2870, 2045, 2380, 3060, 2865, 2835, 2790, 2625, 
    3504, 4354, 3761, 4376, 2264, 3329, 4209, 2408, 2408, 4274, 
    4098, 3988, 4464, 3278, 4997, 2401, 4082, 3336, 2542, 3781, 
    3459, 4440, 3897, 3221, 4215, 3353, 2035, 2164, 4055, 3880, 
    3520, 4165, 2051, 3155, 3425, 2155, 3840, 3605, 2595, 2120, 
    2678, 2725, 2605, 2640, 2395, 2950, 2720, 2130, 2130, 1613, 
    2288, 2379, 1950, 2003, 2545, 1990, 1945, 2815, 2070, 2300, 
    2405, 2020, 2019, 2434, 2110, 2910, 1975, 2615, 2930, 1995, 
    2372, 2228, 1773, 2278, 2506, 2100, 2279, 2807, 1836, 1649, 
    2171, 2702, 2155, 2930, 2265, 2560, 2515, 1968, 2711, 2265, 
    1755, 2050, 2350, 2900, 2245, 2665), acceleration = c(10.5, 
    10, 16, 14, 15.5, 13.5, 12, 14.5, 16.5, 13, 16, 16, 14.5, 
    13, 16.5, 12.5, 18.5, 16.5, 14, 16, 19.5, 13.5, 12, 17, 14.5, 
    12.800000000000001, 17.600000000000001, 21, 13.6, 15, 19, 
    14.5, 15.800000000000001, 14.4, 15.800000000000001, 15.4, 
    11.199999999999999, 17.300000000000001, 13.4, 14.300000000000001, 
    18.100000000000001, 16.199999999999999, 20.699999999999999, 
    17.100000000000001, 16.399999999999999, 14.699999999999999, 
    15.6, 18.600000000000001, 12, 9, 9.5, 15, 15.5, 15.5, 12, 
    19, 19.5, 12, 14, 13, 12, 18, 14, 19.5, 13, 17, 17, 17, 16, 
    14, 18.5, 13.5, 13, 14.5, 22.199999999999999, 22.100000000000001, 
    12, 12.5, 16.399999999999999, 11.4, 17, 18.199999999999999, 
    13.199999999999999, 16.5, 15.4, 15, 11.300000000000001, 15.5, 
    16.5, 12.6, 19.600000000000001, 18.600000000000001, 18, 17.300000000000001, 
    19.399999999999999, 14.5, 14.5, 18, 17, 16.5, 19, 19, 17, 
    17, 16.800000000000001, 14.5, 18.600000000000001, 14.699999999999999, 
    14.9, 19.199999999999999, 16.399999999999999, 15, 19.199999999999999, 
    11.4, 19.399999999999999, 14.800000000000001, 13.800000000000001, 
    16.199999999999999, 15, 14, 19, 15.5, 14.5, 16.5, 19, 13.5, 
    21, 16.5, 16, 13.5, 16.399999999999999, 15.5, 18.199999999999999, 
    14.199999999999999, 14.800000000000001, 18.800000000000001, 
    15.5, 15.199999999999999, 16.899999999999999, 17.300000000000001, 
    16.800000000000001, 12.6, 16.899999999999999, 13.9), year = c(70, 
    70, 70, 70, 71, 71, 71, 71, 72, 72, 72, 72, 73, 73, 73, 73, 
    73, 74, 74, 74, 75, 75, 75, 75, 75, 76, 76, 76, 76, 76, 77, 
    77, 77, 78, 78, 78, 78, 79, 79, 79, 80, 81, 81, 81, 82, 82, 
    82, 82, 70, 70, 70, 70, 71, 71, 71, 71, 72, 72, 72, 73, 73, 
    73, 73, 73, 73, 74, 74, 74, 75, 75, 75, 75, 76, 76, 76, 76, 
    76, 77, 77, 77, 77, 78, 78, 78, 79, 79, 79, 80, 80, 81, 82, 
    82, 82, 82, 82, 70, 71, 71, 72, 73, 74, 74, 75, 76, 77, 77, 
    78, 78, 78, 79, 80, 80, 80, 80, 81, 81, 81, 82, 70, 71, 71, 
    72, 72, 72, 73, 73, 74, 74, 75, 75, 76, 76, 77, 78, 78, 80, 
    80, 80, 81, 81, 81, 81, 82, 82)), class = "data.frame", row.names = c("5", 
"6", "18", "26", "37", "41", "44", "49", "62", "66", "75", "81", 
"89", "93", "101", "105", "113", "131", "137", "140", "156", 
"160", "167", "169", "175", "191", "194", "201", "207", "215", 
"229", "233", "237", "246", "255", "256", "265", "283", "287", 
"291", "315", "352", "353", "366", "374", "389", "393", "396", 
"1", "7", "13", "27", "31", "36", "39", "47", "61", "63", "74", 
"88", "92", "99", "104", "110", "116", "129", "133", "134", "154", 
"158", "162", "166", "188", "193", "196", "197", "214", "222", 
"226", "231", "238", "254", "263", "267", "286", "292", "307", 
"312", "314", "342", "368", "369", "370", "392", "397", "19", 
"30", "55", "82", "111", "130", "146", "174", "205", "221", "242", 
"248", "269", "274", "304", "313", "321", "325", "334", "348", 
"358", "363", "385", "15", "32", "54", "58", "83", "85", "109", 
"124", "132", "145", "168", "172", "206", "211", "236", "268", 
"271", "311", "319", "322", "344", "349", "357", "362", "382", 
"390"))
