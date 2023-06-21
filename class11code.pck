class11code.pck <-
c("class11code.pck", "my.mvknorm")
my.mvknorm <-
function(x0,k,thresh.stop=1e-9,do.scale=F){
if(do.scale){
x<-scale(x0)
}else{
x<-x0
}
n<-length(x[,1])
k1<-length(x[1,])

pmat<-NULL
#minit:
muvec<-NULL
sigvec<-list()
for(i in 1:k){
muvec<-rbind(muvec,x[i,])
sigvec[[i]]<-diag(k1)
}
musubt<-function(vec){
        vec-mu0
        }

minfunx<-function(z){
nz<-length(z)
v1<-NULL
for(i in 1:nz){
	v1<-c(v1,eigen(z[[i]])$val)
}
min(v1)
}
likdif<-1e10
likstart<-(-1e10)
pivec<-rep(1/k,k)
while((abs(likdif)>thresh.stop)&&(minfunx(sigvec)>.0007)){
#Estepn
pmat<-NULL
Nvec<-NULL
for(i in 1:k){
pvec0<-pivec[i]*dmvnorm(x,muvec[i,],sigvec[[i]])
pmat<-rbind(pmat,c(pvec0))
}
piden<-apply(pmat,2,sum)
for(i in 1:k){
pmat[i,]<-pmat[i,]/piden
Nvec<-c(Nvec,sum(pmat[i,]))
}
#print(pmat)
#Mstep
for(i in 1:k){
muvec[i,]<-apply(pmat[i,]%*%x,2,sum)/Nvec[i]
mu0<-muvec[i,]
musubt<-function(vec){
        vec-mu0
        }

z9<-apply(x,1,musubt)
sigvec[[i]]<-((z9)%*%diag(pmat[i,])%*%t(z9))/Nvec[i]
#print(sigvec[[i]])
}
#Evalstep
liknewa<-NULL
for(i in 1:k){
pi0<-sum(pmat[i,])/(sum(Nvec))
liknewa<-rbind(liknewa,(pi0*dmvnorm(x,muvec[i,],sigvec[[i]])))
pivec[i]<-pi0
}
liknew0<-apply(liknewa,2,sum)
liknew<-sum(log(liknew0))
likdif<-liknew-likstart
#print(muvec)
print(likdif)
likstart<-liknew
#print(sigvec)
print(liknew)
}
pmat1<-t(pmat)
z<-rep(0,length(pmat1[,1]))
for(i in 1:k){
	z[pmat1[,i]>.5]<-i
}
x99<-cbind(z,x)
doh<-Fisher.disc.space(x99,1,k)
list(mu=muvec,sig=sigvec,pivec=pivec,lik=likstart,post=pmat,cluster=z)
}
