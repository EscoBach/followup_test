
Gene.data<-function(n,K,rho,para.true,ep.mean,ep.var,cen){
  epsi<-rmvnorm(n,ep.mean,ep.var)
  epsi<-as.vector(t(epsi))
  Z.1<-rnorm(n*K,0,1)
  Z.2<-rbinom(n*K,1,0.5); Z.3<-runif(n*K,0,1)
  Z<-cbind(Z.1,Z.2,Z.3)
  para.true<-as.matrix(para.true)
  log.T<-as.vector(t(para.true)%*%t(Z))+epsi
  T<-exp(log.T)
  Cen<-runif(n*K,0,cen)
  Y<-pmin(T,Cen)
  delta<-I(T<=Cen)*1
  id<-rep(1:n,each=K)
  Data<-data.frame(id,delta,Z,Y)
  return(Data)
}

aftGmm.ext<-function(seed,n=mm,K=KK,rho=rhoo,cen=cenn,
                     ep.mean=ep.meann,ep.var=ep.varr,para.true=para.truee){
  set.seed(seed)
  data1<-Gene.data(n,K,rho,para.true,ep.mean,ep.var,cen)
  
  ###gehan estimation 
  out<-aftgee(Surv(Y, delta)~0+Z.2+Z.3,id = id,data = data1,corstr = "ind",B = 2, binit = "lm")$coef.res
  
  return(round(out,5))
}


pw.Survf<-function(s,para,ZZ,TT,delta,n,K)
{
  e<-as.vector(log(TT)-ZZ%*%para)
  ord<-order(e)
  ei<-e[ord]
  deltai<-delta[ord]
  ei.fai<-unique(ei[deltai==1],fromLast=TRUE)
  riskn<-function(x)
  {
    dimv<-min(which(ei==x)):(K*n)
    riskn<-length(dimv)
    return(riskn)
  }
  riskv<-sapply(ei.fai,FUN = riskn)
  riskd<-riskv/(riskv+1)
  in.ei<-I(ei.fai<=s)*1
  F.pw<-prod(riskd[in.ei==1])
  return(F.pw)
}

phi<-function(u,para,ZZ,TT,delta,Weight,n,K)
{
  e<-as.vector(log(TT)-ZZ%*%para)
  if(Weight=="logrank"){obj<-1}
  if(Weight=="gehan"){obj<-n^(-1)*sum(e>=u)}
  if(Weight=="pw"){obj<-pw.Survf(u,para,ZZ,TT,delta,n,K)}
  return(obj)
}

C_l<-function(M,n,K)
{
  Matri<-matrix(0,K*n,K*n)
  for(i in 1:n)
  {
    Matri[((i-1)*K+1):(i*K),((i-1)*K+1):(i*K)]<-M
  }
  return(Matri)
}


Phi.f<-function(para,ZZ,TT,delta,Weight,n,K){
  
  C.11<-diag(1,K); C.12<-matrix(1,K,K)-diag(1,K)
  
  e<-as.vector(log(TT)-ZZ%*%para)
  # I.e: matrix(I(e>=e.11),,I(e>=e.nK)), the first col is I(e>=e.11). matrix
  I.e<-matrix(unlist(lapply(e,function(u){I(e>=u)*1})), nrow=n*K,ncol=n*K) 
  # S0.e: (S(e.11),,S(e.nK)). vector
  S0.e<-sapply(e,function(u){sum(e>=u)})
  #phi_e: (w(e.11),,w(e.nK)). vector
  phi_e<-sapply(e,function(u){phi(u,para,ZZ,TT,delta,Weight,n,K)})
  
  Eq.1<-apply(t(t(I.e)*(phi_e*delta*S0.e^(-1))),1,sum)
  
  Sn.11<-n^(-1)*(t(ZZ)%*%C_l(C.11,n,K)%*%(phi_e*delta-Eq.1))
  Sn.12<-n^(-1)*(t(ZZ)%*%C_l(C.12,n,K)%*%(phi_e*delta-Eq.1))
  Sn<-rbind(Sn.11,Sn.12)
  
  return(Sn)
}

Covari.f<-function(para,d,ZZ,TT,delta,Weight,n,K){
  
  C.11<-diag(1,K); C.12<-matrix(1,K,K)-diag(1,K)
  
  e<-as.vector(log(TT)-ZZ%*%para)
  # I.e: matrix(I(e>=e.11),,I(e>=e.nK)), the first col is I(e>=e.11). matrix
  I.e<-matrix(unlist(lapply(e,function(u){I(e>=u)*1})), nrow=n*K,ncol=n*K) 
  # S0.e: (S(e.11),,S(e.nK)). vector
  S0.e<-sapply(e,function(u){sum(e>=u)})
  #phi_e: (w(e.11),,w(e.nK)). vector
  phi_e<-sapply(e,function(u){phi(u,para,ZZ,TT,delta,Weight,n,K)})
  
  Eq.1<-apply(t(t(I.e)*(phi_e*delta*S0.e^(-1))),1,sum)
  
  Sn.11<-n^(-1)*(t(ZZ)%*%C_l(C.11,n,K)%*%(phi_e*delta-Eq.1))
  Sn.12<-n^(-1)*(t(ZZ)%*%C_l(C.12,n,K)%*%(phi_e*delta-Eq.1))
  Sn<-rbind(Sn.11,Sn.12)-d
  
  phii<-function(MM){
    
    # a.matrix: (a11(e11),,anK(e11),,a11(enK),,anK(enK)) matrix
    a.matrix<-t(t(matrix(as.vector(t(ZZ)%*%C_l(MM,n,K)),dim(ZZ)[2],n*K*n*K))*rep(phi_e,each=n*K))
    
    I.e.vector<-as.vector(I.e)
    Eq.4<-t(t(a.matrix)*I.e.vector)
    
    # a.bar: (a.bar(e11),,a.bar(enK)) matrix
    a.bar<-NULL
    for(i in 1:(n*K)){
      a.bar<-cbind(a.bar,apply(Eq.4[,((i-1)*n*K+1):(i*n*K)],1,sum)/S0.e[i])
    }
    
    Eq.5<-t(t(a.bar)*delta*S0.e^(-1))
    
    Eq.2<-NULL
    for(i in 1:n){
      
      Eq.3<-t(ZZ)[,((i-1)*K+1):(i*K)]%*%MM%*%
        (phi_e[((i-1)*K+1):(i*K)]*delta[((i-1)*K+1):(i*K)]-Eq.1[((i-1)*K+1):(i*K)])
      
      Eq.6<-apply(t(t(a.bar[,((i-1)*K+1):(i*K)])*delta[((i-1)*K+1):(i*K)]),1,sum)-
        apply(t(t(Eq.5)*apply(I.e[((i-1)*K+1):(i*K),],2,sum)),1,sum)
      
      phi.i<-Eq.3-Eq.6
      Eq.2<-cbind(Eq.2,phi.i)
    }
    
    return(Eq.2)
  }
  
  phiall<-rbind(phii(C.11),phii(C.12))-d
  
  Sig1<-0
  for(i in 1:n){
    Sig1<-Sig1+phiall[,i]%*%t(phiall[,i])
  }
  
  out<-n^(-1)*Sig1-Sn%*%t(Sn)
  
  return(out)
}



Sn.f<-function(para,ZZ,TT,delta,Weight,n,K,MA){
  
  e<-as.vector(log(TT)-ZZ%*%para)
  # I.e: matrix(I(e>=e.11),,I(e>=e.nK)), the first col is I(e>=e.11). matrix
  I.e<-matrix(unlist(lapply(e,function(u){I(e>=u)*1})), nrow=n*K,ncol=n*K) 
  # S0.e: (S(e.11),,S(e.nK)). vector
  S0.e<-sapply(e,function(u){sum(e>=u)})
  #phi_e: (w(e.11),,w(e.nK)). vector
  phi_e<-sapply(e,function(u){phi(u,para,ZZ,TT,delta,Weight,n,K)})
  
  Eq.1<-apply(t(t(I.e)*(phi_e*delta*S0.e^(-1))),1,sum)
  
  out<-n^(-1)*(t(ZZ)%*%C_l(MA,n,K)%*%(phi_e*delta-Eq.1))
  
  return(out)
}

u.f<-function(para,ZZ,TT,delta,Weight,n,K,MA){
  
  e<-as.vector(log(TT)-ZZ%*%para)
  # I.e: matrix(I(e>=e.11),,I(e>=e.nK)), the first col is I(e>=e.11). matrix
  I.e<-matrix(unlist(lapply(e,function(u){I(e>=u)*1})), nrow=n*K,ncol=n*K) 
  # S0.e: (S(e.11),,S(e.nK)). vector
  S0.e<-sapply(e,function(u){sum(e>=u)})
  #phi_e: (w(e.11),,w(e.nK)). vector
  phi_e<-sapply(e,function(u){phi(u,para,ZZ,TT,delta,Weight,n,K)})
  
  Eq.1<-apply(t(t(I.e)*(phi_e*delta*S0.e^(-1))),1,sum)
  
  phii<-function(MM){
    
    # a.matrix: (a11(e11),,anK(e11),,a11(enK),,anK(enK)) matrix
    a.matrix<-t(t(matrix(as.vector(t(ZZ)%*%C_l(MM,n,K)),dim(ZZ)[2],n*K*n*K))*rep(phi_e,each=n*K))
    
    I.e.vector<-as.vector(I.e)
    Eq.4<-t(t(a.matrix)*I.e.vector)
    
    # a.bar: (a.bar(e11),,a.bar(enK)) matrix
    a.bar<-NULL
    for(i in 1:(n*K)){
      a.bar<-cbind(a.bar,apply(Eq.4[,((i-1)*n*K+1):(i*n*K)],1,sum)/S0.e[i])
    }
    
    Eq.5<-t(t(a.bar)*delta*S0.e^(-1))
    
    Eq.2<-NULL
    for(i in 1:n){
      
      Eq.3<-t(ZZ)[,((i-1)*K+1):(i*K)]%*%MM%*%
        (phi_e[((i-1)*K+1):(i*K)]*delta[((i-1)*K+1):(i*K)]-Eq.1[((i-1)*K+1):(i*K)])
      
      Eq.6<-apply(t(t(a.bar[,((i-1)*K+1):(i*K)])*delta[((i-1)*K+1):(i*K)]),1,sum)-
        apply(t(t(Eq.5)*apply(I.e[((i-1)*K+1):(i*K),],2,sum)),1,sum)
      
      phi.i<-Eq.3-Eq.6
      Eq.2<-cbind(Eq.2,phi.i)
    }
    
    return(Eq.2)
  }
  
  out<-phii(MA)
  
  
  return(out)
}

CovariAI.f<-function(para,d,ZZ,ZZAI,TT,delta,Weight,n,K,gamma){
  
  C.11<-diag(1,K); C.12<-matrix(1,K,K)-diag(1,K)
  
  u1<-u.f(para=para,ZZ=ZZ,TT=TT,delta=delta,Weight=Weight,n,K,MA=C.11)
  u2<-u.f(para=para,ZZ=ZZ,TT=TT,delta=delta,Weight=Weight,n,K,MA=C.12)
  u3<-u.f(para=gamma,ZZ=ZZAI,TT=TT,delta=delta,Weight=Weight,n,K,MA=C.11)
  u<-rbind(u1,u2,u3)-d
  U1<-Sn.f(para=para,ZZ=ZZ,TT=TT,delta=delta,Weight=Weight,n,K,MA=C.11)
  U2<-Sn.f(para=para,ZZ=ZZ,TT=TT,delta=delta,Weight=Weight,n,K,MA=C.12)
  U3<-Sn.f(para=gamma,ZZ=ZZAI,TT=TT,delta=delta,Weight=Weight,n,K,MA=C.11)
  U<-rbind(U1,U2,U3)-d
  out<--U%*%t(U)
  for(s in 1:n){
    out<-out+n^(-1)*u[,s]%*%t(u[,s])
  }
  
  return(out)
}



GMM<-function(seed,n=nn,K=KK,rho=rhoo,cen=cenn,
              ep.mean=ep.meann,ep.var=ep.varr,
              Weight=Weightt,para.true=para.truee,gamma.e=gamma.in){
  set.seed(seed*4-1)
  
  data1<-Gene.data(n,K,rho,para.true,ep.mean,ep.var,cen)
  Y<-data1$Y;delta<-data1$delta;id<-data1$id;X<-cbind(data1$Z.1,data1$Z.2,data1$Z.3)
  Z<-cbind(data1$Z.2,data1$Z.3)
  
  
  C.11<-diag(1,K); C.12<-matrix(1,K,K)-diag(1,K)
  
  #############The initial value
  
  para.init<-aftgee(Surv(Y,delta)~0+Z.1+Z.2+Z.3,id =id,data=data1,corstr="ind",B =2,binit = "lm")$coef.res
  
  
  ##############################GIF Li and Yin 2009
  
  ###############Empirical variance
  
  Covari<-Covari.f(para=para.init,d=0,ZZ=X,TT=Y,delta=delta,Weight=Weight,n,K)
  
  QIF<-function(para1){
    Phi<-Phi.f(para=para1,ZZ=X,TT=Y,delta=delta,Weight=Weight,n,K)
    Qn<-t(Phi)%*%solve(Covari)%*%Phi
    return(Qn)
  }
  est.QIF<-optim(par = para.init,fn=QIF,
                 method="Nelder-Mead",control = list(maxit=200,ndeps=0.01))$par
  
  #####Variance estimation
  
  CD<-t(chol(n^(-1)*Covari.f(para=est.QIF,d=0,ZZ=X,TT=Y,delta=delta,Weight=Weight,n,K)))
  
  QIF.tilde<-function(para1,dd,para2){
    Phi<-Phi.f(para=para1,ZZ=X,TT=Y,delta=delta,Weight=Weight,n,K)-dd
    Cova<-Covari.f(para=para2,d=dd,ZZ=X,TT=Y,delta=delta,Weight=Weight,n,K)
    obj.fun<-t(Phi)%*%solve(Cova)%*%Phi
  }
  
  eta<-NULL
  for(s in 1:dim(CD)[1]){
    est.tilde<-optim(par = para.init,
                     fn=function(para){QIF.tilde(para1=para,dd=CD[,s],para2=para.init)},
                     method="Nelder-Mead",control = list(maxit=200,ndeps=0.01))$par
    
    eta<-cbind(eta,est.tilde)
  }
  eta<-eta-est.QIF; Cov.QIF<-eta%*%t(eta)
  se.QIF<-as.vector((diag(Cov.QIF))^(1/2)) #variance for GMM
  
  
  ##################### GIFAI auxiliary information
  
  CovariAI<-CovariAI.f(para=para.init,
                       d=0,ZZ=X,ZZAI=Z,TT=Y,delta=delta,Weight=Weight,n,K,gamma=gamma.e)
  
  QIFAI<-function(para2){
    U1<-Sn.f(para=para2,ZZ=X,TT=Y,delta=delta,Weight=Weight,n,K,MA=C.11)
    U2<-Sn.f(para=para2,ZZ=X,TT=Y,delta=delta,Weight=Weight,n,K,MA=C.12)
    U3<-Sn.f(para=gamma.e,ZZ=Z,TT=Y,delta=delta,Weight=Weight,n,K,MA=C.11)
    U<-rbind(U1,U2,U3)
    out<-t(U)%*%solve(CovariAI)%*%U
    return(out)
  }
  est.AI<-optim(par = para.init,fn=QIFAI,
                method="Nelder-Mead",control = list(maxit=200,ndeps=0.01))$par
  
  #####Variance estimation
  
  data2<-as.vector(t(data1))
  data3<-as.data.frame(matrix(data2,n,K*dim(data1)[2],byrow = TRUE))
  
  training<-NULL
  for(s in 1:50){
    iid<-sample(1:n,n,replace = TRUE)
    data.B<-data3[iid,]
    data4<-as.data.frame(matrix(as.vector(t(data.B)),n*K,dim(data1)[2],byrow=TRUE))
    colnames(data4)<-c("id","delta","Z.1","Z.2","Z.3","Y")
    delta<-data4$delta;Y<-data4$Y;X<-cbind(data4$Z.1,data4$Z.2,data4$Z.3)
    Z<-cbind(data4$Z.2,data4$Z.3)
    est.AIboots<-optim(par = para.init,fn=QIFAI,
                       method="Nelder-Mead",control = list(maxit=200,ndeps=0.01))$par
    training<-rbind(training,est.AIboots)
  }
  
  se.AI<-apply(training,2,sd)
  
  return(round(c(est.QIF,se.QIF,est.AI,se.AI),7))
  
}

install.packages("survival")
library(survival)
install.packages("mvtnorm")
library(mvtnorm)
install.packages("miscTools")
library(miscTools)
install.packages("maxLik")
library(maxLik)
install.packages("aftgee")
library(aftgee)
install.packages("foreign")
library(foreign)
install.packages("foreach")
library(foreach)
install.packages("parallel")
library(parallel)
install.packages("iterators")
library(iterators)
install.packages("doParallel")
library(doParallel)

source("P3_T1_Fun.R")

out.f<-function(data1){
  esti<-data1
  bias<-c(apply(esti[,c(1:3,7:9)],2,mean)-rep(para.truee,2))
  sd<-c(apply(esti[,c(1:3,7:9)],2,sd))
  se<-c(apply(esti[,-c(1:3,7:9)],2,mean))
  mse<-c(apply(t((t(esti[,c(1:3,7:9)])-rep(para.truee,2))^2),2,mean))
  re<-mse/mse[1:3]
  cp<-c(apply(I(abs(t(t(esti[,c(1:3,7:9)])-rep(para.truee,2))/esti[,-c(1:3,7:9)])
                <qnorm(0.975))*1,2,mean))
  output<-round(data.frame(bias,sd,se,mse,re,cp),3)
  colnames(output)<-c("bias","sd","se","mse","re","cp")
  output1<-data.frame(output[1:3,],output[4:6,])
  return(output1)
}


packagelist = c('aftgee','mvtnorm','miscTools', 'boot', 'foreign','foreach','parallel','iterators','doParallel','survival','Matrix','base','maxLik') #将需要的包全部都放进来
#######1
no_cores <- detectCores(logical = FALSE)-2  # 电脑核数
cl<-makeCluster(no_cores)  # 分群环境
registerDoParallel(cl)  # 注册并行




###duizhao
in.time<-proc.time()
Reap<-27
nn<-75;mm<-nn*40;KK<-3;
rhoo<-0.2;cenn<-16;
ep.meann<-rep(0,each=KK)
ep.varr<-diag(1,KK)+rhoo*(matrix(1,KK,KK)-diag(1,KK))
para.truee<-c(-1,1,0.5)
Weightt<-"gehan"

Auxin<-NULL
for(i in 1:50){
  Auxin<-rbind(Auxin,aftGmm.ext(seed=i))
}
gamma.in<-c(apply(Auxin,2,mean))

esti<-NULL
esti= foreach(j = 1:Reap, .combine = rbind, .multicombine = TRUE, .packages = packagelist)  %dopar%  GMM(j)

output1<-out.f(data=esti);output1

fi.time<-proc.time()-in.time
fi.time


######################################################
Out<-NULL

###################################Table 1

#############rho=0.2

###Cen=25
in.time<-proc.time()
Reap<-1000
nn<-75;mm<-nn*40;KK<-3;
rhoo<-0.2;cenn<-16;
ep.meann<-rep(0,each=KK)
ep.varr<-diag(1,KK)+rhoo*(matrix(1,KK,KK)-diag(1,KK))
para.truee<-c(-1,1,0.5)
Weightt<-"gehan"
Auxin<-NULL
for(i in 1:50){
  Auxin<-rbind(Auxin,aftGmm.ext(seed=i))
}
gamma.in<-c(apply(Auxin,2,mean))

esti<-NULL
esti= foreach(j = 1:Reap, .combine = rbind, .multicombine = TRUE, .packages = packagelist)  %dopar%  GMM(j)
write.table( esti,file="save/T1_02_25.txt",sep = ' ',append = TRUE,row.names = TRUE,col.names = FALSE)

output1<-out.f(data=esti);output1
Out<-rbind(Out,output1)

fi.time<-proc.time()-in.time
fi.time

write.table(round(Out,3),file="save/Table1.csv",sep = ",")

#repli：模拟次数
#Solve(i): 自己的估计过程写成的函数
#MyResult: 估计的结果，可以提取结果

stopImplicitCluster()  # 结束并行
