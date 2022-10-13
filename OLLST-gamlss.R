
#-----------------------------ROTINA SEM PACOTES-----------------------------------------------
require(numDeriv)

OLLST <- function (mu.link = "identity", sigma.link="log", nu.link = "identity", tau.link = "log")
{
    mstats <- checklink(   "mu.link", "odd log logistic skew t-student", substitute(mu.link), 
                           c("inverse", "log", "identity", "own"))
    dstats <- checklink("sigma.link", "odd log logistic skew t-student", substitute(sigma.link), 
                           c("inverse", "log", "identity", "own"))
    vstats <- checklink(   "nu.link", "odd log logistic skew t-student", substitute(nu.link),    
                           c("1/nu^2", "log", "identity", "own"))
    tstats <- checklink(  "tau.link", "odd log logistic skew t-student", substitute(tau.link),   
                           c("1/tau^2", "log", "identity", "own")) 
    structure(
          list(family = c("OLLST", "odd log logistic skew t-student"),#T STUDENT
           parameters = list(mu=TRUE, sigma=TRUE, nu=TRUE, tau=TRUE), 
                nopar = 4, 
                 type = "Continuous",
              mu.link = as.character(substitute(mu.link)),  
           sigma.link = as.character(substitute(sigma.link)), 
              nu.link = as.character(substitute(nu.link)), 
             tau.link = as.character(substitute(tau.link)), 
           mu.linkfun = mstats$linkfun, 
        sigma.linkfun = dstats$linkfun, 
           nu.linkfun = vstats$linkfun,
          tau.linkfun = tstats$linkfun,  
           mu.linkinv = mstats$linkinv, 
        sigma.linkinv = dstats$linkinv,
           nu.linkinv = vstats$linkinv,
          tau.linkinv = tstats$linkinv, 
                mu.dr = mstats$mu.eta, 
             sigma.dr = dstats$mu.eta,
                nu.dr = vstats$mu.eta,
               tau.dr = tstats$mu.eta, 

   dldm = function(y,mu,sigma,nu,tau){ #----------------------------------------------------- ok
	lpdf<-function(t,x,sigma,nu,tau){log(dauxiOLLST(t,x,sigma,nu,tau))}
	dldm<-grad(func=lpdf,t=y,x=mu,sigma=sigma,nu=nu,tau=tau,method='simple')
        dldm
           },
   d2ldm2 = function(y,mu,sigma,nu,tau){#----------------------------------------------------- ok
	lpdf<-function(t,x,sigma,nu,tau){log(dauxiOLLST(t,x,sigma,nu,tau))}
	dldm<-grad(func=lpdf,t=y,x=mu,sigma=sigma,nu=nu,tau=tau,method='simple')
        d2ldm2 <- -dldm * dldm
        d2ldm2 <- ifelse(d2ldm2 < -1e-15, d2ldm2,-1e-15) 
        d2ldm2 
           },     
   dldd = function(y,mu,sigma,nu,tau){#----------------------------------------------------- ok  
	lpdf<-function(t,mu,x,nu,tau){log(dauxiOLLST(t,mu,x,nu,tau))}
	dldd<-grad(func=lpdf,t=y,mu=mu,x=sigma,nu=nu,tau=tau,method='simple')
        dldd
           } ,
   d2ldd2 = function(y,mu,sigma,nu,tau){#----------------------------------------------------- ok
	lpdf<-function(t,mu,x,nu,tau){log(dauxiOLLST(t,mu,x,nu,tau))}
	dldd<-grad(func=lpdf,t=y,mu=mu,x=sigma,nu=nu,tau=tau,method='simple')
        d2ldd2 <- -dldd*dldd
     	d2ldd2 <- ifelse(d2ldd2 < -1e-15, d2ldd2,-1e-15)  
    	d2ldd2
           },   
     dldv = function(y,mu,sigma,nu,tau){#----------------------------------------------------- ok
	lpdf<-function(t,mu,sigma,x,tau){log(dauxiOLLST(t,mu,sigma,x,tau))}
	dldv<-grad(func=lpdf,t=y,mu=mu,sigma=sigma,x=nu,tau=tau,method='simple')
	dldv
           },
    d2ldv2 = function(y,mu,sigma,nu,tau){#----------------------------------------------------- ok 
	lpdf<-function(t,mu,sigma,x,tau){log(dauxiOLLST(t,mu,sigma,x,tau))}
	dldv<-grad(func=lpdf,t=y,mu=mu,sigma=sigma,x=nu,tau=tau,method='simple')
	d2ldv2<- -dldv * dldv
        d2ldv2 <- ifelse(d2ldv2 < -1e-15, d2ldv2,-1e-15)                   
        d2ldv2
           },
      dldt = function(y,mu,sigma,nu,tau){#----------------------------------------------------- ok
	lpdf<-function(t,mu,sigma,nu,x){log(dauxiOLLST(t,mu,sigma,nu,x))}
	dldt<-grad(func=lpdf,t=y,mu=mu,sigma=sigma,nu=nu,x=tau)
	dldt
           } ,
      d2ldt2 = function(y,mu,sigma,nu,tau){ #----------------------------------------------------- ok
	lpdf<-function(t,mu,sigma,nu,x){log(dauxiOLLST(t,mu,sigma,nu,x))}
	dldt<-grad(func=lpdf,t=y,mu=mu,sigma=sigma,nu=nu,x=tau,method='simple')
	d2ldt2<- -dldt * dldt
	d2ldt2 <- ifelse(d2ldt2 < -1e-15, d2ldt2,-1e-15)                                    
	d2ldt2
           },
  d2ldmdd = function(y,mu,sigma,nu,tau){#----------------------------------------------------- ok
	lpdf<-function(t,x,sigma,nu,tau){log(dauxiOLLST(t,x,sigma,nu,tau))}
	dldm<-grad(func=lpdf,t=y,x=mu,sigma=sigma,nu=nu,tau=tau,method='simple')
	lpdf<-function(t,mu,x,nu,tau){log(dauxiOLLST(t,mu,x,nu,tau))}
	dldd<-grad(func=lpdf,t=y,mu=mu,x=sigma,nu=nu,tau=tau)
	 d2ldmdd = -(dldm * dldd)
	 d2ldmdd<-ifelse(is.na(d2ldmdd)==TRUE,0,d2ldmdd)
         d2ldmdd                 
           },
  d2ldmdv = function(y,mu,sigma,nu,tau){#----------------------------------------------------- ok
	lpdf<-function(t,x,sigma,nu,tau){log(dauxiOLLST(t,x,sigma,nu,tau))}
	dldm<-grad(func=lpdf,t=y,x=mu,sigma=sigma,nu=nu,tau=tau,method='simple')
	lpdf<-function(t,mu,sigma,x,tau){log(dauxiOLLST(t,mu,sigma,x,tau))}
	dldv<-grad(func=lpdf,t=y,mu=mu,sigma=sigma,x=nu,tau=tau)
       		 d2ldmdv = -(dldm * dldv)
       		 d2ldmdv				
           },

  d2ldmdt = function(y,mu,sigma,nu,tau){#----------------------------------------------------- ok
	lpdf<-function(t,x,sigma,nu,tau){log(dauxiOLLST(t,x,sigma,nu,tau))}
	dldm<-grad(func=lpdf,t=y,x=mu,sigma=sigma,nu=nu,tau=tau,method='simple')
	lpdf<-function(t,mu,sigma,nu,x){log(dauxiOLLST(t,mu,sigma,nu,x))}
	dldt<-grad(func=lpdf,t=y,mu=mu,sigma=sigma,nu=nu,x=tau)
			d2ldmdt <- -(dldm*dldt)
			d2ldmdt
          },

  d2ldddv = function(y,mu,sigma,nu,tau){#----------------------------------------------------- ok
	lpdf<-function(t,mu,x,nu,tau){log(dauxiOLLST(t,mu,x,nu,tau))}
	dldd<-grad(func=lpdf,t=y,mu=mu,x=sigma,nu=nu,tau=tau,method='simple')
	lpdf<-function(t,mu,sigma,x,tau){log(dauxiOLLST(t,mu,sigma,x,tau))}
	dldv<-grad(func=lpdf,t=y,mu=mu,sigma=sigma,x=nu,tau=tau)
			d2ldddv = -(dldd * dldv)
            		d2ldddv	
           },
  d2ldddt = function(y,mu,sigma,nu,tau){#----------------------------------------------------- ok
	lpdf<-function(t,mu,x,nu,tau){log(dauxiOLLST(t,mu,x,nu,tau))}
	dldd<-grad(func=lpdf,t=y,mu=mu,x=sigma,nu=nu,tau=tau,method = 'simple')
	lpdf<-function(t,mu,sigma,nu,x){log(dauxiOLLST(t,mu,sigma,nu,x))}
	dldt<-grad(func=lpdf,t=y,mu=mu,sigma=sigma,nu=nu,x=tau)
			d2ldddt <- -(dldd*dldt) 
			d2ldddt 
           },
  d2ldvdt = function(y,mu,sigma,nu,tau){#----------------------------------------------------- ok
	lpdf<-function(t,mu,sigma,x,tau){log(dauxiOLLST(t,mu,sigma,x,tau))}
	dldv<-grad(func=lpdf,t=y,mu=mu,sigma=sigma,x=nu,tau=tau,method = 'simple')
	lpdf<-function(t,mu,sigma,nu,x){log(dauxiOLLST(t,mu,sigma,nu,x))}
	dldt<-grad(func=lpdf,t=y,mu=mu,sigma=sigma,nu=nu,x=tau)
			d2ldvdt <- -(dldv*dldt) 
			d2ldvdt 
           },
#----------------------------------------------------- ok

 G.dev.incr  = function(y,mu,sigma,nu,tau,...) 
                       { 
           -2*dOLLST(y,mu,sigma,nu,tau,log=TRUE)
                        } ,                     
         rqres = expression(   
               rqres(pfun="pOLLST", type="Continuous", y=y, mu=mu, sigma=sigma, nu=nu, tau=tau)) ,
    mu.initial = expression(mu <- (y+mean(y))/2), 
    sigma.initial = expression(sigma <- rep(sd(y), length(y))), 
    nu.initial = expression(nu <- rep(1, length(y))), 
   tau.initial = expression(tau <-rep(1, length(y))), 
      mu.valid = function(mu) TRUE, 
   sigma.valid = function(sigma) all(sigma > 0),
      nu.valid = function(nu) TRUE, 
     tau.valid = function(tau) all(tau > 0), 
       y.valid = function(y)  TRUE
          ),
            class = c("gamlss.family","family"))
}



#---------------------PDF--------------------------------------------
dOLLST <- function(x, mu = 0, sigma = 1, nu= 2, tau=0.2, log = FALSE){
          if (any(tau < 0))  stop(paste("tau must be positive", "\n", ""))  
 
  theta=300
  z=(x-mu)/sigma 
  
  pdfst <- (2/sigma)*dt(z,theta,0,log = FALSE)*pt(nu*z,theta,0,lower.tail = TRUE,log.p = FALSE)   
  cdfst <-  pt(z,theta,0,lower.tail = TRUE, log.p = FALSE)-2*owen(z,nu)  
  fy1 <- (tau*pdfst*(cdfst*(1-cdfst))^(tau-1))/((cdfst^tau+(1-cdfst)^tau)^2) 
     if(log==FALSE) fy<-fy1 else fy<-log(fy1)
  fy
  }    
#--------------------------CDF---------------------------------------  
pOLLST <- function(q, mu = 0, sigma = 1, nu= 2, tau=0.2, lower.tail = TRUE, log.p = FALSE){  
         if (any(tau < 0))  stop(paste("tau must be positive", "\n", ""))
  
  theta=300
  z=(q-mu)/sigma
  cdfst <-  pt(z,theta,0,lower.tail = TRUE, log.p = FALSE)-2*owen(z,nu)  
  cdf1 <- (cdfst^tau)/(cdfst^tau+(1-cdfst)^tau)
    if(lower.tail==TRUE) cdf<-cdf1 else cdf<- 1-cdf1
    if(log.p==FALSE) cdf<- cdf else cdf<- log(cdf)
cdf
 } 
#---------------------------QUANTILICA--------------------------------------  

qOLLST <-  function(p, mu=0, sigma=1, nu= 2, tau=0.2, lower.tail = TRUE, log.p = FALSE){   
         if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", "")) 
	 if (any(tau < 0))    stop(paste("tau must be positive", "\n", "")) 
         if (any(p < 0)|any(p > 1))  stop(paste("p must be between 0 and 1", "\n", ""))    
  	     if (log.p==TRUE) p <- exp(p) else p <- p
  	     if (lower.tail==TRUE) p <- p else p <- 1-p
	     
  theta=300
  u <- (p^(1/tau))/((1-p)^(1/tau)+p^(1/tau))
  q <- mu + sigma*qst(u,mu, omega = 1,tau, theta) 
  q
}
 
#---------------------------RANDOM--------------------------------------  
rOLLST <- function(n, mu=0, sigma=1, nu= 2, tau=0.2){
    if (any(tau < 0))  stop(paste("tau must be positive", "\n", ""))  
     
    uni<- runif(n = n,0,1)
  	r <- qOLLST(uni,mu =mu, sigma =sigma, nu=nu, tau=tau)
    r
 }

#-------------------------DENSIDADE AUXILIAR----------------------------------------  

dauxiOLLST <- function(t,mu,sigma,nu,tau){ 
  
  theta=300
  z=(t-mu)/sigma
  pdfst <- (2/sigma)*dt(z,theta,0,log = FALSE)*pt(nu*z,theta,0,lower.tail = TRUE,log.p = FALSE)   
  cdfst <-  pt(z,theta,0,lower.tail = TRUE, log.p = FALSE)-2*owen(z,nu)  
  fy1 <- (tau*pdfst*(cdfst*(1-cdfst))^(tau-1))/((cdfst^tau+(1-cdfst)^tau)^2)
  fy1}

#-------------------------OWEN----------------------------------------  
owen<-function(h,a){
  func<-function(x,h)((exp(-0.5*h^2*(1+x^2)))/(1+x^2))*(1/(2*pi))
  temp1<-c()
  if(length(h)>1 & length(a)>1){ for(i in 1:length(h)){
    int<-integrate(f=func,lower =0,upper=a[i],h=h[i])
    temp1<-c(temp1,int$value)
  }}
  if(length(h)>1 & length(a)==1){ for(i in 1:length(h)){
    int<-integrate(f=func,lower =0,upper=a,h=h[i])
    temp1<-c(temp1,int$value)
  }}
  if(length(h)==1 & length(a)==1){ 
    int<-integrate(f=func,lower =0,upper=a,h=h)
    temp1<-c(temp1,int$value)
  }
  return(temp1)}




