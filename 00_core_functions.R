# ================================================================================
# Script Name：00_core_functions.R
# Functional Description：All core functions
# Package：[simsurv、progress]
# ================================================================================

# Load dependency packages
library("simsurv")
library("progress")

# --Function1：thetacal(Program for Solving Key Parameter θ in Censoring Time Distribution) ----------------
# --When both beta1 and beta2 are nonzero, use kernel density estimation to handle covariate effects--------
# --Includes four scenarios: no covariates, one binary covariate, one normally distributed covariate, and multiple covariates---

thetacal <- function(seedtheta,gamma1,lambda1,gamma2,lambda2,mixp,cens.p,size,binomp,
                     beta1,beta2,maxt,mean,std,threshold,censdist,unirootlower,unirootupper)
{
  if(beta1==0 && beta2==0){
    upper <-"maxt"
    if(censdist=="Uniform"||censdist=="Uni"||censdist=="U"){
      theta0 <- maxt
      cenmaxt.P3<-integrate(function(t){
        part2<-dunif(t,min=0,max=theta0)
        part3<-(mixp*exp(-lambda1*(t^gamma1))+(1-mixp)*exp(-lambda2*(t^gamma2)))
        return(part2*part3)
      } , 0,maxt)$value
      
      cenmaxt.P1<-integrate(function(t){
        part2<- dunif(t,min=0,max=theta0)
        part3<-(mixp*exp(-lambda1*(maxt^gamma1))+(1-mixp)*exp(-lambda2*(maxt^gamma2)))
        return(part2*part3)
      } , maxt,Inf)$value
      
      
      if(cenmaxt.P1+cenmaxt.P3-cens.p>0){
        upper  <-  "maxt"
      }else{upper  <-  "theta"}
    }
    
    if(censdist=="Exponential"||censdist=="Expo"||censdist=="E"){ 
      censor.prop<-function(theta){
        cen.P3<-integrate(function(t){
          part2<-dexp(t,theta)
          part3<-(mixp*exp(-lambda1*(t^gamma1))+(1-mixp)*exp(-lambda2*(t^gamma2)))
          return(part2*part3)
        }, 0,maxt)$value
        
        cen.P1<-integrate(function(t){
          part2<-dexp(t,theta)
          part3<-(mixp*exp(-lambda1*(maxt^gamma1))+(1-mixp)*exp(-lambda2*(maxt^gamma2)))
          return(part2*part3)
        } , maxt,Inf)$value
        
        
        return(cen.P3+cen.P1-cens.p)}
      
      theta<-uniroot(censor.prop,c(0,200),tol=0.00000000000000000001)$root
      return(list(theta))
    }
    
    if(censdist=="Uniform"||censdist=="Uni"||censdist=="U"){
      censor.prop<-function(theta){
        if(censdist=="Uniform"||censdist=="Uni"||censdist=="U" && upper=="maxt"){
          cen.P3<-
            integrate(function(t){
              part2<-dunif(t,0,theta)
              part3<-(mixp*exp(-lambda1*(t^gamma1))+(1-mixp)*exp(-lambda2*(t^gamma2)))
              return(part2*part3)
            } , 0,maxt)$value
        }
        
        if(censdist=="Uniform"||censdist=="Uni"||censdist=="U" && upper=="theta"){
          cen.P3<-
            integrate(function(t){
              part2<-dunif(t,0,theta)
              part3<-(mixp*exp(-lambda1*(t^gamma1))+(1-mixp)*exp(-lambda2*(t^gamma2)))
              return(part2*part3)
            } , 0,theta)$value
          
        }
        
        
        if(censdist=="Uniform"||censdist=="Uni"||censdist=="U" && upper=="maxt"){
          cen.P1<-
            integrate(function(t){
              
              part2<-dunif(t,0,theta)
              part3<-(mixp*exp(-lambda1*(maxt^gamma1))+(1-mixp)*exp(-lambda2*(maxt^gamma2)))
              return(part2*part3)
            } , maxt,Inf)$value 
          
        }
        
        if(censdist=="Uniform"||censdist=="Uni"||censdist=="U" && upper=="theta"){
          cen.P1<-
            integrate(function(t){
              
              part2<-dunif(t,0,theta)
              part3<-(mixp*exp(-lambda1*(maxt^gamma1))+(1-mixp)*exp(-lambda2*(maxt^gamma2)))
              return(part2*part3)
            } , maxt,Inf)$value 
          
        }
        return(cen.P1+cen.P3-cens.p)
      }
    }      
    theta<-uniroot(censor.prop,c(unirootlower,unirootupper),tol=0.000000000000000000001)$root
    return(list(theta))
  }
  
 
  if(beta1!=0 && beta2==0){
    upper <-"maxt"
    if(censdist=="Uniform"||censdist=="Uni"||censdist=="U"){
      theta0 <- maxt
      
      cenmaxt.control.P3<-
        integrate(function(t){
          part2<-dunif(t,0,theta0)
          part3<-(mixp*exp(-lambda1*(t^gamma1))+(1-mixp)*exp(-lambda2*(t^gamma2))) 
          return(part2*part3)
        } , 0,maxt)$value
      cenmaxt.control.P1<-
        integrate(function(t){
          part2<-dunif(t,0,theta0)
          part3<-(mixp*exp(-lambda1*(maxt^gamma1))+(1-mixp)*exp(-lambda2*(maxt^gamma2))) 
          return(part2*part3)
        } , maxt,Inf)$value
      
      cenmaxt.exposure.P3<-
        integrate(function(t){
          part2<-dunif(t,0,theta0)
          part3<-(mixp*exp(-lambda1*(t^gamma1))+(1-mixp)*exp(-lambda2*(t^gamma2)))^exp(beta1)
          return(part2*part3)
        } , 0,maxt)$value
      cenmaxt.exposure.P1<-
        integrate(function(t){
          part2<-dunif(t,0,theta0)
          part3<-(mixp*exp(-lambda1*(maxt^gamma1))+(1-mixp)*exp(-lambda2*(maxt^gamma2)))^exp(beta1)
          return(part2*part3)
        } , maxt,Inf)$value
      
      if((cenmaxt.control.P3+cenmaxt.control.P1)*(1-binomp)+(cenmaxt.exposure.P3+cenmaxt.exposure.P1)*(binomp)-cens.p>0){
        upper  <-  "maxt"
      }else{upper  <-  "theta"}
    }
    
    if(censdist=="Exponential"||censdist=="Expo"||censdist=="E"){ 
      censor.prop<-function(theta){
        cen.control.P3<-
          integrate(function(t){
            part2<-dexp(t,theta)
            part3<-(mixp*exp(-lambda1*(t^gamma1))+(1-mixp)*exp(-lambda2*(t^gamma2))) 
            return(part2*part3)
          } , 0,maxt)$value
        cen.control.P1<-
          integrate(function(t){
            part2<-dexp(t,theta)
            part3<-(mixp*exp(-lambda1*(maxt^gamma1))+(1-mixp)*exp(-lambda2*(maxt^gamma2))) 
            return(part2*part3)
          } , maxt,Inf)$value
        
        cen.exposure.P3<-
          integrate(function(t){
            part2<-dexp(t,theta)
            part3<-(mixp*exp(-lambda1*(t^gamma1))+(1-mixp)*exp(-lambda2*(t^gamma2)))^exp(beta1)
            return(part2*part3)
          } , 0,maxt)$value
        cen.exposure.P1<-
          integrate(function(t){
            part2<-dexp(t,theta)
            part3<-(mixp*exp(-lambda1*(maxt^gamma1))+(1-mixp)*exp(-lambda2*(maxt^gamma2)))^exp(beta1)
            return(part2*part3)
          } , maxt,Inf)$value
        
        return((cen.control.P3+cen.control.P1)*(1-binomp)+(cen.exposure.P3+cen.exposure.P1)*(binomp)-cens.p)
      }
      
      theta<-uniroot(censor.prop,c(0.0000000000000000000001,200),tol=0.0000000000001)$root
      return(theta)
    }
    
    if(censdist=="Uniform"||censdist=="Uni"||censdist=="U"){
      censor.prop<-function(theta){
        if(censdist=="Uniform"||censdist=="Uni"||censdist=="U" && upper=="maxt"){
          
          cen.control.P3<-
            integrate(function(t){
              part2<-dunif(t,0,theta)
              part3<-(mixp*exp(-lambda1*(t^gamma1))+(1-mixp)*exp(-lambda2*(t^gamma2))) 
              return(part2*part3)
            } , 0,maxt)$value
          cen.control.P1<-
            integrate(function(t){
              part2<-dunif(t,0,theta)
              part3<-(mixp*exp(-lambda1*(maxt^gamma1))+(1-mixp)*exp(-lambda2*(maxt^gamma2))) 
              return(part2*part3)
            } , maxt,Inf)$value
          
          cen.exposure.P3<-
            integrate(function(t){
              part2<-dunif(t,0,theta)
              part3<-(mixp*exp(-lambda1*(t^gamma1))+(1-mixp)*exp(-lambda2*(t^gamma2)))^exp(beta1)
              return(part2*part3)
            } , 0,maxt)$value
          cen.exposure.P1<-
            integrate(function(t){
              part2<-dunif(t,0,theta)
              part3<-(mixp*exp(-lambda1*(maxt^gamma1))+(1-mixp)*exp(-lambda2*(maxt^gamma2)))^exp(beta1)
              return(part2*part3)
            } , maxt,Inf)$value
          
          return((cen.control.P3+cen.control.P1)*(1-binomp)+(cen.exposure.P3+cen.exposure.P1)*(binomp)-cens.p)
        }
      }
      
      if(censdist=="Uniform"||censdist=="Uni"||censdist=="U" && upper=="theta"){
        censor.prop<-function(theta){
          cen.control.P3<-
            integrate(function(t){
              part2<-dunif(t,0,theta)
              part3<-(mixp*exp(-lambda1*(t^gamma1))+(1-mixp)*exp(-lambda2*(t^gamma2))) 
              return(part2*part3)
            } , 0,theta)$value
          cen.control.P1<-
            integrate(function(t){
              part2<-dunif(t,0,theta)
              part3<-(mixp*exp(-lambda1*(maxt^gamma1))+(1-mixp)*exp(-lambda2*(maxt^gamma2))) 
              return(part2*part3)
            } , maxt,Inf)$value
          
          cen.exposure.P3<-
            integrate(function(t){
              part2<-dunif(t,0,theta)
              part3<-(mixp*exp(-lambda1*(t^gamma1))+(1-mixp)*exp(-lambda2*(t^gamma2)))^exp(beta1)
              return(part2*part3)
            } , 0,theta)$value
          cen.exposure.P1<-
            integrate(function(t){
              part2<-dunif(t,0,theta)
              part3<-(mixp*exp(-lambda1*(maxt^gamma1))+(1-mixp)*exp(-lambda2*(maxt^gamma2)))^exp(beta1)
              return(part2*part3)
            } , maxt,Inf)$value
          
          return((cen.control.P3+cen.control.P1)*(1-binomp)+(cen.exposure.P3+cen.exposure.P1)*(binomp)-cens.p)
        }
      }
      
      
      theta<-uniroot(censor.prop,c(unirootlower,unirootupper),tol=0.0000000000001)$root
      return(theta)
      
    }  
  }
 
  if(beta1==0 && beta2!=0){
    upper <-"maxt"
    if(censdist=="Uniform"||censdist=="Uni"||censdist=="U"){
      theta0 <- maxt
      cenmaxt.P3<-integrate(function(u){
        sapply(u,function(u){
          integrate(function(t){
            part1<-dnorm(u,mean,std)
            part2<-dunif(t,min=0,max=theta0)
            part3<-(mixp*exp(-lambda1*(t^gamma1))+(1-mixp)*exp(-lambda2*(t^gamma2)))^exp(u*beta2)
            return(part1*part2*part3)
          } , 0,maxt)$value
        })
      },-Inf,Inf)$value
      cenmaxt.P1<-integrate(function(u){
        sapply(u,function(u){
          integrate(function(t){
            part1<-dnorm(u,mean,std)
            part2<- dunif(t,min=0,max=theta0)
            part3<-(mixp*exp(-lambda1*(maxt^gamma1))+(1-mixp)*exp(-lambda2*(maxt^gamma2)))^exp(u*beta2)
            return(part1*part2*part3)
          } , maxt,Inf)$value
        })
      },-Inf,Inf)$value
      
      if(cenmaxt.P1+cenmaxt.P3-cens.p>0){
        upper  <-  "maxt"
      }else{upper  <-  "theta"}
    }
    
    if(censdist=="Exponential"||censdist=="Expo"||censdist=="E"){ 
      censor.prop<-function(theta){
        cen.P3<-integrate(function(u){
          sapply(u,function(u){
            integrate(function(t){
              part1<-dnorm(u,mean,std)
              part2<-dexp(t,theta)
              part3<-(mixp*exp(-lambda1*(t^gamma1))+(1-mixp)*exp(-lambda2*(t^gamma2)))^exp(u*beta2)
              return(part1*part2*part3)
            }, 0,maxt)$value
          })
        },-Inf,Inf)$value
        cen.P1<-integrate(function(u){
          sapply(u,function(u){
            integrate(function(t){
              
              part1<-dnorm(u,mean,std)
              part2<-dexp(t,theta)
              part3<-(mixp*exp(-lambda1*(maxt^gamma1))+(1-mixp)*exp(-lambda2*(maxt^gamma2)))^exp(u*beta2)
              return(part1*part2*part3)
            } , maxt,Inf)$value
          })
        },-Inf,Inf)$value
        
        return(cen.P3+cen.P1-cens.p)}
      
      theta<-uniroot(censor.prop,c(0,200),tol=0.00000000000000000001)$root
      return(list(theta))
    }
    
    if(censdist=="Uniform"||censdist=="Uni"||censdist=="U"){
      censor.prop<-function(theta){
        if(censdist=="Uniform"||censdist=="Uni"||censdist=="U" && upper=="maxt"){
          cen.P3<-integrate(function(u){
            sapply(u,function(u){
              integrate(function(t){
                
                part1<-dnorm(u,mean,std)
                part2<-dunif(t,0,theta)
                part3<-(mixp*exp(-lambda1*(t^gamma1))+(1-mixp)*exp(-lambda2*(t^gamma2)))^exp(u*beta2)
                return(part1*part2*part3)
              } , 0,maxt)$value
            })
          },-Inf,Inf)$value
        }
        
        if(censdist=="Uniform"||censdist=="Uni"||censdist=="U" && upper=="theta"){
          cen.P3<-integrate(function(u){
            sapply(u,function(u){
              integrate(function(t){
                
                part1<-dnorm(u,mean,std)
                part2<-dunif(t,0,theta)
                part3<-(mixp*exp(-lambda1*(t^gamma1))+(1-mixp)*exp(-lambda2*(t^gamma2)))^exp(u*beta2)
                return(part1*part2*part3)
              } , 0,theta)$value
            })
          },-Inf,Inf)$value
        }
        
        
        if(censdist=="Uniform"||censdist=="Uni"||censdist=="U" && upper=="maxt"){
          cen.P1<-integrate(function(u){
            sapply(u,function(u){
              integrate(function(t){
                part1<-dnorm(u,mean,std)
                part2<-dunif(t,0,theta)
                part3<-(mixp*exp(-lambda1*(maxt^gamma1))+(1-mixp)*exp(-lambda2*(maxt^gamma2)))^exp(u*beta2)
                return(part1*part2*part3)
              } , maxt,Inf)$value 
            })
          },-Inf,Inf)$value
        }
        
        if(censdist=="Uniform"||censdist=="Uni"||censdist=="U" && upper=="theta"){
          cen.P1<-integrate(function(u){
            sapply(u,function(u){
              integrate(function(t){
                part1<-dnorm(u,mean,std)
                part2<-dunif(t,0,theta)
                part3<-(mixp*exp(-lambda1*(maxt^gamma1))+(1-mixp)*exp(-lambda2*(maxt^gamma2)))^exp(u*beta2)
                return(part1*part2*part3)
              } , maxt,Inf)$value 
            })
          },-Inf,Inf)$value
        }
        return(cen.P1+cen.P3-cens.p)
      }
      theta<-uniroot(censor.prop,c(unirootlower,unirootupper),tol=0.000000000001)$root
      return(list(theta))
    }
  }
 
  if(beta1!=0 && beta2!=0){
    set.seed(seedtheta)
    beta1 <- beta1
    x1 <- rbinom(10000,1L,binomp)
    beta2 <- beta2
    x2 <- rnorm(10000,mean,std)
    m<- exp(x1*beta1+x2*beta2)
    #dens<-density(m,n=2000,bw=0.1,from=0,to=400,na.rm=TRUE)
    #x<-dens$x
    #y<-dens$y
    #y.loess <-loess(y~x,span=0.1)
    ### check estimated density function and get the range of lambda ###
    #plot(limit,lty=1,lwd=2)
    #lines(y.loess,col=2)
    
    dens<-density(m,bw = "nrd0",na.rm=TRUE,n=10000)
    x <- dens$x
    y <- dens$y
    y.loess <-loess(y~x,span=0.01)
    #plot(limit,lty=1,lwd=2)
    #lines(y.loess,col=2)
    densdata <- data.frame(x,y)
    subdens <- subset(densdata,y>=threshold)
    max.m.i<-max(subdens$x)
    min.m.i<-min(subdens$x)
    
    density.fun.lambda<-function(u){
      pred.y <- predict(y.loess, newdata=u)
      return(pred.y)
    }
    upper <-"maxt"
    if(censdist=="Uniform"||censdist=="Uni"||censdist=="U"){
      theta0 <- maxt
      cenmaxt.P3<-integrate(function(u){
        sapply(u,function(u){
          integrate(function(t){
            m.i<-u
            part1<-density.fun.lambda(m.i)
            part2<- dunif(t,min=0,max=theta0)
            part3<-(mixp*exp(-lambda1*(t^gamma1))+(1-mixp)*exp(-lambda2*(t^gamma2)))^m.i
            return(part1*part2*part3)
          } , 0,maxt)$value
        })
      },min.m.i,max.m.i)$value
      cenmaxt.P1<-integrate(function(u){
        sapply(u,function(u){
          integrate(function(t){
            m.i<-u
            part1<-density.fun.lambda(m.i)
            part2<- dunif(t,min=0,max=theta0)
            part3<-(mixp*exp(-lambda1*(maxt^gamma1))+(1-mixp)*exp(-lambda2*(maxt^gamma2)))^m.i
            return(part1*part2*part3)
          } , maxt,Inf)$value
        })
      },min.m.i,max.m.i)$value
      
      if(cenmaxt.P1+cenmaxt.P3-cens.p>0){
        upper  <-  "maxt"
      }else{upper  <-  "theta"}
    }
    
    if(censdist=="Exponential"||censdist=="Expo"||censdist=="E"){ 
      censor.prop<-function(theta){
        cen.P3<-integrate(function(u){
          sapply(u,function(u){
            integrate(function(t){
              m.i<-u
              part1<-density.fun.lambda(m.i)
              part2<- dexp(t,theta)
              part3<-(mixp*exp(-lambda1*(t^gamma1))+(1-mixp)*exp(-lambda2*(t^gamma2)))^m.i
              return(part1*part2*part3)
            },0,maxt)$value
          })
        },min.m.i,max.m.i)$value
        cen.P1<-integrate(function(u){
          sapply(u,function(u){
            integrate(function(t){
              m.i<-u
              part1<-density.fun.lambda(m.i)
              part2<- dexp(t,theta)
              part3<-(mixp*exp(-lambda1*(maxt^gamma1))+(1-mixp)*exp(-lambda2*(maxt^gamma2)))^m.i
              return(part1*part2*part3)
            }, maxt,Inf)$value
          })
        },min.m.i,max.m.i)$value
        
        return(cen.P1+cen.P3-cens.p)}
      
      theta<-uniroot(censor.prop,c(0,200),tol=0.00000000000000000001)$root
      return(list(theta))
    }
    
    if(censdist=="Uniform"||censdist=="Uni"||censdist=="U"){
      censor.prop<-function(theta){
        if(censdist=="Uniform"||censdist=="Uni"||censdist=="U" && upper=="maxt"){
          cen.P3<-integrate(function(u){
            sapply(u,function(u){
              integrate(function(t){
                m.i<-u
                part1<-density.fun.lambda(m.i)
                part2<-dunif(t,0,theta)
                part3<-(mixp*exp(-lambda1*(t^gamma1))+(1-mixp)*exp(-lambda2*(t^gamma2)))^m.i
                return(part1*part2*part3)
              } , 0,maxt)$value
            })
          },min.m.i,max.m.i)$value
        }
        
        if(censdist=="Uniform"||censdist=="Uni"||censdist=="U" && upper=="theta"){
          cen.P3<-integrate(function(u){
            sapply(u,function(u){
              integrate(function(t){
                m.i<-u
                part1<-density.fun.lambda(m.i)
                part2<-dunif(t,0,theta)
                part3<-(mixp*exp(-lambda1*(t^gamma1))+(1-mixp)*exp(-lambda2*(t^gamma2)))^m.i
                return(part1*part2*part3)
              } , 0,theta)$value
            })
          },min.m.i,max.m.i)$value
        }
        
        if(censdist=="Uniform"||censdist=="Uni"||censdist=="U" && upper=="maxt"){
          cen.P1<-(1-(maxt/theta))*integrate(function(u){
            sapply(u,function(u){
              m.i<-u
              part1<-density.fun.lambda(m.i)
              part3<-(mixp*exp(-lambda1*(maxt^gamma1))+(1-mixp)*exp(-lambda2*(maxt^gamma2)))^m.i
              return(part1*part3)
            } 
            )
          },min.m.i,max.m.i)$value
        }
        if(censdist=="Uniform"||censdist=="Uni"||censdist=="U" && upper=="theta"){
          cen.P1<-integrate(function(u){
            sapply(u,function(u){
              integrate(function(t){
                m.i<-u
                part1<-density.fun.lambda(m.i)
                part2<-dunif(t,0,theta)
                part3<-(mixp*exp(-lambda1*(maxt^gamma1))+(1-mixp)*exp(-lambda2*(maxt^gamma2)))^m.i
                return(part1*part2*part3)
              } , maxt,Inf)$value 
            })
          },min.m.i,max.m.i)$value
        }
        return(cen.P1+cen.P3-cens.p)
      }
      theta<-uniroot(censor.prop,c(unirootlower,unirootupper),tol=0.000000000000000000001)$root
      return(list(theta))
    }
  }
}


# --Function2：thetacal_a(Program for Solving Key Parameter θ in Censoring Time Distribution)--------------
# --If either beta1 or beta2 is nonzero, use kernel density estimation to handle covariate effects------------------

thetacal_a <- function(seedtheta,gamma1,lambda1,gamma2,lambda2,mixp,cens.p,size,
                       binomp,beta1,beta2,maxt,mean,std,threshold,censdist,unirootlower,unirootupper)
{
  
  if(beta1!=0 || beta2!=0){
    set.seed(seedtheta)
    beta1 <- beta1
    x1 <- rbinom(10000,1L,binomp)
    beta2 <- beta2
    x2 <- rnorm(10000,mean,std)
    m<- exp(x1*beta1+x2*beta2)
    #dens<-density(m,n=2000,bw=0.1,from=0,to=400,na.rm=TRUE)
    #x<-dens$x
    #y<-dens$y
    #y.loess <-loess(y~x,span=0.1)
    ### check estimated density function and get the range of lambda ###
    #plot(limit,lty=1,lwd=2)
    #lines(y.loess,col=2)
    
    dens<-density(m,bw = "nrd0",na.rm=TRUE,n=10000)
    x <- dens$x
    y <- dens$y
    y.loess <-loess(y~x,span=0.01)
    #plot(limit,lty=1,lwd=2)
    #lines(y.loess,col=2)
    densdata <- data.frame(x,y)
    subdens <- subset(densdata,y>=threshold)
    max.m.i<-max(subdens$x)
    min.m.i<-min(subdens$x)
    
    density.fun.lambda<-function(u){
      pred.y <- predict(y.loess, newdata=u)
      return(pred.y)
    }
    upper <-"maxt"
    if(censdist=="Uniform"||censdist=="Uni"||censdist=="U"){
      theta0 <- maxt
      cenmaxt.P3<-integrate(function(u){
        sapply(u,function(u){
          integrate(function(t){
            m.i<-u
            part1<-density.fun.lambda(m.i)
            part2<- dunif(t,min=0,max=theta0)
            part3<-(mixp*exp(-lambda1*(t^gamma1))+(1-mixp)*exp(-lambda2*(t^gamma2)))^m.i
            return(part1*part2*part3)
          } , 0,maxt)$value
        })
      },min.m.i,max.m.i)$value
      cenmaxt.P1<-integrate(function(u){
        sapply(u,function(u){
          integrate(function(t){
            m.i<-u
            part1<-density.fun.lambda(m.i)
            part2<- dunif(t,min=0,max=theta0)
            part3<-(mixp*exp(-lambda1*(maxt^gamma1))+(1-mixp)*exp(-lambda2*(maxt^gamma2)))^m.i
            return(part1*part2*part3)
          } , maxt,Inf)$value
        })
      },min.m.i,max.m.i)$value
      
      if(cenmaxt.P1+cenmaxt.P3-cens.p>0){
        upper  <-  "maxt"
      }else{upper  <-  "theta"}
    }
    
    if(censdist=="Exponential"||censdist=="Expo"||censdist=="E"){ 
      censor.prop<-function(theta){
        cen.P3<-integrate(function(u){
          sapply(u,function(u){
            integrate(function(t){
              m.i<-u
              part1<-density.fun.lambda(m.i)
              part2<- dexp(t,theta)
              part3<-(mixp*exp(-lambda1*(t^gamma1))+(1-mixp)*exp(-lambda2*(t^gamma2)))^m.i
              return(part1*part2*part3)
            },0,maxt)$value
          })
        },min.m.i,max.m.i)$value
        cen.P1<-integrate(function(u){
          sapply(u,function(u){
            integrate(function(t){
              m.i<-u
              part1<-density.fun.lambda(m.i)
              part2<- dexp(t,theta)
              part3<-(mixp*exp(-lambda1*(maxt^gamma1))+(1-mixp)*exp(-lambda2*(maxt^gamma2)))^m.i
              return(part1*part2*part3)
            }, maxt,Inf)$value
          })
        },min.m.i,max.m.i)$value
        
        return(cen.P1+cen.P3-cens.p)}
      
      theta<-uniroot(censor.prop,c(0,200),tol=0.00000000000000000001)$root
      return(list(theta))
    }
    
    if(censdist=="Uniform"||censdist=="Uni"||censdist=="U"){
      censor.prop<-function(theta){
        if(censdist=="Uniform"||censdist=="Uni"||censdist=="U" && upper=="maxt"){
          cen.P3<-integrate(function(u){
            sapply(u,function(u){
              integrate(function(t){
                m.i<-u
                part1<-density.fun.lambda(m.i)
                part2<-dunif(t,0,theta)
                part3<-(mixp*exp(-lambda1*(t^gamma1))+(1-mixp)*exp(-lambda2*(t^gamma2)))^m.i
                return(part1*part2*part3)
              } , 0,maxt)$value
            })
          },min.m.i,max.m.i)$value
        }
        
        if(censdist=="Uniform"||censdist=="Uni"||censdist=="U" && upper=="theta"){
          cen.P3<-integrate(function(u){
            sapply(u,function(u){
              integrate(function(t){
                m.i<-u
                part1<-density.fun.lambda(m.i)
                part2<-dunif(t,0,theta)
                part3<-(mixp*exp(-lambda1*(t^gamma1))+(1-mixp)*exp(-lambda2*(t^gamma2)))^m.i
                return(part1*part2*part3)
              } , 0,theta)$value
            })
          },min.m.i,max.m.i)$value
        }
        
        if(censdist=="Uniform"||censdist=="Uni"||censdist=="U" && upper=="maxt"){
          cen.P1<-(1-(maxt/theta))*integrate(function(u){
            sapply(u,function(u){
              m.i<-u
              part1<-density.fun.lambda(m.i)
              part3<-(mixp*exp(-lambda1*(maxt^gamma1))+(1-mixp)*exp(-lambda2*(maxt^gamma2)))^m.i
              return(part1*part3)
            } 
            )
          },min.m.i,max.m.i)$value
        }
        if(censdist=="Uniform"||censdist=="Uni"||censdist=="U" && upper=="theta"){
          cen.P1<-integrate(function(u){
            sapply(u,function(u){
              integrate(function(t){
                m.i<-u
                part1<-density.fun.lambda(m.i)
                part2<-dunif(t,0,theta)
                part3<-(mixp*exp(-lambda1*(maxt^gamma1))+(1-mixp)*exp(-lambda2*(maxt^gamma2)))^m.i
                return(part1*part2*part3)
              } , maxt,Inf)$value 
            })
          },min.m.i,max.m.i)$value
        }
        return(cen.P1+cen.P3-cens.p)
      }
      theta<-uniroot(censor.prop,c(unirootlower,unirootupper),tol=0.000000000000000000001)$root
      return(list(theta))
    }
  }
}


# --Function3：datagen(Survival Data Generator with Predefined Censoring Rate)------------
datagen <- function(seednum,gamma1,lambda1,gamma2,lambda2,mixp,cens.p,size,binomp,
                    beta1,beta2,maxt,mean,std,theta,censdist)
{set.seed(seednum)#指定种子数
  #对S(maxt)进行判断，若S(maxt)大于总删失率了，此时就不合理
  Smaxt=mixp*exp(-lambda1*(maxt^gamma1))+(1-mixp)*exp(-lambda2*(maxt^gamma2))
  #alerts1 <-  "Smaxt大于总删失率，Smaxt为"
  #alerts2 <-as.character(round(Smaxt,3)) 
  #try(if (Smaxt>cens.p) stop(paste(alerts1,alerts2)))
  covs <- data.frame(id = 1:size, trt = stats::rbinom(size, 1L, binomp),age=stats::rnorm(size,mean,std))#协变量生成
  
  if(censdist=="Uniform"||censdist=="Uni"||censdist=="U"){
    data<-data.frame(
      covs,
      simsurv(lambdas = c(lambda1, lambda2), gammas = c(gamma1, gamma2),
              mixture = TRUE, pmix = mixp, betas = c(trt = beta1,age=beta2),x = covs),
      C  <-  runif(size,0,theta))
    
  }
  if(censdist=="Exponential"||censdist=="Expo"||censdist=="E"){
    data<-data.frame(
      covs,
      simsurv(lambdas = c(lambda1, lambda2), gammas = c(gamma1, gamma2),
              mixture = TRUE, pmix = mixp, betas = c(trt = beta1,age=beta2),x = covs),
      C  <-  rexp(size,theta))
    
  }
  
  
  cens.data <- data.frame(id=data$id,theta=theta,trt=data$trt,eventtime=data$eventtime,status=data$status,censor3time=data$C,censor1time=maxt,simulationid=seednum,
                          Type3time = with(data,pmin(eventtime,C)),
                          statusT3= ifelse(data$eventtime<=data$C,1,0),
                          Type1time = with(data,pmin(eventtime,maxt)),
                          statusT1 = ifelse(data$eventtime<=maxt,1,0),
                          Type13time = with(data,pmin(eventtime,C,maxt)),
                          statusT13 = ifelse(with(data,pmin(eventtime,C,maxt))==data$eventtime, 1,0),
                          Outcome=ifelse(with(data,pmin(eventtime,C,maxt))==data$eventtime,"Event",ifelse(with(data,pmin(eventtime,C,maxt))==maxt,"Type1RC","Type3RC")),
                          T10=ifelse(data$eventtime<=maxt,1,0),
                          Smaxt=Smaxt)
  #status只有TMAX为0，censor是是只要eventtime<C就为1，status1是只要满足eventtime=10或eventtime>=C即为1，以status1作为结局变量，time做结局时间
  
  return(cens.data)
}

# --Function4：censoringcal(Calculation of the Censoring Rate for 1000 Simulated Datasets)---------

censoringcal <- function(starti,simusize,simusizeplus,seednum,gamma1,lambda1,gamma2,lambda2,
               mixp,cens.p,size,binomp,beta1,beta2,maxt,mean,std,censdist)
{
  errornum <- 0
  Cens1PControl<- c()
  CensPControl<- c()
  temp <- c()
  out <- data.frame()
  abT10 <- c()
  mm <- c()
  pb <- progress_bar$new(
    format="[:bar] Percent of Completed Iteration: :percent   ETA: :eta  ELA: :elapsedfull",
    total=simusize, clear=FALSE, width=80
  )
  for (i in starti:(starti+simusize+simusizeplus-1)){
    #if(nrow(out_13_a)==simusize ) {stop("Simulation for 1000 times reached")}
    tryCatch({
      
      mm <- data.frame(datagen(seednum=i,gamma1=gamma1,lambda1=lambda1,
                               gamma2=gamma2,lambda2=lambda2,mixp=mixp,
                               cens.p=cens.p,size=size,binomp=binomp,
                               beta1=beta1,beta2=beta2,maxt=maxt,mean=mean,std=std,theta=theta,censdist=censdist))
      mytable <- xtabs(~ trt+Outcome, data=mm)
      temp$simid <- i
      temp$CensPAll <- round(1-prop.table(mytable)[1,1]-prop.table(mytable)[2,1],5)
      out <- rbind(out,temp)
      pb$tick()
    },
    error=function(e) {
      cat("Error occur in:",i,"\n")
      cat("Total Error:",errornum <<-  errornum+1,"\n")
    } 
    )
  }
  return(out)
}