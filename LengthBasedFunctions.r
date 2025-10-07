# Babcock functions for length based analysis

##### Froese and Binohlan functions to convert between life history parameters
#Froese, R., Binohlan, C., 2000. Empirical relationships to estimate asymptotic length,length at first maturity and length at maximum yield per recruit in fishes, with a simple method to evaluate length frequency data. J. Fish. Biol. 56, 758-773.
#Each method converts the parameter after the period to the paramteter before
Linf.Lmax <- function(Lmax) 10^(0.044+0.9841*log10(Lmax))
Lopt.Linf <- function(Linf) 10^(1.0421*log10(Linf)-.2742) 
Lopt.Lm <- function(Lm) 10^(1.053*log10(Lm)-0.0565)
Lm.Linf <- function(Linf) 10^(0.898*log10(Linf)-0.0781)
## See the file popdynJFB.xls for the rest of the functions from this paper

##### Functions to calculate Beverton invariants
# Beverton, R.J.H., 1992. Patterns of reproductive strategy parameters in some marine teleost fishes. J. Fish. Biol. 41, 137-160.Lopt.Bev=function(Linf,M,K)  3*Linf/(3+M/K)
Lopt.Bev <- function(Linf,M,K)  3*Linf/(3+M/K)
Lm.Linf.bev <- function(Linf) 0.66*Linf

#### Functions to calculate the Froese length-based indicator
#Froese, R., 2004. Keep it simple: three indicators to deal with overfishing. Fish. Fish 5, 86-91.
PmatFunc <- function(x,Lm) {
  x=x[!is.na(x)]
  temp=length(x[x>=Lm])/length(x)
  if(is.na(Lm)|length(x)<5) temp=NA
  temp
}
PoptFunc <- function(x,Lopt) {
  x=x[!is.na(x)]
  temp=length(x[x>=Lopt*0.9 & x<=Lopt*1.1])/length(x)
  if(is.na(Lopt)|length(x)<5) temp=NA
  temp
}
PmegaFunc=function(x,Lopt) {
  x=x[!is.na(x)]
  temp=length(x[x>=1.1*Lopt])/length(x)
  if(is.na(Lopt)|length(x)<5) temp=NA
  temp
}

## Z from Lbar using either Bevertion and Holt (BH.Z) or Ehrhardt and Ault (EA) 
# Beverton, R.J.H., Holt, S.J., 1957. On the Dynamics of Exploited Fish Populations.Chapman and Hall, London.
# Ehrhardt, N.M., Ault, J.S., 1992. Analysis of two length-based mortality models applied to bounded catch length frequencies. J. Am. Fish Soc. 121, 115-122.
BH.Z <- function(K,Linf,Lave,Lc)     K*(Linf-Lave)/(Lave-Lc)
EA.Z <- function(K,Linf,Lave,Lc,Lup,Zest,Zmin) {
  f1 = function(Z) abs(((Linf-Lup)/(Linf-Lc))^(Z/K)-((Z)*(Lc-Lave)+K*(Linf-Lave))/(Z*(Lup-Lave)+K*(Linf-Lave)))
  nlminb(Zest,f1,lower=Zmin)
  # Note that the iterative function may not converge on the first try
  # Function returns a list containing Z, the objective function etc.
}

# Cope and Punt function to calculate whether population is overfished
#Cope, J.M., Punt, A.E., 2009. Length-based reference points for data-limited situations:applications and restrictions. Mar. Coastal Fish. 1, 1-18.
#Output is a vector containing status (0=not overfished, 1=overfished), selectivity category (1-5), and Px/Ptarget, which is less than 1 for overfished populations
cope <- function(Pmat,Popt,Pmega,Lmat,Lopt) {
  Pobj=Pmat+Popt+Pmega
  Px=Pmat
  if(Pobj<=1) {
    if(Popt+Pmega==0)  {
      selectivity=1 
      Ptarg=0.25
    } else {
      selectivity=2
      if(Lmat>0.825*Lopt) 
        Ptarg=0.25
      else  
        Ptarg=0.4
    }
  } 
  if(Pobj>1 & Pobj<2) {
    selectivity=3
    if(Lmat>0.825*Lopt)
      Ptarg=0.9
    else
      Ptarg=0.95
  }
  if(Pobj==2) {
    if(Popt<1) {
      selectivity=4 
      Px=Popt
      Ptarg=0.65
    } else 
      if(Popt==1) {
        selectivity=5 
        Px=NA
        Ptarg=NA
      } 
  }
  Px.Ptarg=Px/Ptarg
  if(is.na(Px.Ptarg)) status=-1 else status=ifelse(Px.Ptarg>=1,0,1)
  c("status"=status,"selectivity"=selectivity,"Px.Ptarg"=Px.Ptarg)
}


#Find Lower bound function. This function finds the mode of the length frequency distribution
#The Lc used to calculate Z from Lbar should be this number or higher
Lc.func <- function(x) {
  #  x is length sample.Returns Lc
  z=table(x)
  z1=cumsum(z)
  z1=z1/sum(z)
  a=as.numeric(names(z))
  a1=seq(trunc(min(a))+1,max(a),by=1)
  d=loess(z1~a)
  d1=predict(d,newdata=a1)
  d2=d1[2:length(a1)]-d1[1:(length(a1)-1)]
  d3=a1[1:length(d2)][d2==max(d2)]
  d3
}

#Von Bertalanffy growth cruve
vonbert<-function(age,Linf,K,t0) {
  Linf*(1-exp(-K*(age-t0)))
}


# Functions to calculate M from life history parameters
# Hewitt, D.A., Lambert, D.M., Hoenig, J.M., Lipcius, R.N., Bunnell, D.B., Miller, T.J., 2007.Direct and indirect estimates of natural mortality for Chesapeake Bay blue crab.J. Am. Fish. Soc. 136, 1030-1040.
# Jensen, A.L., 1996. Beverton and Holt life history invariants result from optimal trade-off of reproduction and survival. Can. J. Fish. Aquat. Sci. 53, 820-822.
# Methods 1-8 are from Hewitt.
M1func <- function(tm)  round(1.65/tm,3)
M2func <- function(K)   round(1.6*K,3)    
M3func <- function(tmax,K)  round(3*K/(exp(0.38*K*tmax)-1),3)
M4func <- function(tmax) round(exp(1.44-0.982*log(tmax)),3)
M5func <- function(Linf,K,T) round(10^(-0.0066-0.279*log10(Linf)+0.6543*log10(K)+0.4634*log10(T)),3) #Pauly
M6func <- function(Winf,K,T) round(10^(-.2107-0.0824*log10(Winf)+0.6757*log10(k)+0.4627*log10(T)),3) #Pauly
M7func <- function(K,tm) round(3*K/(exp(K*tm)-1),3)
M8func <- function(W)  round(3*W^(-2.88),3)        
M9func <- function(tmax) round(-log(0.01)/tmax,3) # With lower number than Ault
M10func <- function(tmax) round(-log(0.05)/tmax,3) #Alagaraja as used by Ault, J.S., Smith, S.G., Luo, J.G., Monaco, M.E., Appeldoorn, R.S., 2008. Length-based assessment of sustainability benchmarks for coral reef fishes in Puerto Rico. Environ. Conserv. 35, 221-231.
M11func <- function(tmax) round(-log(0.005)/tmax,3) # With lower bound than Ault
M12func <- function(K) round(0.21+1.45*K,3) #Jenson 
M13func <- function(K) round(exp(0.51+0.95*log(K)),3) #Jenson
#From Then et al. 2014, https://doi.org/10.1093/icesjms/fsu136
M.Then <- function(tmax) 4.899*tmax^(-.916)
M.Then.K <- function(Linf,K) 4.118*K^0.73*(Linf*10)^(-0.33) #converted to cm from mm
M.Then.1par <- function(tmax) 5.90/tmax 
#From Dureuil et al. 2021 https://doi.org/10.3354/meps13704
M.Dureuil.tmax<-function(tmax) {
  exp(1.551-1.066*log(tmax))
}


#Function to bootstrap a length frequency sample
boot.sample<-function(x) {
  sample(x,length(x),replace=TRUE)
}


#Function to bootstrap a BH estimate of Z
boot.BH<-function(x) {
  x<-sample(x,length(x),replace=TRUE)
  Lc<-Lc.func(x)
  BH.Z(K,Linf,mean(x[x>Lc]),Lc)
}

#Function to bootstrap a BH estimate of F/M, with parameter uncertainty
#defined by lognormal distribution
boot.MonteCarlo.BH<-function(x,Kmed,Ksd,Linfmed,Linfsd, Mmed,Msd) {
  K<-rlnorm(1,Kmed,Ksd)
  Linf<-rlnorm(1,Linfmed,Linfsd)
  M<-rlnorm(1,Mmed,Msd)
  x<-sample(x,length(x),replace=TRUE)
  Lc<-Lc.func(x)
  Z<-BH.Z(K,Linf,mean(x[x>Lc]),Lc)
  (Z-M)/M
}

#### data summary functions
mostfreqfunc=function(x) {
  y=table(x)
  names(y)[y==max(y)][1]
}
standard.error=function(x) {
  sd(x,na.rm=T)/sqrt(length(x))
}

## Logistic selectivity

sel.func <- function(L,SL50,SL95) {
  1/(1+exp(-log(19)*(L-SL50)/(SL95-SL50)))
}

logit<-function(p) log(p/(1-p))
ilogit<-function(x)  exp(x)/(exp(x)+1)

##lognormal to normal and back functions
lnorm.mean <- function(x1,x1e) {
  #convert lognormal to normal
  exp(x1+0.5*x1e^2)
}
lnorm.se <- function(x1,x1e) {
  #convert lognormal to normal
  ((exp(x1e^2)-1)*exp(2*x1+x1e^2))^0.5
}

#Calculate variance of sum

var.sum<-function(var.x,var.y,cov.xy) var.x+var.y+2*cov(var.xy)
