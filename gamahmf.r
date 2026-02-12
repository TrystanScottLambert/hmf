#
# Plot segim for individual object
#
library(celestial)
library(devtools)
library(mvtnorm)
library(Cairo)
library(sm)
library(Rfits)
library('magicaxis')
library('data.table')
library('plotrix')
library(foreign)
library(MASS)
#
# Define functions
#
mymagtri=function (chains, samples = 1000, thin = 1, samptype = "end", 
          grid = FALSE, do.tick = FALSE, refvals = NULL, bestvals=NULL, lab = NULL, 
          ...) 
{
  chains = as.data.frame(chains)
  chaincolnames = colnames(chains)
  Nsamp = dim(chains)[1]
  Npar = dim(chains)[2]
  if (!is.null(refvals)) {
    if (length(refvals) != Npar) {
      stop("Length of refvales must be equal to number of parameters!")
    }
  }
  if (Npar <= 1) {
    stop("Need 2+ parameters!")
  }
  if (thin > 1) {
    chains = chains[seq(1, Nsamp, by = thin), , drop = FALSE]
    Nsamp = dim(chains)[1]
  }
  if (samples > Nsamp) {
    samples = Nsamp
  }
  layout(matrix(1:Npar^2, Npar, Npar)[Npar:1, ])
  meanvec = {
  }
  sdvec = {
  }
  if (samptype == "end") {
    usesamps = (Nsamp - samples + 1):Nsamp
  }
  if (samptype == "ran") {
    usesamps = sample(Nsamp, samples)
  }
  if (samptype == "thin") {
    usesamps = seq(1, Nsamp, length = samples)
  }
  for (i in 1:Npar) {
    meanvec = c(meanvec, mean(chains[usesamps, i]))
    sdvec = c(sdvec, sd(chains[usesamps, i]))
  }
  par(oma = c(4.1, 4.1, 1.1, 1.1))
  for (i in 1:Npar) {
    for (j in 1:Npar) {
      par(mar = c(0, 0, 0, 0))
      xrange = range(chains[usesamps, i])
      yrange = range(chains[usesamps, j])
      if (xrange[1] == xrange[2]) {
        val = xrange[1]
        xrange[1] = val - 0.05
        xrange[2] = val + 0.05
      }
      if (yrange[1] == yrange[2]) {
        val = yrange[1]
        yrange[1] = val - 0.05
        yrange[2] = val + 0.05
      }
      if(i==1){
        xrange=c(12.25,14.75)
      }else if(i==2){
        xrange=c(-5.5,-2.5)
      }else if(i==3){
        xrange=c(-2,0)
      }else if(i==4){
        xrange=c(0,1)
      }
      if(j==1){
        yrange=c(12.25,14.75)
      }else if(j==2){
        yrange=c(-5.5,-2.5)
      }else if(j==3){
        yrange=c(-2,0)
      }else if(j==4){
        yrange=c(0,1)
      }
      if (i == j) {
        xtemp = chains[usesamps, i]
        if (sd(xtemp) == 0) {
          xtemp = xtemp + rnorm(samples, sd = 0.001)
        }
        plot(density(xtemp), axes = FALSE, main = "", 
             xlim = xrange,col="darkgrey")
        magaxis(1, grid = grid, grid.col = "darkgrey", 
                labels = FALSE, do.tick = do.tick)
        abline(v = bestvals[i], lty = 1, col = "blue")
        abline(v = meanvec[i], lty = 1, col = "red")
        abline(v = meanvec[i] - sdvec[i], lty = 3, col = "red")
        abline(v = meanvec[i] + sdvec[i], lty = 3, col = "red")
        if (!is.null(refvals)) {
          abline(v = refvals[i], lty = 2, col = "black")
        }
        box()
        if (i == 1) {
          plot.window(xlim = xrange, ylim = yrange)
          if (is.null(lab)) {
            magaxis(1, xlab = chaincolnames[i])
          }
          else {
            magaxis(1, xlab = lab[[i]])
          }
          if (is.null(lab)) {
            magaxis(2, ylab = chaincolnames[j])
          }
          else {
            magaxis(2, ylab = lab[[j]])
          }
        }
      }
      else {
        if (i > j) {
          plot.new()
          plot.window(xlim = xrange, ylim = yrange)
          xtemp = chains[usesamps, i]
          ytemp = chains[usesamps, j]
          if (sd(xtemp) == 0) {
            xtemp = xtemp + rnorm(samples, sd = 0.001)
          }
          if (sd(ytemp) == 0) {
            ytemp = ytemp + rnorm(samples, sd = 0.001)
          }
          magaxis(1:2, grid = grid, grid.col = "darkgrey", 
                  labels = FALSE, do.tick = do.tick)
          magcon(xtemp, ytemp, dobar = FALSE, doim = TRUE, 
                 add = TRUE, lty = c(1, 1, 1), xlim = xrange, imcol=c(NA,terrain.colors(1000,alpha=1,rev=TRUE)),
                 ylim = yrange, h = c(diff(xrange), diff(yrange))/50,col="darkgrey",
                 ...)
          points(meanvec[i], meanvec[j], col = "red", 
                 pch = 4, cex = 1)
          points(refvals[i], refvals[j], col = "black", 
                 pch = 1, cex = 1)
          points(bestvals[i], bestvals[j], col = "blue", 
                 pch = 4, cex = 1)
          box()
          abline(v = bestvals[i], lty = 1, col = "blue")
          abline(v = meanvec[i], lty = 1, col = "red")
          abline(v = meanvec[i] - sdvec[i], lty = 3, 
                 col = "red")
          abline(v = meanvec[i] + sdvec[i], lty = 3, 
                 col = "red")
          if (!is.null(refvals)) {
            abline(v = refvals[i], lty = 2, col = "black")
          }
          if (j == 1) {
            if (is.null(lab)) {
              magaxis(1, xlab = chaincolnames[i])
            }
            else {
              magaxis(1, xlab = lab[[i]])
            }
          }
        }
        else {
          plot.new()
          plot.window(xlim = xrange, ylim = yrange)
          magaxis(1:2, grid = grid, grid.col = "darkgrey", 
                  labels = FALSE, do.tick = do.tick)
          points(chains[usesamps, c(i, j)], pch = ".", 
                 col = rgb(169/255,169/255,169/255,0.1))
          points(meanvec[i], meanvec[j], col = "red", 
                 pch = 4, cex = 1)
          points(refvals[i], refvals[j], col = "black", 
                 pch = 1, cex = 1)
          points(bestvals[i], bestvals[j], col = "blue", 
                 pch = 4, cex = 1)
          box()
          if (i == 1) {
            if (is.null(lab)) {
              magaxis(2, ylab = chaincolnames[j])
            }
            else {
              magaxis(2, ylab = lab[[j]])
            }
          }
        }
      }
    }
  }
  output = cbind(mean = meanvec, sd = sdvec)
  rownames(output) = chaincolnames
  return(invisible(output))
}
#
massfn <- function(x){
  mstar <- x[1]
  phi <- x[2]
  alpha <- x[3]
  beta <- x[4]
  allxxx=c(max(allx)+seq(1,10,1)*logbin)
  penalty=2*vlimit*sum(log(10)*beta*exp(-10^(beta*(allxxx-mstar)))*(phi*(10^allxxx/10^mstar)^(alpha+1)))*logbin
  model=log10(beta*log(10)*(exp(-10^(beta*(allx-mstar))))*(phi*((10^allx/10^mstar)^(alpha+1))))
  sum(((ally-model)/(allf/log(10)))^2)+penalty
  }
#
# Define constants
#
ho=67.37
omegam=0.3147
omegamerr=0.0074
omegal=1-omegam
G=6.67408E-11
msol=1.988E30
parsec=3.0857E16
cosvar=function(V,N){((219.7-52.4*log10(V)+3.21*(log10(V))^2)/N^0.5)/100.0}
logbin=0.2
fitbinwid=0.01
rhocrit=3*(1000*ho/(1E6*parsec))^2/(8*pi*G)
RAmid=c(135,180,217,5)
name=c("G09","G12","G15")
zlimit=0.25
dlimit=cosdist(zlimit,OmegaM=omegam,OmegaL=omegal,H0=ho)$CoDist
vlimit=179.92/(360^2/pi)*1E9*cosdist(zlimit,OmegaM=omegam,OmegaL=omegal,H0=ho)$CoVol
vlimitmin=vlimit/1000.0
#
# Input arguments
#
#inputargs=commandArgs(TRUE)
#multi=as.integer(inputargs[1])
#mlimit=as.numeric(inputargs[2])
#myoption=as.character(inputargs[3])
#zmin=as.character(inputargs[4])
magica=13.9
multi=5
mlimit=12.7
myoption="GAMA"
zmin=0.015
#
# Prepared HMFcalc curve
#
betamrp=0.7097976
A=1.727006E-19
mstarmrp=14.42947
alphamrp=-1.864908
mrpx=seq(0,17,0.001)+log10(100/ho)
mrpy=A*betamrp*10^((alphamrp+1)*(mrpx-mstarmrp))*exp(-10^(betamrp*(mrpx-mstarmrp)))*(ho/100)^3
factor=sum(10^mrpx*mrpy)*0.001*msol/(1E6*parsec)^3/(omegam*rhocrit)
phimrp=A/factor
#
# Read GAMA cats
#
g3cx=Rfits_read_table("/Users/sdriver/Drpbx/active/hmf/G3CFoFGroupv10.fits")
g3c=g3cx[Nfof > multi-1 & Zfof < zlimit & Zfof > zmin & MassAfunc > 1E1 & IterCenDec > -3.5,]
g3c$MassAfunc=g3c$MassAfunc*100/ho
g3c$MassA=g3c$MassA*100/ho
g3c$MassWL=(1.E14*100/ho)*(g3c$VelDisp/500)^1.89
#
g3c$mymass=magica*(g3c$VelDisp*1000)^2*g3c$Rad50*parsec*1E6/(G*msol)*(100/ho)
xx=seq(3,22)
yy=c(0.68389355,0.38719116,0.40325591,0.32696735,0.27680685,00.24018684,0.20226682,0.18645475,0.17437005,0.14271506,0.13922450,0.13482418,0.13741619,0.11715141,0.12134983,0.10078830,0.09944761,0.09913166,0.08590223,0.07588408)
if(myoption=="OrigErr"){
  xx=c(3.0,4.3,6.25,8.1,12.1,19.6,31.0,47.5,68.9,81.9,100.0)
  yy=c(0.866,0.763,0.716,0.679,0.572,0.419,0.310,0.246,0.189,0.171,0.126)
}
g3c$log10MassErr=approx(xx,yy,g3c$Nfof)$y
g3c$log10MassErr[is.na(g3c$log10MassErr)]=0.03
g3c$log10MassErr[g3c$log10MassErr<0.1]=0.1
masscorr=c(0.0,0.0,-2.672595e-01,-1.513503e-01,-1.259069e-01,-9.006064e-02,-5.466009e-02,-6.666895e-02,-1.988694e-02,-2.439581e-02,-2.067060e-02,-1.812964e-02,-1.556899e-02,-1.313664e-02,-1.743112e-02,-7.965513e-03,-1.257178e-02,-7.064037e-03,-3.963656e-03,-1.271533e-02,-2.664687e-03,-1.691287e-03)
g3c$masscorr=masscorr[g3c$Nfof]
g3c$masscorr[is.na(g3c$masscorr)]=0.0
g3c$mymasscorr=g3c$mymass/10^g3c$masscorr
#
# Which mass to use?
#
if(myoption=="MassA"){
  g3c$MassAfunc=g3c$MassA
}else if(myoption=="MassAfunc"){
  g3c$MassAfunc=g3c$MassAfunc
}else if(myoption=="MassWL"){
  g3c$MassAfunc=g3c$MassWL
}else if(myoption=="AngDiam"){
  g3c$MassAfunc=g3c$mymass/(1+g3c$Zfof)
}else if(myoption=="GAMA" | myoption=="OrigErr"){
  g3c$MassAfunc=g3c$mymasscorr
}
#
# Determine Vmax values
#
gig=fread("/Users/sdriver/Drpbx/active/hmf/GAMAGalsInGroups.csv")
for (i in 1:length(g3c$GroupID)){
  if(g3c$Nfof[i]==2){
    g3c$zmax[i]=sort(gig[GroupID==g3c$GroupID[i],zmax_19p8],decreasing=TRUE)[2]
  }else{
    g3c$zmax[i]=sort(gig[GroupID==g3c$GroupID[i],zmax_19p8],decreasing=TRUE)[multi]
  }}
g3c$cod=cosdist(g3c$Zfof,OmegaM=omegam,OmegaL=omegal,H0=ho)$CoDist
g3c$zmax=ifelse(g3c$zmax<g3c$Zfof,g3c$Zfof,g3c$zmax)
g3c$zmax=ifelse(g3c$zmax>zlimit,zlimit,g3c$zmax)
g3c$vmax=179.92/(360^2/pi)*1E9*cosdist(g3c$zmax,OmegaM=omegam,OmegaL=omegal,H0=ho)$CoVol
g3c$vmax=g3c$vmax-179.92/(360^2/pi)*1E9*cosdist(zmin,OmegaM=omegam,OmegaL=omegal,H0=ho)$CoVol
g3c$dmax=cosdist(g3c$zmax,OmegaM=omegam,OmegaL=omegal,H0=ho)$CoDist
g3c$weightszlimit=ifelse(g3c$vmax>vlimit,vlimit,g3c$vmax)
g3c$weightszlimit=ifelse(g3c$vmax<vlimitmin,vlimitmin,g3c$vmax)
g3c$cosvar=cosvar(g3c$vmax,1)
g3c$MassAfunc[g3c$GroupID==100622]=1E9
#
# GAMA HMF and cosmic variance per galaxy
#
massx=seq(10.3,16.1,logbin)
gamahmf=maghist(log10(g3c$MassAfunc),breaks=massx,plot=FALSE,verbose=FALSE)
gamahmf2=weighted.hist(log10(g3c$MassAfunc),w=1/g3c$weightszlimit,breaks=massx,plot=FALSE)
#cosvariance=((weighted.hist(log10(g3c$MassAfunc),w=g3c$cosvar,breaks=massx,plot=FALSE)$counts))/gamahmf$counts
#cosvariance[is.na(cosvariance)]=0
cosvariance=cosvar(vlimit/3,3)
#
# Looks at vmax values
#
png(filename=paste0("/Users/sdriver/Drpbx/active/hmf/dlimits.png"),width=18.0,height=12.0,units="cm",res=240)
g3c$dmax=cosdist(g3c$zmax,OmegaM=omegam,OmegaL=omegal,H0=ho)$CoDist
magplot(log10(g3c$MassAfunc),g3c$vmax,xlim=c(10,16),grid=FALSE,pch=".",xlab=expression("Halo Mass (M"["\u0298"]~")"),ylab=expression("Comoving Volume limit (Mpc"^3~")"))
#points(log10(g3c$MassAfunc[g3c$MassAfunc>10^12.7 & g3c$MassAfunc<10^12.9]),g3c$vmax[g3c$MassAfunc>10^12.7 & g3c$MassAfunc<10^12.9],pch=4,col="blue")
medianvmax=0.0
for (i in 2:length(massx)){medianvmax[i-1]=mean(g3c$vmax[g3c$MassAfunc > 10^massx[i-1] & g3c$MassAfunc < 10^massx[i]])}
points(gamahmf$mids,medianvmax,col="green",pch=16,cex=1)
lines(gamahmf$mids,medianvmax,col="green",lwd=1)
g3c$vmax=ifelse(g3c$vmax>vlimit,vlimit,g3c$vmax)
for (i in 2:length(massx)){medianvmax[i-1]=median(g3c$vmax[g3c$MassAfunc > 10^massx[i-1] & g3c$MassAfunc < 10^massx[i]])}
points(gamahmf$mids,medianvmax,col="red",pch=19,cex=1)
lines(gamahmf$mids,medianvmax,col="red",lwd=1)
tmp=dev.off()
gamahmf3=gamahmf$counts/medianvmax
#
#
#
png(filename=paste0("/Users/sdriver/Drpbx/active/hmf/gamahmf",myoption,multi,".png"),width=18.0,height=12.0,units="cm",res=240)
#
par(fig=c(0,1,0.8,1),mar=c(0,3.35,0.5,0.25),oma=c(0,0,0,0))
magplot(gamahmf$mids-0.5*logbin,gamahmf$counts,grid=FALSE,xlim=c(12,16),type="s",ylab="Number",majorn=c(5,2),labels=c(0,1))
text(16.2,0.9*max(gamahmf$counts),paste0("GAMA groups at z < ",zlimit),pos=2,cex=1.0)
text(16.2,0.6*max(gamahmf$counts),paste0("Multiplicity > ",multi-1),pos=2,cex=1.0)
text(16.2,0.3*max(gamahmf$counts),paste0("N groups = ",sum(gamahmf$counts)),pos=2,cex=1.0)
#
par(fig=c(0,1,0,0.795),mar=c(3.0,3.35,0.25,0.25),oma=c(0,0,0,0),new=TRUE)
magplot(0,0,xlim=c(12,16),ylim=c(-8,-2),grid=FALSE,xlab=expression("Halo Mass (M"["\u0298"]~")"),ylab=expression("log"[10]~"(number density) [Mpc"^-3~"dex"^-1~"]"))
#lines(mrpx,log10(mrpy)-log10(factor),lty=2,lwd=2)
#lines(mrpx-0.125,log10(mrpy)-log10(factor)+0.11,lty=2,lwd=2)
lines(mrpx-0.08,log10(mrpy)-log10(factor)+0.08,lty=2,lwd=2)
#
# Monty Carlo for edB
#
meancounts=0.0
mediancounts=0.0
upcounts=0.0
docounts=0.0
mcerr=0.0
mockcounts=matrix(0,nrow=1001,ncol=length(gamahmf$mids))
for (i in 1:1001){
  mockmass=log10(g3c$MassAfunc)+rnorm(length(g3c$log10MassErr),0.0,g3c$log10MassErr)
  mockcounts[i,1:length(gamahmf$mids)]=weighted.hist(mockmass,w=1/g3c$weightszlimit,breaks=massx,plot=FALSE)$counts
#  points(gamahmf$mids+rnorm(1,0,0.01),log10(mockcounts[i,1:length(gamahmf$mids)]/logbin),pch=".",col=rgb(0.5,0.5,0.5,0.1))
  }
# meancounts
for(i in 1:length(gamahmf$mids)){meancounts[i]=mean(mockcounts[,i])}
for(i in 1:length(gamahmf$mids)){mediancounts[i]=quantile(mockcounts[,i],0.5,na.rm=TRUE)}
for(i in 1:length(gamahmf$mids)){upcounts[i]=quantile(mockcounts[,i],0.84,na.rm=TRUE)}
for(i in 1:length(gamahmf$mids)){docounts[i]=quantile(mockcounts[,i],0.16,na.rm=TRUE)}
# edb
edb=meancounts/gamahmf2$counts
edb[is.infinite(edb)]=1.0
edb[is.na(edb)]=1.0
#
# MC err
#
for(i in 1:length(gamahmf$mids)){mcerr[i]=quantile((meancounts[i]-mockcounts[,i])^2,0.66)^0.5/(meancounts[i])}
#
# Poisson Error
#
rootnerr=(1/gamahmf$counts^0.5)
#
# Generate edB corrected HMF and errors
#
gamax=gamahmf2$mids
gamay=gamahmf2$counts/(logbin*edb)
gamaf=((mcerr^2+rootnerr^2)^0.5)
gamaf[is.na(gamaf)]=0.9999
gamaf[is.infinite(gamaf)]=0.0
gamaf=ifelse(gamaf>=1,0.9999,gamaf)
#
allx=c(gamax[gamay > 0 & !is.na(gamay) & gamax>mlimit])
ally=c(log10(gamay[gamay > 0 & !is.na(gamay) & gamax>mlimit]))
allf=c(gamaf[gamay > 0 & !is.na(gamay) & gamax>mlimit])
gamafit <- optim(par=c(mstarmrp,phimrp,alphamrp,betamrp),fn=massfn,control=list(maxit=500,reltol=1e-8,parscale=c(1,1,1,0.5))) # parscale=c(1,1,1,10)
#
xfit=seq(0.0,20,fitbinwid)
yfit=gamafit$par[4]*log(10)*(exp(-10^(gamafit$par[4]*(xfit-gamafit$par[1]))))*(gamafit$par[2]*((10^xfit/10^gamafit$par[1])^(gamafit$par[3]+1)))
gamaomegamatter=sum((yfit*10^xfit)*fitbinwid)*msol/(1E6*parsec)^3/rhocrit
#
mstar=0
phistar=0
alphastar=0
betastar=0
#
for (i in 1:10001){
  cv=gamay*rnorm(length(gamaf),0.0,cosvariance)
  mockgamay=gamay+gamay*rnorm(length(gamaf),0.0,gamaf)+cv
  allx=c(gamax[mockgamay > 0 & !is.na(gamay) & gamax>mlimit])
  ally=c(log10(mockgamay[mockgamay > 0 & !is.na(gamay) & gamax>mlimit]))
  allf=c(gamaf[mockgamay > 0 & !is.na(gamay) & gamax>mlimit])
  montyfit <- optim(par=c(mstarmrp,phimrp,alphamrp,betamrp),fn=massfn,control=list(maxit=500,reltol=1e-8,parscale=c(1,1,1,0.5)))
  ymonty=montyfit$par[4]*log(10)*(exp(-10^(montyfit$par[4]*(xfit-montyfit$par[1]))))*(montyfit$par[2]*((10^xfit/10^montyfit$par[1])^(montyfit$par[3]+1)))
  if (i<1001){lines(xfit,log10(ymonty),col=rgb((100/255),(149/255),(237/255),0.01))}
#    points(gamax,log10(mockgamay),pch=".",col=rgb(0.5,0.5,0.5,0.1))
    mstar[i]=montyfit$par[1]
    phistar[i]=montyfit$par[2]
    alphastar[i]=montyfit$par[3]
    betastar[i]=montyfit$par[4]
}
#
points(gamax,log10(gamahmf2$counts/logbin),pch=5,col="limegreen",cex=0.75)
#points(gamax,log10(gamahmf2$counts/logbin),pch="|",cex=0.75,col="limegreen")
points(gamax[gamax>mlimit],log10(gamay[gamax>mlimit]),pch=16,col="red")
lines(xfit,log10(yfit),col=rgb((100/255),(149/255),(237/255),1),lwd=2)
points(gamax[gamax<mlimit],log10(gamay[gamax<mlimit]),pch=1,col="red",cex=1)
magerr(gamax,log10(gamay),ylo=log10(1-gamaf),yhi=log10(1+gamaf),col="red")
#
lines(c(14.4,14.6),c(-3.7,-3.7),col="black",lty=2,lwd=2)
text(14.5,-3.7,"  LCDM prediction",col="black",pos=4)
lines(c(14.4,14.6),c(-3.2,-3.2),col=rgb((100/255),(149/255),(237/255)),lwd=2)
text(14.5,-3.2,"  Best MRP function fit",col=rgb((100/255),(149/255),(237/255)),pos=4)
points(14.5,-2.7,pch=5,cex=0.75,col="limegreen")
#points(14.5,-2.7,pch="|",col="limegreen",cex=0.75)
text(14.5,-2.7,paste0("  GAMA z<",zlimit," Raw HMF"), col="limegreen",pos=4,cex=1)
points(14.5,-2.2,pch=16,col="red",cex=1)
text(14.5,-2.2,paste0("  GAMA z<",zlimit," Corrected HMF"), col="red",pos=4,cex=1)
if(myoption=="OrigErr"){text(12,-7,"with original errors",pos=4)}
if(myoption=="AngDiam"){text(12,-7,"using angular diameter distance",pos=4)}
if(myoption=="MassAfunc"){text(12,-7,"using MassAfunc",pos=4)}
if(myoption=="MassA"){text(12,-7,"using MassA (with A=10)",pos=4)}
if(myoption=="MassWL"){text(12,-7,"using Weak Lensing masses (Viola et al 2015)",pos=4)}
tmp=dev.off()
#
# write
#
gamatable=paste(rev(gamax),rev(round(gamahmf$counts,3)),rev(round(log10(gamahmf2$counts),3)),rev(round(log10(gamay),3)),rev(round(rootnerr,3)),rev(round(mcerr,3)),rev(round(cosvariance,3)),rev(round(gamaf,3)),sep=" $&$ ")
write(gamatable,file=paste0("/Users/sdriver/Drpbx/active/hmf/gamatable",myoption,multi,".txt"),append=FALSE)
gamatable2=as.data.frame(cbind(rev(gamax),rev(gamahmf$counts),rev(log10(gamahmf2$counts)),rev(log10(gamay)),rev(rootnerr),rev(mcerr),rev(cosvariance),rev(gamaf)))
write.csv(gamatable2,paste0("/Users/sdriver/Drpbx/active/hmf/gamahmf",myoption,multi,".csv"),row.names=FALSE)
#
# Ed bias plot?
#
png(filename=paste0("/Users/sdriver/Drpbx/active/hmf/edbias",multi,".png"),width=8.0,height=8.0,units="cm",res=240)
par(mar=c(3,3,0.25,0.25),oma=c(0,0,0,0))
#
x=c(0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0)
y=c(1.0,1.092,1.47,2.24,3.27,5.37,8.52,11.16,13.74,16.17,19.32)
magplot(x,y,xlab=expression(sigma["log10(Mass)"]),ylab=expression("Eddington Bias (w"[i]~")"),xlim=c(0,0.8),ylim=c(0,10))
lines(x,y)
x=seq(0,1,0.025)
masshist=maghist(g3c$log10MassErr[g3c$Zfof < zlimit & g3c$MassAfunc > 1E13],breaks=x,plot=FALSE,verbose=FALSE)
points(masshist$mids,masshist$counts/150,type="s")
#
text(0.0,9,expression("N"[fof]~">2"),pos=4)
tmp=dev.off()
#
# Covariance plot
#
png(filename=paste0("/Users/sdriver/Drpbx/active/hmf/gamacov",myoption,multi,".png"),width=16.0,height=16.0,units="cm",res=240)
par(mar=c(3,3.5,0.25,0.25),oma=c(0,0,0,0))
params=cbind(mstar[!is.na(mstar)],log10(phistar[!is.na(mstar)]),alphastar[!is.na(mstar)],betastar[!is.na(mstar)])
fitvals=mymagtri(params,samples=length(params[,1]),samptype="ran",
               #               refvals=c(allfit3$par[1],log10(allfit3$par[2]),allfit3$par[3],allfit3$par[4]),
               lab=c("log\u2081\u2080(M*) [M\u2092]","log\u2081\u2080(\u03c6*) [Mpc\u207B\u00b3]","\u03b1","\u03b2"),
               refvals=c(mstarmrp-0.075,log10(phimrp)+0.075,alphamrp,betamrp),
               bestvals=c(gamafit$par[1],log10(gamafit$par[2]),gamafit$par[3],gamafit$par[4]))
cat(paste0(myoption,multi),
    gamafit$par[1],quantile(mstar,0.16,na.rm=TRUE),quantile(mstar,0.84,na.rm=TRUE),
    log10(gamafit$par[2]),log10(quantile(phistar,0.16,na.rm=TRUE)),log10(quantile(phistar,0.84,na.rm=TRUE)),
    gamafit$par[3],quantile(alphastar,0.16,na.rm=TRUE),quantile(alphastar,0.84,na.rm=TRUE),
    gamafit$par[4],quantile(betastar,0.16,na.rm=TRUE),quantile(betastar,0.84,na.rm=TRUE),"\n")
tmp=dev.off()
