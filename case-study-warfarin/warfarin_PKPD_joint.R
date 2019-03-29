## Example courtesy of Tomoo Funaki and Nick Holford
## Some changes from Matthew Fidler and Wenping Wang
## Load libraries
library(nlmixr)
library(xpose.nlmixr) # Needed for GOF
nlmixrver=packageVersion("nlmixr")
nlmixrver

# Read data
homedir="."
data.pkpd <- read.csv(file = paste(homedir,"warfarin_dat.csv",sep="/"))
data.pkpd$cmt  <- factor(data.pkpd$dvid,c(1,2),c("center","effect"))

pk.turnover.emax <- function() {
  ini({
    tktr <- log(0.01)
    tka <- log(1)
    tcl <- log(0.1)
    tv <- log(1)

    eta.ktr ~ 1
    eta.ka ~ 1
    eta.cl ~ 2
    eta.v ~ 1
    prop.err <- 1
    pkadd.err <- 0.00002

    poplogit <- 2
    #temax <- 7.5
    tec50 <- log(0.5)
    tkout <- log(0.05)
    te0 <- log(100)

    eta.emax ~ .5
    eta.ec50  ~ .5
    eta.kout ~ .5
    eta.e0 ~ .5

    pdadd.err <- 4
  })
  model({
    ktr <- exp(tktr + eta.ktr)
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl + eta.cl)
    v <- exp(tv + eta.v)

    #poplogit = log(temax/(1-temax))
    logit=exp(poplogit+eta.emax)
    #logit=temax+eta.emax
    emax = logit/(1+logit)
    ec50 =  exp(tec50 + eta.ec50)
    kout = exp(tkout + eta.kout)
    e0 = exp(te0 + eta.e0)

    DCP = center/v
    PD=1-emax*DCP/(ec50+DCP)

    effect(0) = e0
    kin = e0*kout

    d/dt(depot) = -ktr * depot
    d/dt(gut) =  ktr * depot -ka * gut
    d/dt(center) =  ka * gut - cl / v * center
    f(center)=1
    lag(center)=0
    rate(center)=0
    dur(center)=0
    d/dt(effect) = kin*PD -kout*effect

    cp = center / v
    cp ~ prop(prop.err) + add(pkadd.err) | center
    effect ~ add(pdadd.err) | effect
  })
}

nlmixr(pk.turnover.emax)
fit.TOF <- nlmixr(pk.turnover.emax, data.pkpd, est="saem")

pdf("warfarin-joint-saem.pdf")
plot(fit.TOF)
vpc.ui(fit.TOF, stratify="CMT", show=list(obs_dv=T))
dev.off()


fit.TOF2 <- nlmixr(pk.turnover.emax, data.pkpd, est="focei")
fit.TOF2
pdf("warfarin-joint-focei.pdf");
plot(fit.TOF2)
vpc.ui(fit.TOF2, stratify="CMT", show=list(obs_dv=T))
dev.off();
