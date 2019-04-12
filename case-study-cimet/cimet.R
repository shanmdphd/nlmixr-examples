## Source Tomoo Funaki
## Modified by Matthew Fidler

library(nlmixr)
library(ggplot2)

##str(theo_sd)
cimet_sd <- readRDS("cimet_sd.rds")
#cimet_sd$CMT <- ifelse(cimet_sd$CMT==2, 3, cimet_sd$CMT)
cimet_sd$dose <- as.numeric(300)
is.numeric(cimet_sd$dose)


ggplot(cimet_sd, aes(TIME, DV)) +
  geom_line(aes(group=ID), col="red") +
  scale_x_continuous("Time (h)") +
  scale_y_continuous("Concentration") +
    labs(title="Cimetidine single-dose", subtitle="Concentration vs. time by individual")


## solved system (pred)


one.compartment.ehc <- function() {
  ini({
    tka <- log(0.5)
    tcl <- log(60)
    tv <- log(25)
    ttgap <- fix(log(2))
    trkeb <- log(0.5)

    eta.cl ~ 0.1
    eta.ka ~ 0.1
    eta.v ~ 0.1
    eta.tgap ~ 0.1
    eta.rkeb ~ 0.1

    #prop.err <- 0.1
    add.err <- 0.01
  })
  model({
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl + eta.cl)
    v <- exp(tv + eta.v)
    tgap <- exp(ttgap + eta.tgap)
    rkeb <- exp(trkeb + eta.rkeb)

    #tgap <- 2
    bile <- 1
    if (t < tgap){
      bile <- 0
    }

    ha <- exp(-(cl/v)*tgap)/((cl/v) - ka)
    hb <- exp(-ka*tgap)*(cl/v)/ka/((cl/v) - ka)
    tote <- ka*dose*(1/ka + ha - hb)

    hc <- exp(-(cl/v)*t) - exp(-ka*t)
    timh <- bile*(t - tgap)
    hd <- exp(-(cl/v)*timh) - exp(-ka*timh)

    cp <-dose/v*ka/(ka-(cl/v))*hc + bile*rkeb*tote/v*ka/(ka - (cl/v))*hd

    cp ~ add(add.err) # + prop(prop.err)
  })
}

fit <- nlmixr(one.compartment.ehc, cimet_sd, est="saem")
## Currently CWRES/FOCEi is disabled due to time-based derivatives on conditions.
## fit.one.ehc <- nlmixr(one.compartment.ehc, cimet_sd, est="focei")

AIC(fit)  ## Calculates Objective Function

pdf("cimet-solve.pdf")
plot(fit)
dev.off();


one.compartment.saem <- function() {
  ini({
    tktr <- log(0.1)
    tka <- log(0.5)
    tcl <- log(60)
    tv <- log(25)
    tfgb <- fix(log(0.1))
    tto <- fix(log(1.5))
    ttc <- fix(log(3))

    eta.ktr ~ 1
    eta.ka ~ 1
    eta.cl ~ 2
    eta.v ~ 1

    eta.fgb ~ 0.5
    eta.to ~ 0.5
    eta.tc ~ 0.5

    add.err <- 0.1
    #prop.err <- 0.1
  })
  model({
    ktr <- exp(tktr + eta.ktr)
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl + eta.cl)
    v <- exp(tv + eta.v)

    fgb <- exp(tfgb + eta.fgb)
    to <- exp(tto + eta.to)
    tc <- exp(ttc + eta.tc)

    Q1 <- 0

    #if (t >= 1.5 && t < 3) {
    if (t >= to && t < tc) {
      Q1 <- 1
    }
    else {
      Q1 <- 0
    }

    d/dt(depot) = -ktr * depot
    d/dt(gut) =  ktr * depot -ka * gut
    d/dt(center) = ka * gut - (cl / v) * center + Q1 * gb
    d/dt(gb) = (cl / v * fgb * center)
    cp = center / v
    #cp ~ prop(prop.err) + add(add.err)
    cp ~ add(add.err)
  })
}


# Model fitting
fit.one <- nlmixr(one.compartment.saem, cimet_sd, est="saem", table=list(npde=TRUE))

AIC(fit.one);

print(fit.one);

# Show results

pdf("fit-ode.pdf")
plot(fit.one)
dev.off();



