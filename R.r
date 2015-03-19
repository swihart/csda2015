## Author:  Bruce Swihart
## Date:  2013/02
## Warranty:  None

## This .R file accompanies the paper "Modeling sleep fragmentation in populations of sleep hypnograms".  


## 5 STATE / 20 Transition Types:  Poisson format / loglinear modeling;
## with one possible bootstrapping procedure coding block; 

## analysis of 5598
library("geepack")
library("gee")
library("Matrix")
## matlab library for tic()/toc() timing of function calls.
library("matlab")
library("survival")

## allows us to make nice tables of the linear combinations that ultimately are the effects of interest
linear.comb <- function(fit, cntrst.mtrx1){
	hold <- cbind(exp(esticon(fit, cntrst.mtrx1, joint.test=FALSE)[,c("Estimate", "Lower", "Upper")]),
		      esticon(fit, cntrst.mtrx1, joint.test=FALSE)[,c("Pr(>|X^2|)")])
	hold <- round(hold,2)
	colnames(hold) <- c("RR est", "95% L", "95% U", "p-value")
hold
}


## read in data
d.5598 <- read.csv("LLfull_allgroups_agesexracesmoking.csv", TRUE)

## the order I want these in the table; models, etc
## I put the line breaks where the 3-state analogues would be.
shifttype <- c("1R", "2R", "SR",                   ## 3-state analogue:  NR
               "1W", "2W", "SW",                   ## 3-state analogue:  NW
               "R1", "R2", "RS",                   ## 3-state analogue:  RN
               "RW",                               ## 3-state analogue:  RW
               "W1", "W2", "WS",                   ## 3-state analogue:  WN
               "WR",                               ## 3-state analogue:  WR
               "S1", "S2", "12", "1S", "21", "2S") ## 3-state analogue:  --; these are intra-NREM.

## establish order of types
d.5598$shift <- factor(d.5598$shift,
                       levels=c("1R", "2R", "SR", "1W", "2W", "SW", "R1", "R2", "RS", "RW", "W1", "W2", "WS", "WR", "S1", "S2", "12", "1S", "21", "2S"))

## waves argument for geeglm(); reinforces data row order.
f.wav <- d.5598$type
tic() 
## fit with geeglm()
f.off.y <- geeglm(counts~ offset(I(log(tar))) + I(factor(race)) + sex + I(smokstatus) + age + shift*grouplabel,
		  id=pptid,
		  data=d.5598,
		  family="poisson",
		  corstr="independence",
		  scale.fix=TRUE,
		  wave=f.wav,
		  control=geese.control(epsilon=1e-4, maxit=as.integer(10), trace=TRUE, scale.fix=TRUE))
toc()
## check if error-free
f.off.y$geese["error"]
summary(f.off.y)
## construct contrast matrix, picking up main effects for group-label
## and corresponding interaction term for shift*group-label For 4 SDB
## severity groups and 20 transitions, we have 3*20 = 60 estimates of
## interest.  The 3 non-referent groups each have 20 transition-type
## specific effects.  To construct an overall contrast matrix
## ("contog"), we first make matrices for each group ("con0515",
## "con1530", "con30pp", for the mild, moderate and severe SDB groups,
## respectively).  For each group, the first row has just one 1 to
## pick up the main effect for the group (the interaction term is zero
## because the first transition type is the referent), stacked on top
## of the remaining 19 rows which need to pick up the same main effect
## AND a corresponding interaction term.
p0 <- 28 ## number of leading zeros in contrast matrix.
con0515<-rbind(c(rep(0,p0),1,0,0,rep(0,19),rep(0,19),rep(0,19)),
               cbind(matrix(0,ncol=p0,nrow=19),rep(1,19),rep(0,19), rep(0,19), diag(19)*1, diag(19)*0, diag(19)*0))
con1530<-rbind(c(rep(0,p0),0,1,0,rep(0,19),rep(0,19),rep(0,19)),
               cbind(matrix(0,ncol=p0,nrow=19),rep(0,19),rep(1,19), rep(0,19), diag(19)*0, diag(19)*1, diag(19)*0))
con30pp<-rbind(c(rep(0,p0),0,0,1,rep(0,19),rep(0,19),rep(0,19)),
               cbind(matrix(0,ncol=p0,nrow=19),rep(0,19),rep(0,19), rep(1,19), diag(19)*0, diag(19)*0, diag(19)*1))
contog<-rbind(con0515, con1530, con30pp) ## stack the group specific contrast matrices 
results<-linear.comb(f.off.y, contog)    ## use our function to linearly combine; exp()
cbind(shifttype, results)                ## long table
cbind(shifttype, results[1:20,],results[21:40,],results[41:60,])  ## wide (utilize common rows)



## bootstrapping can be done to "correct" confidence intervals.  Doing so can be computationally intensive.
## This is my BEE function:  "Bootstrap Estimating Equations"
bee <- function(gee.fit, data, id, B=1000) {
	mat <- matrix(nrow=B, ncol=length(coef(gee.fit)))
	unique.ids <- unique(id)
	
	for (b in 1:B) {
		# Boostrap ids
		ids.b <- sample(unique.ids, replace=TRUE)
		dat.b <- data[as.vector(sapply(ids.b, function(x) which(id==x))),]
		dat.b$pptid <- id
		
		# Fit GEE and get coefficients
		mat[b,] <- coef(update(gee.fit, data=dat.b, control=geese.control(maxit=as.integer(1))))
	}

      mat
}


## took          168    sec   for    2 reps
## estimating at   9.33 hours for 1000 reps  (actually took:  79143 seconds = 21.98 hours).  Ouch.
tic()
coeff<-bee(f.off.y, d.5598, d.5598$pptid, B=1000)
toc()

## linearly combine and exp() as after the original geeglm call
ests <- exp(contog%*%coeff[1,])
cbind(shifttype, ests[1:20], ests[21:40], ests[41:60])

## Gives in the 60 rows, 1000 samples of my estimate(s) in the columns
hold<-exp(contog%*%t(coeff))

apply(hold, 1, function(W){quantile(W, c(.025,.975))})
cbind(shifttype,t(round(apply(hold, 1, function(W){quantile(W, c(.025,.975))}),2)))

write.csv(coeff, "geeIND_boot1000_coeff.csv", row.names=FALSE)
write.csv(hold,  "geeIND_boot1000_rr.csv", row.names=FALSE)


    



## 3 STATE /  6 Transition Types:  Poisson format / loglinear modeling;
d.5598 <- read.csv("LLfull_allgroups_agesexracesmoking_0306.csv", TRUE)

##quick calculation for example in paper: 671 has 0 RN.  Whereas 6269 have 0 RS.
table(d.5598$count[d.5598$shift=="RN"])
table(d.5598$count[d.5598$shift=="RN"]<1)

## the order I want these in the table; models, etc
shifttype <- c("NR", "NW", "RN", "RW", "WN", "WR")

## establish order of types
d.5598$shift <- factor(d.5598$shift,
       levels=c("NR", "NW", "RN", "RW", "WN", "WR"))

## waves argument
f.wav <- d.5598$type
tic()
f.off.y <- geeglm(counts~ offset(I(log(tar))) + I(factor(race)) + sex + I(smokstatus) + age + shift*grouplabel,
		  id=pptid,
		  data=d.5598,
		  family="poisson",
		  corstr="independence",
		  scale.fix=TRUE,
		  wave=f.wav,
		  control=geese.control(epsilon=1e-4, maxit=as.integer(10), trace=TRUE, scale.fix=TRUE))
toc()
f.off.y$geese["error"]
summary(f.off.y)
p0 <- 14 ##28 ##32+4
con0515<-rbind(c(rep(0,p0),1,0,0,rep(0,5),rep(0,5),rep(0,5)), cbind(matrix(0,ncol=p0,nrow=5),rep(1,5),rep(0,5), rep(0,5), diag(5)*1, diag(5)*0, diag(5)*0))
con1530<-rbind(c(rep(0,p0),0,1,0,rep(0,5),rep(0,5),rep(0,5)), cbind(matrix(0,ncol=p0,nrow=5),rep(0,5),rep(1,5), rep(0,5), diag(5)*0, diag(5)*1, diag(5)*0))
con30pp<-rbind(c(rep(0,p0),0,0,1,rep(0,5),rep(0,5),rep(0,5)), cbind(matrix(0,ncol=p0,nrow=5),rep(0,5),rep(0,5), rep(1,5), diag(5)*0, diag(5)*0, diag(5)*1))
contog<-rbind(con0515, con1530, con30pp)
results<-linear.comb(f.off.y, contog)
cbind(shifttype, results)  ## long table
cbind(shifttype, results[1:6,],results[7:12,],results[13:18,])  ## wide (utlze common rows)


## bootstrap  ... 
# This is my BEE function
bee <- function(gee.fit, data, id, B=1000) {
	mat <- matrix(nrow=B, ncol=length(coef(gee.fit)))
	unique.ids <- unique(id)
	
	for (b in 1:B) {
		# Boostrap ids
		ids.b <- sample(unique.ids, replace=TRUE)
		dat.b <- data[as.vector(sapply(ids.b, function(x) which(id==x))),]
		dat.b$pptid <- id
		
		# Fit GEE and get coefficients
		mat[b,] <- coef(update(gee.fit, data=dat.b, control=geese.control(maxit=as.integer(1))))
	}
      mat
}

tic()
coeff<-bee(f.off.y, d.5598, d.5598$pptid, B=1000)
toc()

%% tinker
ests <- exp(contog%*%coeff[1,])
cbind(shifttype, ests[1:20], ests[21:40], ests[41:60])


exp(contog%*%t(coeff[1:2,]))


hold<-exp(contog%*%t(coeff))


apply(hold, 1, function(W){quantile(W, c(.025,.975))})
cbind(shifttype,t(round(apply(hold, 1, function(W){quantile(W, c(.025,.975))}),2)))

write.csv(coeff, "geeIND_boot1000_coeff_0306.csv", row.names=FALSE)
write.csv(hold,  "geeIND_boot1000_rr_0306.csv", row.names=FALSE)


## Survival analysis only has interaction terms, so linearly combining
## as in the Poisson examples above is not required.  The naming of the
## predictors gX.tYY.ZZ is for the X group (1-mild, 2-moderate, 3-severe SDB),
## YY transition-type (1-20) (ZZ was the label version of YY).
## The line breaks occur where the analogue sets would be and a space between groups.
## The interaction terms were preprocessed; that is they are found in the accompanying datasets ready-to-go.


## 5 STATE / 20 Transition Types:  multistate survival;
d.5598 <- read.csv("SAfull_allgroups_agesexracesmoking.csv", TRUE)


## "SA 5/20";
tic()
SA520 <- coxph(Surv(time, status) ~ 
                 g1t01.1R + g1t02.2R + g1t03.SR
               + g1t04.1W + g1t05.2W + g1t06.SW
               + g1t07.R1 + g1t08.R2 + g1t09.RS
               + g1t10.RW
               + g1t11.W1 + g1t12.W2 + g1t13.WS
               + g1t14.WR
               + g1t15.S1 + g1t16.S2 + g1t17.12 + g1t18.1S + g1t19.21 + g1t20.2S
               
               + g2t01.1R + g2t02.2R + g2t03.SR
               + g2t04.1W + g2t05.2W + g2t06.SW
               + g2t07.R1 + g2t08.R2 + g2t09.RS
               + g2t10.RW
               + g2t11.W1 + g2t12.W2 + g2t13.WS
               + g2t14.WR
               + g2t15.S1 + g2t16.S2 + g2t17.12 + g2t18.1S + g2t19.21 + g2t20.2S
        
               + g3t01.1R + g3t02.2R + g3t03.SR
               + g3t04.1W + g3t05.2W + g3t06.SW
               + g3t07.R1 + g3t08.R2 + g3t09.RS
               + g3t10.RW
               + g3t11.W1 + g3t12.W2 + g3t13.WS
               + g3t14.WR
               + g3t15.S1 + g3t16.S2 + g3t17.12 + g3t18.1S + g3t19.21 + g3t20.2S
               
               + age + sex + I(factor(race)) + smokstatus
               
               
               + strata(type) + cluster(pptid),
               data=d.5598)
toc()


## output from SA520: 
##   Warning messages:
## 1: In matrix(0, n, nvar) :
##   Reached total allocation of 16296Mb: see help(memory.size)
## 2: In matrix(0, n, nvar) :
##   Reached total allocation of 16296Mb: see help(memory.size)
## > elapsed time is 25881.810000 seconds 
## > 25881/3600
## [1] 7.189167  ## took 7 hours to run on 64 bit machine
## > SA520
## Call:
## coxph(formula = Surv(time, status) ~ g1t01.1R + g1t02.2R + g1t03.SR + 
##     g1t04.1W + g1t05.2W + g1t06.SW + g1t07.R1 + g1t08.R2 + g1t09.RS + 
##     g1t10.RW + g1t11.W1 + g1t12.W2 + g1t13.WS + g1t14.WR + g1t15.S1 + 
##     g1t16.S2 + g1t17.12 + g1t18.1S + g1t19.21 + g1t20.2S + g2t01.1R + 
##     g2t02.2R + g2t03.SR + g2t04.1W + g2t05.2W + g2t06.SW + g2t07.R1 + 
##     g2t08.R2 + g2t09.RS + g2t10.RW + g2t11.W1 + g2t12.W2 + g2t13.WS + 
##     g2t14.WR + g2t15.S1 + g2t16.S2 + g2t17.12 + g2t18.1S + g2t19.21 + 
##     g2t20.2S + g3t01.1R + g3t02.2R + g3t03.SR + g3t04.1W + g3t05.2W + 
##     g3t06.SW + g3t07.R1 + g3t08.R2 + g3t09.RS + g3t10.RW + g3t11.W1 + 
##     g3t12.W2 + g3t13.WS + g3t14.WR + g3t15.S1 + g3t16.S2 + g3t17.12 + 
##     g3t18.1S + g3t19.21 + g3t20.2S + age + sex + I(factor(race)) + 
##     smokstatus + strata(type) + cluster(pptid), data = d.5598)


##                       coef exp(coef) se(coef) robust se         z       p
## g1t01.1R          0.051333     1.053 0.025571  0.047676  1.076704 2.8e-01
## g1t02.2R         -0.008662     0.991 0.013936  0.015198 -0.569921 5.7e-01
## g1t03.SR          0.069485     1.072 0.067863  0.073931  0.939854 3.5e-01
## g1t04.1W          0.185714     1.204 0.014154  0.024106  7.703975 1.3e-14
## g1t05.2W          0.095205     1.100 0.007919  0.014662  6.493418 8.4e-11
## g1t06.SW          0.034416     1.035 0.022924  0.029099  1.182730 2.4e-01
## g1t07.R1          0.227274     1.255 0.028898  0.056279  4.038363 5.4e-05
## g1t08.R2          0.002209     1.002 0.021905  0.030290  0.072930 9.4e-01
## g1t09.RS         -0.089109     0.915 0.258597  0.324726 -0.274414 7.8e-01
## g1t10.RW          0.197191     1.218 0.012814  0.022395  8.805263 0.0e+00
## g1t11.W1          0.000618     1.001 0.007798  0.015284  0.040425 9.7e-01
## g1t12.W2         -0.016763     0.983 0.010507  0.024466 -0.685129 4.9e-01
## g1t13.WS          0.000122     1.000 0.077133  0.134980  0.000905 1.0e+00
## g1t14.WR          0.167531     1.182 0.019698  0.037370  4.483068 7.4e-06
## g1t15.S1         -0.096194     0.908 0.260157  0.261271 -0.368177 7.1e-01
## g1t16.S2          0.028888     1.029 0.006712  0.011333  2.549086 1.1e-02
## g1t17.12         -0.083974     0.919 0.009074  0.016778 -5.004943 5.6e-07
## g1t18.1S         -0.251288     0.778 0.232362  0.239032 -1.051274 2.9e-01
## g1t19.21          0.211747     1.236 0.064400  0.080785  2.621104 8.8e-03
## g1t20.2S         -0.080072     0.923 0.006446  0.017631 -4.541512 5.6e-06
## g2t01.1R          0.084791     1.088 0.034605  0.067002  1.265496 2.1e-01
## g2t02.2R         -0.040070     0.961 0.020062  0.023330 -1.717539 8.6e-02
## g2t03.SR         -0.064428     0.938 0.105998  0.129175 -0.498767 6.2e-01
## g2t04.1W          0.310558     1.364 0.018358  0.033774  9.195062 0.0e+00
## g2t05.2W          0.206118     1.229 0.010684  0.022556  9.138093 0.0e+00
## g2t06.SW         -0.051814     0.950 0.035067  0.041355 -1.252903 2.1e-01
## g2t07.R1          0.295024     1.343 0.040109  0.082452  3.578133 3.5e-04
## g2t08.R2          0.028918     1.029 0.031612  0.042723  0.676868 5.0e-01
## g2t09.RS          0.504464     1.656 0.296996  0.368477  1.369050 1.7e-01
## g2t10.RW          0.316460     1.372 0.017612  0.036406  8.692633 0.0e+00
## g2t11.W1         -0.011131     0.989 0.010680  0.021737 -0.512074 6.1e-01
## g2t12.W2          0.053191     1.055 0.013978  0.032770  1.623170 1.0e-01
## g2t13.WS         -0.185127     0.831 0.112896  0.179751 -1.029905 3.0e-01
## g2t14.WR          0.198719     1.220 0.026087  0.052292  3.800175 1.4e-04
## g2t15.S1          0.059581     1.061 0.362156  0.363195  0.164046 8.7e-01
## g2t16.S2          0.031723     1.032 0.009898  0.016537  1.918269 5.5e-02
## g2t17.12         -0.142443     0.867 0.012824  0.024448 -5.826452 5.7e-09
## g2t18.1S         -0.465896     0.628 0.357861  0.392251 -1.187750 2.3e-01
## g2t19.21          0.209726     1.233 0.089929  0.114378  1.833622 6.7e-02
## g2t20.2S         -0.163465     0.849 0.009536  0.027096 -6.032872 1.6e-09
## g3t01.1R          0.250123     1.284 0.042165  0.095344  2.623382 8.7e-03
## g3t02.2R         -0.219493     0.803 0.028474  0.031382 -6.994323 2.7e-12
## g3t03.SR          0.402372     1.495 0.125769  0.142208  2.829457 4.7e-03
## g3t04.1W          0.442225     1.556 0.022617  0.055627  7.949854 1.9e-15
## g3t05.2W          0.394735     1.484 0.012809  0.042642  9.256964 0.0e+00
## g3t06.SW         -0.106631     0.899 0.053805  0.076468 -1.394450 1.6e-01
## g3t07.R1          0.763347     2.145 0.045443  0.108070  7.063433 1.6e-12
## g3t08.R2          0.064954     1.067 0.043853  0.055104  1.178757 2.4e-01
## g3t09.RS          1.053801     2.869 0.324010  0.558388  1.887218 5.9e-02
## g3t10.RW          0.230134     1.259 0.025348  0.050193  4.585019 4.5e-06
## g3t11.W1         -0.029422     0.971 0.013993  0.037796 -0.778451 4.4e-01
## g3t12.W2          0.346996     1.415 0.016023  0.048640  7.133894 9.8e-13
## g3t13.WS         -0.086156     0.917 0.137445  0.231173 -0.372691 7.1e-01
## g3t14.WR         -0.100069     0.905 0.037097  0.078893 -1.268409 2.0e-01
## g3t15.S1          0.121554     1.129 0.519840  0.513861  0.236550 8.1e-01
## g3t16.S2          0.077289     1.080 0.013962  0.021743  3.554692 3.8e-04
## g3t17.12         -0.183433     0.832 0.017211  0.027394 -6.696129 2.1e-11
## g3t18.1S         -0.643352     0.526 0.516701  0.518080 -1.241800 2.1e-01
## g3t19.21          0.666314     1.947 0.096308  0.138849  4.798850 1.6e-06
## g3t20.2S         -0.332604     0.717 0.013536  0.046203 -7.198731 6.1e-13
## age              -0.000681     0.999 0.000119  0.000246 -2.771197 5.6e-03
## sex               0.001470     1.001 0.002586  0.005069  0.289959 7.7e-01
## I(factor(race))2 -0.016677     0.983 0.004977  0.009598 -1.737470 8.2e-02
## I(factor(race))3 -0.077522     0.925 0.004464  0.009624 -8.055465 7.8e-16
## I(factor(race))4  0.030283     1.031 0.009950  0.019409  1.560305 1.2e-01
## I(factor(race))5 -0.032723     0.968 0.006143  0.011640 -2.811127 4.9e-03
## smokstatusFormer  0.039114     1.040 0.004337  0.008545  4.577186 4.7e-06
## smokstatusNever   0.037947     1.039 0.004313  0.008455  4.488375 7.2e-06

## Likelihood ratio test=4883  on 68 df, p=0  n= 2716188, number of events= 677295 


## 3 STATE / 6 Transition Types:  multistate survival;
d.5598.36 <- read.csv("SAfull_allgroups_agesexracesmoking_0306.csv", TRUE)

## "SA 3/06";
tic()
SA306 <- coxph(Surv(time, status) ~ 
               g1t01.NR
               + g1t02.NW
               + g1t03.RN
               + g1t04.RW
               + g1t05.WN
               + g1t06.WR
               + g2t01.NR
               + g2t02.NW
               + g2t03.RN
               + g2t04.RW
               + g2t05.WN
               + g2t06.WR
               + g3t01.NR
               + g3t02.NW
               + g3t03.RN
               + g3t04.RW
               + g3t05.WN
               + g3t06.WR
               
               + age + sex + I(factor(race)) + smokstatus
               
               
               + strata(type) + cluster(pptid),
               data=d.5598.36)
toc()



## Bonus:  time varying SA 5/20;
## takes a while to run:  possibly double the previous SA520 fit
## 5 STATE / 20 Transition Types:  multistate survival;
d.5598 <- read.csv("SAfull_allgroups_agesexracesmoking.csv", TRUE)

## "time-varying SA 5/20";
tic()
SA520.tt <- coxph(Surv(time, status) ~ 
                 g1t01.1R + g1t02.2R + g1t03.SR
               + g1t04.1W + g1t05.2W + g1t06.SW
               + g1t07.R1 + g1t08.R2 + g1t09.RS
               + g1t10.RW
               + g1t11.W1 + g1t12.W2 + g1t13.WS
               + g1t14.WR
               + g1t15.S1 + g1t16.S2 + g1t17.12 + g1t18.1S + g1t19.21 + g1t20.2S
               
               + g2t01.1R + g2t02.2R + g2t03.SR
               + g2t04.1W + g2t05.2W + g2t06.SW
               + g2t07.R1 + g2t08.R2 + g2t09.RS
               + g2t10.RW
               + g2t11.W1 + g2t12.W2 + g2t13.WS
               + g2t14.WR
               + g2t15.S1 + g2t16.S2 + g2t17.12 + g2t18.1S + g2t19.21 + g2t20.2S
        
               + g3t01.1R + g3t02.2R + g3t03.SR
               + g3t04.1W + g3t05.2W + g3t06.SW
               + g3t07.R1 + g3t08.R2 + g3t09.RS
               + g3t10.RW
               + g3t11.W1 + g3t12.W2 + g3t13.WS
               + g3t14.WR
               + g3t15.S1 + g3t16.S2 + g3t17.12 + g3t18.1S + g3t19.21 + g3t20.2S

               + tt(g1t01.1R) + tt(g1t02.2R) + tt(g1t03.SR)
               + tt(g1t04.1W) + tt(g1t05.2W) + tt(g1t06.SW)
               + tt(g1t07.R1) + tt(g1t08.R2) + tt(g1t09.RS)
               + tt(g1t10.RW)
               + tt(g1t11.W1) + tt(g1t12.W2) + tt(g1t13.WS)
               + tt(g1t14.WR)
               + tt(g1t15.S1) + tt(g1t16.S2) + tt(g1t17.12) + tt(g1t18.1S) + tt(g1t19.21) + tt(g1t20.2S)
               
               + tt(g2t01.1R) + tt(g2t02.2R) + tt(g2t03.SR)
               + tt(g2t04.1W) + tt(g2t05.2W) + tt(g2t06.SW)
               + tt(g2t07.R1) + tt(g2t08.R2) + tt(g2t09.RS)
               + tt(g2t10.RW)
               + tt(g2t11.W1) + tt(g2t12.W2) + tt(g2t13.WS)
               + tt(g2t14.WR
               + tt(g2t15.S1) + tt(g2t16.S2) + tt(g2t17.12) + tt(g2t18.1S) + tt(g2t19.21) + tt(g2t20.2S)
        
               + tt(g3t01.1R) + tt(g3t02.2R) + tt(g3t03.SR)
               + tt(g3t04.1W) + tt(g3t05.2W) + tt(g3t06.SW)
               + tt(g3t07.R1) + tt(g3t08.R2) + tt(g3t09.RS)
               + tt(g3t10.RW)
               + tt(g3t11.W1) + tt(g3t12.W2) + tt(g3t13.WS)
               + tt(g3t14.WR)
               + tt(g3t15.S1) + tt(g3t16.S2) + tt(g3t17.12) + tt(g3t18.1S) + tt(g3t19.21) + tt(g3t20.2S)
               
               + age + sex + I(factor(race)) + smokstatus
               
               
               + strata(type) + cluster(pptid),
               data=d.5598,
               tt=function(x,t, ...){x*log(t)})
toc()
