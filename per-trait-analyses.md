MCMCglmm: analysing traits separatly
================

``` r
load("~/Documents/Analyses/Post-zygotic mechanisms/Model/mglmm1.RData")
library(tidyr)
library(MCMCglmm)
library(ggplot2)
library(plyr)
library(QGglmm)
library(gridExtra)
library(tidybayes)
```

Import dataset

``` r
mod.length <- read.csv("PATH.TO/dataset.per.trait.analyses.csv", h =T)
```

# I. Variation in length across time and cross types

## 1.b. Random regression

``` r
prior1.2<-list(G=list(
  G1=list(V=1,nu=0.002),
  G2=list(V=diag(2),alpha.mu=c(0,0), alpha.V=diag(2)*1000,nu=1)),R=list(V=1,nu=0.002))

iter<-5
rr.log.length<-MCMCglmm(log(length)~ cross2*poly(mc.dd,2) - 1,
                       random =~ fam + us(1 + mc.dd):ind, prior = prior1.2, rcov = ~units, data = mod.length[mod.length$fam != "SBPL8" & complete.cases(mod.length$age.dd),],verbose = T,nitt = 13000*iter,thin = 10*iter, burnin = 3000*iter)
summary(rr.log.length)
```

Predicted growth trajectories and credible intervals

``` r
summary(rr.log.length)
cross2<-c("PL","SB","hybrid")
age<-439:1109
pred.val<-expand.grid(cross2,age)
pred.data<-mod.length[mod.length$fam != "SBPL8" & complete.cases(mod.length$age.dd),]
pred.data<-rbind(pred.data,pred.data,pred.data,pred.data)
pred.data<-pred.data[1:nrow(pred.val),]
pred.data$cross2<-pred.val[,1]
pred.data$mc.dd<-pred.val[,2]
pred.rr<-predict.MCMCglmm(rr.log.length,newdata = pred.data, interval = "confidence")
growth.df<-data.frame(age = pred.data$mc.dd,
                      cross2 =pred.data$cross2,
                      predict = pred.rr[,1], 
                      lwr = pred.rr[,2], 
                      upr = pred.rr[,3])
p.gth<-ggplot(data=growth.df, aes(x=age, y=predict, colour=cross2)) +  geom_line(size = 0.5) +geom_ribbon(aes(ymin=growth.df$lwr, ymax=growth.df$upr,fill= growth.df$cross2),linetype=0, alpha=0.1) + theme_classic() + 
  xlab("Growing Degree Day") + ylab("log(body length in mm)") +
  theme(axis.text = element_text(size = 15))
```

# II. Variations in yolksac resorption

Here we consider the volume of the yolksac as the area of a curve (mm2)
traced on of each individual at two developmental time points (hatching
and hatching +20 days).

``` r
dfa3 <- <- read.csv("~/Documents/Writtings/Hybrid_traits/BMC/dataset.yolksac.csv", h = T)
```

#### Fit a Multi response model

``` r
prior.yk<-list(G=list(
  G1=list(V=1,nu=2),
  G2=list(V=diag(2),nu=2)),
  R=list(V=diag(2),nu=2))
itr <- 300
m.yolk<-MCMCglmm(cbind(log(D1_yolk.area),log(D2_yolk.area))~trait-1 + 
                   at.level(trait,1):scale(log(D1_length)) + # log body length at D1  
                   at.level(trait,2):scale(log(D2_length)) +  # log body length at D2
                   trait:cross2,
               random =~fam + us(trait):ind, prior = prior.yk, rcov = ~us(trait):units, 
               data = dfa3[complete.cases(dfa3$D2_relative.yolk),], family = c("gaussian","gaussian"),
               nitt = 13000*itr,thin = 10*itr, burnin = 3000*itr,verbose = T)
```

# III : First feeding date and Feeding behaviour

##### 4.1 Variation in the date of first feeding

``` r
plot(mod.length$fam,mod.length$degree.day.frist.feeding, main = "Among family variation in age at first feeding", ylab = "Degree days")
```

… recall that SBPL8 was fertillized and hatched 15 days later than the
rest. They however started feeding at about the same date than the
others, hence the difference in age at first feeding.

``` r
pf <-list(G=list(
  G1=list(V=1,nu=0.002)),
  R=list(V=diag(3),nu=0.002))
itr<-20
m.ff<-MCMCglmm(degree.day.frist.feeding~cross2,
               random =~fam, prior = pf, rcov = ~idh(cross2):units, data = mod.length[!duplicated(mod.length$ind), ], verbose = F)
# Without SBPL8
m.ff2<-MCMCglmm(degree.day.frist.feeding~cross2,
               random =~fam, prior = pf, rcov = ~idh(cross2):units, 
               data = mod.length[!duplicated(mod.length$ind) &
                                   mod.length$fam != "SBPL8", ], nitt = 13000*itr,thin = 10*itr, burnin = 3000*itr,verbose = T)
```

###### Residual components

``` r
ff.res<-data.frame(p.mode=posterior.mode(m.ff2$VCV)[-1], 
                     lowr = HPDinterval(m.ff2$VCV)[-1,1],
                     upp =  HPDinterval(m.ff2$VCV)[-1,2],
                     cross = names(posterior.mode(m.ff2$VCV)[-1]))

p.ff.res<-ggplot(ff.res, aes(x=cross, y=p.mode)) +
  geom_hline(yintercept = 0.0, color = "lightgrey") +
  geom_errorbar(aes(ymin=lowr, ymax=upp), colour="cyan4",width=.01, size = 0.8) +
  geom_point(size = 2) + theme_bw() +
  scale_x_discrete(limits = c("cross2PL.units","cross2SB.units","cross2hybrid.units"),
                   labels= c("PL","SB","Hybrid")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.text.x = element_text(size = 10, angle = 90),
        axis.title.y = element_text(size = 10),
        strip.background = element_rect(fill = "white", colour = "white")) + 
  labs(title="Residual components of the model on date at fist feeding") + ylab("Posterior modes") + xlab("")
```

##### 5\. Variation in propensity to feed

``` r
pr<-list(G=list(
  G1=list(V=1,nu=0.002),
  G2=list(V=diag(3),nu=0.002)), 
  R=list(V=1,fix=1))
iter<-350
m.rf<-MCMCglmm(feedingresp~factor(trial)*cross2,
               random =~fam + idh(cross2):ind, prior = pr, rcov = ~units, family = "categorical",data = m.ll,nitt = 13000*iter,thin = 10*iter, burnin = 3000*iter)

# Reduced model for observed R
m.rf.r<-MCMCglmm(feedingresp~1,
               random =~fam + idh(cross2):ind, prior = pr, rcov = ~units, family = "categorical",data = m.ll,nitt = 13000*iter,thin = 10*iter, burnin = 3000*iter)
```

``` r
# Find estimates on the observed scale for each iteration
par.rf <- data.frame(pl.mean.obs = 1:length(m.rf$Sol[,1]),pl.var.obs = 1:length(m.rf$Sol[,1]),pl.var.a.obs = 1:length(m.rf$Sol[,1]),pl.h2.obs = 1:length(m.rf$Sol[,1]),sb.mean.obs = 1:length(m.rf$Sol[,1]),sb.var.obs = 1:length(m.rf$Sol[,1]),sb.var.a.obs = 1:length(m.rf$Sol[,1]),sb.h2.obs = 1:length(m.rf$Sol[,1]),hyb.mean.obs = 1:length(m.rf$Sol[,1]),hyb.var.obs = 1:length(m.rf$Sol[,1]),hyb.var.a.obs = 1:length(m.rf$Sol[,1]),hyb.h2.obs = 1:length(m.rf$Sol[,1])) # virgin dataframe

for (i in 1:length(m.rf$Sol[,1])) {
  param.pl<-QGparams(mu=mean(m.rf$Sol[i,1]),var.a = mean(m.rf$VCV[i,2]),var.p = sum(m.rf$VCV[i,c(1,2)]),model = "binom1.logit")
  par.rf[i,1:4] <- param.pl
  param.sb<-QGparams(mu=mean(m.rf$Sol[i,1])+mean(m.rf$Sol[i,4]),var.a = mean(m.rf$VCV[i,3]),var.p = sum(m.rf$VCV[i,c(1,3)]),model = "binom1.logit")
  par.rf[i,5:8] <- param.sb
param.hyb<-QGparams(mu=mean(m.rf$Sol[i,1])+mean(m.rf$Sol[i,5]),var.a = mean(m.rf$VCV[i,4]),var.p = sum(m.rf$VCV[i,c(1,4)]),model = "binom1.logit")
  par.rf[i,9:12] <- param.hyb
}
par.rf<-mcmc(data = par.rf,start = 1,end = length(par.rf[,1]),thin = 1) # Convert to a MCMC object
HPDinterval(par.rf) # HPD interval
posterior.mode(par.rf) # posterior estimates

# Dataframe of posterior estimates and 95% CIs
prop.feed.df<-data.frame(gp = row.names(HPDinterval(par.rf)), posterior.mode = posterior.mode(par.rf), lower = HPDinterval(par.rf)[,1], upper =  HPDinterval(par.rf)[,2])
prop.feed.df <- prop.feed.df %>% separate(gp, c("morph","estimate"), extra = "merge")
levels(prop.feed.df$estimate <- c(levels(prop.feed.df$estimate),"Mean","Var. total","Var. individual","h2"))
prop.feed.df$estimate[prop.feed.df$estimate == "mean.obs"] <- "Mean"
prop.feed.df$estimate[prop.feed.df$estimate == "var.obs"] <- "Var. total"
prop.feed.df$estimate[prop.feed.df$estimate == "var.a.obs"] <- "Var. individual"
prop.feed.df$estimate[prop.feed.df$estimate == "R.obs"] <- "h2"
prop.feed.df$estimate.f<-factor(prop.feed.df$estimate, levels = c("Mean","Var. individual","Var. total","h2"))
View(prop.feed.df)
# Plot 
ggplot(prop.feed.df[prop.feed.df$estimate == c("Mean","h2"),], aes(x=morph, y=posterior.mode)) +
  facet_wrap(~estimate.f, ncol = 2, scales = "free_y") +
  geom_hline(yintercept = 0.0, color = "lightgrey") +
  geom_errorbar(aes(ymin=lower, ymax=upper), colour="cyan4",width=.01, size = 1) +
  geom_point(size = 2.5) + theme_bw() + 
  scale_x_discrete(limits = c("pl","sb","hyb"),
                   labels = c("PL","SB","HY")) + theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        strip.background = element_rect(fill = "white", colour = "white"),
        strip.text.x = element_text(size = 10, face = "bold")) + 
  labs(title="Mean and variance components propensity to feed") + ylab("Posterior mode")
```

Repeatability estimates, from Nakagawa et al. (2018)

``` r
prop.VarFc<-numeric(1000)
for(i in 1:1000){
  Var<-var(as.vector(m.rf$Sol[i,] %*% t(m.rf$X[,])))
  prop.VarFc[i]<-Var
}

# PL
Vt <- m.rf$VCV[,"fam"] + m.rf$VCV[,"cross2PL.ind"]
pmean <- plogis(as.numeric(prop.VarFc) - 0.5 * Vt * tanh(as.numeric(prop.VarFc) * 
                                                                     (1 + 2 * exp(-0.5 * Vt))/6))
VarOL <- 1/(pmean * (1 - pmean))
R.pl<- m.rf$VCV[,"cross2PL.ind"]/(m.rf$VCV[,"cross2PL.ind"]+m.rf$VCV[,"fam"]+VarOL)

# SB
Vt <- m.rf$VCV[,"fam"] + m.rf$VCV[,"cross2SB.ind"]
pmean <- plogis(as.numeric(prop.VarFc) - 0.5 * Vt * tanh(as.numeric(prop.VarFc) * 
                                                                     (1 + 2 * exp(-0.5 * Vt))/6))
VarOL <- 1/(pmean * (1 - pmean))
R.sb<- m.rf$VCV[,"cross2SB.ind"]/(m.rf$VCV[,"cross2SB.ind"]+m.rf$VCV[,"fam"]+VarOL)


# Hybrid
Vt <- m.rf$VCV[,"fam"] + m.rf$VCV[,"cross2hybrid.ind"]
pmean <- plogis(as.numeric(prop.VarFc) - 0.5 * Vt * tanh(as.numeric(prop.VarFc) * 
                                                                     (1 + 2 * exp(-0.5 * Vt))/6))
VarOL <- 1/(pmean * (1 - pmean))
R.hybrid<- m.rf$VCV[,"cross2hybrid.ind"]/(m.rf$VCV[,"cross2hybrid.ind"]+m.rf$VCV[,"fam"]+VarOL)

prop.feed.df2<- rbind(
  prop.feed.df[,1:5],
  c("pl", "R", posterior.mode(R.pl), HPDinterval(R.pl)),
  c("sb", "R", posterior.mode(R.sb), HPDinterval(R.sb)),
  c("hyb", "R", posterior.mode(R.hybrid), HPDinterval(R.hybrid))
)
prop.feed.df2$posterior.mode<- as.numeric(prop.feed.df2$posterior.mode)
prop.feed.df2$lower<- as.numeric(prop.feed.df2$lower)
prop.feed.df2$upper<- as.numeric(prop.feed.df2$upper)
ggplot(prop.feed.df2[prop.feed.df2$estimate == "Mean" | prop.feed.df2$estimate == "R",], aes(x=morph, y=posterior.mode)) +
  facet_wrap(~estimate, ncol = 2, scales = "free_y") +
  geom_hline(yintercept = 0.0, color = "lightgrey") +
  geom_errorbar(aes(ymin=lower, ymax=upper), colour="cyan4",width=.01, size = 1) +
  geom_point(size = 2.5) + theme_bw() + 
  scale_x_discrete(limits = c("pl","sb","hyb"),
                   labels = c("PL","SB","HY")) + theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        strip.background = element_rect(fill = "white", colour = "white"),
        strip.text.x = element_text(size = 10, face = "bold")) + 
  labs(title="Mean and variance components propensity to feed") + ylab("Posterior mode")
```

##### 6\. Variation in the latency to start feeding

``` r
p.lat<-list(G=list(
  G1=list(V=1,nu=0.002),
  G2=list(V=diag(3),nu=0.002)), 
  R=list(V=diag(3), nu = 0.002))
itr<-200
latency.m<-MCMCglmm(log(time.response)~factor(trial)*cross2 - 1,
               random =~fam + idh(cross2):ind, prior = p.lat, rcov = ~idh(cross2):units,data = m.ll[complete.cases(m.ll$time.response),],nitt = 13000*itr,thin = 10*itr, burnin = 3000*itr, verbose = T)
```

# Repeatability

``` r
rep.pl<-latency.m$VCV[,"cross2PL.ind"]/(latency.m$VCV[,"cross2PL.ind"]+latency.m$VCV[,"cross2PL.units"])
rep.sb<-latency.m$VCV[,"cross2SB.ind"]/(latency.m$VCV[,"cross2SB.ind"]+latency.m$VCV[,"cross2SB.units"])
rep.hyb<-latency.m$VCV[,"cross2hybrid.ind"]/(latency.m$VCV[,"cross2hybrid.ind"]+latency.m$VCV[,"cross2hybrid.units"])

lat.rep<-data.frame(cross = c("PL","SB","HY"), 
                     estimate = c(posterior.mode(rep.pl),posterior.mode(rep.sb),posterior.mode(rep.hyb)),
                     lower = c(HPDinterval(rep.pl)[,1],HPDinterval(rep.sb)[,1],HPDinterval(rep.hyb)[,1]),
                     upper = c(HPDinterval(rep.pl)[,2],HPDinterval(rep.sb)[,2],HPDinterval(rep.hyb)[,2]))

ggplot(lat.rep, aes(x=cross, y=estimate)) +
  geom_hline(yintercept = 0.0, color = "lightgrey") +
  geom_errorbar(aes(ymin=lower, ymax=upper), colour="cyan4",width=.01, size = 0.8) +
  geom_point(size = 2) + theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.text.x = element_text(size = 10),
        axis.title.y = element_text(size = 10)) + 
  labs(title="Repeatability latency to start first feeding") + ylab("Posterior mode")
```

##### 7\. Variation in the number of feeding attempts

``` r
hist(m.ll$nb.bouts,breaks = 10)
p.bout<-list(G=list(
  G1=list(V=1,nu=0.002),
  G2=list(V=diag(3),nu=0.002)), 
  R=list(V=diag(3), nu = 0.002))
itr<-4500
bouts.m<-MCMCglmm(nb.bouts~factor(trial)*cross2,
               random =~fam + idh(cross2):ind, prior = p.bout, rcov = ~idh(cross2):units, family = "poisson",data = m.ll,nitt = 13000*itr,thin = 10*itr, burnin = 3000*itr)
```

``` r
# Find estimates on the observed scale for each iteration
par.bouts <- data.frame(pl.mean.obs = 1:length(bouts.m$Sol[,1]),pl.var.obs = 1:length(bouts.m$Sol[,1]),pl.var.a.obs = 1:length(bouts.m$Sol[,1]),pl.h2.obs = 1:length(bouts.m$Sol[,1]),sb.mean.obs = 1:length(bouts.m$Sol[,1]),sb.var.obs = 1:length(bouts.m$Sol[,1]),sb.var.a.obs = 1:length(bouts.m$Sol[,1]),sb.h2.obs = 1:length(bouts.m$Sol[,1]),hyb.mean.obs = 1:length(bouts.m$Sol[,1]),hyb.var.obs = 1:length(bouts.m$Sol[,1]),hyb.var.a.obs = 1:length(bouts.m$Sol[,1]),hyb.h2.obs = 1:length(bouts.m$Sol[,1])) # Virgin dataframe
for (i in 1:length(bouts.m$Sol[,1])) {
  parbpl <- QGparams(mu=bouts.m$Sol[i,1], var.a = bouts.m$VCV[i,2], var.p = sum(bouts.m$VCV[i,c(1,2,5)]), model = "Poisson.log")
  par.bouts [i,1:4] <- parbpl
  parbsb <- QGparams(mu=bouts.m$Sol[i,1]+bouts.m$Sol[i,4], var.a = bouts.m$VCV[i,3], var.p = sum(bouts.m$VCV[i,c(1,3,6)]), model = "Poisson.log")
  par.bouts[i,5:8] <- parbsb
  parbhyb<-QGparams(mu=bouts.m$Sol[i,1]+bouts.m$Sol[i,5],var.a = bouts.m$VCV[i,4],var.p = sum(bouts.m$VCV[i,c(1,4,7)]),model = "Poisson.log")
  par.bouts[i,9:12] <- parbhyb
}
par.bouts<-mcmc(data = par.bouts,start = 1,end = length(par.bouts[,1]),thin = 1) # Convert to a MCMC object
HPDinterval(par.bouts) # HPD interval
posterior.mode(par.bouts) # posterior estimates

# Dataframe of posterior estimates and 95% CIs
bouts<-data.frame(gp = row.names(HPDinterval(par.bouts)), posterior.mode = posterior.mode(par.bouts), mean.estimate = colMeans(par.bouts), lower = HPDinterval(par.bouts)[,1], upper =  HPDinterval(par.bouts)[,2])
bouts <- bouts %>% separate(gp, c("morph","estimate"), extra = "merge")

# Distribution h2 estimates
par(mfrow = c(2,2))
hist(par.bouts[,4],main = "h2 PL")
hist(par.bouts[,8],main = "h2 SL")
hist(par.bouts[,12],main = "h2 hyb")
```

``` r
tem<-theme(title = element_text(size=9),legend.position = "none")
means.bouts <- ggplot(bouts[bouts$estimate == "mean.obs",], aes(x=morph, y=posterior.mode, colour= morph,group = morph)) + 
  geom_point(position=pd,size = 2) + 
  geom_errorbar(aes(ymin=lower, ymax=upper), colour="black",width=.1, position=pd) + 
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  labs(title="n bouts: posterior mode of 
       phenotypic mean +/- 95% CI") + ylab("n bouts") + xlab("cross2") + tem
var.bouts<-ggplot(bouts[bouts$estimate == "var.obs",], aes(x=morph, y=posterior.mode, colour= morph,group = morph)) + 
  geom_point(position=pd,size = 2) + 
  geom_errorbar(aes(ymin=lower, ymax=upper), colour="black",width=.1, position=pd) + 
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  labs(title="n bouts: posterior mode of
       variance estimates +/- 95% CI") + ylab("n bouts") + xlab("cross2") + tem
h2.bouts<-ggplot(bouts[bouts$estimate == "h2.obs",], aes(x=morph, y=posterior.mode, colour= morph,group = morph)) + 
  geom_point(position=pd,size = 2) + 
  geom_errorbar(aes(ymin=lower, ymax=upper), colour="black",width=.1, position=pd) + 
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  labs(title="n bouts: posterior mode of
       heritability estimates +/- 95% CI") + ylab("n bouts") + xlab("cross2") + theme(title = element_text(size=9))
grid.arrange(means.bouts,var.bouts,h2.bouts,nrow = 2, ncol = 2)
#
```

##### Estimates on the Observed scale, following Nakagawa et al. 2018

``` r
p.bout<-list(G=list(
  G1=list(V=1,nu=0.002),
  G2=list(V=diag(3),nu=0.002)), 
  R=list(V=diag(3), nu = 0.002))
itr<-4500
bouts.m.r<-MCMCglmm(nb.bouts~1,
               random =~fam + idh(cross2):ind, prior = p.bout, rcov = ~idh(cross2):units, family = "poisson",data = m.ll,nitt = 13000*itr,thin = 10*itr, burnin = 3000*itr)

nb.VarFc<-numeric(1000)
for(i in 1:1000){
  Var<-var(as.vector(bouts.m$Sol[i,] %*% t(bouts.m$X[,])))
  nb.VarFc[i]<-Var
}

o.nb.pl <- as.numeric(bouts.m.r$VCV[,"cross2PL.units"])
lambda.pl <- (bouts.m.r$Sol[,"(Intercept)"]) + 0.5 * (bouts.m.r$VCV[,"fam"]+bouts.m.r$VCV[,"cross2PL.ind"])
tri.nb.pl <- trigamma(lambda.pl/o.nb.pl)
```

##### 8\. Variation in the tendency to forage on the bottom / column / surface

``` r
itr<-800
bott.m<-MCMCglmm(bottom~nb.bouts + cross2,
               random =~fam + idh(cross2):ind, prior = p.bout, rcov = ~idh(cross2):units, family = "poisson",data = m.ll[complete.cases(m.ll$bottom),],nitt = 13000*itr,thin = 10*itr, burnin = 3000*itr)

col.m<-MCMCglmm(collumn~nb.bouts + cross2,
               random =~fam + idh(cross2):ind, prior = p.bout, rcov = ~idh(cross2):units, family = "poisson",data = m.ll[complete.cases(m.ll$collumn),],nitt = 13000*itr,thin = 10*itr, burnin = 3000*itr)

surf.m<-MCMCglmm(surface~nb.bouts + cross2,
               random =~fam + idh(cross2):ind, prior = p.bout, rcov = ~idh(cross2):units, family = "poisson",data = m.ll[complete.cases(m.ll$surface),],nitt = 13000*itr,thin = 10*itr, burnin = 3000*itr)
```

Bottom

``` r
# Find estimates on the observed scale for each iteration
par.bott <- data.frame(pl.mean.obs = 1:length(bott.m$Sol[,1]),pl.var.obs = 1:length(bott.m$Sol[,1]),pl.var.a.obs = 1:length(bott.m$Sol[,1]),pl.h2.obs = 1:length(bott.m$Sol[,1]),sb.mean.obs = 1:length(bott.m$Sol[,1]),sb.var.obs = 1:length(bott.m$Sol[,1]),sb.var.a.obs = 1:length(bott.m$Sol[,1]),sb.h2.obs = 1:length(bott.m$Sol[,1]),hyb.mean.obs = 1:length(bott.m$Sol[,1]),hyb.var.obs = 1:length(bott.m$Sol[,1]),hyb.var.a.obs = 1:length(bott.m$Sol[,1]),hyb.h2.obs = 1:length(bott.m$Sol[,1])) # virgin dataframe
for (i in 1:length(bott.m$Sol[,1])) {
  parbpl <- QGparams(mu=bott.m$Sol[i,1], var.a = bott.m$VCV[i,2], var.p = sum(bott.m$VCV[i,c(1,2,5)]), model = "Poisson.log")
  par.bott [i,1:4] <- parbpl
  parbsb <- QGparams(mu=bott.m$Sol[i,1]+bott.m$Sol[i,3], var.a = bott.m$VCV[i,3], var.p = sum(bott.m$VCV[i,c(1,3,6)]), model = "Poisson.log")
  par.bott[i,5:8] <- parbsb
  parbhyb<-QGparams(mu=bott.m$Sol[i,1]+bott.m$Sol[i,4],var.a = bott.m$VCV[i,4],var.p = sum(bott.m$VCV[i,c(1,4,7)]),model = "Poisson.log")
  par.bott[i,9:12] <- parbhyb
}
par.bott<-mcmc(data = par.bott,start = 1,end = length(par.bott[,1]),thin = 1) # Convert to a MCMC object
HPDinterval(par.bott) # HPD interval
posterior.mode(par.bott) # posterior estimates

bott<-data.frame(gp = row.names(HPDinterval(par.bott)), posterior.mode = posterior.mode(par.bott), mean.estimate = colMeans(par.bott), lower = HPDinterval(par.bott)[,1], upper =  HPDinterval(par.bott)[,2])
bott <- bott %>% separate(gp, c("morph","estimate"), extra = "merge")

ggplot(bott, aes(x=morph, y=posterior.mode)) +
  facet_grid(cols = vars(estimate)) +
  geom_hline(yintercept = 0.0, color = "lightgrey") +
  geom_errorbar(aes(ymin=lower, ymax=upper), colour="cyan4",width=.01, size = 0.8) +
  geom_point(size = 2) + theme_bw() 
  scale_x_discrete(limits = c("PL","SB","HB"),labels = c("PL","SB","HY")) + theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.text.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        strip.background = element_rect(fill = "white", colour = "white")) + 
  labs(title="Personality trait covariance") + ylab("Covariance") + xlab("")
```

Collumn

``` r
# Find estimates on the observed scale for each iteration
par.col <- data.frame(pl.mean.obs = 1:length(col.m$Sol[,1]),pl.var.obs = 1:length(col.m$Sol[,1]),pl.var.a.obs = 1:length(col.m$Sol[,1]),pl.h2.obs = 1:length(col.m$Sol[,1]),sb.mean.obs = 1:length(col.m$Sol[,1]),sb.var.obs = 1:length(col.m$Sol[,1]),sb.var.a.obs = 1:length(col.m$Sol[,1]),sb.h2.obs = 1:length(col.m$Sol[,1]),hyb.mean.obs = 1:length(col.m$Sol[,1]),hyb.var.obs = 1:length(col.m$Sol[,1]),hyb.var.a.obs = 1:length(col.m$Sol[,1]),hyb.h2.obs = 1:length(col.m$Sol[,1])) # virgin dataframe
for (i in 1:length(col.m$Sol[,1])) {
  parbpl <- QGparams(mu=col.m$Sol[i,1], var.a = col.m$VCV[i,2], var.p = sum(col.m$VCV[i,c(1,2,5)]), model = "Poisson.log")
  par.col [i,1:4] <- parbpl
  parbsb <- QGparams(mu=col.m$Sol[i,1]+col.m$Sol[i,3], var.a = col.m$VCV[i,3], var.p = sum(col.m$VCV[i,c(1,3,6)]), model = "Poisson.log")
  par.col[i,5:8] <- parbsb
  parbhyb<-QGparams(mu=col.m$Sol[i,1]+col.m$Sol[i,4],var.a = col.m$VCV[i,4],var.p = sum(col.m$VCV[i,c(1,4,7)]),model = "Poisson.log")
  par.col[i,9:12] <- parbhyb
}
par.col<-mcmc(data = par.col,start = 1,end = length(par.col[,1]),thin = 1) # Convert to a MCMC object
HPDinterval(par.col) # HPD interval
posterior.mode(par.col) # posterior estimates

col<-data.frame(gp = row.names(HPDinterval(par.col)), posterior.mode = posterior.mode(par.col), mean.estimate = colMeans(par.col), lower = HPDinterval(par.col)[,1], upper =  HPDinterval(par.col)[,2])
col <- col %>% separate(gp, c("morph","estimate"), extra = "merge")

ggplot(col, aes(x=morph, y=posterior.mode)) +
  facet_grid(cols = vars(estimate)) +
  geom_hline(yintercept = 0.0, color = "lightgrey") +
  geom_errorbar(aes(ymin=lower, ymax=upper), colour="cyan4",width=.01, size = 0.8) +
  geom_point(size = 2) + theme_bw() 
  scale_x_discrete(limits = c("PL","SB","HB"),labels = c("PL","SB","HY")) + theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.text.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        strip.background = element_rect(fill = "white", colour = "white")) + 
  labs(title="Personality trait covariance") + ylab("Covariance") + xlab("")
```

Surface

``` r
# Find estimates on the observed scale for each iteration
par.surf <- data.frame(pl.mean.obs = 1:length(surf.m$Sol[,1]),pl.var.obs = 1:length(surf.m$Sol[,1]),pl.var.a.obs = 1:length(surf.m$Sol[,1]),pl.h2.obs = 1:length(surf.m$Sol[,1]),sb.mean.obs = 1:length(surf.m$Sol[,1]),sb.var.obs = 1:length(surf.m$Sol[,1]),sb.var.a.obs = 1:length(surf.m$Sol[,1]),sb.h2.obs = 1:length(surf.m$Sol[,1]),hyb.mean.obs = 1:length(surf.m$Sol[,1]),hyb.var.obs = 1:length(surf.m$Sol[,1]),hyb.var.a.obs = 1:length(surf.m$Sol[,1]),hyb.h2.obs = 1:length(surf.m$Sol[,1])) # virgin dataframe
for (i in 1:length(surf.m$Sol[,1])) {
  parbpl <- QGparams(mu=surf.m$Sol[i,1], var.a = surf.m$VCV[i,2], var.p = sum(surf.m$VCV[i,c(1,2,5)]), model = "Poisson.log")
  par.surf [i,1:4] <- parbpl
  parbsb <- QGparams(mu=surf.m$Sol[i,1]+surf.m$Sol[i,3], var.a = surf.m$VCV[i,3], var.p = sum(surf.m$VCV[i,c(1,3,6)]), model = "Poisson.log")
  par.surf[i,5:8] <- parbsb
  parbhyb<-QGparams(mu=surf.m$Sol[i,1]+surf.m$Sol[i,4],var.a = surf.m$VCV[i,4],var.p = sum(surf.m$VCV[i,c(1,4,7)]),model = "Poisson.log")
  par.surf[i,9:12] <- parbhyb
}
par.surf<-mcmc(data = par.surf,start = 1,end = length(par.surf[,1]),thin = 1) # Convert to a MCMC object
HPDinterval(par.surf) # HPD interval
posterior.mode(par.surf) # posterior estimates

surf<-data.frame(gp = row.names(HPDinterval(par.surf)), posterior.mode = posterior.mode(par.surf), mean.estimate = colMeans(par.surf), lower = HPDinterval(par.surf)[,1], upper =  HPDinterval(par.surf)[,2])
surf <- surf %>% separate(gp, c("morph","estimate"), extra = "merge")

ggplot(surf[surf$estimate == "mean.obs",],aes(x=morph, y=posterior.mode)) +
  facet_grid(cols = vars(estimate),scales = "free") +
  geom_hline(yintercept = 0.0, color = "lightgrey") +
  geom_errorbar(aes(ymin=lower, ymax=upper), colour="cyan4",width=.01, size = 0.8) +
  geom_point(size = 2) + theme_bw() 
  scale_x_discrete(limits = c("PL","SB","HB"),labels = c("PL","SB","HY")) + theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.text.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        strip.background = element_rect(fill = "white", surfour = "white")) + 
  labs(title="Trait covariance") + ylab("Covariance") + xlab("")
```
