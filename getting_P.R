library(tidyr)
library(MCMCglmm)
library(ggplot2)
library(plyr)
library(gridExtra)
library(tidybayes)

### Set priors
p.big.t<-list(G=list(
  G1=list(V=diag(7),nu=7)),
  R=list(V=diag(7),nu=7))
itr<-500

# Run MCMC GLMM for PL charr
m.big.t.pl2<-MCMCglmm(cbind(D1/mean(na.omit(D1)),D3/mean(na.omit(D3)),
                            growth.d1.d2/mean(na.omit(growth.d1.d2)),
                            growth.d3.d4/mean(na.omit(growth.d3.d4)),
                            r.yolk.D1/mean(na.omit(r.yolk.D1)),y.conversion/mean(na.omit(y.conversion)), 
                            log(time.response)/mean(complete.cases(log(time.response))))~ 
                    trait - 1 + fam,
                  random =~ us(trait):ind, prior = p.big.t, 
                  rcov = ~us(trait):units, data = df.big.test2[df.big.test2$cross2 == "PL" & complete.cases(df.big.test2$fam),], 
                  family = c("gaussian","gaussian","gaussian","gaussian","gaussian",
                             "gaussian","gaussian"),
                  nitt = 13000*itr,thin = 10*itr, burnin = 3000*itr,verbose = T)

# Run MCMC GLMM for SB charr
m.big.t.sb2<-MCMCglmm(cbind(D1/mean(na.omit(D1)),D3/mean(na.omit(D3)),
                            growth.d1.d2/mean(na.omit(growth.d1.d2)),
                            growth.d3.d4/mean(na.omit(growth.d3.d4)),
                            r.yolk.D1/mean(na.omit(r.yolk.D1)),y.conversion/mean(na.omit(y.conversion)), 
                            log(time.response)/mean(complete.cases(log(time.response))))~
                       trait - 1 + fam,
                     random =~ us(trait):ind, prior = p.big.t, 
                     rcov = ~us(trait):units, data = df.big.test2[df.big.test2$cross2 == "SB" & complete.cases(df.big.test2$fam),], 
                     family = c("gaussian","gaussian","gaussian","gaussian","gaussian",
                                 "gaussian","gaussian"),
                     nitt = 13000*itr,thin = 10*itr, burnin = 3000*itr,verbose = T)

# Run MCMC GLMM for hybtrids
m.big.t.hyb2<-MCMCglmm(cbind(D1/mean(na.omit(D1)),D3/mean(na.omit(D3)),
                             growth.d1.d2/mean(na.omit(growth.d1.d2)),
                             growth.d3.d4/mean(na.omit(growth.d3.d4)),
                             r.yolk.D1/mean(na.omit(r.yolk.D1)),y.conversion/mean(na.omit(y.conversion)),
                             log(time.response)/mean(complete.cases(log(time.response))))~
                        trait - 1 + fam,
                      random =~ us(trait):ind, prior = p.big.t, 
                      rcov = ~us(trait):units, data = df.big.test2[df.big.test2$cross2 == "hybrid" &
                                                                     complete.cases(df.big.test2$fam),], 
                      family = c("gaussian","gaussian","gaussian","gaussian","gaussian",
                                 "gaussian","gaussian"),
                      nitt = 13000*itr,thin = 10*itr, burnin = 3000*itr,verbose = T)


#### Second set of model

data.big.t2<-data.frame(
  cross = c(rep("PL",49),rep("SB",49),rep("F1",49)),
  post.mod = c(posterior.mode(m.big.t.pl2$VCV)[1:49],posterior.mode(m.big.t.sb2$VCV)[1:49],posterior.mode(m.big.t.hyb2$VCV)[1:49]),
  upper = c(HPDinterval(m.big.t.pl2$VCV)[1:49,1],HPDinterval(m.big.t.sb2$VCV)[1:49,1],HPDinterval(m.big.t.hyb2$VCV)[1:49,1]),
  lower = c(HPDinterval(m.big.t.pl2$VCV)[1:49,2],HPDinterval(m.big.t.sb2$VCV)[1:49,2],HPDinterval(m.big.t.hyb2$VCV)[1:49,2]),
  trait = rep(names(posterior.mode(m.big.t.pl2$VCV)[1:49]),3))

# Cross differences in the posteriors estimates of the elements of P
ggplot(data.big.t2, aes(x=cross, y=post.mod, color = cross)) +
  facet_wrap(vars(trait),nrow = 7,ncol = 7, scales = "free") +
  geom_hline(yintercept = 0.0, color = "grey") +
  geom_errorbar(aes(ymin=lower, ymax=upper), colour="grey",width=.01, size = 0.8) +
  geom_point(size = 1.5) + theme_bw() +
  scale_color_manual(values = c("#3333FF","black","#009933")) +
  scale_x_discrete(limits = c("PL","SB","F1")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #strip.text = element_blank(), # remove here to see labels
        axis.text.x =element_blank(),
        axis.title.y = element_text(size = 5),
        axis.text.y = element_text(size = 5),
        strip.background = element_blank()) +
  labs(title="") + ylab("") + xlab("")



##### Randomized datasets

vperm<-sample(levels(df.big.test2$ind))
plperm<-vperm[1:length(vperm[grep("^PL.._",vperm)])]
sbperm<-vperm[(1+length(vperm[grep("^PL.._",vperm)])):(length(vperm[grep("^PL.._",vperm)])+length(vperm[grep("^SB.._",vperm)]))]
hybperm<-vperm[(1+length(vperm[grep("^PL.._",vperm)])+length(vperm[grep("^SB.._",vperm)])):length(vperm)]
r.pl<-df.big.test2[df.big.test2$ind == vperm[1:5],]

m.big.t.plperm<-MCMCglmm(cbind(D1/mean(na.omit(D1)),D3/mean(na.omit(D3)),
                            growth.d1.d2/mean(na.omit(growth.d1.d2)),
                            growth.d3.d4/mean(na.omit(growth.d3.d4)),
                            r.yolk.D1/mean(na.omit(r.yolk.D1)),y.conversion/mean(na.omit(y.conversion)), 
                            log(time.response)/mean(complete.cases(log(time.response))))~ 
                        trait - 1 + fam,
                      random =~ us(trait):ind, prior = p.big.t, 
                      rcov = ~us(trait):units, data = df.big.test2[df.big.test2$ind %in% plperm & complete.cases(df.big.test2$fam),], 
                      family = c("gaussian","gaussian","gaussian","gaussian","gaussian",
                                 "gaussian","gaussian"),
                      nitt = 13000*itr,thin = 10*itr, burnin = 3000*itr,verbose = T)

m.big.t.sbperm<-MCMCglmm(cbind(D1/mean(na.omit(D1)),D3/mean(na.omit(D3)),
                               growth.d1.d2/mean(na.omit(growth.d1.d2)),
                               growth.d3.d4/mean(na.omit(growth.d3.d4)),
                               r.yolk.D1/mean(na.omit(r.yolk.D1)),y.conversion/mean(na.omit(y.conversion)), 
                               log(time.response)/mean(complete.cases(log(time.response))))~ 
                           trait - 1 + fam,
                         random =~ us(trait):ind, prior = p.big.t, 
                         rcov = ~us(trait):units, data = df.big.test2[df.big.test2$ind %in% sbperm & complete.cases(df.big.test2$fam),], 
                         family = c("gaussian","gaussian","gaussian","gaussian","gaussian",
                                    "gaussian","gaussian"),
                         nitt = 13000*itr,thin = 10*itr, burnin = 3000*itr,verbose = T)

m.big.t.hybperm<-MCMCglmm(cbind(D1/mean(na.omit(D1)),D3/mean(na.omit(D3)),
                               growth.d1.d2/mean(na.omit(growth.d1.d2)),
                               growth.d3.d4/mean(na.omit(growth.d3.d4)),
                               r.yolk.D1/mean(na.omit(r.yolk.D1)),y.conversion/mean(na.omit(y.conversion)), 
                               log(time.response)/mean(complete.cases(log(time.response))))~ 
                           trait - 1 + fam,
                         random =~ us(trait):ind, prior = p.big.t, 
                         rcov = ~us(trait):units, data = df.big.test2[df.big.test2$ind %in% hybperm & complete.cases(df.big.test2$fam),], 
                         family = c("gaussian","gaussian","gaussian","gaussian","gaussian",
                                    "gaussian","gaussian"),
                         nitt = 13000*itr,thin = 10*itr, burnin = 3000*itr,verbose = T)



##########################################################################################
########################## Extract P matrices (code adapted from Aguirre et al. 2014) ####
##########################################################################################

#number of MCMC samples
MCMCsamp <- 1000 
#number of traits 
n <- 7
#number of matrices to compare
m <- 3
#number of random effects specified in the model. 
r <- 2

traitnames <- c("log_size_age_1","log_size_age_3","growth_D1_to_D2","growth_D3_to_D4","yolk_sac_age_1",
                "yolk.conversion.rate","feeding_latency")
Gnames <- c("PL","SB","F1")
#matrix labels 

############ Make an empty array and fill with MCMC output
MCMCarray <- array( ,c(MCMCsamp,(n^2)*r,m)) 
#empty array 
MCMCarray[,,1] <- as.matrix(m.big.t.pl2$VCV)
#G1 stored as the 1st element of dim[3] 
MCMCarray[,,2] <- as.matrix(m.big.t.sb2$VCV)
#G2 stored as the 2nd element of dim[3]
MCMCarray[,,3] <- as.matrix((m.big.t.hyb2$VCV))


Garray <- array(,c(n,n,m,MCMCsamp))
dimnames(Garray) <- list(traitnames,traitnames,Gnames)
dimnames(Parray) <- list(traitnames,traitnames,Gnames)
for (i in 1:m){
  for (j in 1:MCMCsamp){
    G <- matrix(MCMCarray[j,1:(n^2),i],ncol= n)
    # CE <- matrix(MCMCarray[j,((n^2)+1):((n^2)*2),i],ncol= n)
    R <- matrix(MCMCarray[j,((n^2)+1):((n^2)*2),i],ncol= n)
    Garray[,,i,j] <- G
  }
}

###### P of randomized populations 

############ Make an empty array and fill with MCMC output
MCMCarray.r2 <- array( ,c(MCMCsamp,(n^2)*r,m)) 
#empty array 
MCMCarray.r2[,,1] <- as.matrix(m.big.t.plperm$VCV)
#G1 stored as the 1st element of dim[3] 
MCMCarray.r2[,,2] <- as.matrix(m.big.t.sbperm$VCV)
#G2 stored as the 2nd element of dim[3]
MCMCarray.r2[,,3] <- as.matrix((m.big.t.hybperm$VCV))


Garray.r2 <- array(,c(n,n,m,MCMCsamp))
dimnames(Garray.r2) <- list(traitnames,traitnames,Gnames)
for (i in 1:m){
  for (j in 1:MCMCsamp){
    G <- matrix(MCMCarray.r2[j,1:(n^2),i],ncol= n)
    # CE <- matrix(MCMCarray[j,((n^2)+1):((n^2)*2),i],ncol= n)
    R <- matrix(MCMCarray[j,((n^2)+1):((n^2)*2),i],ncol= n)
    Garray.r2[,,i,j] <- G
  }
}
