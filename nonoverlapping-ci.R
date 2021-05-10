
pair.cov<- matrix(,ncol(m.big.t.pl2$VCV),4) # correct length to only .ind
colnames(pair.cov) <- c("Trait","PLxPL","SBxSB","F1")
Vnonoverlap <- function(CrI) {
for (i in 1:ncol(m.big.t.pl2$VCV)) {
   if (HPDinterval(m.big.t.pl2$VCV[,i], prob = CrI)[,1] > HPDinterval(m.big.t.sb2$VCV[,i], prob = CrI)[,2] ||
       HPDinterval(m.big.t.sb2$VCV[,i], prob = CrI)[,1] > HPDinterval(m.big.t.pl2$VCV[,i], prob = CrI)[,2] ||
       
       HPDinterval(m.big.t.pl2$VCV[,i], prob = CrI)[,1] > HPDinterval(m.big.t.hyb2$VCV[,i], prob = CrI)[,2] ||
       HPDinterval(m.big.t.hyb2$VCV[,i], prob = CrI)[,1] > HPDinterval(m.big.t.pl2$VCV[,i], prob = CrI)[,2] ||
       
       HPDinterval(m.big.t.sb2$VCV[,i], prob = CrI)[,1] > HPDinterval(m.big.t.hyb2$VCV[,i], prob = CrI)[,2] ||
       HPDinterval(m.big.t.hyb2$VCV[,i], prob = CrI)[,1] > HPDinterval(m.big.t.sb2$VCV[,i], prob = CrI)[,2] )
   { pair.cov[i,1] <- colnames(m.big.t.pl2$VCV)[i]
   pair.cov[i,"PLxPL"] <- paste(round(posterior.mode(m.big.t.pl2$VCV[,i]),2)," [", round(HPDinterval(m.big.t.pl2$VCV[,i])[,1],2),
                                "; ",round(HPDinterval(m.big.t.pl2$VCV[,i])[,2],2),"]",sep = "")
   pair.cov[i,"SBxSB"] <- paste(round(posterior.mode(m.big.t.sb2$VCV[,i]),2)," [", round(HPDinterval(m.big.t.sb2$VCV[,i])[,1],2),
                                "; ",round(HPDinterval(m.big.t.sb2$VCV[,i])[,2],2),"]",sep = "")
   pair.cov[i,"F1"] <- paste(round(posterior.mode(m.big.t.hyb2$VCV[,i]),2)," [", round(HPDinterval(m.big.t.hyb2$VCV[,i])[,1],2),
                                "; ",round(HPDinterval(m.big.t.hyb2$VCV[,i])[,2],2),"]",sep = "")
   } 
   else {
     pair.cov[i,1] <- NA
   }
}
  return(pair.cov[complete.cases(pair.cov[,1]),])
}
Vnonoverlap(.90)
rm(Vn)


########## For correlations 

corr.PL<- posterior.cor(m.big.t.pl2$VCV[,1:n^2])
corr.PL[,c(1,9,17,25,33,41,49)] <- m.big.t.pl2$VCV[,c(1,9,17,25,33,41,49)]
colnames(corr.PL)<- colnames(m.big.t.pl2$VCV[,1:n^2])

corr.SB<- posterior.cor(m.big.t.sb2$VCV[,1:n^2])
corr.SB[,c(1,9,17,25,33,41,49)] <- m.big.t.sb2$VCV[,c(1,9,17,25,33,41,49)]
colnames(corr.SB)<- colnames(m.big.t.sb2$VCV[,1:n^2])

corr.hyb<- posterior.cor(m.big.t.hyb2$VCV[,1:n^2])
corr.hyb[,c(1,9,17,25,33,41,49)] <- m.big.t.hyb2$VCV[,c(1,9,17,25,33,41,49)]
colnames(corr.hyb)<- colnames(m.big.t.hyb2$VCV[,1:n^2])

pair.corr<- matrix(,ncol(corr.PL),4) # correct length to only .ind
colnames(pair.corr) <- c("Trait","PLxPL","SBxSB","F1")
Vnonoverlap.corr <- function(CrI) {
  for (i in 1:ncol(corr.PL)) {
    if (HPDinterval(corr.PL[,i], prob = CrI)[,1] > HPDinterval(corr.SB[,i], prob = CrI)[,2] ||
        HPDinterval(corr.SB[,i], prob = CrI)[,1] > HPDinterval(corr.PL[,i], prob = CrI)[,2] ||
        
        HPDinterval(corr.PL[,i], prob = CrI)[,1] > HPDinterval(corr.hyb[,i], prob = CrI)[,2] ||
        HPDinterval(corr.hyb[,i], prob = CrI)[,1] > HPDinterval(corr.PL[,i], prob = CrI)[,2] ||
        
        HPDinterval(corr.SB[,i], prob = CrI)[,1] > HPDinterval(corr.hyb[,i], prob = CrI)[,2] ||
        HPDinterval(corr.hyb[,i], prob = CrI)[,1] > HPDinterval(corr.SB[,i], prob = CrI)[,2] )
    { pair.corr[i,1] <- colnames(corr.PL)[i]
    pair.corr[i,"PLxPL"] <- paste(round(posterior.mode(corr.PL[,i]),2)," [", round(HPDinterval(corr.PL[,i], prob = CrI)[,1],2),
                                 "; ",round(HPDinterval(corr.PL[,i], prob = CrI)[,2],2),"]",sep = "")
    pair.corr[i,"SBxSB"] <- paste(round(posterior.mode(corr.SB[,i]),2)," [", round(HPDinterval(corr.SB[,i], prob = CrI)[,1],2),
                                 "; ",round(HPDinterval(corr.SB[,i], prob = CrI)[,2],2),"]",sep = "")
    pair.corr[i,"F1"] <- paste(round(posterior.mode(corr.hyb[,i]),2)," [", round(HPDinterval(corr.hyb[,i], prob = CrI)[,1],2),
                              "; ",round(HPDinterval(corr.hyb[,i], prob = CrI)[,2],2),"]",sep = "")
    } 
    else {
      pair.corr[i,1] <- NA
    }
  }
  return(pair.corr[complete.cases(pair.corr[,1]),])
}

nno95 <- Vnonoverlap.corr(.95)
nno95 <- as.data.frame(nno95[!duplicated(nno95[,2]) & !duplicated(nno95[,3]),])

nno90 <- Vnonoverlap.corr(.90)
nno90 <- as.data.frame(nno90[!duplicated(nno90[,2]) & !duplicated(nno90[,3]),])

nno80 <-Vnonoverlap.corr(.80)
nno80 <- as.data.frame(nno80[!duplicated(nno80[,2]) & !duplicated(nno80[,3]),])

write.csv(rbind(data.frame(target.prob = "95% CrI",t(nno95)),
                data.frame(target.prob = rep("90% CrI", nrow(nno90)),nno90),
                data.frame(target.prob = rep("90% CrI", nrow(nno80)),nno80)),"Desktop/nonoverlapping-corr.txt")




####### Plot correlations

data.big.t2<-data.frame(
  cross = c(rep("PL",49),rep("SB",49),rep("F1",49)),
  post.mod = c(posterior.mode(corr.PL),posterior.mode(corr.SB),posterior.mode(corr.hyb)),
  upper = c(HPDinterval(corr.PL)[,1],HPDinterval(corr.SB)[,1],HPDinterval(corr.hyb)[,1]),
  lower = c(HPDinterval(corr.PL)[,2],HPDinterval(corr.SB)[,2],HPDinterval(corr.hyb)[,2]),
  trait = rep(names(posterior.mode(m.big.t.pl2$VCV)[1:49]),3))

#data.big.t1<- data.big.t[data.big.t$cross != "SB",] # dataset without SB
#data.big.t1<-droplevels(data.big.t1)

ggplot(data.big.t2, aes(x=cross, y=post.mod, color = cross)) +
  facet_wrap(vars(trait),nrow = 10,ncol = 10, scales = "free") +
  geom_hline(yintercept = 0.0, color = "grey") +
  geom_errorbar(aes(ymin=lower, ymax=upper), colour="grey",width=.01, size = 0.8) +
  geom_point(size = 1.5) + theme_bw() +
  scale_color_manual(values = c("#3333FF","orange","#009933")) +
  scale_x_discrete(limits = c("PL","SB","F1")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_blank(), # remove here to see labels
        axis.text.x =element_blank(),
        axis.title.y = element_text(size = 5),
        axis.text.y = element_text(size = 7),
        strip.background = element_blank()) +
  labs(title="") + ylab("") + xlab("")
