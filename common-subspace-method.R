####### P matrices comparisions using Krzanowski's common subspaces

#START
kr.subspace <- function(Gs, vec){
  if (dim(Gs)[[1]] != dim(Gs)[[2]]){
    stop("G array must be of order n x n x m x MCMCsamp")
  }
  if (is.na(dim(Gs)[4])) {
    stop("There are no MCMCsamples")
  }
  n <- dim(Gs)[[1]]
  m <- dim(Gs)[[3]]
  MCMCsamp <- dim(Gs)[[4]] 
  if(length(vec) != m){stop("vec must have length = m")}
  h <- function (g, v){
    AA <- array(, c(n, n, m))  
    for (k in 1:m){
      g.vec <- eigen(g[,,k])$vectors[,1:(v[k])] 
      AA[,,k] <- g.vec %*% t(g.vec)
    }
    H <- apply(AA, 1:2, sum)
    list(H = H, AA = AA)
  }
  #internal function to calculate AA and H
  MCMC.H <- array(, c(n, n, MCMCsamp))
  dimnames(MCMC.H) <- list(dimnames(Gs)[[1]], dimnames(Gs)[[1]], dimnames(Gs)[[4]])      
  MCMC.AA <- array(, c(n, n, m, MCMCsamp))
  dimnames(MCMC.AA) <- list(dimnames(Gs)[[1]], dimnames(Gs)[[1]], dimnames(Gs)[[3]], dimnames(Gs)[[4]])
  for (i in 1:MCMCsamp){
    kr <- h(Gs[,,,i], v = vec)
    MCMC.H[,,i] <- kr$H
    MCMC.AA[,,,i] <- kr$AA
  }	
  #calculate AA and H for the ith MCMC sample of the G array		
  avH <- apply(MCMC.H, 1:2, mean)
  rownames(avH) <- dimnames(Gs)[[1]]
  colnames(avH) <- dimnames(Gs)[[1]]
  #calculate the posterior mean H
  avAA <- apply(MCMC.AA, 1:3, mean)
  dimnames(avAA) <- list(dimnames(Gs)[[1]], dimnames(Gs)[[1]], dimnames(Gs)[[3]])
  #calculate the posterior mean AA
  avH.vec <- eigen(avH)$vectors
  #eigenanalysis of posterior mean H	
  proj<- function(a, b) t(b) %*% a %*% b
  #internal function to do projection
  avH.theta <- matrix(, n, m)
  for (i in 1:n){
    for (i in 1:n){
      avH.theta[i,] <- acos(sqrt(apply(avAA, 3, proj, b = avH.vec[,i]))) * (180/pi)
    }
  }
  #angles between the eigenvectors posterior mean H and the posterior mean subspaces of each population
  MCMC.H.val <- matrix(, MCMCsamp, n)
  colnames(MCMC.H.val) <- paste("h", 1:n, sep="")
  for (i in 1:n){
    MCMC.H.val[,i] <- apply(MCMC.H, 3, proj, b = avH.vec[,i])
  }
  #posterior distribution of the genetic variance for the eigenvectors of posterior mean H 
  MCMC.H.theta <- array(, c(n, m, MCMCsamp))
  rownames(MCMC.H.theta) <- paste("h", 1:n, sep="")
  colnames(MCMC.H.theta) <- dimnames(Gs)[[3]]
  for(i in 1:n){
    for(j in 1:MCMCsamp){
      MCMC.H.theta[i,,j] <- acos(sqrt(apply(MCMC.AA[,,,j], 3, proj, b = avH.vec[,i]))) * (180/pi)
    }
  }
  #posterior distribution of the angles between the eigenvectors of posterior mean H and the MCMC samples of the subspaces of each population
  list(avAA = avAA, avH = avH, MCMC.AA = MCMC.AA, MCMC.H = MCMC.H, MCMC.H.val = MCMC.H.val, MCMC.H.theta = MCMC.H.theta)
}
#END

val <- matrix(, n, m)
for (i in 1:m){
  avG <- apply(Garray, 1:3, mean)
  val[,i] <- round(cumsum(t(eigen(avG[,,i])$values))/sum(eigen(avG[,,i])$values)*100)
}

n.vec <- apply(ifelse(round(val,1) < 90, 1, 0), 2, sum)+1
# Common subspace of observed P matrices
MCMCG.kr <- kr.subspace(Garray, vec = n.vec)


# Common subspace of randomized P matrices
MCMCG.kr.rand <- kr.subspace(Garray.r2, vec = n.vec)
HPD.H.val <- cbind(HPDinterval(as.mcmc(MCMCG.kr$MCMC.H.val)), HPDinterval(as.mcmc(MCMCG.kr.rand$MCMC.H.val)))


### Plot
dfsub<- rbind(pivot_longer(as.data.frame(MCMCG.kr$MCMC.H.val),cols = 1:7, names_to = c("vec")),
              pivot_longer(as.data.frame(MCMCG.kr.rand$MCMC.H.val),cols = 1:7, names_to = c("vec")))
dfsub$mat <- c(rep("Observed",nrow(dfsub)/2),rep("Random",nrow(dfsub)/2))

# dfsub<- merge(pivot_longer(as.data.frame(MCMCG.kr$MCMC.H.val),cols = 1:7, names_to = c("vec")),
#               pivot_longer(as.data.frame(MCMCG.kr.rand$MCMC.H.val),cols = 1:7, names_to = c("vec")),
#               by = "vec",suffixes = c("Observed","Random"))
# dfsub$mat <- c(rep("Observed",nrow(dfsub)/2),rep("Random",nrow(dfsub)/2))

psub <-ggplot(dfsub, aes(y = value , x = vec, colour = mat)) +  
  stat_pointinterval(aes(y = value, colour = mat), 
                     point_interval = mode_hdi, position = position_dodge(width = .3)) +
  scale_color_manual(values = c("darkblue","grey"))+
  ggtitle("(g)") +
  theme_classic() +
  xlab(expression(paste("Eigenvectors of ",bold(H)))) + 
  ylab(expression(paste("Eigenvalues of ",bold(H)))) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        plot.title = element_text(size = 17, face = "bold"),
        axis.text.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14, margin = margin(t = 10)),
        axis.text.y = element_text(size = 15),
        legend.position = "none")

grid.arrange(pVtot,ptheta,p.omega,psub,ncol = 2,nrow= 2)

# Trait combinations that underline common subspace
round(eigen(MCMCG.kr$avH)$vectors, 2)

# Assess how close are populations to the eigenvector of H : Examine angle between H and AA^t of each population.
round(apply(MCMCG.kr$MCMC.H.theta, 1:2, posterior.mode), 1)