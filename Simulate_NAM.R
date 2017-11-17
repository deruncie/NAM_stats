# simulate a locus bounded by two typed markers (expanded to several consecutive loci)
# With a known percent recombination between markers
# simulate random SNPs across locus in parents that can be imputed in lines
# simulate a single causal locus (with potentially multiple alleles across founders) in locus

nPop = 7
nLines = 100
r = 0.1  # percent recombination between markers
nM = 10 # this many regions bounded by typed markers
nSNP = 100 # per region
allele_freq = .2 # Average frequency in founders at all SNPs

# parent_genotypes = t(sapply(1:2^nPop,function(x) as.integer(intToBits(x)))[1:nPop,])
# topi = diff(cumsum(c(2^(0:7))))-1
# parent_genotypes = rbind(parent_genotypes[topi,],parent_genotypes[-topi,])
# nSNP = nrow(parent_genotypes)

p = dpois(1:nPop,allele_freq*nPop);p=p/sum(p)
parent_genotypes = matrix(0,nrow = nSNP*nM,ncol = nPop)
for(i in 1:nrow(parent_genotypes)){
  freq = sample(1:nPop,1,prob=p)
  parent_genotypes[i,sample(1:nPop,freq)] = 1
}

n = nPop*nLines
marker_genotypes = expand.grid(Rep=1:(nLines/2),M0 = c(0,1),Pop = factor(1:nPop))
marker_genotypes = marker_genotypes[,c('Pop','Rep','M0')]
prev_marker = 'M0'
for(marker in 1:nM){
  current_marker = paste0('M',marker)
  marker_genotypes[[current_marker]] = marker_genotypes[[prev_marker]]
  for(pop in 1:nPop){
    i = which(marker_genotypes$Pop == pop)
    i = sample(i,r*length(i))
    marker_genotypes[[current_marker]][i] = (marker_genotypes[[current_marker]][i]+1)%%2
  }
  prev_marker = current_marker
}
Image(marker_genotypes[,-c(1:2)])

# parent_genotypes = matrix(sample(c(0,1),nPop*nSNP,replace=T),nc=nPop)
# marker_genotypes = data.frame(Pop = factor(rep(1:nPop,each = nLines)),M1 = sample(c(0,1),n,replace=T))
# marker_genotypes$M2 = marker_genotypes$M1
# recombinations = sample(1:n,r*n)
# marker_genotypes$M2[recombinations] = (marker_genotypes$M2[recombinations] + 1) %% 2
#
# marker_genotypes = marker_genotypes[order(marker_genotypes$Pop,marker_genotypes$M1,marker_genotypes$M2),]

imputed_genotypes = matrix(0,nSNP*nM,n)
imputed_pops = array(0,dim=c(nSNP*nM,n,nPop))

dist = seq(0,1,length=nSNP)#rep(1,nSNP)#

prev_marker = 'M0'
for(marker in 1:nM){
  current_marker = paste0('M',marker)
  snps = nSNP*(marker-1) + 1:nSNP
  for(pop in 1:nPop){
    i = marker_genotypes$Pop == pop
    imputed_genotypes[snps,i & marker_genotypes[[prev_marker]] == 1 & marker_genotypes[[current_marker]] == 1] = parent_genotypes[snps,pop]
    imputed_genotypes[snps,i & marker_genotypes[[prev_marker]] == 0 & marker_genotypes[[current_marker]] == 1] = parent_genotypes[snps,pop] * dist
    imputed_genotypes[snps,i & marker_genotypes[[prev_marker]] == 1 & marker_genotypes[[current_marker]] == 0] = parent_genotypes[snps,pop] * (1-dist)

    imputed_pops[snps,i & marker_genotypes[[prev_marker]] == 1 & marker_genotypes[[current_marker]] == 1,pop] = 1
    imputed_pops[snps,i & marker_genotypes[[prev_marker]] == 0 & marker_genotypes[[current_marker]] == 1,pop] = dist
    imputed_pops[snps,i & marker_genotypes[[prev_marker]] == 1 & marker_genotypes[[current_marker]] == 0,pop] = 1-dist
  }
  prev_marker = current_marker
}


data = marker_genotypes
s2_pop = 1+0*seq(1,4,length=nPop)/4
mean_pop = rnorm(nPop,0,4)
data$s2 = s2_pop[data$Pop]
data$mean = mean_pop[data$Pop]
# data$base = with(data,mean + rnorm(nLines,0,sqrt(s2)))
data$base = with(data,mean + rnorm(n,0,sqrt(s2)))

causal = nSNP*floor(nM/2) + sample(1:nSNP,1)
# causal = which(rowSums(parent_genotypes)==7)[1]

effect = 1
# data$effect = imputed_genotypes[causal,] * effect
data$effect = NA
data$effect[imputed_genotypes[causal,] == 0] = 0
data$effect[imputed_genotypes[causal,] == 1] = effect
i = which(imputed_genotypes[causal,]%%1 != 0)
if(length(i)>0) data$effect[i] = sapply(i,function(x) sample(c(0,effect),size = 1,prob = c(1-imputed_genotypes[causal,x],imputed_genotypes[causal,x])))

# effect = rnorm(4,1,.002);effect = c(effect,rep(0,nPop-length(effect)))
# effect = sample(effect)
# data$effect = imputed_pops[90,,] %*% effect

data$y = data$base + data$effect

library(lmerTest)
res1=sapply(1:nrow(imputed_genotypes),function(i) {
  data$SNP = imputed_genotypes[i,]
  # lm1 = lmer(y~SNP+(1|Pop),data)
  lm1 = lm(y~Pop+SNP,data)
  anova(lm1)$P[2]
})
# plot(sort(-log10(res1)))
# plot(-log10(res1));abline(v=causal)

# library(meta)
# res2 = sapply(1:nrow(imputed_genotypes),function(i) {
#   data$SNP = imputed_genotypes[i,]
#   res = sapply(1:nPop,function(p) {
#     sub_data = subset(data,Pop==p)
#     if(sum(sub_data$SNP) == 0) return(c(NA,NA))
#     lm1 = lm(y~SNP,sub_data)
#     c(summary(lm1)$coef[2,1:2])
#   })
#   d = data.frame(est= res[1,],se = res[2,])
#   z=metagen(d$est,d$se,comb.random=F)
#   z$pval.fixed
#   # summary(z)
#   # weighted.mean(d[1,],1/d[2,]^2)
#   # summary(lm(est~1,weights = 1/se^2,d))$coef[1,4]
# })
#
res3=sapply(1:nrow(imputed_genotypes),function(i) {
  X = imputed_pops[i,,]
  # lm1 = lmer(y~X+(1|Pop),data)
  lm1 = lm(y~Pop+X,data)
  anova(lm1)$P[2]
})


plot(-log10(res1));abline(v=causal,col=2);abline(v=seq(0,nM)*nSNP)
# plot(-log10(res2));abline(v=causal,col=2);abline(v=seq(0,nM)*nSNP)
plot(-log10(res3));abline(v=causal,col=2);abline(v=seq(0,nM)*nSNP)

# plot(-log10(res1),abs(seq(1,nM*nSNP)-causal))
d = data.frame(y = -log10(res1),dist = abs(seq(1,nM*nSNP)-causal))[-causal,]
ggplot(d,aes(x=y,y=dist)) + geom_point() + geom_smooth()
