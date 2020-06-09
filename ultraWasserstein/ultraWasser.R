

desired_length=10
list=vector(mode = "list", length = desired_length)
dist.vec=c(sort(unique(dist.mat[1,]), decreasing = T),0)
set = numeric(length(dist.vec))
ultraWasser = 0
for(i in 1:length(dist.vec)){
  set=index[dist.mat[1,]<dist.vec[i]]
  list[position]=distmat[1,dist.mat[1,]>dist.vec[i]]
  value = value + abs(mu[set]-nu[set])*((dist.vec[i]/2)^p-(dist.vec[i-1]/2)^p)
  dist.mat =dist.mat
}