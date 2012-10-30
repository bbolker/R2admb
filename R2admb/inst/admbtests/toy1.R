library(lme4)
gm1 <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
                family = binomial, data = cbpp)
  
X <- model.matrix(~period,data=cbpp)
Zherd <- model.matrix(~herd-1,data=cbpp)
library(R2admb)
setup_admb()

tmpdat <- list(X=X,Zherd=Zherd,
                 incidence=cbpp$incidence,size=cbpp$size,
                 nobs=nrow(cbpp))
d1 <- do_admb("toy1",
              data=tmpdat,
              re=TRUE,
              params=list(beta=rep(0,ncol(X)),sigma_herd=0.1),
              bounds=list(sigma_herd=c(0.0001,20)),
              re_vectors=c(u_herd=ncol(Zherd)),
              checkdata="write",checkparam="write",
              mcmc=TRUE,
              mcmcpars=c("beta","sigma_herd"),
              clean=FALSE)

save("d1",file="toy1_runs.RData")
unlink(c("toy1","toy1_gen.tpl"))
