library(R2admb)
if ((s <- setup_admb())!="") {
    ## run only if we can
    file.copy(system.file("tplfiles","ReedfrogSizepred0.tpl",package="R2admb"),"tadpole.tpl")
    tadpoledat <-  data.frame(TBL = rep(c(9,12,21,25,37),each=3),
                              Kill = c(0,2,1,3,4,5,0,0,0,0,1,0,0,0,0L),
                              nexposed=rep(10,15))
    par1 <- list(c=0.45,d=13,g=1)
    ## getvals(par1,"value",valsOK=TRUE)
    m1 <- do_admb("tadpole",
                  data=c(list(nobs=15),tadpoledat),
                  params=par1,
                  bounds=list(c=c(0,1),d=c(0,50),g=c(-1,25)),
                  run.opts=run.control(checkparam="write",
                  checkdata="write",clean="all"))
    par2 <- list(c=list(0.45,bounds=c(0,1)),
                 d=list(13,bounds=c(0,50)),
                 g=list(1,bounds=c(-1,25)))
    ## getvals(par2)
    ## getvals(par2,"value",TRUE)
    m2 <- do_admb("tadpole",
                  data=c(list(nobs=15),tadpoledat),
                  params=par2,
                  run.opts=run.control(checkparam="write",
                  checkdata="write",clean="all"))
    stopifnot(all.equal(m1,m2))
    unlink("tadpole.tpl")
}
