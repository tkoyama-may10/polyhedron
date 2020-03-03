 x <- read.table("tukey.txt")

 postscript(file="./tukey.eps",width=6,height=6,horizontal=FALSE,
  onefile=FALSE,paper="special",family="Helvetica")
 N <- length(x$V1)
 xl <- c(min(x$V1),max(x$V1))
 yl <- c(0,1.1)
 plot(x$V1,x$V2, xlim=xl,ylim=yl,xlab="t",ylab="F(t)",lty=1,type="b")
 par(new=T)
 plot(x$V1,x$V3,,xlim=xl,ylim=yl,xlab="",ylab="",lty=2,type="b",axes=F)
 par(new=T)
 plot(x$V1,x$V4, xlim=xl,ylim=yl,xlab="",ylab="",lty=3,type="b",axes=F)
 par(new=T)
 plot(x$V1,x$V5, xlim=xl,ylim=yl,xlab="",ylab="",lty=4,type="b",axes=F)
 par(new=T)
 legend("bottomright",
        legend = c("dim=2", "dim=5", "dim=10", "dim=20"),
        col = c("black", "black", "black", "black"),
        pch = c(1,1,1,1),
        lty = c(1,2,3,4)
       )
 dev.off()

