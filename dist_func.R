##
## This R script draws the distribution function P(c)
##

tab <- read.table("raw_data4.txt")
x <- tab$V2 
n <- length(x) / 13
mat <- matrix(x, ncol = n)

points <- seq(0,6, by=0.5)
i <- 1:13


postscript(file="./fig1.eps",width=5,height=5,horizontal=FALSE,onefile=FALSE,paper="special",family="Helvetica")
par(new=F)
plot(points, mat[i,1], ylim=c(0.0,1.0), xlab="c", ylab="P(c)", lty=1,type="b")
par(new=T)
plot(points, mat[i,2], ylim=c(0.0,1.0), xlab="", ylab="",lty=2, type="b",xaxt="n",yaxt="n")
par(new=T)
plot(points, mat[i,3], ylim=c(0.0,1.0), xlab="", ylab="",lty=3, type="b",xaxt="n",yaxt="n")

labels <- c("dim=2","dim=5","dim=10")
cols <-   c("black","black","black")
pchs <-   c(1,1,1)
ltys <-   c(1,2,3)
legend("bottomright", legend = labels, col = cols, pch = pchs, lty=ltys)
dev.off()