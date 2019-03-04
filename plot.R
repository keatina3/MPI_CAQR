dat=read.table("a.out", sep=",",header=TRUE)

xvals=log10(dat$block)

plot(xvals[1:5],dat$t[1:5],pch="|",col='red',cex=2,
        ylim=c(min(dat$t)-5,max(dat$t)+5), xaxt='n',xlab = "Panel Size", ylab = "Time (s)")
points(xvals[6:10],dat$t[6:10],pch="|",col='blue',cex=2)
