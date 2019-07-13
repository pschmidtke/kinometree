manning=read.table("3cs9_manning_full_sa_blosum62_manning.nwk_enrichment.csv",header=T)
ours=read.table("3cs9_blosum90.nwk_enrichment.csv",header=T)

svg(filename="3cs9_ef.svg", 
    width=5, 
    height=4, 
    pointsize=8)

ylim=c(0,1.0)
par(las=2)
plot(ours[,"meanActives"],ty="l",ylim=ylim,col="blue",lwd=3,ylab="relative Enrichment",xlab="Clade Size")

arrows(ours[,"cladeSize"], ours[,"meanActives"]-ours[,"stdActives"], ours[,"cladeSize"], ours[,"meanActives"]+ours[,"stdActives"], length=0.05, angle=90, code=3,col="#1111FF15")

par(new=T)

plot(manning[,"meanActives"],ty="l",ylim=ylim,col="orange",lwd=3,ylab="",xlab="")
arrows(manning[,"cladeSize"], manning[,"meanActives"]-manning[,"stdActives"], manning[,"cladeSize"], manning[,"meanActives"]+manning[,"stdActives"], length=0.05, angle=90, code=3,col="#ffc60055")

par(new=T)
plot(manning[,"meanInactives"],ty="l",ylim=ylim,col="red",lwd=3,ylab="",xlab="")
arrows(manning[,"cladeSize"], manning[,"meanInactives"]-manning[,"stdInactives"], manning[,"cladeSize"], manning[,"meanInactives"]+manning[,"stdInactives"], length=0.05, angle=90, code=3,col="#ff010115")


par(new=T)
plot(ours[,"meanInactives"],ty="l",ylim=ylim,col="grey",lwd=3,ylab="",xlab="")
arrows(ours[,"cladeSize"], ours[,"meanInactives"]-ours[,"stdInactives"], ours[,"cladeSize"], ours[,"meanInactives"]+ours[,"stdInactives"], length=0.05, angle=90, code=3,col="#10101015")

legend("bottomright",c("EF actives - local Sequence","EF inactives - local Sequence","EF actives - Manning","EF inactives - Manning"),fill=c("blue","grey","orange","red",pch=c(22,22,22,22)))

dev.off()