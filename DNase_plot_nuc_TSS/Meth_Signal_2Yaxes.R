args<-commandArgs(T)

data_meth<-read.table(args[1],header=F)
data_mnase<-read.table(args[2],header=F)

out_pdf<-paste("/datd/lilin/Project/NOM-Seq/mouse_NOME_seq_select_brok-reRun/05.per_plot/GCH_MNaseSignal/",args[3],".GCH_MNaseSignal.pdf",sep="")
pdf(out_pdf,width=6, height=6)

bin_num <- 100
data_meth <- data_meth[,-1]	
data_mnase <- data_mnase[,-1]  

left<-seq(from=1,by=1,length.out=bin_num/2) 
right<-seq(from=max(left+1),to=max(left+1)+bin_num/2 - 1,length.out=bin_num/2)
x<-c(left,right)
labelx<-c(seq(from=1,by=1,length.out=bin_num/2),seq(from=max(left+1),to=max(left+1)+bin_num/2 - 1,length.out=bin_num/2))
   
par(mar=c(2,0,0,1.5),oma=c(5,5,5,5))

plot(x=x,y=colMeans(data_meth,na.rm=T),type="l",xaxs="i",yaxs="i",ylim=c(0,0.7),axes=FALSE,xlab="",ylab="",col="red",lwd=2)
box(bty="l")
axis(side=2,at=pretty(c(0,0.7),4),tcl=-0.2,las=2)
mtext(side=2,"Chromatin accessibility(GCH) ",cex=0.8,line=2.5)
labeltext<-c("-1kb","0%","1kb")
mtext(side=1,at=labelx[c(1,50,100)],c("-1kb","0","+1kb"),cex=0.7,line=1)

par(new=TRUE)  #no new picture
plot(x=x,y=colMeans(data_mnase,na.rm=T),type="l",xaxs="i",yaxs="i",axes=FALSE,xlab="",ylab="",ylim=c(2,7),col="blue",lwd=2)
axis(side=4,at=pretty(c(2,7),4),tcl=-0.2,las=2)
mtext(side=4,"Signal of MNase-Seq ",cex=0.8,line=2.5)

legend("topright",legend=c("Methylation","Signal"),col=c("red","blue"),lwd=1,cex=0.6)
dev.off()


