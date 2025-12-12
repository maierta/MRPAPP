# 3D FS plot
require(misc3d)
require(rgl)
options(rgl.printRglwidget = TRUE)
open3d()
read.table("ek_test1.txt", sep=",", header = F) -> b
nkx = 32
nky = 32
nkz = 8
nbands = 17
col=c("red","blue","green","orange","yellow","brown","purple","grey","pink","beige")
add <- F
for (i in 4:(nbands+3)) {
	FS <- range(b[,i])[1]*range(b[,i])[2]
	if (FS < 0) {
		contour3d(array(b[,i],dim=c(nkx,nky,nkz)),level=0,color=col[i-3],aspect=c(1,10),frames=T,perspective=T,distance=2,smooth=T,lwd=2,add=add,ticktype="")
		add <- T
	}
}
aspect3d(1,1,1.44592)
# axes3d(edges = "bbox",labels=FALSE,tick=FALSE)
box3d()
my.ticks=c(1,16,32)
# axis3d('x',at=my.ticks,labels=c(-pi,0,pi))
# axis3d('y',at=my.ticks,labels=c(-pi,0,pi))
# axis3d('z',at=my.ticks,labels=c(-pi,0,pi))
title3d("Fermi surface","","kx","ky","kz")

