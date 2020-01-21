library(tidyr)
library(tidyverse)
library(rgl)

setwd("~/Documents/Github/ai/ecoli")
landscape<-read.csv("matrix3D.txt", sep="\t")

r3dDefaults$windowRect <- c(0,0, 4000, 4000)

landscape<-landscape %>% distinct(genomeSize, err_rate, .keep_all = TRUE)

m_QV<-pivot_wider(landscape, id_cols = genomeSize, names_from = err_rate, values_from = QV)
m_QV<- m_QV[order(m_QV$genomeSize),]

m_cost<-pivot_wider(landscape, id_cols = genomeSize, names_from = err_rate, values_from = cost)
m_cost<- m_cost[order(m_cost$genomeSize),]
z_cost<-data.matrix(unname(m_cost[,-1]))

m_ctgn<-pivot_wider(landscape, id_cols = genomeSize, names_from = err_rate, values_from = ctgn)
m_ctgn<- m_ctgn[order(m_ctgn$genomeSize),]
z_ctgn<-data.matrix(unname(m_ctgn[,-1]))

m_gsize<-pivot_wider(landscape, id_cols = genomeSize, names_from = err_rate, values_from = gsize)
m_gsize<-m_gsize[order(m_gsize$genomeSize),]

z_gsize<-data.matrix(unname(m_gsize[,-1]))

y<-na.omit(as.numeric(colnames(m_QV)))
x<-deframe(tibble(m_QV$genomeSize))
z_QV<-data.matrix(unname(m_QV[,-1]))

nbcol = 4
color = rev(heat.colors(nbcol))
zcol  = cut(z_QV, nbcol)

draw_lines <- function(file,color,x,y,z,z_labels) {

  path<-read.csv(file, sep="\t")
  seg_x <- c(path$x)
  seg_y <- c(path$y)
  seg_z <- c(path$cost)
  
  draw_landscape(x,y,z,z_labels)
  
  #peak<-filter(landscape, QV == max(na.omit(z_QV)))
  #arrow3d(cbind(peak[,1],peak[,2],peak[,3]+0.2),cbind(peak[,1],peak[,2],peak[,3]+0.03), col = "red", size=20)
  
  
  start_time=0
  dur=3

  for(i in 1:dim(path)[1])
  {
    if (i==1) {
      points3d(cbind(seg_x[1],seg_y[1],seg_z[1]), col = color, size=6)
      movie3d(spin3d(axis = c(0, 0, 1),rpm = 0.5), startTime = start_time, duration = dur,dir = "figures/", convert = FALSE)  
      start_time=start_time+dur
      dur=1

      
    } else {

      M <- par3d("userMatrix")
      
      points3d(cbind(seg_x[i],seg_y[i],seg_z[i]), col = "black", size=3)
      lines3d(cbind(rbind(seg_x[i],seg_x[i-1]),rbind(seg_y[i],seg_y[i-1]),rbind(seg_z[i],seg_z[i-1])), add=TRUE, lwd = 3, col = color)

      movie3d( par3dinterp(time = c(start_time,start_time+dur), userMatrix = list(M,rotate3d(M, 0.05, 0, 0, 1) ) ), startTime = start_time, duration = start_time+dur,dir = "figures/", convert = FALSE)

      start_time=start_time+dur
    }
  }
  
  movie3d( par3dinterp(time = c(start_time,dur), userMatrix = list(M,M) ), startTime = start_time, duration = start_time+1,dir = "figures/", convert = FALSE)
    
}

draw_landscape <- function(x,y,z,z_labels) {

par3d(zoom=1.5)
persp3d(x,sort(y),z, col=color[zcol], xlab="",ylab="",zlab="", axes = FALSE, box = TRUE, alpha=0.8)
#title3d('Assembly size and QV by Canu parameters Genome size and Error rate')
axis3d(edge= 'x', at=seq(min(x), max(x),by=200000), labels = round(seq(min(x), max(x),by=200000)/1000000, digits=4))
axis3d(edge= 'y-', at=seq(min(y), max(y),by=0.01))
axis3d(edge= 'z++', at=seq(min(z), max(z),by=z_labels), labels = round(seq(min(z), max(z),by=z_labels), digits=1))

mtext3d("Canu error rate", edge= 'y-', line = 4, at = NULL, pos = NA)
mtext3d("Canu genome size (Gbp)", edge= 'x', line = 4, at = NULL, pos = NA)
mtext3d("Cost (dQV,dGenome_size,dContig_N)", edge= 'z-+', line = 4, at = NULL, pos = NA)

}

draw_landscape(x,y,z_QV,0.1)
peak<-filter(landscape, QV == max(na.omit(z_QV)))
arrow3d(cbind(peak[,1],peak[,2],peak[,3]+0.2),cbind(peak[,1],peak[,2],peak[,3]+0.03), col = "red", size=20)
movie3d(spin3d(axis = c(0, 0, 1),rpm = 0.5), startTime = 0, duration = 40,dir = "figures/", convert = TRUE)

draw_lines("test4.txt", "purple")
draw_lines("test15.txt", "green",x,y,z_QV,0.1)
draw_lines("test21.txt", "blue")
draw_lines("test16.txt", "aquamarine")
draw_lines("test19.txt", "darkslateblue")
draw_lines("test25.txt", "darkolivegreen1",x,y,z_QV,0.1)

draw_landscape(z_QV)
rgl.quads( x = c(3915000,3915000,3915000,3915000), y = c(min(y), max(y), max(y), min(y)),
           z = c(min(na.omit(z_QV)),min(na.omit(z_QV)),max(na.omit(z_QV)),max(na.omit(z_QV))), alpha=0.6)

rgl.quads( x = c(4395000,4395000,4395000,4395000), y = c(min(y), max(y), max(y), min(y)),
           z = c(min(na.omit(z_QV)),min(na.omit(z_QV)),max(na.omit(z_QV)),max(na.omit(z_QV))), alpha=0.6)
movie3d(spin3d(axis = c(0, 0, 1),rpm = 0.5), startTime = 0, duration = 40,dir = "figures/", convert = TRUE)

draw_landscape(z_QV)
arrow3d(cbind(3550000,0.23,36.7),cbind(3550000,0.24,36.7), col = "red", s = 1/3, theta = pi/12, type="flat")
arrow3d(cbind(3500000,0.227,36.7),cbind(3900000,0.227,36.7), col = "red", s = 1.5/3, theta = pi/3, type="flat")


draw_landscape(z_cost,1)
draw_lines("test25.txt", "darkolivegreen1",x,y,z_cost,0.5)

fil_landscape<-filter(landscape, between(genomeSize, 4000000, 4200000))
fil_m_cost<-pivot_wider(fil_landscape, id_cols = genomeSize, names_from = err_rate, values_from = cost)
fil_m_cost<- fil_m_cost[order(fil_m_cost$genomeSize),]
fil_z_cost<-data.matrix(unname(fil_m_cost[,-1]))
fil_y<-na.omit(as.numeric(colnames(fil_m_cost)))
fil_x<-deframe(tibble(fil_m_cost$genomeSize))


draw_lines("test25.txt", "darkolivegreen1",fil_x,fil_y,fil_z_cost,1)


