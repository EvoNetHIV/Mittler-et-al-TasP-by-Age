# Process data structures created by CollateEvonetData.R
# rand, cd4_rand, u35_rand, u30_rand, u25_rand, spvl_rand
x_axis <- 10*c(0,1,2,3,4,5,6,7,8,9,10)
jit <- rep(0.3,reps=length(x_axis))
x_axis_label <- "TasP Target [Percent Treated]" # "Growth Rate [Percent per year]" # 
par(mfrow=c(2,2))

treatments <- c("CD4","rand","spvl","u25","u30","u35")
#treatments <- c("CD4","rand","spvl","u25","u30")

# Limits for y-axes
max_incid = 1.6
max_AIDS = 550
max_p_yrs_tx <- 3500
max_outside <- 110
incid_comp_line <- 10.05

# Plot incidence
plot(x_axis+jit,   as.numeric(rand[4,]), xlab=x_axis_label,xaxt = "n",type="b",pch=1,lwd=2,ylim=c(0,max_incid),
                                                        xlim=c(0,max(x_axis)),ylab="Incidence Last 5 Years",col="black")
lines(x_axis,rep(incid_comp_line,length(x_axis)),col="grey",lwd=0.5)
axis(1,x_axis)
arrows(x_axis+jit, as.numeric(rand[5,]), x_axis+jit, as.numeric(rand[6,]), code=3, length=0.02, angle = 90, col="black")

lines( x_axis, u25_rand[4,], col="darkred",pch=5, lwd=2, type="b")
arrows(x_axis, as.numeric(u25_rand[5,]), x_axis, as.numeric(u25_rand[6,]), code=3, length=0.02, angle = 90, col="darkred")

if ("u30" %in% treatments) {
  lines( x_axis-2*jit, u30_rand[4,], col="red",pch=2, lwd=2, type="b")
  arrows(x_axis-2*jit, as.numeric(u30_rand[5,]), x_axis-2*jit, as.numeric(u30_rand[6,]), code=3, length=0.02, angle = 90, col="red")
}
if ("u35" %in% treatments) {
  lines( x_axis-3*jit, u35_rand[4,], pch=4, col="orange",lwd=2,type="b")
  arrows(x_axis-3*jit, as.numeric(u35_rand[5,]), x_axis-3*jit, as.numeric(u35_rand[6,]), code=3, length=0.02, angle = 90, col="orange")
}
if ("spvl" %in% treatments) {
  lines( x_axis-jit, spvl_rand[4,], pch =6, col="blue",lwd=2,type="b")
  arrows(x_axis-jit, as.numeric(spvl_rand[5,]), x_axis-jit, as.numeric(spvl_rand[6,]), code=3, length=0.02, angle = 90,col="blue")
}
if ("CD4" %in% treatments) {
  lines( x_axis+2*jit, cd4_rand[4,], col="green",pch=0, lwd=2,type="b")
  arrows(x_axis+2*jit, as.numeric(cd4_rand[5,]), x_axis+2*jit, as.numeric(cd4_rand[6,]), code=3, length=0.02, angle = 90, col="green")
}
if ("m23w27" %in% treatments) {
  lines( x_axis, u23_u27_rand[4,], col="purple",pch=4, lwd=2,type="b")
  arrows(x_axis, as.numeric(u23_u27_rand[5,]), x_axis, as.numeric(u23_u27_rand[6,]), code=3, length=0.02, angle = 90, col="purple")
}
#pch=c(0,1,6,5,2),
#       col=c("green","black","blue","orange","red"), 
legend( x="topright", 
      legend=c("CD4<500","Random","SPVL","Under age 25","Under age 30","Under age 35"),
#        legend=c("CD4<500","Random","Under age 25"),
       col=  c("green","black","blue","darkred","red","orange"), 
#        col=  c("green","black","darkred"), 
text.col= c("green","black","blue","darkred","red","orange"), # green, black, blue, brown, red, darkred
#text.col= c("green","black","darkred"), # green, black, blue, brown, red, darkred
         lwd=2.0,bty="n",cex = 1,pt.cex=1.2,
       pch=c(0,1,6,5,2, 4)) # 0 1 6 5
#        pch=c(0,1,5)) # 0 1 6 5

#AIDS 16, 17 , 18
#Pills 10, 11, 12
#non-Pr 13, 14, 15

# Plot AIDs Deaths 16, 17, and 18 instead of 4, 5, and 6
plot(x_axis+jit,   as.numeric(rand[16,]), xlab=x_axis_label,xaxt = "n",type="b",pch=1,lwd=2,
     ylim=c(0,max_AIDS),
     xlim=c(0,max(x_axis)),
     ylab="AIDS Deaths",col="black")
axis(1,x_axis)
arrows(x_axis+jit, as.numeric(rand[17,]), x_axis+jit, as.numeric(rand[18,]), code=3, length=0.02, angle = 90, col="black")

lines( x_axis, u25_rand[16,], col="darkred",pch=5, lwd=2, type="b")
arrows(x_axis, as.numeric(u25_rand[17,]), x_axis, as.numeric(u25_rand[18,]), code=3, length=0.02, angle = 90, col="darkred")

if ("u30" %in% treatments) {
  lines( x_axis-2*jit, u30_rand[16,], col="red",pch=2, lwd=2, type="b")
  arrows(x_axis-2*jit, as.numeric(u30_rand[17,]), x_axis-2*jit, as.numeric(u30_rand[18,]), code=3, length=0.02, angle = 90, col="red")
}
if ("u35" %in% treatments) {
  lines( x_axis-3*jit, u35_rand[16,], pch=4, col="orange",lwd=2,type="b")
  arrows(x_axis-3*jit, as.numeric(u35_rand[17,]), x_axis-3*jit, as.numeric(u35_rand[18,]), code=3, length=0.02, angle = 90, col="orange")
}
if ("spvl" %in% treatments) {
  lines( x_axis-jit, spvl_rand[16,], pch =6, col="blue",lwd=2,type="b")
  arrows(x_axis-jit, as.numeric(spvl_rand[17,]), x_axis-jit, as.numeric(spvl_rand[18,]), code=3, length=0.02, angle = 90,col="blue")
}
if ("CD4" %in% treatments) {
  lines( x_axis+2*jit, cd4_rand[16,], col="green",pch=0, lwd=2,type="b")
  arrows(x_axis+2*jit, as.numeric(cd4_rand[17,]), x_axis+2*jit, as.numeric(cd4_rand[18,]), code=3, length=0.02, angle = 90, col="green")
}

# Plot Pills Taken 10, 11, 12 instead of 4-6
plot(x_axis+jit,   as.numeric(rand[10,]), xlab=x_axis_label,xaxt = "n",type="b",pch=1,lwd=2,
     ylim=c(0,max_p_yrs_tx),
     xlim=c(0,max(x_axis)),ylab="Person-Years of Therapy",col="black")
axis(1,x_axis)
arrows(x_axis+jit, as.numeric(rand[11,]), x_axis+jit, as.numeric(rand[12,]), code=3, length=0.02, angle = 90, col="black")

lines( x_axis, u25_rand[10,], col="darkred",pch=5, lwd=2, type="b")
arrows(x_axis, as.numeric(u25_rand[11,]), x_axis, as.numeric(u25_rand[12,]), code=3, length=0.02, angle = 90, col="darkred")

if ("u30" %in% treatments) {
  lines( x_axis-2*jit, u30_rand[10,], col="red",pch=2, lwd=2, type="b")
  arrows(x_axis-2*jit, as.numeric(u30_rand[11,]), x_axis-2*jit, as.numeric(u30_rand[12,]), code=3, length=0.02, angle = 90, col="red")
}
if ("u35" %in% treatments) {
  lines( x_axis-3*jit, u35_rand[10,], pch=4, col="orange",lwd=2,type="b")
  arrows(x_axis-3*jit, as.numeric(u35_rand[11,]), x_axis-3*jit, as.numeric(u35_rand[12,]), code=3, length=0.02, angle = 90, col="orange")
}
if ("spvl" %in% treatments) {
  lines( x_axis-jit, spvl_rand[10,], pch =6, col="blue",lwd=2,type="b")
  arrows(x_axis-jit, as.numeric(spvl_rand[11,]), x_axis-jit, as.numeric(spvl_rand[12,]), code=3, length=0.02, angle = 90,col="blue")
}
if ("CD4" %in% treatments) {
  lines( x_axis+2*jit, cd4_rand[10,], col="green",pch=0, lwd=2,type="b")
  arrows(x_axis+2*jit, as.numeric(cd4_rand[11,]), x_axis+2*jit, as.numeric(cd4_rand[12,]), code=3, length=0.02, angle = 90, col="green")
}


# Plot Treated non-targeted 13, 14, 15 instead of 4, 5, and 6
plot(x_axis+jit,   (unlist(u25_rand[13,])*100/unlist(rand[13,])), xlab=x_axis_label,xaxt = "n",type="b",pch=5,lwd=2,
     ylim=c(0,max_outside),
     xlim=c(0,max(x_axis)),ylab="Percent Treated Outside Target Group",col="darkred")
axis(1,x_axis)
arrows(x_axis+jit, (unlist(u25_rand[14,])*100/unlist(rand[13,])), x_axis+jit, (unlist(u25_rand[15,])*100/unlist(rand[13,])), code=3, length=0.02, angle = 90, col="darkred")

if ("u30" %in% treatments) {
  lines( x_axis-2*jit, (unlist(u30_rand[13,])*100/unlist(rand[13,])), col="red",pch=2, lwd=2, type="b")
  arrows(x_axis-2*jit, (unlist(u30_rand[14,])*100/unlist(rand[13,])), x_axis-2*jit, (unlist(u30_rand[15,])*100/unlist(rand[13,])), code=3, length=0.02, angle = 90, col="red")
}
if ("u35" %in% treatments) {
  lines( x_axis-3*jit, (unlist(u35_rand[13,])*100/unlist(rand[13,])), pch=4, col="orange",lwd=2,type="b")
  arrows(x_axis-3*jit, (unlist(u35_rand[14,])*100/unlist(rand[13,])), x_axis-3*jit, (unlist(u35_rand[15,])*100/unlist(rand[13,])), code=3, length=0.02, angle = 90, col="orange")
}
if ("spvl" %in% treatments) {
  lines( x_axis-jit, (unlist(spvl_rand[13,])*100/unlist(rand[13,])), pch =6, col="blue",lwd=2,type="b")
  arrows(x_axis-jit, (unlist(spvl_rand[14,])*100/unlist(rand[13,])), x_axis-jit, (unlist(spvl_rand[15,])*100/unlist(rand[13,])), code=3, length=0.02, angle = 90,col="blue")
}
if ("CD4" %in% treatments) {
  lines( x_axis+2*jit, (unlist(cd4_rand[13,])*100/unlist(rand[13,])), col="green",pch=0, lwd=2,type="b")
  arrows(x_axis+2*jit, (unlist(cd4_rand[14,])*100/unlist(rand[13,])), x_axis+2*jit, (unlist(cd4_rand[15,])*100/unlist(rand[13,])), code=3, length=0.02, angle = 90, col="green")
}







