require(fifer)
auto.layout(16)
for (i in seq(from=1, to=3, length.out=4)){
  for (j in seq(from=1, to=3, length.out=4)){
    par1()
    hist(rbeta(1000, i,j)*2 - 1, main=paste0("Shape1 = ", round(i, digits=1), ", Shape2 = ", round(j, digits=1)))
  }
}

hist((rbeta(1000, 4, 4)*2-1))
