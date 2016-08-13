library(plyr)



d0 = list()


r = 1

for( n in c(120,240) ){
  for (m in c(80, 160) ){
    for (C1 in c(0.25, 0.5, 0.75, 1) ){
      tit = paste( "DGP_linearIV_rho_0.6_B_corr_n_", n, "_m_", m, "_Rep_500_C1_", C1, "_seed_301.csv", sep = ""  )
      d0[[r]] = read.csv(file = tit, header = T, sep = ",")    
      r = r + 1
    }
  }
}

d1 = ldply(d0)

dREL = read.csv(file = "IV_REL_MSE.csv", header = T, sep = ",")
dREL1 = dREL[, -(1:2)]
dREL1 = t(as.matrix( dREL1)) 
dREL1 = as.vector(dREL1)
dREL1

d1 = cbind( dREL1, d1)
write.csv(d1, file = "sense_all_IV.csv")
