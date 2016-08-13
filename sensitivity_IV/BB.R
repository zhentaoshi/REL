
library(ggplot2)
library(grid)

d0 <- read.csv("BB_n_100_m_123_Rep_50_C1_0.5_seed_667BB.csv", sep = ",", header = T)
Rep <- nrow(d0)
d0 <- d0[1:Rep,3:6]
names(d0) <- c("B1.x", "B1.y", "B2.x", "B2.y")

R = nrow(d0)
# alpha.v = 1/R^(0.25)
alpha.v = 1
size.v = 6
center = c(1,1)

# d0 = d0[ sample(1:R, 15), ]
#########################
rm(p)

p <- ggplot(data = d0) + geom_point(aes(x = B1.x, y = B1.y), size = size.v, alpha = alpha.v )
p <- p + geom_point(aes(x = B2.x, y = B2.y), col = "red", size = size.v, alpha = alpha.v)
p <- p + geom_vline(xintercept = 1) + geom_hline(yintercept = 1)
p <- p + labs(x = expression(beta_1), y = expression(beta_2))

# a = 0.02
# 
#   x.lim = max( abs( quantile(c(d0$B1.x, d0$B2.x), c(a, 1-a) )  - center[1] ) )
#   y.lim = max( abs( quantile(c(d0$B1.y, d0$B2.y), c(a, 1-a) )  - center[2] ) )
#   xy.lim = max(x.lim, y.lim)

# p <- p +  xlim(c(-xy.lim, xy.lim) + center[1])  + ylim( c(-xy.lim, xy.lim) + center[2]) 
p

#############################333
p <- p + geom_segment(data = d0, aes(x = B1.x, y = B1.y, xend = B2.x, yend = B2.y ), arrow = arrow(), col = "red") 
p



