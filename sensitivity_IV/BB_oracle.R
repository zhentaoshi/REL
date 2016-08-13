
library(ggplot2)
library(grid)

d0 <- read.csv("BB.csv", sep = ",", header = T)
d1 <- d0[1:100, ]

xbase <- seq(-4, 4, by = 0.01)
pdfN <- dnorm(xbase)

plot(xbase, pdfN, type = "l", lwd = 3)
lines(density(d0[,1]), col = "red",)
lines(density(d1[,3]), col = "blue")
lines(density(d1[,5]), col = "brown", lwd = 2)





plot(xbase, pdfN, type = "l", lwd = 3)
lines(density(d0[,2]), col = "red")
lines(density(d1[,4]), col = "blue")
lines(density(d1[,6]), col = "brown", lwd = 2)


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



