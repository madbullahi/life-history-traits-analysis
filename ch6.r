1 + 1
2 + 5 - 3
?sqrt()


require(stats)
require(graphics)
xx <- -9:9
plot(xx, sqrt(xx), col="red")
lines(spline(xx, sqrt(abs(xx)), n=101), col="pink")