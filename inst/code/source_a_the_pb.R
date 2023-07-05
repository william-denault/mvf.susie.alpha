rm(list=ls())
devtools::load_all(".")

source("D:/Document/Serieux/Travail/Package/mvf.susie.alpha/R/computational_routine.R")


load("D:/Document/Serieux/Travail/Package/mvf.susie.alpha/inst/pb_input.RData")
X <- pb_input[[3]]$X
Y <- pb_input[[3]]$Y
true_pos  <- pb_input[[3]]$true_pos

source("D:/Document/Serieux/Travail/Package/mvf.susie.alpha/R/operation_on_multfsusie_obj.R")
Y <-list(Y_f= list(pb_input[[3]]$Y$Y_f[[1]]),
         Y_u=NULL)
m1 <- multfsusie(Y=Y,
                 X=X,
                 L=11)
m1$cs

Y <-list(Y_f= list(pb_input[[3]]$Y$Y_f[[1]],
                   pb_input[[3]]$Y$Y_f[[2]]),
         Y_u=NULL)
m2 <- multfsusie(Y=Y,
                 X=X,
                 L=11)



Y <-list(Y_f= list(
                   pb_input[[3]]$Y$Y_f[[2]]),
         Y_u= pb_input[[3]]$Y$Y_u)
m3 <- multfsusie(Y=Y,
                 X=X,
                 L=11)

m3$cs
m1$cs
m2$cs
plot(m1$fitted_func[[1]][[1]])
lines(m2$fitted_func[[3]][[1]])
plot(m1$fitted_func[[2]][[1]])
lines(m2$fitted_func[[1]][[1]])
plot(m1$fitted_func[[4]][[1]])
lines(m2$fitted_func[[4]][[1]])


plot(m1$fitted_func[[5]][[1]])
lines(m2$fitted_func[[2]][[1]])
plot(m1$fitted_func[[3]][[1]], type = "l", ylim=c(-8,4), ylab = "estimated effect")
lines(m2$fitted_func[[2]][[1]], col="red")
lines( m1$fitted_func[[3]][[1]], col="blue")
lines( m1$fitted_func[[5]][[1]], col="orange")
lines( m1$fitted_func[[3]][[1]] + m1$fitted_func[[5]][[1]], col="green")
legend(x=1, y=-4,
       legend= c("fsusie CS 1",
                 "fsusie CS 2",
                 "estimated effect CS 1+CS2 fsusie",
                 "CS 1 mvfsusie"),
       col= c("blue", "orange", "green", "red"),
       lty= rep(1,4))


plot(m3$fitted_func[[2]][[1]])
