set.seed(1)
r1<- rnorm(n= 50, mean = 0, sd= 1)

# variance changes dramatically
set.seed(2)
r2<- rnorm(n= 100, mean = 0, sd= 10)

# variance changes less dramatically
set.seed(3)
r3<- rnorm(n= 100, mean = 0, sd= 5)

# change in mean
set.seed(4)
r4<- rnorm(n= 100, mean = 10, sd= 5)

# smaller change in mean
set.seed(5)
r5<- rnorm(n= 100, mean = 5, sd= 5)

# change in mean and variance
set.seed(6)
r6<- rnorm(n= 100, mean = 15, sd= 3)

exsim<- c(r1,r2,r3, r4, r5, r6)

plot(exsim)

exocpd<- onlineCPD(exsim, getR = TRUE, hazard_func = function(x, lambda){const_hazard(x, lambda=10)})

plot(exocpd, data=exsim, main_title = "Example Simulation")
