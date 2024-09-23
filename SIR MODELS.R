#Write down the equations of a frequency-dependent SIR model with no births or
#deaths and transmission rate β and recovery rate γ

# dS/dt= −βSI/N

# dI/dt= βSI/N − γI

# dR/dt=γI

#Show that for this model the population size remains constant, i.e.dN/dt=0
#substitute each function with its parameters
#dN/dt=ds/dt+dI/dt+dR/dt
#     =(−βSI/N)+(βSI/N − γI)+(γI)
#     =0

#Coding the model
# SIR Epidemic model - Practical 1
library(deSolve)

SIR_dyn <- function(t,var,par) {
  # Rename the variables and parameters
  S <- var[1]
  I <- var[2]
  R <- var[3]
  N <- S+I+R
  beta <- par[1]
  gamma <- par[2]
  
  # Derivatives
  dS <- -beta*S*I/N
  dI <- beta*S*I/N-gamma*I
  dR <- gamma*I
  
  # Return the 3 values
  list(c(dS,dI,dR))
}

beta <- 1
gamma <- 0.25
SIR_par <- c(beta,gamma)
SIR_init <- c(99,1,0)
SIR_t <- seq(0,30,by=0.1)

# The numerical solution is given by
SIR_sol <- lsoda(SIR_init,
                 SIR_t,
                 SIR_dyn,
                 SIR_par)

#renaming variables
TIME <- SIR_sol[,1]
S <- SIR_sol[,2]
I <- SIR_sol[,3]
R <- SIR_sol[,4]
N <- S + I + R

#Plot S(t),I(t)and R(t)together on a graph, using the commands plot(), lines() and legend().

plot(S/N, I/N, type = "l", lwd = 3, col = "darkgreen")

#It can be difficult to read this graph at first. Find where TIME=1
#(at the start S/N≈1 and I≈0) and label it on the graph

index <- which(TIME == 1)

text(S[index]/N[index],I[index]/N[index],"t=1",pos=2)

#Also label t=2,t=5,t=10 &t=15. What is happening to time in this plot?

index <- which(TIME == 5)

text(S[index]/N[index],I[index]/N[index],"t=5",pos=2)


index <- which(TIME == 10)

text(S[index]/N[index],I[index]/N[index],"t=10",pos=2)


index <- which(TIME == 15)

text(S[index]/N[index],I[index]/N[index],"t=15",pos=2)

# set up a 2 by 3 grid for the plots:
par(mfrow=c(2,3), xaxs='i', yaxs='i')
# define the parameters:
infperiod <- c(1,1,28,5,7,14)
Rzero <- c(1.2,1.5,2,5,17,17)
for (i in 1:6) {
  gamma = signif(1/infperiod[i],2)
  beta = signif(Rzero[i]*gamma,2)
  SIR_par <- c(beta,gamma)
  SIR_sol<-lsoda(SIR_init,
                 SIR_t,
                 SIR_dyn,
                 SIR_par)
  plot(....., ylim=c(0,1),type="l", main=paste0("beta = ",beta,"gamma = ",gamma))
  lines(.....) # I(t)
  lines(.....) # R(t)
}
par(mfrow = c(1, 1))

# set up a 2 by 3 grid for the plots:
par(mfrow=c(2,3), xaxs='i', yaxs='i')
# define the parameters:
infperiod <- c(1,1,28,5,7,14)
Rzero <- c(1.2,1.5,2,5,17,17)
SIR.t =seq(0,500, by=0.1)
for (i in 1:6) {
  gamma = signif(1/infperiod[i],2)
  beta = signif(Rzero[i]*gamma,2)
  SIR_par <- c(beta,gamma)
  SIR.sol <- lsoda(SIR_init,SIR_t,SIR_dyn,SIR_par)
  TIME <- SIR.sol[,1]
  S <- SIR.sol[,2]
  I <- SIR.sol[,3]
  R <- SIR.sol[,4]
  N <- S + I + R
  plot(TIME,S/N,col='blue', type='l', ylim=c(0,1), 
       main=paste("beta = ",beta,"gamma = ",gamma))
  lines(TIME,I/N,col="red")
  lines(TIME,R/N,col="green")
  legend("right", col = c("blue","red","green"),legend = c("S","I","R"), lty = 2)
  print(paste("peak size = ",max(I/N)))
  print(paste("peak time = ",TIME[which.max(I/N)]) )
  print(paste("Epidemic size = ",max(R/N)))
}
par(mfrow = c(1, 1))


#Which of the last 3 quantities appear linked to R0?
#ans-Peak size and epidemic size.

#Plot log(I(t))for t in[0,4] for the same parameter combinations above 
#and check that the initial gradient is equal to β−γ

#checking initial gradient
par(mfrow=c(2,3))
# define the parameters:
infperiod <- c(1,1,28,5,7,14)
Rzero <- c(1.2,1.5,2,5,17,17)
SIR.t =seq(0,500, by=0.1)
for (i in 1:6) {
  gamma = signif(1/infperiod[i],2)
  beta = signif(Rzero[i]*gamma,2)
  SIR_par <- c(beta,gamma)
  SIR.sol <- lsoda(SIR_init,SIR_t,SIR_dyn,SIR_par)
  TIME <- SIR.sol[,1]
  S <- SIR.sol[,2]
  I <- SIR.sol[,3]
  R <- SIR.sol[,4]
  N <- S + I + R
  plot(TIME[1:41],log( I[1:41]/N[1:41] ),col='red', type='l', 
       main=paste("beta = ",beta,"gamma = ",gamma))
  abline(a = log( I[1]/N[1] ), b = beta-gamma, col = "black")
}
par(mfrow=c(1,1))

#What dynamics do you expect with R0=1? Plot the graph forβ=1,γ=1
#Why does I(t)decrease?

#Ans-As our initial conditions are with  99susceptible people, 
#and R0=1 this renders the effective reproduction number less than 1
#so I(t)decreases.

#Explore the behaviour of the system around R0≈1
#using the loop you wrote before. Try changing 0.9<R0<1.1
#and altering the initial number of individuals infected.

### looking around R0 = 1
par(mfrow=c(2,3))
# define the parameters:
SIR.init <- c(99,1,0)
infperiod <- rep(2,6)
Rzero <- c(0.9,0.95,1,1.05,1.1, 1.15)
SIR.t =seq(0,100, by=0.1)
for (i in 1:6) {
  gamma<- signif(1/infperiod[i],2)
  beta <- signif(Rzero[i]*gamma,2)
  SIR_par <- c(beta,gamma)
  SIR.sol <- lsoda(SIR_init,SIR_t,SIR_dyn,SIR_par)
  TIME <- SIR.sol[,1]
  S <- SIR.sol[,2]
  I <- SIR.sol[,3]
  R <- SIR.sol[,4]
  N <- S + I + R
  plot(TIME, I/N,col='red', type='l', ylim=c(0,0.03),
       main=paste("R0 = ",Rzero[i]))
}
par(mfrow=c(1,1))

#Calculating the infectious period

#If the recovery rate of a disease is γ=0.2days^−1
#what is the average infectious period?


# set up a vector of times:
x=seq(0,30,0.1)

# define the recovery rate, gamma:
gamma=0.2

# define the exponential function:
fx=exp(-gamma*x)

# plot:
plot(x,fx,type="l")

#calculatung the average infectious period using the fnx defined above
print(sum(fx*x)/sum(fx))


#How do your two answers compare?
#ans-Not quite the same but perhaps will have a better approximation as you use more points


#extending  SIR to other compartments


#Include loss of immunity, so that Recovered individuals become Susceptible again after a certain period of time.
#Include a latent state post-infection but before an individual becomes infectious (this type of model is often called an
#SEIRmodel, where  E stands for the Exposed class).
#Include treatment in the model, allowing Infected individuals to recover more quickly.
#Sketch a diagram of your model, showing the compartments and the flows between them.
#Write down the equations that govern your model
#Write a function, such as SEIR_dyn() that calculates the derivatives of your model for use with lsoda().
#Plot the epidemic curve, I(t). How does this compare to an equivalent  
#SIR model?

#### Extra exercise
# I am putting everything in at once

SEITRS.dyn <- function(t, var, par) {
  # Rename the variables and parameters
  S <- var[1]
  E <- var[2]
  I <- var[3]
  tr <- var[4]
  R <- var[5]
  N <- S + E + I + tr + R
  beta <- par[1]
  gamma <-par[2]
  epsilon <- par[3] #rate of latency loss
  tau <- par[4] #treatment rate
  omega <- par[5] #recovery rate for treated people
  
  # Derivatives
  dS <- -beta*S*I/N
  dE <- beta*S*I/N - epsilon*E
  dI <- epsilon*E - gamma*I - tau*I
  dT <- tau*I - omega*tr
  dR <- gamma*I + omega*tr
  # Return the 5 values
  list(c(dS, dE, dI, dT, dR))
}

beta <- 2
gamma <- 0.25
epsilon <- 1
tau <- 1
omega <- gamma*2

SEITRS.par <- c(beta,gamma,epsilon,tau,omega)
SEITRS.init <- c(99,1,0,0,0)
SEITRS.t <- seq(0,30,by=0.1)

SEITRS.sol <- lsoda(SEITRS.init,SEITRS.t,SEITRS.dyn,SEITRS.par)

TIME <- SEITRS.sol[,1]
S <- SEITRS.sol[,2]
E <- SEITRS.sol[,3]
I <- SEITRS.sol[,4]
tr <- SEITRS.sol[,5]
R <- SEITRS.sol[,6]

N <- S + E + I + tr + R

plot(TIME,S/N,col='blue', type='l', ylim=c(0,1))
lines(TIME,E/N,col="orange")
lines(TIME,I/N,col="red")
lines(TIME,tr/N,col="cyan")
lines(TIME,R/N,col="green")
legend("right", col = c("blue","orange", "red","cyan","green"), legend = c("S","E","I","T","R"), lty = 2)

#estimate the epidemic size for any value of R0

my_fun <- function(x,R0) {
  (x+log(1-x)/R0)ˆ2
}
optimize(my_fun, c(0,1), R0=2)
optimize(my_fun, c(0,1), R0=2)
