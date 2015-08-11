#For computer initial values of beta

mmmeth <-function(ti)
{
  S=function(ti)
  {
    n <- length(ti)
    return((1/n)*sum(ti))
  }

  R=function(ti)
  {
    n <- length(ti)
    return(((1/n)*sum(ti^(-1)))^(-1))
  }

  beta0ini  <- (S(ti)*R(ti))^0.5
  alpha0ini <- sqrt(2)*((S(ti)/R(ti))^0.5 - 1)^0.5

  result   <- list(beta0init = beta0ini,alpha0ini=alpha0ini, n = length(ti))
  return(result)
}

#mmmeth(ti)

#beta0<-mmmeth(ti)$beta0init
#alpha0<-mmmeth(ti)$alpha0ini
#shape0<-0
#delta0<-shape0/sqrt(1+shape0^2)
