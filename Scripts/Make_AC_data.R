Make_AC_data<- function(N, sd, AC){
  numbers <- rep(NA, N)
  numbers[1] <- rnorm(1,0)  # Mean=0
  
  for(i in 2:N) { 
    numbers[i]<-numbers[i-1]* AC + rnorm(1,0)
  }
  return(  scale(numbers)  *sd)
}