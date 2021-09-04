library(quantmod)
library(data.table)
library(bizdays)
setwd("C:/Grad School/FE 620 Pricing and Hedging/Project")



########################################################################################################
#   PART 1: Data Gathering
########################################################################################################

#Read in Option Data from Bloomberg Excel File
Option.Data <- read.csv(file = "Option_Data.csv")
Option.Data



########################################################################################################
#   PART 2: Data Cleaning
########################################################################################################

# Load in previously generated datasets



#Short Term Interest Rates
# Risk free interest rate taken as 13 week T-Bill rate ^IRX on 4/7/20
int <- 0.128/100



#Calculate asset volatilities from 1 yr of Historical closing prices
#getSymbols(Symbols = c("^GSPC","UAL","AMZN", "BRK-B", "FB"), from = "2019-07-15", to = "2020-07-15")
getSymbols(Symbols = c("^GSPC","UAL","AMZN", "BRK-B", "FB"), from = "2019-07-15", to = "2020-04-15")
BRK <- `BRK-B`
SPX <- GSPC
colnames(SPX) <- c("SPX.Open", "SPX.High", "SPX.Low", "SPX.Close", "SPX.Volume", "SPX.Adjusted")


Historic <- as.data.frame(cbind(SPX$SPX.Adjusted, UAL$UAL.Adjusted, AMZN$AMZN.Adjusted, BRK$`BRK-B.Adjusted`, FB$FB.Adjusted))
Historic <- cbind(rownames(Historic), Historic)
colnames(Historic) <- c("Date", "SPX.Adjusted", "UAL.Adjusted", "AMZN.Adjusted", "BRK.B.Adjusted", "FB.Adjusted")
rownames(Historic) <- NULL

Historic$SPXrets <- c(NA,log(Historic$SPX.Adjusted[-1]/Historic$SPX.Adjusted[-length(Historic$SPX.Adjusted)]))
Historic$UALrets <- c(NA,log(Historic$UAL.Adjusted[-1]/Historic$UAL.Adjusted[-length(Historic$UAL.Adjusted)]))
Historic$AMZNrets <- c(NA,log(Historic$AMZN.Adjusted[-1]/Historic$AMZN.Adjusted[-length(Historic$AMZN.Adjusted)]))
Historic$BRKrets <- c(NA,log(Historic$BRK.B.Adjusted[-1]/Historic$BRK.B.Adjusted[-length(Historic$BRK.B.Adjusted)]))
Historic$FBrets <- c(NA,log(Historic$FB.Adjusted[-1]/Historic$FB.Adjusted[-length(Historic$FB.Adjusted)]))

SPX.vol <- sd(Historic$SPXrets, na.rm = TRUE)/sqrt(1/252)
UAL.vol <- sd(Historic$UALrets, na.rm = TRUE)/sqrt(1/252)
AMZN.vol <- sd(Historic$AMZNrets, na.rm = TRUE)/sqrt(1/252)
BRK.vol <- sd(Historic$BRKrets, na.rm = TRUE)/sqrt(1/252)
FB.vol <- sd(Historic$FBrets, na.rm = TRUE)/sqrt(1/252)



# Add Current Stock price (at time of option data) and volatility to datasets
for(i in 1:nrow(Option.Data)){
  if(Option.Data[i, "Underlying"] == "AMZN"){
      Option.Data[i, "Volatility"] = AMZN.vol
      temp = Historic[Historic$Date == as.character(as.Date(Option.Data[i, "Price.Date"], format = "%m/%d/%Y")), "AMZN.Adjusted"]
      Option.Data[i, "Stock.Price"] = temp
      
  }else if(Option.Data[i, "Underlying"] == "SPX"){
    Option.Data[i, "Volatility"] = SPX.vol
    temp = Historic[Historic$Date == as.character(as.Date(Option.Data[i, "Price.Date"], format = "%m/%d/%Y")), "SPX.Adjusted"]
    Option.Data[i, "Stock.Price"] = temp
    
  }else if(Option.Data[i, "Underlying"] == "UAL"){
    Option.Data[i, "Volatility"] = UAL.vol
    temp = Historic[Historic$Date == as.character(as.Date(Option.Data[i, "Price.Date"], format = "%m/%d/%Y")), "UAL.Adjusted"]
    Option.Data[i, "Stock.Price"] = temp
    
  }else if(Option.Data[i, "Underlying"] == "BRK"){
    Option.Data[i, "Volatility"] = BRK.vol
    temp = Historic[Historic$Date == as.character(as.Date(Option.Data[i, "Price.Date"], format = "%m/%d/%Y")), "BRK.B.Adjusted"]
    Option.Data[i, "Stock.Price"] = temp
    
  }else{
    Option.Data[i, "Volatility"] = FB.vol
    temp = Historic[Historic$Date == as.character(as.Date(Option.Data[i, "Price.Date"], format = "%m/%d/%Y")), "FB.Adjusted"]
    Option.Data[i, "Stock.Price"] = temp
  }
  
}



#Add Time to Maturity for each option
working_calendar <-  create.calendar(name = "mycal", weekdays=c("saturday", "sunday"))
Option.Data$TTM = bizdays(from = as.Date(Option.Data[,"Price.Date"], format = "%m/%d/%Y"), to = as.Date(Option.Data[,"Maturity"], format = "%m/%d/%Y"), working_calendar)/252



########################################################################################################
#   PART 3: Definition of Option Pricing Models
########################################################################################################

#Define a function to output an additive Binomial Tree of potential X's (log returns)
#The change up / down was held constant. Xd = Xu (Trigeorgis Tree Construction)
Binom.Tree <- function(S0, K, TTM, r, d, sigma, N){
  del.t <- TTM/N
  del.x <- sqrt((r-d-sigma^2/2)^2*del.t^2+sigma^2*del.t)
  pu <- 1/2 + 1/2*(r-d-sigma^2/2)*del.t/del.x
  pd <- 1 - pu
  
  Xtree <- list()
  
  for(i in 1:N){
    jvec <- vector("numeric", length = i+1)
    for(j in 0:i){
      jvec[j+1] <- log(S0) + (i-j)*del.x - j*del.x
    }
    Xtree[[i]] <- jvec 
  }
  
  return(Xtree)
}
#Function Testing
#Binom.Tree (S0 = 30, K = 30, TTM = 1, r = 5/100, d=0, sigma = 0.25, N = 5)


#Define Function for pricing American Calls via the Binomial Tree
American.Call <- function(S0, K, TTM, r, d, sigma, N){
  
  Xtree <- Binom.Tree (S0=S0, K=K, TTM=TTM, r=r, d=d, sigma=sigma, N=N)
  del.t <- TTM/N
  del.x <- sqrt((r-d-sigma^2/2)^2*del.t^2+sigma^2*del.t)
  pu <- 1/2 + 1/2*(r-d-sigma^2/2)*del.t/del.x
  pd <- 1 - pu
  
  temp <- Xtree[[length(Xtree)]]
  terminal <- exp(temp)-K
  terminal <- mapply(FUN = max, terminal, 0)
  
  #Note loop terminates at time 2, due to lookback from immediate excersise condition
  for(i in length(Xtree):2){
    temp1 <- vector("numeric", length = i)
    temp2 <- vector("numeric", length = i)
    temp <- vector("numeric", length = i)
    for(j in 1:i){
      
      #temp1 vector for expected option value working backwards
      temp1[j] <- exp(-r*del.t)*(pu*terminal[j]+pd*terminal[j+1])
      
      #temp2 vector for immediate excersise condition
      temp2[j] <- exp(Xtree[[i-1]][j])-K
      
      #temp vector to do the optimum of holding option or excersising at that time
      temp <- mapply(FUN = max, temp1, temp2)
    }
    terminal <- temp
  }
  
  #Complete last iteration
  temp1 <- vector("numeric")
  temp2 <- vector("numeric")
  temp <- vector("numeric")
  
  temp1 <- exp(-r*del.t)*(pu*terminal[1]+pd*terminal[2])
  temp2 <- S0-K
  option.price <- mapply(FUN = max, temp1, temp2)
  
  return(option.price)
}

#Testing Function
#American.Call(S0 = 30, K = 30, TTM = 1, r = 5/100, d=0, sigma = 0.25, N = 100)


#Define Function for pricing European Calls via the Binomial Tree
European.Call <- function(S0, K, TTM, r, d, sigma, N){
  
  Xtree <- Binom.Tree (S0=S0, K=K, TTM=TTM, r=r, d=d, sigma=sigma, N=N)
  del.t <- TTM/N
  del.x <- sqrt((r-d-sigma^2/2)^2*del.t^2+sigma^2*del.t)
  pu <- 1/2 + 1/2*(r-d-sigma^2/2)*del.t/del.x
  pd <- 1 - pu
  
  temp <- Xtree[[length(Xtree)]]
  terminal <- exp(temp)-K
  terminal <- mapply(FUN = max, terminal, 0)
  
  for(i in length(Xtree):1){
    temp <- vector("numeric", length = i)
    for(j in 1:i){
      temp[j] <- exp(-r*del.t)*(pu*terminal[j]+pd*terminal[j+1])
    }
    
    terminal <- temp
  }
  
  option.price <- terminal
  return(option.price)
}
#Testing Function
#European.Call(S0 = 30, K = 30, TTM = 1, r = 5/100, d=0, sigma = 0.25, N = 100)


#Define Function for pricing American Puts via the Binomial Tree
American.Put <- function(S0, K, TTM, r, d, sigma, N){
  
  Xtree <- Binom.Tree (S0=S0, K=K, TTM=TTM, r=r, d=d, sigma=sigma, N=N)
  del.t <- TTM/N
  del.x <- sqrt((r-d-sigma^2/2)^2*del.t^2+sigma^2*del.t)
  pu <- 1/2 + 1/2*(r-d-sigma^2/2)*del.t/del.x
  pd <- 1 - pu
  
  temp <- Xtree[[length(Xtree)]]
  terminal <- K-exp(temp)
  terminal <- mapply(FUN = max, terminal, 0)
  
  #Note loop terminates at time 2, due to lookback from immediate excersise condition
  for(i in length(Xtree):2){
    temp1 <- vector("numeric", length = i)
    temp2 <- vector("numeric", length = i)
    temp <- vector("numeric", length = i)
    for(j in 1:i){
      
      #temp1 vector for expected option value working backwards
      temp1[j] <- exp(-r*del.t)*(pu*terminal[j]+pd*terminal[j+1])
      
      #temp2 vector for immediate excersise condition
      temp2[j] <- K-exp(Xtree[[i-1]][j])
      
      #temp vector to do the optimum of holding option or excersising at that time
      temp <- mapply(FUN = max, temp1, temp2)
    }
    terminal <- temp
  }
  
  #Complete last iteration
  temp1 <- vector("numeric")
  temp2 <- vector("numeric")
  temp <- vector("numeric")
  
  temp1 <- exp(-r*del.t)*(pu*terminal[1]+pd*terminal[2])
  temp2 <- K-S0
  option.price <- mapply(FUN = max, temp1, temp2)
  
  return(option.price)
}
#Testing Function
#American.Put(S0 = 30, K = 30, TTM = 1, r = 5/100, d=0, sigma = 0.25, N = 100)


#Define Function for pricing European Puts via the Binomial Tree
European.Put <- function(S0, K, TTM, r, d, sigma, N){
  
  Xtree <- Binom.Tree (S0=S0, K=K, TTM=TTM, r=r, d=d, sigma=sigma, N=N)
  del.t <- TTM/N
  del.x <- sqrt((r-d-sigma^2/2)^2*del.t^2+sigma^2*del.t)
  pu <- 1/2 + 1/2*(r-d-sigma^2/2)*del.t/del.x
  pd <- 1 - pu
  
  temp <- Xtree[[length(Xtree)]]
  terminal <- K-exp(temp)
  terminal <- mapply(FUN = max, terminal, 0)
  
  for(i in length(Xtree):1){
    temp <- vector("numeric", length = i)
    for(j in 1:i){
      temp[j] <- exp(-r*del.t)*(pu*terminal[j]+pd*terminal[j+1])
    }
    
    terminal <- temp
  }
  
  option.price <- terminal
  return(option.price)
}
#Testing Function
#European.Put(S0 = 30, K = 30, TTM = 1, r = 5/100, d=0, sigma = 0.25, N = 100)


#Define a function to output an additive Trinomial Tree of potential X's (log returns)
#The change up / down was held constant. Xd = Xu (Trigeorgis Tree Construction)
Trinom.Tree <- function(S0, K, TTM, r, d, sigma, N){
  del.t <- TTM/N
  del.x <- sigma*sqrt(3*del.t)
  D <- r-d-sigma^2/2
  pu <- 1/2*((sigma^2*del.t+D^2*del.t^2)/(del.x^2)+D*del.t/del.x)
  pm <- 1 - (sigma^2*del.t+D^2*del.t^2)/(del.x^2)
  pd <- 1/2*((sigma^2*del.t+D^2*del.t^2)/(del.x^2)-D*del.t/del.x)
  
  Xtree <- list()
  
  for(i in 1:N){
    jvec <- vector("numeric", length = 2*i+1)
    for(j in -i:i){
      jvec[j+i+1] <- log(S0) - j*del.x
    }
    Xtree[[i]] <- jvec 
  }
  
  return(Xtree)
}
#Function Testing
#Trinom.Tree (S0 = 30, K = 30, TTM = 1, r = 5/100, d=0, sigma = 0.25, N = 5)


#Define Function for pricing American Calls via the Trinomial Tree
Tri.American.Call <- function(S0, K, TTM, r, d, sigma, N){
  
  Xtree <- Trinom.Tree (S0=S0, K=K, TTM=TTM, r=r, d=d, sigma=sigma, N=N)
  del.t <- TTM/N
  del.x <- sigma*sqrt(3*del.t)
  D <- r-d-sigma^2/2
  pu <- 1/2*((sigma^2*del.t+D^2*del.t^2)/(del.x^2)+D*del.t/del.x)
  pm <- 1 - (sigma^2*del.t+D^2*del.t^2)/(del.x^2)
  pd <- 1/2*((sigma^2*del.t+D^2*del.t^2)/(del.x^2)-D*del.t/del.x)
  
  temp <- Xtree[[length(Xtree)]]
  terminal <- exp(temp)-K
  terminal <- mapply(FUN = max, terminal, 0)
  
  #Note loop terminates at time 2, due to lookback from immediate excersise condition
  for(i in length(Xtree):2){
    temp1 <- vector("numeric", length = 2*i-1)
    temp2 <- vector("numeric", length = 2*i-1)
    temp <- vector("numeric", length = 2*i-1)
    for(j in 1:(2*i-1)){
      
      #temp1 vector for expected option value working backwards
      temp1[j] <- exp(-r*del.t)*(pu*terminal[j]+pm*terminal[j+1]+pd*terminal[j+2])
      
      #temp2 vector for immediate exercise condition
      temp2[j] <- exp(Xtree[[i-1]][j])-K
      
      #temp vector to do the optimum of holding option or excersising at that time
      temp <- mapply(FUN = max, temp1, temp2)
    }
    terminal <- temp
  }
  
  #Complete last iteration
  temp1 <- vector("numeric")
  temp2 <- vector("numeric")
  temp <- vector("numeric")
  
  temp1 <- exp(-r*del.t)*(pu*terminal[1]+pm*terminal[2]+pd*terminal[3])
  temp2 <- S0-K
  option.price <- mapply(FUN = max, temp1, temp2)
  
  return(option.price)
}
#Testing Function
#Tri.American.Call(S0 = 30, K = 30, TTM = 1, r = 5/100, d=0, sigma = 0.25, N = 100)



#Define Function for pricing European Calls via the Trinomial Tree
Tri.European.Call <- function(S0, K, TTM, r, d, sigma, N){
  
  Xtree <- Trinom.Tree (S0=S0, K=K, TTM=TTM, r=r, d=d, sigma=sigma, N=N)
  del.t <- TTM/N
  del.x <- sigma*sqrt(3*del.t)
  D <- r-d-sigma^2/2
  pu <- 1/2*((sigma^2*del.t+D^2*del.t^2)/(del.x^2)+D*del.t/del.x)
  pm <- 1 - (sigma^2*del.t+D^2*del.t^2)/(del.x^2)
  pd <- 1/2*((sigma^2*del.t+D^2*del.t^2)/(del.x^2)-D*del.t/del.x)
  
  temp <- Xtree[[length(Xtree)]]
  terminal <- exp(temp)-K
  terminal <- mapply(FUN = max, terminal, 0)
  
  for(i in length(Xtree):1){
    temp <- vector("numeric", length = 2*i-1)
    for(j in 1:(2*i-1)){
      temp[j] <- exp(-r*del.t)*(pu*terminal[j]+pm*terminal[j+1]+pd*terminal[j+2])
    }
    
    terminal <- temp
  }
  
  option.price <- terminal
  return(option.price)
}
#Function Testing
#Tri.European.Call(S0 = 30, K = 30, TTM = 1, r = 5/100, d=0, sigma = 0.25, N = 100)



#Define Function for pricing American Puts via the Trinomial Tree
Tri.American.Put <- function(S0, K, TTM, r, d, sigma, N){
  
  Xtree <- Trinom.Tree (S0=S0, K=K, TTM=TTM, r=r, d=d, sigma=sigma, N=N)
  del.t <- TTM/N
  del.x <- sigma*sqrt(3*del.t)
  D <- r-d-sigma^2/2
  pu <- 1/2*((sigma^2*del.t+D^2*del.t^2)/(del.x^2)+D*del.t/del.x)
  pm <- 1 - (sigma^2*del.t+D^2*del.t^2)/(del.x^2)
  pd <- 1/2*((sigma^2*del.t+D^2*del.t^2)/(del.x^2)-D*del.t/del.x)
  
  temp <- Xtree[[length(Xtree)]]
  terminal <- K-exp(temp)
  terminal <- mapply(FUN = max, terminal, 0)
  
  #Note loop terminates at time 2, due to lookback from immediate excersise condition
  for(i in length(Xtree):2){
    temp1 <- vector("numeric", length = 2*i-1)
    temp2 <- vector("numeric", length = 2*i-1)
    temp <- vector("numeric", length = 2*i-1)
    for(j in 1:(2*i-1)){
      
      #temp1 vector for expected option value working backwards
      temp1[j] <- exp(-r*del.t)*(pu*terminal[j]+pm*terminal[j+1]+pd*terminal[j+2])
      
      #temp2 vector for immediate exercise condition
      temp2[j] <- K-exp(Xtree[[i-1]][j])
      
      #temp vector to do the optimum of holding option or excersising at that time
      temp <- mapply(FUN = max, temp1, temp2)
    }
    terminal <- temp
  }
  
  #Complete last iteration
  temp1 <- vector("numeric")
  temp2 <- vector("numeric")
  temp <- vector("numeric")
  
  temp1 <- exp(-r*del.t)*(pu*terminal[1]+pm*terminal[2]+pd*terminal[3])
  temp2 <- K-S0
  option.price <- mapply(FUN = max, temp1, temp2)
  
  return(option.price)
}
#Testing Function
#Tri.American.Put(S0 = 30, K = 30, TTM = 1, r = 5/100, d=0, sigma = 0.25, N = 100)


#Define Function for pricing European Puts via the Trinomial Tree
Tri.European.Put <- function(S0, K, TTM, r, d, sigma, N){
  
  Xtree <- Trinom.Tree (S0=S0, K=K, TTM=TTM, r=r, d=d, sigma=sigma, N=N)
  del.t <- TTM/N
  del.x <- sigma*sqrt(3*del.t)
  D <- r-d-sigma^2/2
  pu <- 1/2*((sigma^2*del.t+D^2*del.t^2)/(del.x^2)+D*del.t/del.x)
  pm <- 1 - (sigma^2*del.t+D^2*del.t^2)/(del.x^2)
  pd <- 1/2*((sigma^2*del.t+D^2*del.t^2)/(del.x^2)-D*del.t/del.x)
  
  temp <- Xtree[[length(Xtree)]]
  terminal <- K-exp(temp)
  terminal <- mapply(FUN = max, terminal, 0)
  
  for(i in length(Xtree):1){
    temp <- vector("numeric", length = 2*i-1)
    for(j in 1:(2*i-1)){
      temp[j] <- exp(-r*del.t)*(pu*terminal[j]+pm*terminal[j+1]+pd*terminal[j+2])
    }
    
    terminal <- temp
  }
  
  option.price <- terminal
  return(option.price)
}
#Function Testing
#Tri.European.Put(S0 = 30, K = 30, TTM = 1, r = 5/100, d=0, sigma = 0.25, N = 100)



########################################################################################################
#   PART 4: Run Pricing models on real option data
########################################################################################################

# First Look at Convergence of Binomial and Trinomial Trees on 3Mo option to determine number of tree nodes
Nvec <- c(1:20*5)
Sens.Binomial <- mapply(FUN = American.Call,
                        S0 = 100,
                        K = 100,
                        TTM = 0.25,
                        r = int,
                        d = 0,
                        sigma = 0.3,
                        N = Nvec
)

Sens.Trinomial <- mapply(FUN = Tri.American.Call,
                        S0 = 100,
                        K = 100,
                        TTM = 0.25,
                        r = int,
                        d = 0,
                        sigma = 0.3,
                        N = Nvec
)

Sens.Binomial.Put <- mapply(FUN = American.Put,
                        S0 = 100,
                        K = 100,
                        TTM = 0.25,
                        r = int,
                        d = 0,
                        sigma = 0.3,
                        N = Nvec
)

Sens.Trinomial.Put <- mapply(FUN = Tri.American.Put,
                         S0 = 100,
                         K = 100,
                         TTM = 0.25,
                         r = int,
                         d = 0,
                         sigma = 0.3,
                         N = Nvec
)

Sensitivity <- cbind(Nvec, Sens.Binomial, Sens.Trinomial, Sens.Binomial.Put, Sens.Trinomial.Put)
#write.csv(Sensitivity, "Sensitivity.csv")

#Price American Call Options under Binomial Tree
Option.Data$Binomial[Option.Data$Type=="Call" & Option.Data$Underlying!="SPX"] <- mapply(FUN = American.Call,
                                                         S0 = Option.Data$Stock.Price[Option.Data$Type=="Call" & Option.Data$Underlying!="SPX"],
                                                         K = Option.Data$Strike[Option.Data$Type=="Call" & Option.Data$Underlying!="SPX"],
                                                         TTM = Option.Data$TTM[Option.Data$Type=="Call" & Option.Data$Underlying!="SPX"],
                                                         r = int,
                                                         d = 0,
                                                         sigma = Option.Data$Volatility[Option.Data$Type=="Call" & Option.Data$Underlying!="SPX"],
                                                         N = 50
                                                          )
                                                         

#Price American Put Options under Binomial Tree
Option.Data$Binomial[Option.Data$Type=="Put" & Option.Data$Underlying!="SPX"] <- mapply(FUN = American.Put,
                                                         S0 = Option.Data$Stock.Price[Option.Data$Type=="Put" & Option.Data$Underlying!="SPX"],
                                                         K = Option.Data$Strike[Option.Data$Type=="Put" & Option.Data$Underlying!="SPX"],
                                                         TTM = Option.Data$TTM[Option.Data$Type=="Put" & Option.Data$Underlying!="SPX"],
                                                         r = int,
                                                         d = 0,
                                                         sigma = Option.Data$Volatility[Option.Data$Type=="Put" & Option.Data$Underlying!="SPX"],
                                                         N = 50
)

#Price American Call Options under Trinomial Tree
Option.Data$Trinomial[Option.Data$Type=="Call" & Option.Data$Underlying!="SPX"]<- mapply(FUN = Tri.American.Call,
                                                         S0 = Option.Data$Stock.Price[Option.Data$Type=="Call" & Option.Data$Underlying!="SPX"],
                                                         K = Option.Data$Strike[Option.Data$Type=="Call" & Option.Data$Underlying!="SPX"],
                                                         TTM = Option.Data$TTM[Option.Data$Type=="Call" & Option.Data$Underlying!="SPX"],
                                                         r = int,
                                                         d = 0,
                                                         sigma = Option.Data$Volatility[Option.Data$Type=="Call" & Option.Data$Underlying!="SPX"],
                                                         N = 50
)

#Price American Put Options under Binomial Tree
Option.Data$Trinomial[Option.Data$Type=="Put" & Option.Data$Underlying!="SPX"] <- mapply(FUN = Tri.American.Put,
                                                        S0 = Option.Data$Stock.Price[Option.Data$Type=="Put" & Option.Data$Underlying!="SPX"],
                                                        K = Option.Data$Strike[Option.Data$Type=="Put" & Option.Data$Underlying!="SPX"],
                                                        TTM = Option.Data$TTM[Option.Data$Type=="Put" & Option.Data$Underlying!="SPX"],
                                                        r = int,
                                                        d = 0,
                                                        sigma = Option.Data$Volatility[Option.Data$Type=="Put" & Option.Data$Underlying!="SPX"],
                                                        N = 50
)


###################################################################################
##### Price SPX As European

#Price European Call Options under Binomial Tree
Option.Data$Binomial[Option.Data$Type=="Call" & Option.Data$Underlying=="SPX"] <- mapply(FUN = European.Call,
                                                       S0 = Option.Data$Stock.Price[Option.Data$Type=="Call" & Option.Data$Underlying=="SPX"],
                                                       K = Option.Data$Strike[Option.Data$Type=="Call" & Option.Data$Underlying=="SPX"],
                                                       TTM = Option.Data$TTM[Option.Data$Type=="Call" & Option.Data$Underlying=="SPX"],
                                                       r = int,
                                                       d = 0,
                                                       sigma = Option.Data$Volatility[Option.Data$Type=="Call" & Option.Data$Underlying=="SPX"],
                                                       N = 50
)


#Price American Put Options under Binomial Tree
Option.Data$Binomial[Option.Data$Type=="Put" & Option.Data$Underlying=="SPX"] <- mapply(FUN = European.Put,
                                                                                          S0 = Option.Data$Stock.Price[Option.Data$Type=="Put" & Option.Data$Underlying=="SPX"],
                                                                                          K = Option.Data$Strike[Option.Data$Type=="Put" & Option.Data$Underlying=="SPX"],
                                                                                          TTM = Option.Data$TTM[Option.Data$Type=="Put" & Option.Data$Underlying=="SPX"],
                                                                                          r = int,
                                                                                          d = 0,
                                                                                          sigma = Option.Data$Volatility[Option.Data$Type=="Put" & Option.Data$Underlying=="SPX"],
                                                                                          N = 50
)

#Price American Call Options under Trinomial Tree
Option.Data$Trinomial[Option.Data$Type=="Call" & Option.Data$Underlying=="SPX"]<- mapply(FUN = Tri.European.Call,
                                                                                           S0 = Option.Data$Stock.Price[Option.Data$Type=="Call" & Option.Data$Underlying=="SPX"],
                                                                                           K = Option.Data$Strike[Option.Data$Type=="Call" & Option.Data$Underlying=="SPX"],
                                                                                           TTM = Option.Data$TTM[Option.Data$Type=="Call" & Option.Data$Underlying=="SPX"],
                                                                                           r = int,
                                                                                           d = 0,
                                                                                           sigma = Option.Data$Volatility[Option.Data$Type=="Call" & Option.Data$Underlying=="SPX"],
                                                                                           N = 50
)

#Price American Put Options under Binomial Tree
Option.Data$Trinomial[Option.Data$Type=="Put" & Option.Data$Underlying=="SPX"] <- mapply(FUN = Tri.European.Put,
                                                                                           S0 = Option.Data$Stock.Price[Option.Data$Type=="Put" & Option.Data$Underlying=="SPX"],
                                                                                           K = Option.Data$Strike[Option.Data$Type=="Put" & Option.Data$Underlying=="SPX"],
                                                                                           TTM = Option.Data$TTM[Option.Data$Type=="Put" & Option.Data$Underlying=="SPX"],
                                                                                           r = int,
                                                                                           d = 0,
                                                                                           sigma = Option.Data$Volatility[Option.Data$Type=="Put" & Option.Data$Underlying=="SPX"],
                                                                                           N = 50
)





#Final Data Table
write.csv(Option.Data, "Results_Output.csv")

#Test Using Bloomberg IVM
#Option.Data$Test.Option.Value[Option.Data$Type=="Call"] <- mapply(FUN = Tri.American.Call,
#                                                        S0 = Option.Data$Stock.Price[Option.Data$Type=="Call"],
#                                                        K = Option.Data$Strike[Option.Data$Type=="Call"],
#                                                        TTM = Option.Data$TTM[Option.Data$Type=="Call"],
#                                                        r = int,
#                                                        d = 0,
#                                                        sigma = Option.Data$IVM[Option.Data$Type=="Call"]/100,
#                                                        N = 50
#)


#Option.Data$Test.Option.Value[Option.Data$Type=="Put"] <- mapply(FUN = Tri.American.Put,
#                                                         S0 = Option.Data$Stock.Price[Option.Data$Type=="Put"],
#                                                         K = Option.Data$Strike[Option.Data$Type=="Put"],
#                                                         TTM = Option.Data$TTM[Option.Data$Type=="Put"],
#                                                         r = int,
#                                                         d = 0,
#                                                         sigma = Option.Data$IVM[Option.Data$Type=="Put"]/100,
#                                                         N = 50
#)





########################################################################################################
#   2.) Testing Functions against HULL Example
########################################################################################################
#Test the pricer on the Example 21.1 in Hull. The tree should reproduce the numbers
#in the tree on page 21.3 for an American put option with parameters
#S0 = K = 50 , ?? = 0.4 , r = 0.1 , q = 0

#Defining a modification on Binomial Tree function to return the tree instead of the option price, to match HULL example
American.Put.Tree <- function(S0, K, TTM, r, d, sigma, N){
  
  Xtree <- Binom.Tree (S0=S0, K=K, TTM=TTM, r=r, d=d, sigma=sigma, N=N)
  del.t <- TTM/N
  del.x <- sqrt((r-d-sigma^2/2)^2*del.t^2+sigma^2*del.t)
  pu <- 1/2 + 1/2*(r-d-sigma^2/2)*del.t/del.x
  pd <- 1 - pu
  
  Value_tree <- list()
  
  temp <- Xtree[[length(Xtree)]]
  terminal <- K-exp(temp)
  terminal <- mapply(FUN = max, terminal, 0)
  Value_tree[[length(Xtree)+1]] <- terminal
  
  #Note loop terminates at time 2, due to lookback from immediate excersise condition
  for(i in length(Xtree):2){
    temp1 <- vector("numeric", length = i)
    temp2 <- vector("numeric", length = i)
    temp <- vector("numeric", length = i)
    for(j in 1:i){
      
      #temp1 vector for expected option value working backwards
      temp1[j] <- exp(-r*del.t)*(pu*terminal[j]+pd*terminal[j+1])
      
      #temp2 vector for immediate excersise condition
      temp2[j] <- K-exp(Xtree[[i-1]][j])
      
      #temp vector to do the optimum of holding option or excersising at that time
      temp <- mapply(FUN = max, temp1, temp2)
    }
    terminal <- temp
    Value_tree[[i]] <- terminal
  }
  
  #Complete last iteration
  temp1 <- vector("numeric")
  temp2 <- vector("numeric")
  temp <- vector("numeric")
  
  temp1 <- exp(-r*del.t)*(pu*terminal[1]+pd*terminal[2])
  temp2 <- K-S0
  option.price <- mapply(FUN = max, temp1, temp2)
  Value_tree[[1]] <- option.price
  
  return(Value_tree)
}


# Display Stock Price Tree from Figure 21.3
lapply(FUN = exp, Binom.Tree(S0=50, K=50, TTM=5/12, r=0.1, d=0, sigma = 0.4, N=5))

# Display Option Value Tree from Figure 21.3
American.Put.Tree(S0 = 50, K = 50, TTM = 5/12, r = 0.1, d=0, sigma = 0.4, N = 5)

#Run test at 30, 50, 100, 500 timesteps, to match HULL

# N=30
American.Put(S0 = 50, K = 50, TTM = 5/12, r = 0.1, d=0, sigma = 0.4, N = 30)

# N=50
American.Put(S0 = 50, K = 50, TTM = 5/12, r = 0.1, d=0, sigma = 0.4, N = 50)

# N=100
American.Put(S0 = 50, K = 50, TTM = 5/12, r = 0.1, d=0, sigma = 0.4, N = 100)

# N=500
American.Put(S0 = 50, K = 50, TTM = 5/12, r = 0.1, d=0, sigma = 0.4, N = 500)


#Show that as n-> infinity, the price of an American Call approaches the BSM value

#Define Analytic BSM functions for testing results
BSM_call.div <- function(S0, K, TTM, r, d, sigma){
  dplus <- 1/sigma/sqrt(TTM)*(log(S0/K)+(r-d+sigma^2/2)*TTM)
  dmin <- dplus - sigma*sqrt(TTM)
  calc.price <- S0*exp(-d*TTM)*pnorm(q=dplus)-K*exp(-r*TTM)*pnorm(q=dmin)
  return(calc.price)
}

#Testing with N = 5, 10, 50, 100, 500 
nvec <- c(5, 10, 50, 100, 500)
am_call_limit <- mapply(FUN = American.Call, S0 = 50, K = 50, TTM = 5/12, r = 0.1, d=0, sigma = 0.4, N = nvec)
bsm_call_results <- BSM_call.div(S0=50, K=50, TTM=5/12, r=0.1, d=0, sigma=0.4)
bsm_call_results <- rep(bsm_call_results, 5)

limit.df <- cbind(nvec, am_call_limit, bsm_call_results)
colnames(limit.df) <- c("N", "Binom Tree Call", "BSM Call")
limit.df



########################################################################################################
#  Compute Greeks
########################################################################################################

#Define a Function for Computing Delta, Calling the Binomial Tree
Binom_Delta <- function(S0, K, TTM, r, d, sigma, N, type){
  if (type == "Call"){
    Y1 = American.Call(S0=S0, K=K, TTM=TTM, r=r, d=d, sigma=sigma, N=N)
    Y2 = American.Call(S0=S0+0.1, K=K, TTM=TTM, r=r, d=d, sigma=sigma, N=N)
  } else{
    Y1 = American.Put(S0=S0, K=K, TTM=TTM, r=r, d=d, sigma=sigma, N=N)
    Y2 = American.Put(S0=S0+0.1, K=K, TTM=TTM, r=r, d=d, sigma=sigma, N=N)
  }
  Delta = (Y2-Y1)/0.1
  
  #Multiplying by 100% to be consistent with Bloomberg Reporting
  Delta = Delta*100
  return(Delta)
}
# Function Testing
#Binom_Delta(S0=100, K=100, TTM=1, r=0.05, d=0, sigma=0.2, N=50, type="Call")


#Define a Function for Computing Gamma, Calling the Binomial Tree
Binom_Gamma <- function(S0, K, TTM, r, d, sigma, N, type){
  if (type == "Call"){
    Y1 = American.Call(S0=S0, K=K, TTM=TTM, r=r, d=d, sigma=sigma, N=N)
    Y2 = American.Call(S0=S0+0.1, K=K, TTM=TTM, r=r, d=d, sigma=sigma, N=N)
    Y3 = American.Call(S0=S0+0.2, K=K, TTM=TTM, r=r, d=d, sigma=sigma, N=N)
  } else{
    Y1 = American.Put(S0=S0, K=K, TTM=TTM, r=r, d=d, sigma=sigma, N=N)
    Y2 = American.Put(S0=S0+0.1, K=K, TTM=TTM, r=r, d=d, sigma=sigma, N=N)
    Y3 = American.Put(S0=S0+0.2, K=K, TTM=TTM, r=r, d=d, sigma=sigma, N=N)
  }
  Delta1 = (Y2-Y1)/0.1
  Delta2 = (Y3-Y2)/0.1
  Gamma = (Delta2 - Delta1)/0.1
  
  #Multiplying by 100% to be consistent with Bloomberg Reporting
  Gamma = Gamma*100
  
  return(Gamma)
}
# Function Testing
#Binom_Gamma(S0=100, K=100, TTM=1, r=0.05, d=0, sigma=0.2, N=50, type="Put")



#Define a Function for Computing Vega, Calling the Binomial Tree
Binom_Vega <- function(S0, K, TTM, r, d, sigma, N, type){
  if (type == "Call"){
    Y1 = American.Call(S0=S0, K=K, TTM=TTM, r=r, d=d, sigma=sigma, N=N)
    Y2 = American.Call(S0=S0, K=K, TTM=TTM, r=r, d=d, sigma=sigma+0.005, N=N)
  } else{
    Y1 = American.Put(S0=S0, K=K, TTM=TTM, r=r, d=d, sigma=sigma, N=N)
    Y2 = American.Put(S0=S0, K=K, TTM=TTM, r=r, d=d, sigma=sigma+0.005, N=N)
  }
  Vega = (Y2-Y1)/0.005
  
  # Dividing by 100 % to be consistent with Bloomberg Reporting
  Vega = Vega/100
  return(Vega)
}
# Function Testing
#Binom_Vega(S0=100, K=100, TTM=1, r=0.05, d=0, sigma=0.2, N=50, type="Put")



#Define a Function for Computing Theta, Calling the Binomial Tree
Binom_Theta <- function(S0, K, TTM, r, d, sigma, N, type){
  if (type == "Call"){
    Y1 = American.Call(S0=S0, K=K, TTM=TTM, r=r, d=d, sigma=sigma, N=N)
    Y2 = American.Call(S0=S0, K=K, TTM=TTM-0.001, r=r, d=d, sigma=sigma, N=N)
  } else{
    Y1 = American.Put(S0=S0, K=K, TTM=TTM, r=r, d=d, sigma=sigma, N=N)
    Y2 = American.Put(S0=S0, K=K, TTM=TTM-0.001, r=r, d=d, sigma=sigma, N=N)
  }
  Theta = (Y2-Y1)/0.001
  
  #Diving by 365 days to be consistent with Bloomberg Reporting
  Theta = Theta/365
  return(Theta)
}
# Function Testing
#Binom_Theta(S0=100, K=100, TTM=1, r=0.05, d=0, sigma=0.2, N=50, type="Call")


#Define a Function for Computing Rho, Calling the Binomial Tree
Binom_Rho <- function(S0, K, TTM, r, d, sigma, N, type){
  if (type == "Call"){
    Y1 = American.Call(S0=S0, K=K, TTM=TTM, r=r, d=d, sigma=sigma, N=N)
    Y2 = American.Call(S0=S0, K=K, TTM=TTM, r=r+0.0001, d=d, sigma=sigma, N=N)
  } else{
    Y1 = American.Put(S0=S0, K=K, TTM=TTM, r=r, d=d, sigma=sigma, N=N)
    Y2 = American.Put(S0=S0, K=K, TTM=TTM, r=r+0.0001, d=d, sigma=sigma, N=N)
  }
  Rho = (Y2-Y1)/0.0001
  
  #Diving by 100% to be consistent with Bloomberg Reporting
  Rho = Rho/100
  return(Rho)
}
# Function Testing
#Binom_Rho(S0=100, K=100, TTM=1, r=0.05, d=0, sigma=0.2, N=50, type="Call")


# Read in Option Data from Bloomberg for Greek Comparison
Greeks.Data <- read.csv(file = "Greeks_Data.csv")


#Calculate Delta From Binomial Tree
Greeks.Data$Delta <- mapply(FUN = Binom_Delta,
                            S0 = Greeks.Data$Stock.Price,
                            K = Greeks.Data$Strike,
                            TTM = Greeks.Data$TTM,
                            r = int,
                            d = 0,
                            sigma = Greeks.Data$Vol,
                            N = 50,
                            type = Greeks.Data$Type)

#Calculate Gamma From Binomial Tree
Greeks.Data$Gamma <- mapply(FUN = Binom_Gamma,
                            S0 = Greeks.Data$Stock.Price,
                            K = Greeks.Data$Strike,
                            TTM = Greeks.Data$TTM,
                            r = int,
                            d = 0,
                            sigma = Greeks.Data$Vol,
                            N = 50,
                            type = Greeks.Data$Type)

#Calculate Vega From Binomial Tree
Greeks.Data$Vega<- mapply(FUN = Binom_Vega,
                            S0 = Greeks.Data$Stock.Price,
                            K = Greeks.Data$Strike,
                            TTM = Greeks.Data$TTM,
                            r = int,
                            d = 0,
                            sigma = Greeks.Data$Vol,
                            N = 50,
                            type = Greeks.Data$Type)

#Calculate Theta From Binomial Tree
Greeks.Data$Theta <- mapply(FUN = Binom_Theta,
                          S0 = Greeks.Data$Stock.Price,
                          K = Greeks.Data$Strike,
                          TTM = Greeks.Data$TTM,
                          r = int,
                          d = 0,
                          sigma = Greeks.Data$Vol,
                          N = 50,
                          type = Greeks.Data$Type)

#Calculate Rho From Binomial Tree
Greeks.Data$Rho <- mapply(FUN = Binom_Rho,
                            S0 = Greeks.Data$Stock.Price,
                            K = Greeks.Data$Strike,
                            TTM = Greeks.Data$TTM,
                            r = int,
                            d = 0,
                            sigma = Greeks.Data$Vol,
                            N = 50,
                            type = Greeks.Data$Type)

#Final Data Table
write.csv(Greeks.Data, "Greeks_Output.csv")

########################################################################################################
#   5.) Hedging Exercise
########################################################################################################
# Generate a sample of e.g. 10 daily stock prices Si following from the Black-Scholes model. 
# For each day, compute: American put option prices P(ti), Delta of the option ???(ti).

#Plot the option prices P(ti) vs ti, and compare with the plot of the prices of the hedged portfolio P(ti) ??? ???(ti)Si

#For this trial we consider a stock with S0 = 100, vol = 40%, K = 100, TTM = 1M, r = 5%, u = 8%

#Define a function to output a realization of the stock price by the BSM model, using Eulers Method
BSM_Stock_Path <- function(n, S0, r, d, TTM, vol){
  dW <- rnorm(n, mean = 0, sd = 1)
  dT <- TTM/n
  Xnext <- log(S0)
  Xpath <- vector("numeric")
  Xpath[1] <- Xnext
  for(i in 1:n){
    Xnext = (r-d-1/2*vol^2)*dT+vol*dW[i]*sqrt(dT)+Xnext
    Xpath[i+1] <- Xnext
  }
  return(exp(Xpath))
}

# Generate 10 daily stock price realizations via the BSM model
stk_path <- BSM_Stock_Path(n=10, S0=100, r=0.05, d=0, TTM=1/12, vol = 0.4)

# Use American.Put Function to price options based on the Binomial Tree
opt_prices <- mapply(FUN = American.Put, S0 = stk_path, K = 100, TTM = 1/12, r = 0.05, d=0, sigma = 0.4, N = 50)

# Compute the option delta using the Binomial Tree, and a finite difference dS of 0.1
inc_stk_path <- stk_path + 0.1
inc_opt_prices <-  mapply(FUN = American.Put, S0 = inc_stk_path, K = 100, TTM = 1/12, r = 0.05, d=0, sigma = 0.4, N = 50)
opt_deltas <- (inc_opt_prices - opt_prices)/0.1

#Create a vector of delta changes for computing hedge requirement
#shares <- opt_deltas[-1] - opt_deltas[-length(opt_deltas)]
#shares <- c(opt_deltas[1], shares)


# Create a vector of delta stock movements for hedged portfolio comparison
#delta_stk <- stk_path[-1] - stk_path[-length(stk_path)]
#delta_stk <- c(0, delta_stk)

# Calculate value of delta hedged portfolio
hedged_port = opt_prices - opt_deltas*stk_path
opt_pct_change = opt_prices / opt_prices[1] * 100 - 100
port_pct_change = hedged_port / hedged_port[1] * 100 - 100

# Result DataFrame
df <- cbind(1:length(stk_path)-1, stk_path, opt_prices, opt_deltas, hedged_port, opt_pct_change, port_pct_change)
colnames(df) <- c("ti", "Si", "Pi", "Deltas", "Hedged Port", "Pi % Change", "Hedged % Change")
df
