### A simple model for disease spread with isolation probability after infection
# net: the square matrix indicating the amount of contact between two network members
# time50: the amount fo contact that results in a 50% contagion probability
# roundsBeforeRmoval: How long a network member is contagious before it is removed
# periods: how many rounds the simulation runs
# pIsolation: The probability for a network member to immediately isoalte after being infected
# is.directed: is the network directed or undirected

contagionModel <- function(net, time50=300, roundsBeforeRemoval=1, periods=20, pIsolation=0.3, is.directed=T){

  if(is.directed){
    # symmetrise the network
    network <- net + t(net)
  } else {
    network <- net
  }
  
  # chose a random network member to be infected first
  init <- rep(0, nrow(network))
  init[floor(runif(1,1,(nrow(network) + 1)))] <- 1
  
  # create a variable whether any network member is sick
  # and make the initator one sick (state 1)
  is.sick <- init
  
  # decide whether the first infected network member isolates,
  # if yes, move it to removed state (-1)
  isolates <- floor(runif(1) + pIsolation)
  if(isolates == 1){
    is.sick[is.sick==1] <- -1
  }
  
  # create varaible that will be used remove network members 
  # after they were infectious for the specfied number of rounds
  rounds.sick <- rep(0, nrow(network))
  
  # generate help variable
  all <- 1:nrow(network)
  
  # create output variable as list
  sick <- list()
  
  # loop through all periods
  for(i in 1:periods){
    # write current state in output
    sick[[i]] <- is.sick
    
    # get the sick ones
    who.is.sick <- all[all*(is.sick) > 0]
    
    # loop through all sick ones
    for(sick.one in who.is.sick){
      
      # find the contacts they have that can be infected
      at.risk <- all[network[,sick.one] > 0]
      
      # loop though each one in the risk set and decide whether they will be infected
      for(alter in at.risk){
        
        # see if the alter is already infected or removed
        if(is.sick[alter] == 0){
          
          # determine whether the alter gets infected using the time50
          is.sick[alter] <- floor(runif(1) + (((network[alter,sick.one])/ (time50)) / 
                                                (1 + (network[alter,sick.one])/ (time50))))
          
          # if the new other is infected, decide whether it isolates
          # in that case, move it to state removed (-1)
          if(is.sick[alter] == 1){
            isolates <- floor(runif(1) + pIsolation)
            if(isolates == 1){
              is.sick[alter] <- -1
            }
          }
        }
      }
      # find out how long the infector mouse is sick and remove it if for longer than specified
      if(rounds.sick[sick.one] >= roundsBeforeRemoval - 1){
        is.sick[sick.one] <- -1
      } else {
        rounds.sick[sick.one] <- rounds.sick[sick.one] + 1
      }
    }
  }
  
  return(sick)
}

# example
# generate random network
randomNetwork <- matrix(floor(runif(225)*100), 15, 15)

# run simulation
results <- contagionModel(net = randomNetwork, time50 = 75, pIsolation = 0.2)

# results are interpreted as follows:
# each element of the results list refers to the status of all network members at the end of that round
# 1 means infected, -1 means removed (i.e. infected in the past but not contagious anymore), 0 means susceptible
