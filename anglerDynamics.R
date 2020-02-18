
Years = 100

######### Angler fish/fishing choice #######
# First we set up expected catch 
expCa <-  matrix(NA, 1, Years)
expCb <- matrix(NA, 1, Years)
expCa[,1] = 3
expCb[,1] = 2
# we need actual catch to do this. 
actCa <-  matrix(NA, 1, Years)
actCb <- matrix(NA, 1, Years)

# 'Utility' - not sure this is the right thing to call this
Ua <- matrix(0.7, 1, Years) 
Ub <- matrix(0.3, 1, Years) 
Usit = 5 # this should be rescaled based on a few simulations - we want some anglers to sit out. 

#proportions
ProAnga <- matrix(NA, 1, Years)  #starting distribution of anglers
ProAngb <- matrix(NA, 1, Years)  
ProAngsit <- matrix(NA, 1, Years)
# preference determination:
pi <- .7 #preference for species A (.5 = equal preference, 1 = only species A)
# memory term for expected cacth
lambda = 0.7

#####################
for (i in 1:Years){
  actCa[,i] = qa * Ea * Na[,i] # fishing mortality 
  actCa[,i] = qb * Eb * Nb[,i] 
  actCa[,i] = rpois(length(actCa[,i]),actCa[,i]) # actual catch 
  actCb[,i] = rpois(length(actCb[,i]),actCb[,i])
  expCa[,(i+1)] <- lambda*actCa[,i] + (1-lambda)*expCa[,i] #updating expectation
  expCb[,(i+1)] <- lambda*actCb[,i] + (1-lambda)*expCb[,i]
  Ua[,i] = pi*expCa[i] #weight for species 1
  Ub[,i] = (1-pi)*expCb[i] #weight for species 2
  ProAnga[,i] = Ua[,i]/(Ua[,i]+Ub[,i]+Usit) #how anglers distribute
  ProAngb[,i] = Ub[,i]/(Ua[,i]+Ub[,i]+Usit)
  ProAngsit[,i] = Usit/(Ua[,i]+Ub[,i]+Usit)
}

matplot(ProAngb,type='l',lty=1, ylim=c(0,1))
matplot(ProAnga, type='l',lty=1, ylim=c(0,1))
matplot(ProAngsit, type='l',lty=1, ylim=c(0,1))

#plug into model through C = qX*(E*ProAnga) -> ProAng is just a modifier of Effort -
# determines what percent of total effort is allocated to each species
      #basically: Aharv and Bharv
#The utility function s just so we can measure how happy anglers are - some percentage will 
  #likely not fish - which is expected. 
# utility is a function of exp catch - determines the choice of fishing and for what species
#slight concern: current effort model is literally just the angler dynamics
          # do we want a more complex model of effort (in this, effort is farily static, it
           #   just switches between species A, species B, and sitting)


