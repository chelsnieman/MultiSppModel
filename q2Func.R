#script to hold the q2 function along so we can just source it in Rmarkdown 
#parameter definitions
# qE - harvest, species specific, this is what is controlled by 'regulations'
# s - juvenile overwinter survial, species specific
# m - adult natural mortality rate, species specific
# cJA - effect of adults of a given species on juveniles of a given species (cover cannibalism or interspecific predation, both happen in foraging arena)
# cJJ - effect of juveniles of one species on juveniles of the other (can be predation or competition)
# h - rate at which juveniles leave foraging arena for refuge, species specific
# v - rate at which juveniles enter foraging arena from refuge, species specific
# f - fecundity, species specific
# stock1 - annual stocked num spp1
# stock2 - annual stocked num spp2

#base model, closest to Biggs'
# simBiggsQ2<-function(t,y,params){
#   A1<-y[1]
#   A2<-y[2]
#   J1<-y[3]
#   J2<-y[4]
#   with(as.list(params),{
#     dA1dt=-qE1Fun(t)*A1-m1*A1+s1*J1
#     dA2dt=-qE2Fun(t)*A2-m2*A2+s2*J2
#     dJ1dt=-cJ1J2*J2*J1-cJ1A1*J1*A1-((cJ1A2*v1*A2*J1)/(h1Fun(t)+v1+cJ1A2*A2))+f1*A1+st1Fun(t)
#     dJ2dt=-cJ2J1*J1*J2-cJ2A2*J2*A2-((cJ2A1*v2*A1*J2)/(h2Fun(t)+v2+cJ2A1*A1))+f2*A2+st2Fun(t)
#     return(list(c(dA1dt,dA2dt,dJ1dt,dJ2dt)))
#   })
# }

#version with juveniles subrtracted
# simBiggsQ2<-function(t,y,params){
#   A1<-y[1]
#   A2<-y[2]
#   J1<-y[3]
#   J2<-y[4]
#   with(as.list(params),{
#     dA1dt=-qE1Fun(t)*A1-m1*A1+s1*J1
#     dA2dt=-qE2Fun(t)*A2-m2*A2+s2*J2
#     dJ1dt=-cJ1J2*J2*J1-cJ1A1*J1*A1-((cJ1A2*v1*A2*J1)/(h1Fun(t)+v1+cJ1A2*A2))-s1*J1+f1*A1+st1Fun(t)
#     dJ2dt=-cJ2J1*J1*J2-cJ2A2*J2*A2-((cJ2A1*v2*A1*J2)/(h2Fun(t)+v2+cJ2A1*A1))-s2*J2+f2*A2+st2Fun(t)
#     return(list(c(dA1dt,dA2dt,dJ1dt,dJ2dt)))
#   })
# }

#cannnibalism is refuge dependent in this one
# simBiggsQ2<-function(t,y,params){
#   A1<-y[1]
#   A2<-y[2]
#   J1<-y[3]
#   J2<-y[4]
#   with(as.list(params),{
#     dA1dt=-qE1Fun(t)*A1-m1*A1+s1*J1
#     dA2dt=-qE2Fun(t)*A2-m2*A2+s2*J2
#     dJ1dt=-cJ1J2*J2*J1-((cJ1A2*v1*A2*J1)/(h1Fun(t)+v1+cJ1A2*A2))-((cJ1A1*v1*A1*J1)/(h1Fun(t)+v1+cJ1A1*A1))+f1*A1+st1Fun(t)
#     dJ2dt=-cJ2J1*J1*J2-((cJ2A1*v2*A1*J2)/(h2Fun(t)+v2+cJ2A1*A1))-((cJ2A2*v2*A2*J2)/(h2Fun(t)+v2+cJ2A2*A2))+f2*A2+st2Fun(t)
#     return(list(c(dA1dt,dA2dt,dJ1dt,dJ2dt)))
#   })
# }

# non-linear stock-recruitment relationships & refuge dependent cannibalism
#R.a*N*exp(-R.b*N)
simBiggsQ2<-function(t,y,params){
  A1<-y[1]
  A2<-y[2]
  J1<-y[3]
  J2<-y[4]
  with(as.list(params),{
    dA1dt=-qE1Fun(t)*A1-m1*A1+s1*J1
    dA2dt=-qE2Fun(t)*A2-m2*A2+s2*J2
    dJ1dt=-cJ1J2*J2*J1-((cJ1A2*v1*A2*J1)/(h1Fun(t)+v1+cJ1A2*A2))-((cJ1A1*v1*A1*J1)/(h1Fun(t)+v1+cJ1A1*A1))-s1*J1+(5000*A1/(500+A1))+st1Fun(t)
    dJ2dt=-cJ2J1*J1*J2-((cJ2A1*v2*A1*J2)/(h2Fun(t)+v2+cJ2A1*A1))-((cJ2A2*v2*A2*J2)/(h2Fun(t)+v2+cJ2A2*A2))-s2*J2+(5000*A2/(500+A2))+st2Fun(t)
    return(list(c(dA1dt,dA2dt,dJ1dt,dJ2dt)))
  })
}
