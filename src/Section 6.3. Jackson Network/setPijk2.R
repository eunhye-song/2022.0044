#setPijk2.R
# set Pijk parameters - to make incoming rate at every node less than or equal to 1, 
# given interarrival rate of 1
# however, random input samples may give unstable queues so still need capacitated (say 100)
# modified from setPijk to add delta to routing probs for first m/2 to first m/2
# and subtract delta from routing probs for first m/2 to second m/2
# updated 1.28.20 to include nomArate as passed parameter
# PijkRoute is a matrix of dimensions of: (# layers - 1), (numnodes layer 1,)    
#  (numnodes layer 2), where PijkRoute[i,j,k] = routing probability from 
#  jth node in layer i to kth node in layer i+1
#  NOTE: when i = 2, k takes only three values, not numArv

setPijk2 = function(PijkRoute,delta)
{
  imax = dim(PijkRoute)[1]
  jmax = dim(PijkRoute)[2]
  kmax = dim(PijkRoute)[3]
# NOTE: kmax = jmax = m
  klower = (kmax/2) - 1
  khigher = (kmax/2) + 2
  kmid1 = klower + 1
  kmid2 = klower + 2

  i = 1
# stage 1 nodes 1:m/2 unequal between 1:m/2 and 1+m/2:m 
for (j in 1:(jmax/2)) {
      for (k in 1:(kmax/2)) {
          PijkRoute[i,j,k] = (1.0/kmax) + delta
      }
      for (k in ((kmax/2)+1):kmax) {
          PijkRoute[i,j,k] = (1.0/kmax) - delta
      }
  }
# stage 1 nodes 1+m/2:m equally likely to all stage 2 nodes
  for (j in (1+(jmax/2)):jmax) {
      for (k in 1:kmax) {
          PijkRoute[i,j,k] = 1.0/kmax
      }
  }


  i = 2
# lower index stage 2 nodes go either to stage 3 node 1 or to exit (node 3)
# probabilities set to give sum of routing probabilities 3 nodes = 1 
# and average loading to match average layer 2 loading, so need to have
# routing probabilities to 1 and 2 reduced since m/2 + 1 nodes feeding each.
# so PijkRoute to a stage 3 node must average 1/(1.0+(kmax/2.0)) 
#   note kmax = m

  for (j in 1:klower) {
      PijkRoute[i,j,1] = 1/(1.0+(kmax/2.0))            # to upper stage 3 node
      PijkRoute[i,j,3] = 1.0-(1/(1.0+(kmax/2.0)))      # to exit
  }
# middle two nodes in second stage go to both stage 3 final nodes and to exit
  PijkRoute[i,kmid1,1] = 1/(1.0+(kmax/2.0))            # to upper stage 3 node
  PijkRoute[i,kmid1,2] = 1/(1.0+(kmax/2.0))            # to lower stage 3 node
  PijkRoute[i,kmid1,3] = 1.0-2.0*(1/(1.0+(kmax/2.0)))  # to exit
  PijkRoute[i,kmid2,1] = 1/(1.0+(kmax/2.0))            # to upper stage 3 node
  PijkRoute[i,kmid2,2] = 1/(1.0+(kmax/2.0))            # to lower stage 3 node
  PijkRoute[i,kmid2,3] = 1.0-2.0*(1/(1.0+(kmax/2.0)))  # to exit
# higher stage 2 nodes go either to stage 3 node 2 or to exit (node 3)
  for (j in khigher:kmax) {
      PijkRoute[i,j,2] = 1/(1.0+(kmax/2.0))            # to lower stage 3 node
      PijkRoute[i,j,3] = 1.0-(1/(1.0+(kmax/2.0)))      # to exit
  }

return(PijkRoute)
}
