#setPijk.R
# set Pijk parameters - to make incoming rate at every node less than or equal to 1, 
# given interarrival rate of 1
# however, random input samples may give unstable queues so still need capacitated (say 50)

setPijk = function(PijkRoute)
{
  imax = dim(PijkRoute)[1]
  jmax = dim(PijkRoute)[2]
  kmax = dim(PijkRoute)[3]
  klower = (kmax/2) - 1
  khigher = (kmax/2) + 2
  kmid1 = klower + 1
  kmid2 = klower + 2



  i = 1
# all stage 1 nodes equally likely to allstage 2 nodes
  for (j in 1:jmax) {
      for (k in 1:kmax) {
          PijkRoute[i,j,k] = 1.0/kmax
      }
  }

  i = 2
# lower index stage 2 nodes go either to stage 3 node 1 or to exit (node 3)
  for (j in 1:klower) {
      PijkRoute[i,j,1] = .2  # to stage 3 node 1
      PijkRoute[i,j,3] = .8  # to exit
  }
# middle two nodes in second stage go to both stage 3 final nodes and to exit
  PijkRoute[i,kmid1,1] = .20      # to upper stage 3 node
  PijkRoute[i,kmid1,2] = .20      # to lower stage 3 node
  PijkRoute[i,kmid1,3] = .6       # to exit
  PijkRoute[i,kmid2,1] = .20      # to upper stage 3 node
  PijkRoute[i,kmid2,2] = .20      # to lower stage 3 node
  PijkRoute[i,kmid2,3] = .6       # to exit
# higher stage 2 nodes go either to stage 3 node 2 or to exit (node 3)
  for (j in khigher:kmax) {
      PijkRoute[i,j,2] = .2  # to lower stage 3 node
      PijkRoute[i,j,3] = .8  # to exit
  }

return(PijkRoute)
}
