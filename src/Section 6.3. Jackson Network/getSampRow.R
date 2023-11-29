# getSampRow.R

# Purpose: assemble sample interarrival time means and coded sample PIJK to add 
#          to elliptical DOE (include observed parameter estimates in design)
#

# Inputs: 
#   numArv - number of arrival nodes
#   PijkDim -  vector of dimensions of: (# layers - 1), (numnodes layer 1,)    
              # (numnodes layer 2), where PijkRoute = routing probability from 
              #  jth node in layer i to kth node in layer i+1
#   PijkRoute - true routing probabilities, use in constructing reduced parm vector
#   numVecs - number of parameter vectors to generate (e.g. B0 or other)
#   samples - the bootstrapped (BSsamples or newSamples)data
#   probMethod - how to transform probabilities: 
#                'allBut1','pseudoP','Cartesian' 

# Locals:
#   arvMat - a single row matrix with sample interarrival time means (mArvTime)
#   oneVec - argument to data2parms stating only one row vector to xform
#   sampProb - PijkHat placed in format for data2parms
#   origSamp - list of interarrival means and Pijk values passed to data2parms
#   MMDOEsamp - list returned from data2parms (number of MM parameters, parameter matrix) 

getSampRow = function(numArv,mArvTime,PijkDim,PijkHat,probMethod) 
{
arvMat = matrix(0,1,numArv)
arvMat[1,] = mArvTime               # put sample mean interarrival time in a one row matrix
oneVec = 1
sampProb = array(0, dim = c(oneVec,PijkDim))
sampProb[1,,,] = PijkHat            # put sample routing probs in array
origSamp = list(arvMat,sampProb)    # create list like BSsamples for data2parms
MMDOEsamp = data2parms(numArv,PijkDim,PijkRoute,probMethod,oneVec,origSamp)
sampParms = MMDOEsamp[[2]]          # pull out mean arv time and coded prob. values

return(sampParms)
}
