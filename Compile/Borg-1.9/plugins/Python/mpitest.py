# Prerequisites for this example:
#   1. Run on a system with MPI
#   2. Compile the master-slave version of Borg and copy its DLL / shared
#      object into this folder
#   3. Create a submission script for the MPI job.  An example, mpitest.pbs,
#      is provided but may require modifications to work.
#   4. Submit the job (e.g., qsub mpitest.pbs)

from borg import *

Configuration.startMPI()

borg = Borg(2, 2, 0, lambda x,y : [x**2 + y**2, (x-2)**2 + y**2])
borg.setBounds([-50, 50], [-50, 50])
borg.setEpsilons(0.01, 0.01)

result = borg.solveMPI(maxTime=0.1)

# the result will only be returned from one node
if result:
	result.display()

Configuration.stopMPI()
