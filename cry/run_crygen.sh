
# GET ARGUMENTS
CRYOUTPUT=$1
JOBID=$2

# SETUP RAT
tar zxvf rat.tar.gz
source ratpac-kpipe/tier2scripts/env_condornode.sh

# GO TO EXECUTABLE
cd ratpac-kpipe/cry

# GEN CRY FILE
./gen_cosmics_kpipe setup.file ../../$CRYOUTPUT 1000 $JOBID

# GO HOME
#cd ../../