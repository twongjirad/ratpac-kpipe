
# GET JOBID FROM ARGUMENTS
JOBID=$1
CRY_INPUT=$2
KPIPE_CRY_OUT=$3

# SETUP LOCAL PYTHON
tar zxvf python2.6.9.tar >> /dev/null

# SETUP RAT
tar zxvf rat.tar.gz >> /dev/null
#source ratpac-kpipe/tier2scripts/env_condornode.sh
source ratpac-kpipe/tier2scripts/env_condornode_localpython.sh

# GEN MACRO
CRYMAC=cry_job"$JOBID".mac
python ratpac-kpipe/cry/gen_macro.py $CRY_INPUT $CRYMAC $JOBID 1000
cat $CRYMAC

# RUN JOB
rat $CRYMAC -o $KPIPE_CRY_OUT

