# GET JOBID FROM ARGUMENTS
JOBID=$1
KPIPE_KDAR_OUT=$2

# SETUP LOCAL PYTHON
tar zxvf python2.6.9.tar > /dev/null

# SETUP RAT
tar zxvf rat.tar.gz > /dev/null
source ratpac-kpipe/tier2scripts/env_condornode_localpython.sh

# UNPACK KADR SOURCE
bzip2 -d eventsout_nuwroxsec_numu_kpipe_3_25_2015.root.bz2

# GEN MACRO
RATMAC=kpipe_kdar_job"$JOBID".mac
python ratpac-kpipe/tier2scripts/gen_macro.py $RATMAC $JOBID ../../eventsout_nuwroxsec_numu_kpipe_3_25_2015.root ../data/kpipe/kpipe_nosipms.gdml
cat $RATMAC

cd ratpac-kpipe/kdar_muons

# RUN JOB
rat ../../$RATMAC -o ../../$KPIPE_KDAR_OUT

cd ../../

