
# GET ARGUMENTS
INPUT=$1
OUTPUT=$2
PMTINFO=$3

tar zxvf rat.tar.gz

if test -z "${_CONDOR_SCRATCH_DIR}"; then
export _CONDOR_SCRATCH_DIR=/tmp/$USERNAME/mytmpdir/
fi
mkdir -p ${_CONDOR_SCRATCH_DIR}

# SETUP ENVIONMENT
source ratpac-kpipe/tier2scripts/env_condornode.sh

# GO TO EXECUTABLE
cd ratpac-kpipe/analysis

# RUN RAT
./analyze_data ${INPUT} ${OUTPUT}  ${PMTINFO}

