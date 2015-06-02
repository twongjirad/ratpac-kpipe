
# GET ARGUMENTS
INPUT=$1
OUTPUT=$2
PMTINFO=$3

tar zxvf rat.tar.gz

# SETUP ENVIONMENT
source ratpac-kpipe/tier2scripts/env_condornode.sh

# GO TO EXECUTABLE
cd ratpac-kpipe/analysis

# RUN RAT
./analyze_data ../../${INPUT} ../../${OUTPUT}  ../../${PMTINFO}

cd ../../

