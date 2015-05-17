#!/bin/sh
# setup CRY
export SOFTWAREHOME=/cvmfs/cvmfs.cmsaf.mit.edu/t2srv0008/export/app
export CRYHOME=${SOFTWAREHOME}/d-Chooz/Software/cry_v1.7
export CRYDATAPATH=${CRYHOME}/data
PATH=${CRYHOME}/test:${PATH}
LD_LIBRARY_PATH=${CRYHOME}/test:${CRYHOME}/lib:${LD_LIBRARY_PATH}
# setup RAT
RATROOT=ratpac-kpipe
PATH=${RATROOT}/bin:${PATH}
LD_LIBRARY_PATH=${RATROOT}/lib:${LD_LIBRARY_PATH}
GLG4DATA=${RATROOT}/data
PYTHONPATH=${RATROOT}/python:${PYTHONPATH}
# push to environment
export RATROOT PATH LD_LIBRARY_PATH DYLD_LIBRARY_PATH GLG4DATA PYTHONPATH
