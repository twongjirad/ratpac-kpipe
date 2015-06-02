#!/bin/sh

export TLOCAL=/net/hisrv0001/home/taritree/.local
export EDITOR="emacs -nw"
export DISPLAY=hisrv0001.cmsaf.mit.edu:0.0

# SETUP COMMON SOFTWARE
export SOFTWAREHOME=/cvmfs/cvmfs.cmsaf.mit.edu/t2srv0008/export/app
export PATH=${SOFTWAREHOME}/d-Chooz/Software/bin:${PATH}
export LD_LIBRARY_PATH=${SOFTWAREHOME}/d-Chooz/Software/lib:${LD_LIBRARY_PATH}

# SETUP GIT
export PATH=${SOFTWAREHOME}/d-Chooz/Software/git-1.8.5.6:${PATH}

# SETUP LOCAL PYTHON
export PYTHONHOME=export/app/d-Chooz/Software/PYTHON/python2.6.9
export PATH=${PYTHONHOME}/bin:${PATH}
export LD_LIBRARY_PATH=${PYTHONHOME}/lib:${LD_LIBRARY_PATH}

# SETUP ROOT 
source ${SOFTWAREHOME}/d-Chooz/Software/root/root_v5.34.28/bin/thisroot.sh

# SETUP GEANT4
export PATH=${SOFTWAREHOME}/d-Chooz/Software/geant4.9.6.p04/bin:${PATH}
export LD_LIBRARY_PATH=${SOFTWAREHOME}/d-Chooz/Software/geant4.9.6.p04/lib64:${LD_LIBRARY_PATH}
source ${SOFTWAREHOME}/d-Chooz/Software/geant4.9.6.p04/bin/geant4.sh

# setup CRY
export CRYHOME=${SOFTWAREHOME}/d-Chooz/Software/cry_v1.7
export CRYDATAPATH=${CRYHOME}/data
PATH=${CRYHOME}/test:${PATH}
LD_LIBRARY_PATH=${CRYHOME}/test:${CRYHOME}/lib:${LD_LIBRARY_PATH}
# setup RAT
RATROOT=${PWD}/ratpac-kpipe
PATH=${RATROOT}/bin:${PATH}
LD_LIBRARY_PATH=${RATROOT}/lib:${LD_LIBRARY_PATH}
GLG4DATA=${RATROOT}/data
PYTHONPATH=${RATROOT}/python:${PYTHONPATH}
# push to environment
export RATROOT PATH LD_LIBRARY_PATH DYLD_LIBRARY_PATH GLG4DATA PYTHONPATH
