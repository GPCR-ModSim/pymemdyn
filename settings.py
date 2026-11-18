"""This module handles the local settings for pymemdyn on your machine. 
The settings are mostly paths and run settings.
"""

import os

# This is the folder where pymemdyn git repo has been deployed,
# or to be more specific, where the settings.py file is located at.
ROOT_DIR = os.path.dirname(os.path.realpath(__file__))


# The following variable can be given as a full path to the place where
# the templates reside or it can also use a relative path to the previously
# defined ROOT_DIR.
TEMPLATES_DIR = os.path.join(ROOT_DIR, "templates")
# Or using the full path like so:
#TEMPLATES_DIR = "/path/to/your/templates"


# Define a path to the GROMACS (version 5.0>) binaries. See examples of GROMACS_PATHs below:
GROMACS_PATH = "/home2/apps/easybuild/software/GROMACS/2023.4-foss-2023a/bin/" 
#GROMACS cuda
#GROMACS_PATH = "/home2/apps/easybuild/software/GROMACS/2023.3-foss-2023a-CUDA-12.1.1-PLUMED-2.9.0/bin/" 
#GROMACS_PATH = "/usr/bin/gmx"
#GROMACS_PATH = "/Users/hugo/homebrew/bin/gmx"
#GROMACS_PATH = "/usr/local/gromacs/bin/gmx"
#GROMACS_PATH = "/home/hugo/bin/gmx_mpi"

# example GROMACS_PATHs:
#GROMACS_PATH = "/home/apps/apps/.spack/sandybridge/gcc-10.2.0/gromacs-2021-lx52hldze4odq56frfkrf4vdhifqb7do/bin/"                 
                                                                     #csb.bmc.uu.se

# The paths below are outdated, but can be used as reference to get your 
# GROMACS_PATH:
#GROMACS_PATH = "/opt/applications/gromacs/4.0.5/gnu/ib/bin/"
#GROMACS_PATH = "/opt/applications/gromacs/4.0.5/gnu/gige/bin/"
#GROMACS_PATH = "/opt/cesga/gromacs-4.0.7/bin/"
#GROMACS_PATH = "/opt/gromacs405/bin/"                               #cuelebre.inv.usc.es
#GROMACS_PATH = "/software/apps/gromacs/4.6.3/g472/bin/"             #Triolith
#GROMACS_PATH = "/sw/bin/"                                           #Standalone in Mac Fink
#GROMACS_PATH = "/Users/esguerra/software/gromacs-4.6.7/bin/"        #Standalone in Mac
#GROMACS_PATH = "/c3se/apps/Glenn/gromacs/4.6.3-p20130821-gcc48/bin" #Glenn at Chalmers
#GROMACS_PATH = "/c3se/apps/Glenn/gromacs/5.0.4-gcc48-cuda/bin/"     #Glenn GPU at Chalmers
#GROMACS_PATH = "/c3se/NOBACKUP/apps/Hebbe/EB/software/GROMACS/4.6.7-intel-2015b-hybrid.wip/GROMACS/4.6.7-intel-2015b-hybrid/bin" 
                                                                     #Hebbe at Chalmers
#GROMACS_PATH = "/sw/apps/gromacs/4.6.3/tintin/bin"                  #Tintin
#GROMACS_PATH = "/lap/gromacs/4.6.5/bin"                             #Abisko


# Define a path to the clustalw binary.
CLUSTAL_BIN = os.path.join(ROOT_DIR, ".bin/clustalw_linux")
#CLUSTAL_BIN = os.path.join(ROOT_DIR, ".bin/clustalw_mac")
#CLUSTAL_BIN = os.path.join(ROOT_DIR, ".bin/clustalw_mac64")


# Choose how many nodes to use in parallel
QUEUE_NUM_NODES = 1

# Choose how many processor to use in parallel
QUEUE_NUM_PROCS = 32

# Choose the maximum alloted time for your run.
QUEUE_MAX_TIME = "47:59:59"

QUEUE_NUM_TASK = 32
QUEUE_NTS_NODE = 32
