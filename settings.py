#import pymoldyn

#SCRIPT_DIR = os.path.dirname(os.path.realpath(pymoldyn.__file__))
#ROOT_DIR = os.path.dirname(os.path.realpath(__file__))

REPO_DIR = "templates"

#GROMACS_PATH = "/opt/applications/gromacs/4.0.5/gnu/ib/bin/"
#GROMACS_PATH = "/opt/applications/gromacs/4.0.5/gnu/gige/bin/"
GROMACS_PATH = "/opt/gromacs405/bin/"

#Set the queue to use. Look inside queue.py.
#QUEUE = ""
QUEUE = "slurm"
#QUEUE = "pbs"
#QUEUE = "pbs_ib"
