import os

import settings

# TODO: Queuing system requires additional look. Currently PyMemDyn 1.5 works 
#       when queued (in Slurm/CBS) externally (bash-shell calling pymemdyn)
#       Unknown if the same is the case for the other queuing systems.

class Queue(object):
    def __init__(self, *args, **kwargs):
        # Default number of processors, nodes and time alloted in cluster.
        self.num_proc   = getattr(settings, "QUEUE_NUM_PROCS") or 16
        self.num_node   = getattr(settings, "QUEUE_NUM_NODES") or 1
        self.max_time   = getattr(settings, "QUEUE_MAX_TIME") or "47:59:59"
        self.ntasks     = getattr(settings, "QUEUE_NUM_TASK") or 16
        self.ntaskpern  = getattr(settings, "QUEUE_NTS_NODE") or 16
        self.sh = "./mdrun.sh"

    def set_mdrun(self, value):
        """
        Set the md_run command
        """
        self._mdrun = value

    def get_mdrun(self):
        return self._mdrun
    mdrun = property(get_mdrun, set_mdrun)


class NoQueue(Queue):
    """
    Dummy queue when no queue is selected
    """
    def __init__(self, *args, **kwargs):
        super(NoQueue, self).__init__(self, *args, **kwargs)
        self.command = [self.sh]

        self._mdrun = os.path.join(settings.GROMACS_PATH, "gmx mdrun")
#        self._mdrun = os.path.join(settings.GROMACS_PATH, "gmx_mpi mdrun")
#        self._mdrun = os.path.join(settings.GROMACS_PATH, "gmx mdrun_mpi") #For triolith

    def make_script(self, workdir, options):
        """
        workdir is the path to the binary executable
        options is a list with all the options
        """
        sh = open(self.sh, "w")
        sh.write("#!/bin/bash\n")
        sh.write("cd %s\n" % workdir)
        ncpus = str(os.environ.get("SLURM_CPUS_PER_TASK",self.num_proc)) # Get the CPUS from the sbatch SLURM script, otherwise use the settings
        sh.write("%s %s" % (self.mdrun, " ".join(options)))
        sh.write(" -nt %s\n" % ncpus)
        sh.close()
        os.chmod(self.sh, 0o755)

        return True

class Slurm(Queue):
    """
    Queue for SLURM systems
    """
    def __init__(self, *args, **kwargs):
        super(Slurm, self).__init__(self, *args, **kwargs)
        self.command = ["srun",
#            "--ntasks=%s" % str(self.ntasks),
#            "--ntasks-per-node=%s" % str(self.ntaskpern),
#        self.command = ["srun",
#            "-n", str(self.num_node),
            "-c", str(self.num_proc),
            "-t", self.max_time,
            self.sh]
#
#        self._mdrun = os.path.join(settings.GROMACS_PATH, "gmx_mpi mdrun")
        self._mdrun = os.path.join(settings.GROMACS_PATH, "mdrun_slurm") #FOR CUELEBRE
#        self._mdrun = os.path.join(settings.GROMACS_PATH, "gmx mdrun") # FOR CSB

#        self._mdrun = os.path.join(settings.GROMACS_PATH, "mdrun_mpi") # FOR triolith

    def make_script(self, workdir, options):
        """
        workdir is the path to the binary executable
        options is a list with all options
        """
        sh = open(self.sh, "w")
        sh.write("#!/bin/bash\n")
#        sh.write("source /home/apps/gromacs-4.6.5/bin/GMXRC\n")
#        sh.write("source /home/apps/bin/apps.sh\n")
#        sh.write("module load openmpi-x86_64\n")
        sh.write("cd %s  \n" % workdir)
#        sh.write("%s -ntmpi 16 -ntomp 1  %s -v&> mdrun.log\n" % (self.mdrun, " ".join(options)))
#        sh.write("%s -nt 8 %s -v&> mdrun.log\n" % (self.mdrun, " ".join(options)))
#        sh.write("%s %s\n" % (self.mdrun, " ".join(options)))
        sh.write("%s %s -v&> mdrun.log\n" % (self.mdrun, " ".join(options)))
#        sh.write("mpprun %s %s -v&> mdrun.log\n" % (self.mdrun, " ".join(options))) # Triolith needs mpprun
        sh.close()
        os.chmod(self.sh, 0o755)

        return True


class PBS(Queue):
    """
    Queue for PBS systems
    """
    def __init__(self, *args, **kwargs):
        super(PBS, self).__init__(self, *args, **kwargs)
        '''Setting the command to run mdrun in pbs queue with mpi'''
        # These values are here for reference, doesn't do NOTHING      #
        # Calling file run.sh should resemble this lines               #
        self.num_nodes = 1                                             #
        self.proc_per_node = 8                                         #
        self.max_time = getattr(settings, "QUEUE_MAX_TIME") or "72:00:00"#
        self.max_cpu_time = "72:00:00"                                 #
        self.max_mem = "12gb"                                          #
        self.command = ["qsub",                                        #
            "-nodes=%d:ppn=%d" % (self.num_nodes, self.proc_per_node), #
            "-walltime=%s" % self.max_time,                            #
            "-cput=%s" % self.max_cpu_time,                            #
            "-mem=%s" % self.max_mem,                                  #
            self.sh]                                                   #
                                                                       #
        # XXX ##########################################################

#        self._mdrun=os.path.join(settings.GROMACS_PATH, "mdrun_mpi")
        self._mdrun=os.path.join(settings.GROMACS_PATH, "gmx mdrun_")
        self.command = [self.sh]

    def make_script(self, workdir, options):
        """
        PBS must load some modules in each node by shell scripts
        options is a list with all the options
        """
        sh = open(self.sh, "w")
        sh.write("#!/bin/bash\n")
        sh.write("cd %s\n" % os.path.join(os.getcwd(), workdir))
        sh.write("module load gromacs/4.0.5-gige\n")
        sh.write("mpirun %s %s -v &>mdrun.log\n" \
            % (self.mdrun, " ".join(options)))
        sh.close()
        os.chmod(self.sh, 0o755)

        return True


class PBS_IB(Queue):
    def __init__(self, *args, **kwargs):
        super(PBS, self).__init__(self, *args, **kwargs)
        # USELESS, see class PBS for explanation                            #
        self.num_nodes = 10                                                 #
        self.proc_per_node = 4                                              #
        self.max_time = getattr(settings, "QUEUE_MAX_TIME") or "08:00:00"   #
        self.max_cpu_time = "320:00:00"                                     #
        self.max_mem = "12gb"                                               #
                                                                            #
        self.command = [                                                    #
            "qsub",                                                         #
            "-l nodes=%d:ppn=%d:ib" % (self.num_nodes, self.proc_per_node), #
            "-l walltime=%s" % self.max_time,                               #
            "-l cput=%s" % self.max_cpu_time,                               #
            "-l mem=%s" % self.max_mem,                                     #
            self.sh]                                                        #
        #####################################################################

        self._mdrun=os.path.join(settings.GROMACS_PATH, "gmx mdrun_mpi")
        self.command = [self.sh]

    def make_script(self, workdir, options):
        """
        PBS must load some modules in each node by shell scripts
        options is a list with all the options
        """
        sh = open(self.sh, "w")
        sh.write("#!/bin/bash\n")
        sh.write("cd %s\n" % os.path.join(os.getcwd(), workdir))
        sh.write("module load gromacs\n")
        sh.write("mpirun %s %s -v &>mdrun.log\n" \
            % (self.mdrun, " ".join(options)))
        sh.close()
        os.chmod(self.sh, 0o755)

        return True


class Svgd(Queue):
    """
    Queue for the PBS system at svgd.cesga.es
    """
    def __init__(self, *args, **kwargs):
        super(Svgd, self).__init__(self, *args, **kwargs)
        '''Setting the command to run mdrun in pbs queue with mpi'''
        self._mdrun=os.path.join(settings.GROMACS_PATH, "gmx mdrun")
        self.command = [self.sh]

    def make_script(self, workdir, options):
        """
        PBS must load some modules in each node by shell scripts
        options is a list with all the options
        """
        sh = open(self.sh, "w")
        sh.write("#!/bin/bash\n")
        sh.write("cd %s\n" % os.path.join(os.getcwd(), workdir))
        # Somehow impi is loaded, and conflicts with the (see down) tricky way
        # of SVGD of dealing with parallel runnings of mdrun through mpich2
        sh.write("module unload impi\n")
        sh.write("module load mpich2\n")
        sh.write("module load gromacs/4.0.7\n")
        # CESGA SVGD has its own tweaks. The mdrun binary "jumps" to mpiexec,
        # and call back mdrun with as much nslots (~cores) as reserved in the
        # command line calling the whole pipeline, this way:
        #  qsub -l num_proc=1,s_rt=01:00:00,s_vmem=1G,h_fsize=1G,arch=amd \
        #    -pe mpi 4 run.sh
        # That "-pe mpi 4" tells the queue system to allocate 4 cores to run.sh
        # and it's responsability of run.sh (aka the next line) to pass $NSLOTS
        # Note that this SVGD uses THE SAME queue system than PBS and PBS_IB,
        # but lacks of mpirun executable. Always a pleasure to adjust a queue.
        sh.write("%s -np $NSLOTS %s -v &>mdrun.log\n" \
            % (self.mdrun, " ".join(options)))
        sh.close()
        os.chmod(self.sh, 0o755)

        return True
