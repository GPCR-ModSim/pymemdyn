PyMemDyn Version 2.0
====================

PyMemDyn is a standalone *python* package to setup membrane molecular
dynamics calculations using the **GROMACS** set of programs. The package
can be used either in a desktop environment, or in a cluster with
popular queuing systems such as Torque/PBS or Slurm.

**PyMemDyn** is hosted in github at:

https://github.com/GPCR-ModSim/pymemdyn

You can download any version of **PyMemDyn** by cloning the repository
to your local machine using git.

You will need to create a free personal account at github and send and
e-mail to: gpcruser@gmail.com requesting access to the code. After
request processing from us you will be given access to the free
repository.

Dependencies
------------

| **GROMACS**
| Pymemdyn is dependent on GROMACS. Download GROMACS
  `here <https://manual.gromacs.org/current/download.html>`__.
  Instructions for installation are
  `here <https://manual.gromacs.org/current/install-guide/index.html>`__.

| **LigParGen**
| In order to automatically generate .itp files for ligands and
  allosterics, the program ligpargen is used. Install using their
  instructions: https://github.com/Isra3l/ligpargen. Do not forget to
  activate the conda environment in which you installed ligpargen
  (``conda install py37`` if you followed the instructions) before
  running pymemdyn. In case you are using a bash script, this should be
  done inside the script. See also "ligpargen_example" in the folder
  examples.

Testing was done using LigParGen v2.1 using BOSS5.0.

Pymemdyn can also be used without ligpargen installation, but then .itp
files containing the parameters for the ligand(s) should be provided.

| **MODELLER**
| For replacing missing loops or missing side chains in the protein,
  MODELLER is used. The download and installation guid for MODELLER
  are available 
  `here <https://salilab.org/modeller/download_installation.html>`__.

| **Queueing system**
| A queuing system: although not strictly required, this is highly
  advisable since an MD simulation of 2.5 nanoseconds will be
  performed. However, if only membrane insertion and energy
  minimization is requested, this requirement can be avoided.
  Currently, the queuing systems supported include Slurm and PBS.

Installation
------------

To install **PyMemDyn** follow these steps:

1. Clone **PyMemDyn** for python 3.7:

   ::

      git clone https://username@github.com/GPCR-ModSim/pymemdyn.git

   Make sure to change *username* to the one you have created at github.

2. The previous command will create a *pymemdyn* directory. Now you have
   to tell your operating system how to find that folder. You achieve
   this by declaring the location of the directory in a .bashrc file
   .cshrc or .zshrc file in your home folder. An example of what you
   will have to include in your .bashrc file follows:

   ::

      export PYMEMDYN=/home/username/software/pymemdyn
      export PATH=$PYMEMDYN:$PATH

   or if your shell is csh then in your .cshrc file you can add:

   ::

      setenv PYMEMDYN /home/username/software/pymemdyn
      set path = ($path $PYMEMDYN)

   Notice that I have cloned *pymemdyn* in the software folder in my
   home folder, you will have to adapt this to wherever it is that you
   downloaded your *pymemdyn* to.

   After including the route to your *pymemdyn* directory in your
   .bashrc file make sure to issue the command:

   ::

      source .bashrc

   or open a new terminal.

   To check if you have defined the route to the *pymemdyn* directory
   correctly try to run the main program called pymemdyn in a terminal:

   ::

      pymemdyn --help

   You should obtain the following help output:

   ::

    usage: pymemdyn [-h] [-v] [-b OWN_DIR] [-r REPO_DIR] -p PDB [-l LIGAND]
                 [--lc LIGAND_CHARGE] [-w WATERS] [-i IONS] [--res RESTRAINT]
                 [-f LOOP_FILL] [-q QUEUE] [-d] [--debugFast]
                    
    == Setup Molecular Dynamics for Membrane Proteins given a PDB. ==
    
    optional arguments:
      -h, --help            show this help message and exit
      -v, --version         show program's version number and exit
      -b OWN_DIR            Working dir if different from actual dir
      -r REPO_DIR           Path to templates of fixed files. If not provided,
                            take the value from settings.TEMPLATES_DIR.
      -p PDB                Name of the PDB file to insert into membrane for MD
                            (mandatory). Use the .pdb extension. (e.g. -p
                            myprot.pdb)
      -l LIGAND, --lig LIGAND
                            Ligand identifiers of ligands present within the PDB
                            file. If multiple ligands are present, give a comma-
                            delimited list.
      --lc LIGAND_CHARGE    Charge of ligands for ligpargen (when itp file should
                            be generated). If multiple ligands are present, give a
                            comma-delimited list.
      -w WATERS, --waters WATERS
                            Water identifiers of crystalized water molecules
                            present within the PDB file.
      -i IONS, --ions IONS  Ion identifiers of crystalized ions present within the
                            PDB file.
      --res RESTRAINT       Position restraints during MD production run. Options:
                            bw (Ballesteros-Weinstein Restrained Relaxation -
                            default), ca (C-Alpha Restrained Relaxation)
      -f LOOP_FILL, --loop_fill LOOP_FILL
                            Amount of Å per AA to fill cut loops. The total
                            distance is calculated from the coordinates of the
                            remaining residues. The AA contour length is 3.4-4.0
                            Å, To allow for flexibility in the loop, 2.0 Å/AA
                            (default) is suggested. (example: -f 2.0)
      -q QUEUE, --queue QUEUE
                            Queueing system to use (slurm, pbs, pbs_ib and svgd
                            supported)
      -d, --debug
      --debugFast           run pymemdyn in debug mode with less min and eq steps.
                            Do not use for simulation results!

                            

3. Updates are very easy thanks to the git versioning system. Once
   **PyMemDyn** has been downloaded (cloned) into its own *pymemdyn*
   folder you just have to move to it and pull the newest changes:

   ::

      cd /home/username/software/pymemdyn
      git pull   

4. You can also clone older stable versions of **PyMemDyn**. For example
   the stable version 1.4 which works well and has been tested
   extensively again GROMACS version 4.6.7 can be cloned with:

   ::

      git clone https://username@github.com/GPCR-ModSim/pymemdyn.git \
      --branch stable/1.4 --single-branch pymemdyn-1.4

   Now you will have to change your .bashrc or .cshrc files in your home
   folder accordingly.

5. To make sure that your GROMACS installation is understood by
   **PyMemDyn** you will need to specify the path to where GROMACS is
   installed in your system. To do this you will need to edit the
   settings.py file with any text editor (vi and emacs are common
   options in the unix environment). Make sure that only one line is
   uncommented, looking like: GROMACS_PATH = /opt/gromacs-2021/bin
   Provided that in your case gromacs is installed in /opt. The program
   will prepend this line to the binaries names, so calling
   /opt/gromacs-2021/bin/gmx should point to that binary.

Modules
-------


Modeling Modules
~~~~~~~~~~~~~~~~

The following modules define the objects to be modeled.

-  **protein.py**. This module defines the ProteinComplex, Protein,
   Monomer, Dimer, Compound, Ligand, CrystalWaters, Ions, Cholesterol,
   Lipids, and Alosteric objects. These objects are started with the
   required files, and can then be passed to other objects.
-  **membrane.py**. Defines the cellular membrane.
-  **complex.py**. Defines the full complex, protein + membrane.
   It can include any of the previous objects.



Auxiliary Modules
~~~~~~~~~~~~~~~~~

- **checks.py**. Checks continuity of the protein and composition of the 
  residues.
- **aminoAcids.py**. Contains the amino acids class that defines the 1-letter
  and three letter codes, along with the number of different atoms per residue.
- **queue.py**.   Queue  manager.  That  is,  it  receives  objects to  be
  executed.   
- **recipes.py**.   Applies  step by  step instructions  for  carrying a 
  modeling  step.
- **bw4posres.py**. Creates a set of distance restraints based on 
  Ballesteros-Weinstein identities which are gathered by alignment to a 
  multiple-sequence alignment using clustalw.
- **utils.py**.  Puts the  functions done by the previous objects on demand.
  For example, manipulate files, copy  folders, call functions or classes from 
  standalone modules like bw4posres.py, etc.
- **settings.py** This modules sets up the main environment variables needed
  to run the calculation, for example, the path to the gromacs binaries.


Execution Modules
~~~~~~~~~~~~~~~~~

-  **gromacs.py**. Defines the Gromacs and Wrapper objects. \* Gromacs
   will load the objects to be modeled, the modeling recipe, and run it.
   \* Wrapper is a proxy for gromacs commands. When a recipe entry is
   sent to it this returns the command to be run.

Executable
~~~~~~~~~~

-  **pymemdyn** The main program to call which sends the run to a
   cluster.

 

More information about all modules can be found in the Modules chapter.



