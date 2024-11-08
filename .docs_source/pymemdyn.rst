Pymemdyn 
========

This is the main script for the pymemdyn commandline tool.
In this script the following things are accomplished:

        1. Command line arguments are parsed.
        2. (If necessary) a working directory is created.
        3. (If necessary) previous Run files are removed.
        4. A run is done. 


Usage
-----

.. code-block:: bash

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
