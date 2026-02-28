CONGRATULATIONS!!
=================
You have  performed an  MD equilibration  of your  receptor, including
lipids, water molecules and counterions. For more  details on the methods  
followed please take some time to read
reference [1].


The performed equilibration includes the following stages:
----------------------------------------------------------

|   STAGE    | RESTRAINED ATOMS        | FORCE CONSTANT       | TIME           |
|:----------:|:-----------------------:|:--------------------:|:--------------:|
|  -         |   -                     |kJ/(mol�nm^2)        | ns             |
|Minimization|   -                     | -                    |(Max. 500 steps)|
|Equil. 1    |Protein Heavy Atoms      | 1000                 | 0.5            |
|Equil. 2    |Protein Heavy Atoms      | 800                  | 0.5            |
|Equil. 3    |Protein Heavy Atoms      | 600                  | 0.5            |
|Equil. 4    |Protein Heavy Atoms      | 400                  | 0.5            |
|Equil. 5    |Protein Heavy Atoms      | 200                  | 0.5            |
|Equil. 6    |Venkatakrishnan Pairs /  | 200 /                | 2.5 /          |
|            |C-alpha Atoms            | 200                  | 2.5            |

In this folder you will find several files related to this simulation:


INPUT:
------
    - topol.tpr             # Input for the first equilibration stage
    - processed.top         # Topology of the system
    - index.ndx             # Index file with appropriate groups for GROMACS
    - eqCA.mdp              # Example of a parameter file to configure a production run with c-alpha position restraints (see TIPS)
    - dres.mdp              # Example of a parameter file to configure a production run with BW distance restraints (see TIPS)


STRUCTURES:
-----------
    - hexagon.pdb           # Initial structure of the system, with the receptor centered in the box 
    - confout200.pdb        # Pre C-alpha / BW-relaxed structure of the system
    - confout.pdb           # Final structure of the system
    - confout.gro           # Final structure of the system (gromacs file, see TIPS)
    - load_gpcr.pml         # Loads the initial structure and the trajectory in pymol

In the "logs" subfolder, you will find the log files of mdrun:  
    - confout{force_constant}.gro  # Final structure of the system after {force_constant} equilibration step

TRAJECTORY FILES:
-----------------
    - traj_pymol.xtc        # Trajectory of the whole system for visualization in pymol. 1 snapshot/100 ps
    - traj_EQ.xtc           # Trajectory of the whole system in .xtc format: 1 snapshot/50 ps 
    - ener_EQ.edr           # Energy file of the trajectory
    - load_gpcr.pml         # Script to load the equilibration trajectory in pymol.


REPORTS:
--------
In the "reports" subfolder, you will find the following files:  
    - tot_ener.xvg, tot_ener.log    # System total energy plot and log
    - temp.xvg, temp.log            # System temperature plot and log
    - pressure.xvg, pressure.log    # System pressure plot and log
    - volume.xvg, volume.log        # System volume plot and log
    - rmsd-all-atom-vs-start	    # All atoms RMSD plot
    - rmsd-backbone-vs-start.xvg    # Backbone RMSD plot
    - rmsd-calpha-vs-start.xvg      # C-Alpha RMSD plot
    - rmsf-per-residue.xvg          # Residue RMSF plot

LOGS:
-----
    - log.log                       # log of the PyMemDyn run

In the "logs" subfolder, you will find the log files of mdrun:  
    - eq_{force_constant}.log       # log of stages with restrained heavy atoms of the target
    - md_eqCA.log                   # log of the stage with restrained C-alpha atoms of the target
    - md_eqBW.log                   # log of the stage with distance restrained BW pairs of the target

**NOTE ON GROMACS METHODS**
To integrate  the equations of  motion we have selected  the leap-frog
integrator with  a 2 femtosecond timestep.   Longe-range electrostatic
interactions  in periodic  boundary  conditions are  treated with  the
particle mesh  Ewald method.  We  use a Nose-Hoover thermostat  with a
tau_t of 0.5 picoseconds and  a Parinello-Rahman barostat with a tau_p
of 2.0.   The pressure  coupling is  semiisotropic, meaning  that it's
isotropic in the x and y  directions but different in the z direction.
Since  we are  using  pressure coupling  we are  working  with an  NPT
ensemble. This  is done both in  the all-atom restrained steps  and in
the alpha-carbon atom restrained part.   All of these details are more
explicitly stated in the Rodriguez et al. [1] publication.


**TIPS**  

NOTE: these tips work for GROMACS version 2021. For later versions, 
adjustments may be required, but the principle remains the same.

- If you want to configure a .tpr input file for a **production** run, you
can use the template 'prod.mdp' file by introducing the number of 
steps (nsteps) and/or the time between steps (dt), you want to run.  

After that, you just have to type:  
    gmx grompp -f prod.mdp -c confout.gro -p prod.top -n index.ndx -o topol_prod.tpr

    gmx mdrun -s topol_prod.tpr -o traj.trr -e ener.edr -c confout.gro -g prod.log -x traj_prod.xtc

- If you did not execute the **full_relax** protocol to relaxate your system you can run the following commands. If you wish to adjust the relaxation protocol, you can modify the template 'eqCA.mdp' (with C-alpha position restraints) or 'dres.mdp' (with Venkatakrishnan pairs distance restraints). It is highly recommended to first do a (modified) full_relax protocol prior to a **production** run.

With C-alpha position restraints:
    gmx grompp -f eqCA.mdp -c logs/confout200.gro -r logs/confout200.gro -p processed.top -n index.ndx -o topol_relax.tpr -maxwarn 1
With BW distance restraints:
    gmx grompp -f dres.mdp -c logs/confout200.gro -p processed.top -n index.ndx -o topol_relax.tpr -maxwarn 1

    gmx mdrun -s topol_relax.tpr -o traj.trr -e ener.edr -c logs/confout200.gro -g relax.log -x traj_relax.xtc

- If  you  want  to  create  a  PDB file  of  your  system  after  the
equilibration, with the receptor centered in the box, type:  

    echo 1 0 | gmx trjconv -pbc mol -center -ur compact -f confout.gro -s topol.tpr -o confout.pdb

- If you want to create an xmgrace graph of the root mean square
  deviation for c-alpha atoms in the 5.0 ns of simulation you can use:  

    echo 3 3 | gmx rms -f traj_EQ.xtc -s topol.tpr -o rmsd-calpha-vs-start.xvg

- You may want to get a pdb file of your last frame. You can first
check the total time of your trajectory and then use this time to
request the last frame with:

    gmx check -f traj_pymol.xtc
    echo 1 | gmx trjconv -b 5000 -e 5000 -f traj_pymol.xtc -o last51.pdb


References
----------

[1] Rodr�guez D., Pi�eiro �. and Guti�rrez-de-Ter�n H.   
Molecular Dynamics Simulations Reveal Insights into Key Structural Elements of Adenosine Receptors   
Biochemistry (2011), 50, 4194-208.   
