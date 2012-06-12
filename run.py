#!/usr/bin/env python
import argparse
import logging
import os
import shutil
import sys
import textwrap

import complex
import gromacs
import membrane
import protein
import queue
import recipes
import settings

class Run(object):
    #This is a dummy
    def __init__(self, *args, **kwargs):
        self.own_dir = kwargs["own_dir"] or ""
        self.repo_dir = kwargs["repo_dir"] or ""
        self.pdb = kwargs["pdb"] or ""
        self.ligand = kwargs["ligand"] or ""
        self.alosteric = kwargs["alosteric"] or ""
        self.waters = kwargs["waters"] or False
        self.ions = kwargs["ions"] or False
        self.cho = kwargs["cho"] or False
        self.queue = kwargs["queue"] or ""
        self.debug = kwargs["debug"] or False

        if self.pdb:
            self.pdb = protein.Monomer(pdb = self.pdb)
        if self.ligand:
            self.ligand = protein.Ligand(pdb = self.ligand + ".pdb",
                                         itp = self.ligand + ".itp",
                                         ff = self.ligand + ".ff")

        if self.alosteric:
            self.alosteric = protein.Alosteric(pdb = self.alosteric + ".pdb",
                                               itp = self.alosteric + ".itp",
                                               ff = self.alosteric + ".ff")

        if self.waters:
            self.waters = protein.CrystalWaters()

        if self.ions:
            self.ions = protein.Ions()

        if self.cho:
            self.cho = protein.Cholesterol()

        self.membr = membrane.Membrane()

        prot_complex = protein.ProteinComplex(
            monomer = self.pdb,
            ligand = self.ligand or None,
            alosteric = self.alosteric or None,
            waters = self.waters or None,
            ions = self.ions or None,
            cho = self.cho or None)
        full_complex = complex.MembraneComplex()
        full_complex.complex = prot_complex
        full_complex.membrane = self.membr

        self.g = gromacs.Gromacs(membrane_complex = full_complex)

        # Note that if not provided in command line, self.queue is set in
        # settings.py
        if self.queue:
            if self.queue == "slurm":
                my_queue = queue.Slurm()
            elif self.queue == "pbs":
                my_queue = queue.PBS()
        else:
            my_queue = queue.NoQueue()

        self.g.queue = my_queue

    def clean(self):
        '''Removes all previously generated files'''
        to_unlink = ["#index.ndx.1#", "#ligand_ha.ndx.1#", "#mdout.mdp.1#",
            "#mdout.mdp.2#", "#mdout.mdp.3#", "#mdout.mdp.4#", "#mdout.mdp.5#",
            "#mdout.mdp.6#", "#mdout.mdp.7#", "#mdout.mdp.8#", "#mdout.mdp.9#",
            "#output.pdb.1#", "#output.tpr.1#", "#popc.pdb.1#",
            "#posre.itp.1#", "#proteinopls.pdb.1#", "#proteinopls.pdb.2#",
            "#proteinopls.pdb.3#", "#proteinopls.pdb.4#","#protein.top.1#",
            "#protpopc.pdb.1#", "#protpopc.pdb.2#", "#tmp.pdb.1#",
            "#topol.top.1#", "#topol.top.2#", "#topol.top.3#",
            "#topol.tpr.1#", "#topol.tpr.2#", "#topol.tpr.3#",
            "#topol.tpr.4#", "#water.pdb.1#", "ener_EQ.edr", 
            "ffoplsaabon_mod.itp", "ffoplsaa_mod.itp", "ffoplsaanb_mod.itp",
            "genion.log", "hexagon.pdb", "index.ndx",
            "ligand_ha.ndx", "mdout.mdp", "min.pdb",
            "output.pdb", "output.tpr", "popc.pdb", "popc.itp", "posre.itp",
            "posre_lig.itp", "protein.itp", "protein.top",
            "protein_ca200.itp", "proteinopls.pdb", "proteinopls-ligand.pdb",
            "protpopc.pdb", "steep.mdp", "traj.xtc", "traj_EQ.xtc", "tmp.pdb",
            "topol.top", "topol.tpr", "tmp_proteinopls.pdb", "water.pdb"]

        dirs_to_unlink = ["Rmin", "eq", "eqCA"]

        for target in to_unlink:
            if os.path.isfile(target): os.unlink(target)

        for target in dirs_to_unlink:
            if os.path.isdir(target): shutil.rmtree(target)
    
        return True

    def moldyn(self):
        '''Runs all the dynamics'''

        steps = ["Init", "Minimization", "Equilibration", "Relax"]

        for step in steps:
            self.g.select_recipe(stage = step)
            self.g.run_recipe(debug = self.debug)
        
        self.g.recipe = recipes.CAEquilibrate(debug = self.debug)
        self.g.run_recipe()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description = textwrap.dedent('''\
    == This script runs a Molecular Dynamic with a PDB. ==
    '''))

    parser.add_argument('-b',
                        dest = "own_dir",
                        help = "Working dir if different from actual dir",
                        default = os.getcwd())
    parser.add_argument('-r',
                        dest = "repo_dir",
                        help = "Path to templates of fixed files. If not \
                                provided, take the value from \
                                settings.REPO_DIR.",
                        default = settings.REPO_DIR)
    parser.add_argument('-p',
                        dest = "pdb",
                        required = True,
                        help = "Name of the pdb to insert into MD (mandatory)")
    parser.add_argument('-l',
        dest = "ligand",
        help = "Name of the ligand, without extension. Three \
                files must be present along with the molecule \
                pdb: the ligand, its itp and its force field.")
    parser.add_argument("--alo",
        dest = "alosteric",
        help = "Name of the alosteric interaction, without extension. \
                Three files must be present along with the molecule \
                pdb: the alosteric, its itp and its force field.")
    parser.add_argument('--waters',
                        action="store_true",
                        help = "Crystalized water molecules hoh.pdb file \
                                must exist.")
    parser.add_argument('--ions',
                        action="store_true",
                        help = "Crystalized ions ions_local.pdb and \
                                ions_loca.itp file must exist.")
    parser.add_argument('--cho',
                        action="store_true",
                        help = "Crystalized cholesterol molecules cho.pdb \
                                file must exist.")
    parser.add_argument('-q',
                        dest = "queue",
                        help = "Queue system to use (slurm and pbs supported)",
                        default = settings.QUEUE)
    parser.add_argument('--debug',
                        action="store_true")
    args = parser.parse_args()

    run = Run(own_dir = args.own_dir,
              repo_dir = args.repo_dir,
              pdb = args.pdb,
              ligand = args.ligand,
              alosteric = args.alosteric,
              waters = args.waters,
              ions = args.ions,
              cho = args.cho,
              queue = args.queue,
              debug = args.debug)
    run.clean()

    f = open("GROMACS.log", "w")
    f.close()

    if args.debug:
        logging.basicConfig(filename='GROMACS.log', level=logging.DEBUG)
    else:
        logging.basicConfig(filename='GROMACS.log',
                            format='%(asctime)s %(message)s',
                            datefmt='%m/%d/%Y %I:%M:%S',
                            level=logging.DEBUG)
    #
    run.moldyn()

