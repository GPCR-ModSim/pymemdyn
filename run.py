import shutil
import os
import numpy as np

import complex
import gromacs
import membrane
import protein
import queue
import settings
import forcefield_defaults

import logging
run_logger = logging.getLogger('pymemdyn.run')

class Run(object):
    def __init__(self, pdb, *args, **kwargs):
        """
        A molecular dynamics *Run()*  MUST be given a *pdb* file.

        This class tries to initialize a full complex to send to simulation.
        Given a set of molecules (protein, ligand, other ligand, waters, ...),
        this class would try to build a full embedded-in-membrane complex.

        The complex is stored in self.g (a *Gromacs* object), and thus
        can be **run** through g.recipe and g.run_recipe procedure. See
        gromacs.py for more information.

        The queueing system is also created here to be used in certain steps.
        """
        self.pdb = pdb
        self.own_dir = kwargs.get("own_dir") or ""
        self.repo_dir = kwargs.get("repo_dir") or settings.TEMPLATES_DIR
        self.ligand = kwargs.get("ligand") or ""
        self.ligand_charge = kwargs.get("ligand_charge") or ""
        self.protein = 'protein.pdb'
        self.waters = kwargs.get("waters") or ""
        self.ions = kwargs.get("ions") or ""
        self.full_relax = kwargs.get("full_relax")
        self.restraint = kwargs.get("restraint")
        self.loop_fill = kwargs.get("loop_fill")
        self.queue = kwargs.get("queue") or ""
        self.debug = kwargs.get("debug") or False
        self.debugFast = kwargs.get("debugFast") # obsolete
        self.proteinforcefield = kwargs.get("proteinforcefield")
        self.waterforcefield = kwargs.get("waterforcefield")
        self.lipidforcefield = kwargs.get("lipidforcefield")
        self.ligandforcefield = kwargs.get("ligandforcefield")
        
        self.logger = logging.getLogger('pymemdyn.run.Run')
  
        self.logger.debug('Run arguments initialized')
        self.logger.debug(f'self.pdb = {str(self.pdb)}')
        self.logger.debug(f'self.ligand = {str(self.ligand)}')
        self.logger.debug(f'self.ligand_charge = {str(self.ligand_charge)}')
        self.logger.debug(f'self.waters = {str(self.waters)}')
        self.logger.debug(f'self.ions = {str(self.ions)}')
        self.logger.debug(f'self.full_relax = {str(self.full_relax)}')
        self.logger.debug(f'self.restraint = {str(self.restraint)}')

        self.clean()

        # Prepare System
        self.logger.debug(f'{self.pdb} dectected.')                            
        self.logger.debug('Splitting pdb file.')
        protein.System(pdb=self.pdb).split_system(ligand=self.ligand,
                                                    waters=self.waters,
                                                    ions=self.ions)

        # Prepare protein
        self.logger.debug('Checking # of chains of '+ str(self.protein))
        self.proteins = protein.Protein(pdb=self.protein, owndir=self.own_dir, loopfill=self.loop_fill, proteinforcefield=self.proteinforcefield).check_number_of_chains()
        self.logger.debug(f'Protein is a(n) {type(self.proteins)} with {self.proteins.chains} chains')
                
        try:
            self.protein_center = protein.Protein(pdb=self.protein).calculate_center()
            self.logger.info(f'Center of protein at {self.protein_center}')
        except:
            self.logger.warning("Cannot calculate center of protein. Please check alignment manually.")

        # Prepare ligand(s)
        protein.CalculateLigandParameters.__init__(self)

        # Prepare cofactor(s)
        self.cofactors = self.ligand.split(',') + ['HOH' if self.waters else "", self.ions]
        self.cofactors = [value for value in self.cofactors if value] 
        self.logger.info(f'cofactors: {self.cofactors}')
        for index, cofactor in enumerate(self.cofactors):
            ID = cofactor
            if not cofactor:
                continue
            elif cofactor in self.ligand.split(','):
                ID = 'L'+str(index+1).zfill(2)
                cofactor_type = 'Ligand'
            elif cofactor == 'HOH':
                cofactor_type = 'CrystalWaters'
            elif cofactor == self.ions:
                cofactor_type = 'Ions'
 
            setattr(self, cofactor, getattr(protein, cofactor_type)(
                        name=cofactor,
                        ID=ID,
                        pdb=f'{cofactor}.pdb',
                        itp=f'{cofactor}.itp',
                        ff=f'{cofactor}.ff'
                    ))
            
            self.logger.info(f'Checking distance between {cofactor} and protein')

            try:
                center = protein.Compound.calculate_center(f'{cofactor}.pdb')
                self.check_dist(center, self.protein_center)
            except:
                self.logger.warning(f"Cannot check distance between protein and {cofactor}. Please check alignment manually.")
        
        # Prepare membrane
        self.logger = logging.getLogger('pymemdyn.run.Run')
        self.membr = membrane.Membrane(repo_dir=self.repo_dir,
                                        lipidforcefield=self.lipidforcefield)

        # Prepare complex
        self.objects = self.cofactors + ['proteins']
        self.logger.debug(f'objects: {self.objects}')

        prot_complex = protein.ProteinComplex(
            objects= [getattr(self, object) for object in self.objects]
            )

        full_complex = complex.MembraneComplex()
        full_complex.complex = prot_complex
        full_complex.complex.cofactors = self.cofactors
        
        # Check for correct storing of Ligands
        self.logger.debug(f'attributes MembraneComplex.complex: {vars(full_complex.complex)}')
        lig_present = any(isinstance(var, protein.Ligand) for var in vars(full_complex.complex).values())
        self.logger.debug(f'found protein.Ligand: {lig_present}')

        full_complex.membrane = self.membr

        self.g = gromacs.Gromacs(membrane_complex=full_complex,
                                  repo_dir=self.repo_dir,
                                  proteinforcefield=self.proteinforcefield,
                                  waterforcefield=self.waterforcefield,
                                  lipidforcefield=self.lipidforcefield,
                                  ligandforcefield=self.ligandforcefield)

        # NOTE: If not provided in command line, self.queue is set to
        # NoQueue
        if self.queue:
            if self.queue == "slurm":
                my_queue = queue.Slurm()
            elif self.queue == "pbs":
                my_queue = queue.PBS()
            elif self.queue == "pbs_ib":
                my_queue = queue.PBS_IB()
            elif self.queue == "svgd":
                my_queue = queue.Svgd()
        else:
            my_queue = queue.NoQueue()

        self.g.queue = my_queue

    def clean(self):
        """
        Removes all previously generated files
        """
        dirs_to_unlink = []
        files_to_unlink = []

        self.logger.debug(f'current files/dirs: {os.listdir(self.own_dir)}')


        template_files = os.listdir(self.repo_dir)
        generated_files = ["ener_EQ.edr", "ffoplsaanb_mod.itp", "GROMACS.log", "hexagon.pdb", "index.ndx", 
                            "ligand_ha.ndx", "MD_output.tgz", "mdout.mdp", "mdrun.sh", "min.pdb", 
                            "output.pdb", "popc.gro", "popc.pdb", "posre.itp", "posre_lig.itp", 
                            "posre_alo.itp","pressure.log", "pressure.xvg", "pressure2.log", 
                            "pressure2.xvg", "protein.pdb", "protein.itp", "protein.top", 
                            "protein-his.pdb", "protein-modeller-his.pdb", "protein-modeller.pdb",
                            "protein_ca200.itp", "protein_gmx.fasta", "protein_gmx.pdb", 
                            "protein_gmx_bw.aln", "protein_gmx_CA.pdb", "protpopc.pdb",
                            "rmsd-all-atom-vs-start.xvg", "rmsd-backbone-vs-start.xvg", 
                            "rmsd-calpha-vs-start.xvg", "rmsf-per-residue.xvg", "temp.log", "temp.xvg", 
                            "temp2.log", "temp2.xvg", "tmp.pdb", "tmp_protein_gmx.pdb",  "topol.tpr", 
                            "tot_ener.log", "tot_ener.xvg", "tot_ener2.log", "tot_ener2.xvg", 
                            "traj_EQ.xtc", "traj_pymol.xtc", "volume.log", "volume.xvg", 
                            "volume2.log", "volume2.xvg", "water.gro", "water.pdb"]

        for entry in os.listdir(self.own_dir):
            # Safeguard MD_output.tgz
            if entry == "MD_output.tgz":
                counter = 0
                backup_file = f'MD_output_backup_{counter}.tgz'
                while os.path.exists(backup_file):
                    counter += 1
                    backup_file = f'MD_output_backup_{counter}.tgz'
                os.rename(entry, backup_file)
                self.logger.debug(f'Backup: {entry} to {backup_file}')

            if entry.startswith("ligpargenInput_") or entry.startswith("ligpargenOutput_"):
                dirs_to_unlink.append(entry)

            if entry in ["Rmin", "eq", "eqProd", "finalOutput"]:
                dirs_to_unlink.append(entry)

            if entry.startswith("#") and entry.endswith("#"):
                files_to_unlink.append(entry)

            if entry in template_files:
                files_to_unlink.append(entry)

            if entry in generated_files:
                files_to_unlink.append(entry)

            for lig in self.ligand.split(','):
                if entry in [f'{lig}_backup.itp', f'{lig}_backup.pdb',
                             #f'{lig}.ff', f'{lig}.itp', f'{lig}.pdb',         # keep {lig}.ff, {lig}.itp, {lig}.pdb to save time with re-running LigParGen
                             f'{lig}_lpg.pdb', f'ligand_{lig}_ha.ndx', f'posre_{lig}.itp']:
                    files_to_unlink.append(entry)

            if entry in [f'{self.waters}.pdb', f'posre_{self.waters}.itp']:
                files_to_unlink.append(entry)

            if entry in [f'{self.ions}.pdb', f'posre_{self.ions}.itp']:
                files_to_unlink.append(entry)
            
            if (entry.startswith("alignment_") and entry.endswith(".pir")) or \
               (entry.startswith("posre_Protein_chain_") and entry.endswith(".itp")) or \
               (entry.startswith("protein_") and entry.endswith(".pdb")) or \
               (entry.startswith("protein_Protein_chain_") and entry.endswith(".itp")) or \
               (entry.startswith("refined_") and entry.endswith(".B99990001.pdb")) or \
               (entry.startswith("refined_") and entry.endswith(".D00000001")) or \
               (entry.startswith("refined_") and entry.endswith(".ini")) or \
               (entry.startswith("refined_") and entry.endswith(".rsr")) or \
               (entry.startswith("refined_") and entry.endswith(".sch")) or \
               (entry.startswith("refined_") and entry.endswith(".V99990001")) or \
               (entry.startswith("seq_original_") and entry.endswith(".fasta")):  
               files_to_unlink.append(entry)

            # Safeguard self.pdb
            if self.pdb in files_to_unlink:
                files_to_unlink.remove(self.pdb)

        self.logger.debug(f'cleaned dirs: {sorted(dirs_to_unlink)}')
        self.logger.debug(f'cleaned files: {sorted(files_to_unlink)}')

        for target in dirs_to_unlink:
            if os.path.isdir(target): shutil.rmtree(target)

        for target in files_to_unlink:
            if os.path.isfile(target): os.unlink(target)

        return True

    def moldyn(self):
        """
        Run all steps in a molecular dynamics simulation of a membrane protein
        """

        # pdb2gmx looks for <ff_name>.ff in the run folder first before checking on gromacs/share/top
        # used for oplsaam forcefield (oplsaam.ff) and amber14sb (amber14sb.ff), which are not present in the gromacs2023 version
        ff_name = f"{self.proteinforcefield}.ff"
        ff_dir = os.path.join(forcefield_defaults.get_forcefield_defaults(self.proteinforcefield)["templates_dir"], ff_name)
        if os.path.isdir(ff_dir) and not os.path.isdir(ff_name):
            shutil.copytree(ff_dir, ff_name)

        if self.restraint == "bw":
            steps = ["Init", "Minimization", "Equilibration", "Relax", 
                    "BWRelax", "BWCollectResults"]
        elif self.restraint == "ca":
            steps = ["Init", "Minimization", "Equilibration", "Relax", 
                    "CARelax", "CACollectResults"]

        self.logger.debug(f'steps: {steps}')

        for step in steps:
            self.logger.info('\n\n[{}/{}]: {}\n'.format(steps.index(step)+1, len(steps), step))
            self.g.select_recipe(stage=step, full_relax=self.full_relax, debugFast=self.debugFast)
            self.g.run_recipe(debugFast=self.debugFast)

    def check_dist(self, vector1, vector2):
        """Check distance between protein and all possible ligands (if any).
        Raise warning is dist > 50
        """
        d = np.linalg.norm(vector1 - vector2)
        if d > 50:
            self.logger.warning(f'Center of cofactor and protein are unusually far apart ({d}). Did you correctly align both?')
        else:
            self.logger.info(f'Distance between cofactor and protein center is {d}')
                                
        return True

