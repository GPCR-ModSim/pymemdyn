"""This module handles the protein and all submitted molecules around it.
"""

import os
import shutil
import logging

import numpy as np

import checks

protein_logger = logging.getLogger('pymemdyn.protein')

try:
    import ligpargen.ligpargen
except:
    protein_logger.warning("""!! WARNING !! : No installation of ligpargen was found. 
          itp files for ligands cannot be automatically generated and must be provided 
          manually.""")

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from rdkit.Chem import PDBWriter
except:
    protein_logger.warning("""!! WARNING !! : No installation of rdkit was found. 
          itp files for ligands cannot be automatically generated and must be provided 
          manually.""")

class System(object):
    def __init__(self, **kwargs):
        protein_res_names = ["ALA", "ARG", "ASN", "ASP", "ASH", "CYS", "CYX", "GLN", "GLU", "GLH", "GLY", "HIS", 
                             "HIE", "HID", "HIP", "HISE", "HISD", "HISH", "ILE", "LEU", "LYS", 
                             "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",]

        self.pdb = kwargs["pdb"]

        self.logger = logging.getLogger('pymemdyn.protein.System')
        self.logger.info('Initializing system') 

        self.protein_res = set()
        self.cofactor_res = set()

        with open(self.pdb, "r") as src:
            for line in src:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    res_name = line[17:21].strip()
                    if res_name in protein_res_names:
                        self.protein_res.add(res_name)
                    else:
                        self.cofactor_res.add(res_name)

        self.logger.info(f'The PDB file contains the following cofactors: {self.cofactor_res}')
        self.logger.info(f'The PDB file contains the following residues: {self.protein_res}')

    def split_system(self, **kwargs):
        self.ligand = kwargs["ligand"]
        self.waters = kwargs["waters"]
        self.ions = kwargs["ions"]

        # Extract protein.pdb
        tgt_prot = open('protein.pdb', "w")
        
        with open(self.pdb, "r") as src:
            chain = 0
            for line in src:
                if line.startswith("TER"):
                    tgt_prot.write(line)
                    chain += 1    
                      
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    # add protein chain IDs
                    if line[21] == ' ':
                        line = line[:21] + (chr(ord('A') + chain) + line[22:])
                                            
                    # only save protein residues    
                    res_name = line[17:21].strip()
                    if res_name in self.protein_res:
                        tgt_prot.write(line) 
        
        self.logger.info(f'PDB file created containing protein(s): protein.pdb')

        # Extract cofactor pdb(s)
        cofactors = []
        if self.ligand:
            cofactors += self.ligand.split(',')
        if self.waters:
            cofactors += self.waters.split(',')
        if self.ions:
            cofactors += self.ions.split(',')

        for cofactor in cofactors:

            if cofactor == self.waters:
                filename = 'HOH'
            else:
                filename = cofactor

            if not os.path.exists(f'{filename}.pdb'):
                cf_pdb = open(f'{filename}.pdb', "w")
                with open(self.pdb, "r") as src:
                    for line in src:
                        if line.startswith("ATOM") or line.startswith("HETATM"):
                            res_name = line[17:21].strip()
                            if res_name == cofactor:
                                cf_pdb.write(line)
                
                self.logger.info(f'PDB file created containing cofactor: {filename}.pdb')
            else:
                self.logger.info(f'PDB file already exists containing cofactor: {filename}.pdb')


class ProteinComplex(object):
    def __init__(self, *args, **kwargs):
        self.cres = 0  # Central residue
        self.trans = [0, 0, self.cres]  # Module for translating complex
        self.logger = logging.getLogger('pymemdyn.protein.ProteinComplex')
        self.logger.info('Initializing protein complex.')

        for object in kwargs["objects"]:
            self.setObjects(object)
            self.logger.debug(f'set object: {object}')

    def setObjects(self, object):
        """
        Sets an object.
        """
        setattr(self, object.name, object)
    
    def getObjects(self, object):
        return getattr(self, object.name)
    property(getObjects, setObjects)
  
    def set_nanom(self):
        """
        Convert dimension measurements to nanometers for GROMACS
        """
        nanometer = 10
        self.gmx_prot_xy = self.prot_xy / nanometer
        self.gmx_prot_z = self.prot_z / nanometer
        if self.gmx_prot_z <= 15.565:
            self.gmx_emb_z = self.gmx_prot_z
        else:
            self.gmx_emb_z = 15.565


class Protein(object):
    def __init__(self, *args, **kwargs):
        """
        This is a proxy to determine if a protein is a Monomer or a Dimer
        """
        self.pdb = kwargs.get("pdb")
        self.dir = kwargs.get('owndir') or ''
        self.loop_fill = kwargs.get('loopfill') or ''

        self.logger_prot = logging.getLogger('pymemdyn.protein.Protein')
                
        if not os.path.isfile(self.pdb):
            raise IOError("File '{0}' missing".format(self.pdb))

    def check_number_of_chains(self):
        """
        Determine if a PDB is a Monomer or a Oligomer
        """
        chains = []
        with open(self.pdb, "r") as pdb_fp:
            for line in pdb_fp:
                if (len(line) > 21) and (
                        line.startswith(("ATOM", "TER", "HETATM"))):
                    if (line[21] != " ") and (line[21] not in chains):
                        chains.append(line[21])
    
        if len(chains) < 2:
            return Monomer(pdb = self.pdb, chains = chains, dir = self.dir, loopfill = self.loop_fill)
        elif len(chains) >= 2:
            return Oligomer(pdb = self.pdb, chains = chains, dir = self.dir, loopfill=self.loop_fill)

    def calculate_center(self):
        """
        Determine center of the coords in the self.pdb.
        """
        with open(self.pdb, "r") as inf:
            lines = inf.readlines()
        
        n = []
        for line in lines:
            m =  line.split()
            if m[0] in ('ATOM', 'HETATM'):
                n.append(m[:9])
        matrix = np.array(n)
        coord = matrix[:, [6, 7, 8]]
        coord = coord.astype(float)
        mean_coord = np.mean(coord, axis=0)

        return mean_coord
        

class Monomer(object):
    def __init__(self, *args, **kwargs):
        self.name = "proteins"
        self.type = "proteins"
        self.pdb = kwargs["pdb"]
        
        self.logger_monomer = logging.getLogger('pymemdyn.protein.Monomer')
        if not os.path.isfile(self.pdb):
            raise IOError("File '{0}' missing".format(self.pdb))

        self.group = "protlig"

        self.own_dir = kwargs['dir']
        self.loop_fill = kwargs['loopfill']
        self.chains = kwargs['chains']

        self.check_protein = checks.CheckProtein(
                pdb=self.pdb,
                chains=self.chains, 
                tgt='missingLoops.txt', 
                loop_fill = self.loop_fill
                )
        
        self.broken_chains = self.check_protein.make_ml_pir(work_dir=self.own_dir)
        
        if self.broken_chains != []:
            self.logger_monomer.info('Broken chains: {}'.format(self.broken_chains))
            system_file = self.pdb.replace('.pdb', '-modeller.pdb')
            system_pdbs = []
            for chain in self.chains:
                if chain in self.broken_chains:
                    self.logger_monomer.debug('MODELLER input: {}'.format(self.pdb.replace('.pdb', f'_{chain}.pdb')))
                    system_pdbs.append(self.check_protein.refine_protein(knowns=self.pdb.replace('.pdb', f'_{chain}.pdb'), chain=chain))
                else:
                    system_pdbs.append(self.pdb.replace('.pdb', f'_{chain}.pdb'))
 
            with open(system_file, 'w') as output_file:
                for file_idx, pdb_file in enumerate(system_pdbs):
                    with open(pdb_file, 'r') as file:
                        for line in file:
                            output_file.write(line)

            self.pdb = system_file

            # overwrite protein-his.pdb file
            #os.remove(self.pdb_hist)
            #os.rename(system_file, self.pdb_hist)

        self.pdb_hist = self._setRes() 

        return 

    def _setRes(self):
        """
        Change Histidines and Cysteins in pdb to the format preferred by gromacs.
        """
        tgt = open(self.pdb.replace(".pdb", "-his.pdb"), "w")
        self.pdb_his = tgt.name

        for line in open(self.pdb, "r"):
            if len(line.split()) > 3:
                if line.split()[3] == "HIE":
                    tgt.write(line.replace('HIE ','HISE'))
                elif line.split()[3] == "HID":
                    tgt.write(line.replace('HID ','HISD'))
                elif line.split()[3] == "HIP":
                    tgt.write(line.replace('HIP ','HISH'))
                elif line.split()[3] == "CYX":
                    tgt.write(line.replace('CYX ','CYS '))
                else:
                    tgt.write(line)
            else:
                tgt.write(line)
            
        tgt.close()

        
        # os.rename(self.pdb, self.pdb.replace(".pdb", "-oldHIS.pdb"))
        # os.rename(tgt.name, tgt.name.replace("-his.pdb", ".pdb"))

        return tgt.name

class Oligomer(Monomer):
    def __init__(self, *args, **kwargs):
        super(Oligomer, self).__init__(self, *args, **kwargs)

        self.chains = kwargs.get("chains")
        self.points = dict.fromkeys(self.chains, [])
        

class CalculateLigandParameters(object):
    def __init__(self):
        """
        Prepare cofactors for pymemdyn run. 

        .itp files are generated by ligpargen (if not provided).
        """
        self.logger = logging.getLogger('pymemdyn.protein.CalculateLigandParameters')
        self.logger.info('initialization of CalculateLigandParameters started')

        if self.ligand:
            ligands = self.ligand.split(',')
            self.logger.debug(f'List of ligands: {ligands}')
            
            charges = self.ligand_charge.split(',')
            if len(charges) != len(ligands) or charges == ['']:
                charges = [0] * len(ligands)
            self.logger.debug(f'List of charges: {charges}')

            index = 0
            for lig in ligands:
                charge = int(charges[index])
                index += 1
                if os.path.exists(lig + ".ff") == True:
                    self.logger.info('All ligand parameter (.itp and .ff) files present.')
                    pass

                else:
                   self.logger.info(f'Ligand parameter file {lig}.ff not found.')
                   if os.path.exists(lig + ".itp") == False:
                       self.logger.info(f'Ligand parameter file {lig}.itp not found.')
                       self.logger.info('Ligand parameter files will be generated with Ligpargen.')
                       CalculateLigandParameters.create_itp(self, lig, charge)
                CalculateLigandParameters.lpg2pmd(self, lig, index)

    def create_itp(self, ligand: str, charge: int) -> None:
        """Call ligpargen to create gromacs itp file and corresponding openmm
        pdb file. Note that original pdb file will be replaced by opnemm pdb
        file.

        :param pdbfile: string containing local path to pdb of molecule. In commandline -i.
        :param charge: interger charge of molecule. In commandline -c.
        :param numberOfOptimizations: number of optimizations done by ligpargen. In cmdline -o.

        :returns: None

        Writes itp file and new pdf file to current dir. old pdb is saved in dir ligpargenInput. 
        unneccessary ligpargen output is saved in dir ligpargenOutput.
        """
        workdir = 'ligpargenOutput_' + ligand
        if os.path.exists(workdir):
            shutil.rmtree(workdir)
            self.logger.debug(f'Removed {workdir}')        
        
        inputdir = 'ligpargenInput_' + ligand
        if os.path.exists(inputdir):
            shutil.rmtree(inputdir)
            self.logger.debug(f'Removed {inputdir}')  
      
        pdb_lpg = CalculateLigandParameters.process_pdb(self, ligand, charge)

        try: 
            import ligpargen.ligpargen as lpg
            self.logger.info(f'Calculating ligand parameters for {ligand} using LigParGen. ifile: {pdb_lpg}')
            mol = lpg.LigParGen(
                ifile=pdb_lpg, 
                molname=f'{ligand}', 
                workdir=workdir, 
                resname=ligand, 
                charge=charge, 
                numberOfOptimizations=0     # Increasing the optimizations causes the ligand to shift
                )
            self.logger.info(f'Calculated ligand parameters for {ligand} using LigParGen.')
        except:
            self.logger.warning(f'LigParGen unable to calculate ligand parameters for {ligand}.')

        mol.writeAllOuputs()

        # Only save necessary files
        os.mkdir(inputdir)
        os.rename(f'{ligand}.pdb', inputdir+'/'+f'{ligand}.pdb')
        os.rename(workdir+'/'+ligand+'.gmx.itp', ligand+'.itp')
        os.rename(workdir+'/'+ligand+'.openmm.pdb', ligand+'.pdb')

        return
   
    def process_pdb(self, ligand, charge):
        pdbfile_lpg = f'{ligand}_lpg.pdb'
        ref_molecule = Chem.MolFromPDBFile(f'{ligand}.pdb', removeHs=False)
        
        has_hydrogens = CalculateLigandParameters.molecule_has_hydrogens(ref_molecule)
        self.logger.info(f"Molecule has explicit hydrogens: {has_hydrogens}")

        if has_hydrogens == True:
            writer = PDBWriter(pdbfile_lpg)
            writer.write(ref_molecule)
            writer.close()

            return pdbfile_lpg
        
        else:
            self.logger.info(f'Charge of {ligand} was defined by the user as {charge}.')
            self.logger.warning('Hydrogens within Aromatic rings are likely to be incorrectly added. \
                                For compounds with Aromatic rings rings we recommend using explicit hydrogens. ')
            # TODO: protocol for AddHs of aromatic compounds
            # TODO: protocol for AddHs of charged compounds
            if charge == 0:
                self.logger.info(f'Adding hydrogens to {ligand}.')
                pdbfile_tmp = f'{ligand}_tmp.pdb'
                tgt_molecule = Chem.AddHs(ref_molecule)
                AllChem.ConstrainedEmbed(tgt_molecule, ref_molecule)
                writer = PDBWriter(pdbfile_tmp)
                writer.write(tgt_molecule)
                writer.close()

                pdbfile_lpg = f'{ligand}_lpg.pdb'
                pdbfile_lpg_file = open(pdbfile_lpg, "w")

                with open(pdbfile_tmp, "r") as src:
                    count=0              
                    for line in src:
                        if line.startswith("ATOM") or line.startswith("HETATM"):       
                            if count == 0:
                                count += 1
                                res_nr = line[20:27]
                            line = line[:17] + ligand + res_nr + line[27:]
                        pdbfile_lpg_file.write(line)
                pdbfile_lpg_file.close()

            else:
                self.logger.warning('PyMemDyn is not able to add hydrogens for charged molecules. Please provide a molecule file with hydrogens')

            return pdbfile_lpg

    def molecule_has_hydrogens(mol):
        """Check if the molecule has hydrogen atoms."""
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 1:  # Hydrogen has atomic number 1
                return True
        return False

    def lpg2pmd(self, cofactor, index, *args, **kwargs):
        """
        Converts LigParGen structure files to PyMemDyn input files.

        Original files are stored as something_backup.pdb or something_backup.itp.
        """
        # Safeguard for deleting files
        if not os.path.isfile(self.own_dir + "/" + cofactor + ".ff"):    
            shutil.copy(self.own_dir + "/" + cofactor + ".itp", self.own_dir + "/" + cofactor + "_backup.itp")
            shutil.copy(self.own_dir + "/" + cofactor + ".pdb", self.own_dir + "/" + cofactor + "_backup.pdb")
    
            old_itp = open(self.own_dir + "/" + cofactor + ".itp", "r") 
            old_pdb = open(self.own_dir + "/" + cofactor + ".pdb", "r")
            
            lines_itp = old_itp.readlines()
            lines_pdb = old_pdb.readlines()
            old_itp.close()
            old_pdb.close()
                    
            new_itp = open(self.own_dir + "/" + cofactor + ".itp", "w")
            new_ff = open(self.own_dir + "/" + cofactor + ".ff", "w")
            new_pdb = open(self.own_dir + "/" + cofactor + ".pdb", "w")
    
            self.lig_ID = 'L'+str(str(index).zfill(2))
            
            split = False
            count = -1
            tmp_ff = []
            tmp_itp = []
    
            for line in lines_itp:
                line = line.replace("opls_", f"{self.lig_ID}__")
                if "[ defaults ]" in line:
                    # find location of defaults content
                    loc_content_of_defaults = lines_itp.index(line) + 2
                    # also remove content
                    lines_itp.pop(loc_content_of_defaults)
                    # now also do nothing else on this line
                    continue             
                if "[ moleculetype ]" in line:
                    split = True    
                if split == False: 
                    if f"{self.lig_ID}__" not in line:
                        new_ff.write(line)
                    else:
                        tmp_ff.append(line.split())         
                        
                if split == True:
                    count += 1
                    if count == 2:# Added lstrip() to not take starting whitespace into account
                        if line.lstrip()[0:3] != self.lig_ID: 
                            line = line.replace(line.lstrip()[0:3], self.lig_ID)
                        
                    if f"{self.lig_ID}__" in line:
                        tmp_itp.append(line.split())
                        if line[28:31] != self.lig_ID:
                            line = line.replace(line[28:31], self.lig_ID)
    
                    new_itp.write(line)
    
            for i in tmp_itp:
                for j in tmp_ff:
                    if i[1] == j[0]:
                        j[1] = i[4]
                        j.insert(2, i[2])
                        new_ff.write("\t".join(j) + "\n")
    
            for line in lines_pdb:
                cols = line.split()
                if cols[0] == "HETATM":     # Maestro does this.
                    line = line.replace("HETATM", "ATOM  ")
                    cols[0] = "ATOM"
                if "ATOM" in cols[0]: 
                    if cols[3] != self.lig_ID:
                        line = line.replace(cols[3], self.lig_ID)
                if "REMARK" in cols[0]:
                    continue # Do not write remarks in new pdb file
                        
                new_pdb.write(line)
            
            new_itp.close()
            new_ff.close()
            new_pdb.close()


class Compound(object):
    """
    This is a super-class to provide common functions to added compounds
    """
    def __init__(self, *args, **kwargs):
        self.name = kwargs["name"]
        self.ID = kwargs["ID"]
        self.check_files(self.pdb, 
                         # self.itp
                         )

    def check_files(self, *files):
        """
        Check if files passed as *args exist
        """
        for src in files:
            if not os.path.isfile(src):
                raise IOError("File {0} missing".format(src))
            
    def calculate_center(self):
        """
        Determine center of the coords in the self.pdb using column-based parsing.
        """
        try:
            with open(self.pdb, "r") as inf:
                lines = inf.readlines()
        except:
            with open(self, "r") as inf:
                lines = inf.readlines()

        coords = []
        for line in lines:
            record_type = line[:6].strip()
            if record_type in ('ATOM', 'HETATM'):
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                coords.append([x, y, z])

        matrix = np.array(coords)
        mean_coord = np.mean(matrix, axis=0)
        return mean_coord

    def correct_resid(self, pdb, resid):
        """
        Correct the residue id to the specified residue id.
        """
        with open(pdb, 'r') as file:
            lines = file.readlines()

        # Modify the lines with new residue identifier
        with open(pdb, 'w') as file:
            for line in lines:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    # PDB format specifies residue identifier in columns 23-26 (inclusive)
                    modified_line = line[:17] + resid + line[20:]
                    file.write(modified_line)
                else:
                    # Write the line as is if it's not an atom or heteroatom record
                    file.write(line)

class Ligand(Compound):
    def __init__(self, *args, **kwargs):
        self.logger_lig = logging.getLogger('pymemdyn.protein.Ligand')
        self.type = "ligand"
        self.pdb = kwargs["pdb"]
        self.itp = kwargs["itp"]
        super(Ligand, self).__init__(self, *args, **kwargs)

        self.group = "protlig"

        self.force_field = kwargs["ff"]

        self.check_forces()
        try:
            self.center = self.calculate_center()
            self.logger_lig.debug(f'Center of {self.pdb} at {self.center}')
        except:
            self.logger_lig.warning('Could not calculate center of ligand. Please \
                                    check ligand alignment manually.')


    def check_forces(self):
        """
        A force field must give a set of forces which match every atom in
        the pdb file. This showed particularly important to the ligands, as they
        may vary along a very broad range of atoms
        """
        # The itp matches each residue in the ligand pdb with the force field
        atoms_def = False
        molecules = {}
        for line in open(self.itp, "r"):
            if "[ atoms ]" in line:
                atoms_def = True
            if "[ bonds ]" in line: # Assuming here that '[ bonds ]' immediately 
                atoms_def = False   # follows '[ atoms ]'.
            if atoms_def and not line.startswith(";"):
                data = line.split()
                if len(data) > 6:
                    if data[3] not in molecules.keys(): molecules[data[3]] = {}
                    #{"LIG": {"C1": "TC1"},}
                    molecules[data[3]][data[4]] = data[1]

        atoms = {}
        # The force field matches each atom in the pdb with one line
        for line in open(self.force_field, "r"):
            if not line.startswith(";"):
                if (len(line.split()) > 6):
                    #{"TC1": "C1"}
                    atoms[line.split()[0]] = line.split()[1]

        # The pdb has the name of the atom in the third position.
        # Here we cross-check all three files to match their harvested values
        for line in open(self.pdb, "r"):
            data = line.split()
            if len(data) > 6:
                if molecules[data[3]][data[2]] not in atoms.keys():
                    # Some atoms in the pdb have no definition in the parameters
                    # file lig.ff
                    # TODO : Maybe add a guessing function, although it might
                    # just be better to give a better error message stating
                    # to check consistency between the pdb file and the
                    # .ff (parameters) file
                    print ("Atom {0} has no field definition".format(data[1]))

                if atoms[molecules[data[3]][data[2]]] not in\
                    molecules[data[3]].keys():
                    print ("Atom {0} has a wrong field definition. Check .pdb \
                    and .ff files consistency".format(
                        data[1]))
                    print ("Atom names in lig.pdb")
                    print (molecules[data[3]].keys())
                    print ("Atom name in lig.ff")
                    print (atoms[molecules[data[3]][data[2]]])

        return True
    

class CrystalWaters(Compound):
    def __init__(self, *args, **kwargs):
        self.logger_cw = logging.getLogger('pymemdyn.protein.CrystalWaters')
        self.type = "waters"
        self.pdb = kwargs["pdb"]
        self.correct_resid(self.pdb, 'HOH')
        self.itp = kwargs["itp"]
        super(CrystalWaters, self).__init__(self, *args, **kwargs)

        self.group = "wation"
        self.posre_itp = f"posre_{self.ID}.itp"
        self._setITP()
        self._n_waters = self.count_waters()
        try:
            self.center = self.calculate_center()
            self.logger_cw.debug(f'Center of {self.pdb} at {self.center}')
        except:
            self.logger_cw.warning('Could not calculate center of crystal waters. Please check crystal waters alignment manually.')        
                
    def setWaters(self, value):
        """
        Set crystal waters
        """
        self._n_waters = value

    def getWaters(self):
        """
        Get the crystal waters
        """
        return self._n_waters
    number = property(getWaters, setWaters)

    def count_waters(self):
       """
       Count and set the number of crystal waters in the pdb
       """
       return int(len([x for x in open(self.pdb, "r") if self.ID in x])/3)

    def _setITP(self):
        """
        Create the itp to this structure
        """
        s = "\n".join([
            "; position restraints for crystallographic waters (resn HOH)",
            "[ position_restraints ]",
            ";  i funct       fcx        fcy        fcz",
            "   1    1       1000       1000       1000"])

        tgt = open(self.posre_itp, "w")
        tgt.writelines(s)
        tgt.close()


class Ions(Compound):
    def __init__(self, *args, **kwargs):      
        self.logger_ions = logging.getLogger('pymemdyn.protein.Ions')
        self.type = "ions"
        self.pdb = kwargs["pdb"]
        self.itp = kwargs["itp"]
        super(Ions, self).__init__(self, *args, **kwargs)
        
        self.group = "wation"
        self.posre_itp = f"posre_{self.ID}.itp"
        self._setITP()
        self._n_ions = self.count_ions()
        try:
            self.center = self.calculate_center()
            self.logger_ions.debug(f'Center of {self.pdb} at {self.center}')
        except:
            self.logger_ions.warning('Could not calculate center of ions. Please \
                                    check ions alignment manually.')

    def setIons(self, value):
        """
        Sets the crystal ions
        """
        self._n_ions = value

    def getIons(self):
        """
        Get the crystal ions
        """
        return self._n_ions
    number = property(getIons, setIons)

    def count_ions(self):
        """
        Count and set the number of ions in the pdb
        """
        ions = ["NA", "NA+", "SOD", "K", "CA", "MG", "CL", "CL-", 
                "CHL", "LI", "RB", "CS", "F", "BR", "I"]
        ion_count = 0
        for line in open(self.pdb, "r"):
           if len(line.split()) > 2:
               if line.split()[2] in ions:
                   ion_count += 1
        if ion_count == 0:
            self.logger_ions.warning(f'Ion identifier ({self.ID}) not in {ions}. Please make sure the ID links correctly to ions.itp.')
        
        return ion_count

    def _setITP(self):
        """
        Create an itp file for this structure
        """
        s = "\n".join([
            "; position restraints for ions (resn NA, K, CA, MG, CL, ZN)",
            "[ position_restraints ]",
            ";  i funct       fcx        fcy        fcz",
            "   1    1       1000       1000       1000"])

        tgt = open(self.posre_itp, "w")
        tgt.writelines(s)
        tgt.close()
