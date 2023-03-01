import os
import shutil
from utils import create_itp


class ProteinComplex(object):
    def __init__(self, *args, **kwargs):
        self.cres = 0  # Central residue
        self.trans = [0, 0, self.cres]  # Module for translating complex
        self.n_wats = 0  # Number of experimental waters

        if "monomer" in kwargs.keys():
            self.setMonomer(kwargs["monomer"])
        if "ligand" in kwargs.keys():
            self.setLigand(kwargs["ligand"])
        if "waters" in kwargs.keys():
            self.setWaters(kwargs["waters"])
        if "ions" in kwargs.keys():
            self.setIons(kwargs["ions"])
        if "cho" in kwargs.keys():
            self.setCho(kwargs["cho"])
        if "alosteric" in kwargs.keys():
            self.setAlosteric(kwargs["alosteric"])

    def setMonomer(self, value):
        """
        Sets the monomer object.
        """
        self.monomer = value

    def getMonomer(self):
        return self.monomer
    property(getMonomer, setMonomer)

    def setLigand(self, value):
        """
        Sets the ligand object
        """
        self.ligand = value

    def getLigand(self):
        return self.ligand
    property(getLigand, setLigand)

    def setWaters(self, value):
        """
        Sets the crystal waters object
        """
        self.waters = value

    def getWaters(self):
        return self.waters
    property(getWaters, setWaters)

    def setIons(self, value):
        """
        Sets the ions object
        """
        self.ions = value

    def getIons(self):
        return self.ions
    property(getIons, setIons)

    def setCho(self, value):
        """
        Sets the cholesterol object
        """
        self.cho = value

    def getCho(self):
        return self.cho
    property(getCho, setCho)

    def setAlosteric(self, value):
        """
        Sets the alosteric object
        """
        self.alosteric = value

    def getAlosteric(self):
        return self.alosteric
    property(getAlosteric, setAlosteric)

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
        self.pdb = kwargs["pdb"]
        if not os.path.isfile(self.pdb):
            raise IOError("File '{0}' missing".format(self.pdb))

    def check_number_of_chains(self):
        """
        Determine if a PDB is a Monomer or a Dimer
        """
        chains = []
        with open(self.pdb, "r") as pdb_fp:
            for line in pdb_fp:
                if (len(line) > 21) and (
                        line.startswith(("ATOM", "TER", "HETATM"))):
                    if (line[21] != " ") and (line[21] not in chains):
                        chains.append(line[21])
    
        if len(chains) < 2:
            return Monomer(pdb = self.pdb)
        elif len(chains) >= 2:
            return Oligomer(pdb = self.pdb, chains = chains)
        

class Monomer(object):
    def __init__(self, *args, **kwargs):
        self.pdb = kwargs["pdb"]
        if not os.path.isfile(self.pdb):
            raise IOError("File '{0}' missing".format(self.pdb))

        self.group = "protlig"
        self.delete_chain()
        self.chains = []
        self._setHist()

    def delete_chain(self):
        """
        PDBs which have a chain column mess up with pdb2gmx, creating
        an unsuitable protein.itp file by naming the protein ie "Protein_A".
        Here we remove the chain value

        According to http://www.wwpdb.org/documentation/format33/sect9.html,
        the chain value is in column 22
        """
        shutil.move(self.pdb, self.pdb + "~")
        pdb = open(self.pdb + "~", "r")
        pdb_out = open(self.pdb, "w")

        replacing = False
        for line in pdb:
            new_line = line
            if len(line.split()) > 2:
                #Remove chain id
                if line[21] != " ":
                    replacing = True
                    new_line = list(line) #Transform the line into a list...
                    new_line[21] = " " 
                    new_line = "".join(new_line)
            pdb_out.write(new_line)

        if replacing: print ("Removed chain id from your protein pdb!")
        pdb.close()
        pdb_out.close()
 
        return True

    def _setHist(self):
        """
        Change Histidines in pdb to the format preferred by gromacs
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
                else:
                    tgt.write(line)
            else:
                tgt.write(line)
        tgt.close()

        return True

class Oligomer(Monomer):
    def __init__(self, *args, **kwargs):
        super(Oligomer, self).__init__(self, *args, **kwargs)

        self.chains = kwargs.get("chains")
        self.points = dict.fromkeys(self.chains, [])

    def delete_chain(self):
        """
        Overload the delete_chain method from Monomer
        """
        return True


class Sugar_prep(object):
    def __init__(self, *args, **kwargs):
        print('initialization of Sugar_prep started')
        sugars = []
        if self.ligand:
            sugars.append(self.ligand)
        if self.alosteric:
            sugars.append(self.alosteric)
            
        for sugar in sugars:
            if os.path.exists(sugar + ".ff") == True:
                pass # all files exist, so no files need to be generated
            else:
                if os.path.exists(sugar + ".itp") == False:
                    print('no .itp found, it will be generated for '+str(sugar))
                    print('current workdir= '+str(os.getcwd()))
                    create_itp(sugar + ".pdb", 0, 3)
                    # pass # TODO: LigParGen communication can be set here
                    
                Sugar_prep.lpg2pmd(self, sugar)     

   

    def lpg2pmd(self, sugar, *args, **kwargs):
        """
        Converts LigParGen structure files to PyMemDyn input files
        """
        # Safeguard for deleting files
        if not os.path.isfile(self.own_dir + "/" + sugar + ".ff"):    
            shutil.copy(self.own_dir + "/" + sugar + ".itp", self.own_dir + "/" + sugar + "_backup.itp")
            shutil.copy(self.own_dir + "/" + sugar + ".pdb", self.own_dir + "/" + sugar + "_backup.pdb")
    
            old_itp = open(self.own_dir + "/" + sugar + ".itp", "r") 
            old_pdb = open(self.own_dir + "/" + sugar + ".pdb", "r")
            
            lines_itp = old_itp.readlines()
            lines_pdb = old_pdb.readlines()
            old_itp.close()
            old_pdb.close()
                    
            new_itp = open(self.own_dir + "/" + sugar + ".itp", "w")
            new_ff = open(self.own_dir + "/" + sugar + ".ff", "w")
            new_pdb = open(self.own_dir + "/" + sugar + ".pdb", "w")
    
            split = False
            count = -1
            tmp_ff = []
            tmp_itp = []
    
            for line in lines_itp:
                if "[ defaults ]" in line:
                    # find location of defaults content
                    loc_content_of_defaults = lines_itp.index(line) + 2
                    # also remove content
                    lines_itp.pop(loc_content_of_defaults)
                    # now also do nothing else on this line
                    continue
                if sugar == self.alosteric:
                    line = line.replace("opls_8", "opls_a")                
                if "[ moleculetype ]" in line:
                    split = True    
                if split == False: 
                    if line[2:6] != "opls": 
                        new_ff.write(line)
                    else:
                        tmp_ff.append(line.split())         
                        
                if split == True:
                    count += 1
                    if count == 2:# Added lstrip() to not take starting whitespace into account
                        if sugar == self.ligand and line.lstrip()[0:3] != "LIG": 
                            line = line.replace(line.lstrip()[0:3], "LIG")
                        if sugar == self.alosteric and line.lstrip()[0:3] != "ALO":
                            line = line.replace(line.lstrip()[0:3], "ALO")
                        
                    if line[9:13] == "opls":
                        tmp_itp.append(line.split())
                        if sugar == self.ligand and line[28:31] != "LIG":
                            line = line.replace(line[28:31], "LIG")
                        if sugar == self.alosteric and line[28:31] != "ALO":
                            line = line.replace(line[28:31], "ALO")
    
                    new_itp.write(line)
    
            for i in tmp_itp:
                for j in tmp_ff:
                    if i[1] == j[0]:
                        j[1] = i[4]
                        j.insert(2, i[2])
                        new_ff.write("\t".join(j) + "\n")
    
            for line in lines_pdb:
                if line[0:4] == "ATOM":
                    if sugar == self.ligand and line[17:20] != "LIG":
                        line = line.replace(line[17:20], "LIG")
                    if sugar == self.alosteric and line[17:20] != "ALO":
                        line = line.replace(line[17:20], "ALO")
                        
                new_pdb.write(line)
            
            new_itp.close()
            new_ff.close()
            new_pdb.close()


class Compound(object):
    """
    This is a super-class to provide common functions to added compounds
    """
    def __init__(self, *args, **kwargs):
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


class Ligand(Compound):
    def __init__(self, *args, **kwargs):
        self.pdb = kwargs["pdb"]
        self.itp = kwargs["itp"]
        super(Ligand, self).__init__(self, *args, **kwargs)

        self.group = "protlig"

        self.force_field = kwargs["ff"]

        self.check_forces()

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
        self.pdb = kwargs["pdb"]
        self.itp = kwargs["itp"]
        super(CrystalWaters, self).__init__(self, *args, **kwargs)

        self.group = "wation"
        self.posre_itp = "posre_hoh.itp"
        self._setITP()
        self._n_wats = self.count_waters()
        
    def setWaters(self, value):
        """
        Set crystal waters
        """
        self._n_wats = value

    def getWaters(self):
        """
        Get the crystal waters
        """
        return self._n_wats
    number = property(getWaters, setWaters)

    def count_waters(self):
       """
       Count and set the number of crystal waters in the pdb
       """
       return len([x for x in open(self.pdb, "r") if "HOH" in x])/3

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
        self.pdb = kwargs["pdb"]
        self.itp = kwargs["itp"]
        super(Ions, self).__init__(self, *args, **kwargs)
        
        self.group = "wation"
        self.posre_itp = "posre_ion.itp"
        self._setITP()
        self._n_ions = self.count_ions()

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
       ions = ["NA", "CA", "MG", "CL", "ZN"]
       ion_count = 0
       for line in open(self.pdb, "r"):
           if len(line.split()) > 2:
               if line.split()[2] in ions:
                   ion_count += 1
       return ion_count

    def _setITP(self):
        """
        Create an itp file for this structure
        """
        s = "\n".join([
            "; position restraints for ions (resn NA, CA, MG, CL, ZN)",
            "[ position_restraints ]",
            ";  i funct       fcx        fcy        fcz",
            "   1    1       1000       1000       1000"])

        tgt = open(self.posre_itp, "w")
        tgt.writelines(s)
        tgt.close()


class Cholesterol(Compound):
    def __init__(self, *args, **kwargs):
        self.pdb = kwargs["pdb"]
        self.itp = kwargs["itp"]
        super(Cholesterol, self).__init__(self, *args, **kwargs)
        
        self.group = "membr"
        self.posre_itp = "posre_cho.itp"
        self._setITP()
        self._n_cho = self.count_cho()

    def setCho(self, value):
        """
        Sets the crystal cholesterol
        """
        self._n_cho = value

    def getCho(self):
        """
        Get the crystal cholesterols
        """
        return self._n_cho
    number = property(getCho, setCho)

    def check_pdb(self):
       """
       Check the cholesterol file meets some standards
       """
       shutil.move(self.pdb, self.pdb + "~")
       pdb = open(self.pdb + "~", "r")
       pdb_out = open(self.pdb, "w")

       replacing = False
       for line in pdb:
           new_line = line
           if len(line.split()) > 2:
               #Ensure the cholesterol is labeled as CHO
               if line.split()[3] != "CHO":
                   replacing = True
                   new_line = new_line.replace(line.split()[3], "CHO")
           pdb_out.write(new_line)

       if replacing: print ("Made some CHO replacements in cho.pdb!")
       pdb.close()
       pdb_out.close()

       return True

    def count_cho(self):
       """
       Count and set the number of cho in the pdb
       """
       cho_count = 0
       for line in open(self.pdb, "r"):
           print(line.split())
           if len(line.split()) > 3:
               if line.split()[3] in ["CHO", "CLR"]:
                   cho_count += 1
       return cho_count / 74 #Each CHO has 74 atoms

    def _setITP(self):
        """
        Create the itp to this structure
        """
        s = "\n".join([
            "; position restraints for ions (resn CHO)",
            "[ position_restraints ]",
            ";  i funct       fcx        fcy        fcz",
            "   1    1       1000       1000       1000"])

        tgt = open(self.posre_itp, "w")
        tgt.writelines(s)
        tgt.close()


class Alosteric(Compound):
    """
    This is a compound that goes as a ligand but it's placed in an alosteric
    site rather than an orthosteric one.
    """
    def __init__(self, *args, **kwargs):
        self.pdb = kwargs["pdb"]
        self.itp = kwargs["itp"]
        super(Alosteric, self).__init__(self, *args, **kwargs)

        self.check_pdb()

        self.force_field = kwargs["ff"]
        self.check_itp()

        self.group = "protlig"

    def check_pdb(self):
       """
       Check the alosteric file meets some standards
       """
       shutil.move(self.pdb, self.pdb + "~")
       pdb = open(self.pdb + "~", "r")
       pdb_out = open(self.pdb, "w")

       replacing = False
       for line in pdb:
           new_line = line          
           line = line.split()
           try:
               line[3]
               #Ensure the alosteric compound is labeled as ALO
               if line[3] != "ALO":
                   replacing = True
                   new_line = new_line.replace(line[3], "ALO")
           except: IndexError
            
           pdb_out.write(new_line)

       if replacing: print ("Made some ALO replacements in %s!" % self.pdb)
       pdb.close()
       pdb_out.close()

       return True

    def check_itp(self):
        """
        Check the force field is correct
        """
        print('alosteric itp being checked')
        shutil.move(self.itp, self.itp + "~")
        itp = open(self.itp + "~", "r")
        itp_out = open(self.itp, "w")
 
        molecule_type = atoms = False
        for line in itp:
            new_line = line
            if line.startswith("[ moleculetype ]"): molecule_type = True
            if molecule_type:
                if not line.startswith(";"): #Not a comment
                    if len(line.split()) == 2:
                        #Change the user name to "alo"
                        new_line = line.replace(line.split()[0], "alo")
                        molecule_type = False

            if line.startswith("[ "): #Next section (after atoms) reached
                atoms = False
            if line.startswith("[ atoms ]"): atoms = True
            if atoms:
                if not line.startswith(";"):
                    if len(line.split()) > 4:
                        #Change the name of the compound to "ALO"
                        new_line = line.replace(line.split()[3], "ALO")
            itp_out.write(new_line)

        itp.close()
        itp_out.close()
  
        return True
