import os
import shutil
from string import Template

import bw4posres
import protein
import settings
import forcefield_defaults
import logging
logger_utils = logging.getLogger('pymemdyn.utils')


def clean_pdb(src = [], tgt = []):
    """
    Remove incorrectly allocated atom identifiers in pdb file
    """
    source = open(src, "r")
    lines_src = source.readlines() 
    source.close()
    
    target = open(tgt, "w")
    
    count = -1
    for line in lines_src:
        if count <= 3:
            count += 1
            target.write(line)
        else:
            if len(line) >= 8:
                line = line.replace(line[-3:], "\n")
                target.write(line)
            else:
                target.write(line)
    
    target.close()

    return True    
    
def clean_topol(src = [], tgt = []):
    """
    Clean the src topol of path specifics, and paste results in target
    """
    source = open(src, "r")
    target = open(tgt, "w")

    for line in source:
        newline = line
        if line.startswith("#include"):
            newline = line.split()[0] + ' "'
            newline += os.path.split(line.split()[1][1:-1])[1]
            newline += '"\n'
        target.write(newline)

    target.close()
    source.close()

    return True

def concat(**kwargs):
    """
    Make a whole pdb file with all the pdb provided
    """
    for cofactor in kwargs["tgt"].cofactors:
        logger_utils.debug(f'concatting {cofactor} to {kwargs["src"]}')
        _file_append(kwargs["src"], getattr(kwargs["tgt"], cofactor).pdb)


def getbw(**kwargs):
    """
    Call the Ballesteros-Weistein based pair-distance restraint
    module.
    """
    logger_getBW = logging.getLogger('pymemdyn.utils.getbw')
    bw4posres.Run(kwargs["src"]).pdb2fas()
    logger_getBW.debug('fasta file generated')
    bw4posres.Run(kwargs["src"]).clustalalign()
    logger_getBW.debug('clustalalign done')
    bw4posres.Run(kwargs["src"]).getcalphas()
    logger_getBW.debug('calphas gotten')
    bw4posres.Run(kwargs["src"]).makedisre()
    logger_getBW.debug('disre.itp written')

def _file_append(f_src, f2a):
    """
    Add (concatenate) a f2a pdb file to another src pdb file
    """
    src = open(f_src, "r")
    f2a = open(f2a, "r")
    tgt = open("tmp_" + f_src, "w")

    for line in src:
        if ("TER" or "ENDMDL") not in line:
            tgt.write(line)
        else:
            for line_2_add in f2a:
                tgt.write(line_2_add)
            break
    tgt.write("TER\nENDMDL\n")
    tgt.close()
    f2a.close()
    src.close()

    shutil.copy(tgt.name, f_src)

    return True

def make_cat(dir1, dir2, name):
    """
    Very tight function to make a list of files to inject
    in some GROMACS suite programs
    """
    traj_src = [os.path.join(dir1, name)]
    traj_src.extend([os.path.join(dir1, "{0}", name).format(x)
                     for x in range(800, 0, -200)])
    if os.path.isdir(dir2):
        traj_src.extend([os.path.join(dir2, name)])

    return traj_src

def make_ff(complex = None, proteinff = None):
    """
    Join all atomtypes needed to run the simulation.
    """
    base = forcefield_defaults.get_forcefield_defaults(proteinff)["ff_nb"]
    logger_utils.debug(f'base ff: {base}')

    to_concat = [base]
    for key, value in vars(complex).items():
        if isinstance(value, protein.Ligand):
            to_concat.extend([getattr(complex, key).force_field])

    output = "[ atomtypes ]\n"
    logger_utils.debug(f'files to concat: {to_concat}')
    for ff_i in to_concat:
        output += open(ff_i).read()
        if not output.endswith("\n"): output += "\n"

    open("ffnonbonded.itp", "w").write(output)
    return True


########################################################################
'''
    base = "{0}base.itp".format(ff) # This is the ff for proteins and other
    lip = "{0}lip.itp".format(ff)   # This for the lipids
    cho = "{0}cho.itp".format(ff)   # This for cholesterol

    to_concat = []
    if hasattr(complex, "ligand"):
        if hasattr(complex.ligand, "force_field"):
            to_concat.append(complex.ligand.force_field)
    if hasattr(complex, "allosteric"):
        if hasattr(complex.allosteric, "force_field"):
            to_concat.append(complex.allosteric.force_field)
    if hasattr(complex, "cho"):
        to_concat.append(cho)

    to_concat.extend([lip, base])

    output = "[ atomtypes ]\n"
    for ff_i in to_concat:
        output += open(ff_i).read()
        if not output.endswith("\n"): output += "\n"

    open("ffoplsaanb_mod.itp", "w").write(output)

    return True
'''

def make_topol(template_dir = None,
    target_dir = "",  # Dir where topol.top should land
    working_dir = "", # Dir where script is working
    complex = None):  # The MembraneComplex object to deal
    if template_dir is None:
        template_dir = settings.TEMPLATES_DIR
    """
    Make the topol starting from our topol.top template
    """
    if working_dir: working_dir += "/" # Root dir doesn't need to be slashed

    src = open(os.path.join(template_dir, "topol.top"), "r")
    tgt = open(os.path.join(target_dir, "topol.top"), "w")

    t = Template("".join(src.readlines()))
    src.close()
    
    itp_include = []  
    mol_include = []
    for key, value in vars(complex).items():
        # define protein topologies first
        if isinstance(value, protein.Oligomer):
            chainID = getattr(complex, key).chains
            for ID in chainID:
                itp_include.extend([f'#include "protein_Protein_chain_{ID}.itp"',
                                    '#ifdef POSRES',
                                    f'#include "posre_Protein_chain_{ID}.itp"',
                                    '#endif'])   
                mol_include.extend([f'Protein_chain_{ID} 1'])     

        elif isinstance(value, protein.Monomer):
            chainID = getattr(complex, key).chains
            for ID in chainID:
                itp_include.extend(['#include "protein.itp"',
                                    '#ifdef POSRES',
                                    '#include "posre.itp"',
                                    '#endif',
                                    '#ifdef DISRE',             # DISRE (BW) only applicable for monomers
                                    '#include "disre.itp"',
                                    '#endif'])
                mol_include.extend([f'Protein_chain_{ID} 1'])  

    for key, value in vars(complex).items():
        # define cofactor restraints
        if isinstance(value, protein.Ligand):
            itp_include.extend([f'#include "{key}.itp"',
                                '#ifdef POSRES',
                                f'#include "posre_{key}.itp"',
                                '#endif'])
            mol_include.extend([f'{getattr(complex, key).ID} 1'])  

        # TODO: posre ITP doesn't link to correct molecules for Ions and CrystalWaters
        if  isinstance(value, protein.Ions):
            itp_include.extend(['; posre_{key}.itp currently cannot be included"'])
            #itp_include.extend(['; Include Position restraint file',
            #                    f'#include "posre_{key}.itp"'])  
            mol_include.extend([f'{key} {getattr(complex, key)._n_ions}'])

        if isinstance(value, protein.CrystalWaters) :
            itp_include.extend(['; posre_{key}.itp currently cannot be included"'])
            #itp_include.extend(['; Include Position restraint file',
            #                    f'#include "posre_{key}.itp"'])  
            mol_include.extend([f'{key} {getattr(complex, key)._n_waters}'])

    tgt.write(t.substitute(working_dir = working_dir,
                           itp_includes = "\n".join(itp_include),
                           mol_includes = "\n".join(mol_include)))

    tgt.close()

    return True

def tar_out(src_dir = [], tgt = []):
    """
    Tar everything in a src_dir to the tar_file
    """
    import tarfile

    t_f = tarfile.open(tgt, mode="w:gz")
    base_dir = os.getcwd()
    os.chdir(src_dir) # To avoid the include of all parent dirs
    for to_tar in os.listdir(os.path.join(base_dir, src_dir)):
        t_f.add(to_tar)
    t_f.close()
    os.chdir(base_dir)

########################################################################

# AMBER - CHARMM

########################################################################

def get_mapping(csv_path):
    amber_to_charmm = {}
    charmm_to_amber = {}
    with open(csv_path, 'r') as f:
        import csv
        reader = csv.reader(f)
        for row in reader:
            if len(row) > 3 and row[0] == "POPC":
                charmm_str = row[2].strip()
                amber_str = row[3].strip()
                c_parts = charmm_str.split()
                a_parts = amber_str.split()
                if len(c_parts) == 2 and len(a_parts) == 2:
                    c_atom, c_res = c_parts[0], c_parts[1]
                    a_atom, a_res = a_parts[0], a_parts[1]
                    amber_to_charmm[(a_res, a_atom)] = (c_res, c_atom)
                    charmm_to_amber[(c_res, c_atom)] = (a_res, a_atom)
    return amber_to_charmm, charmm_to_amber

def amber2charmm(src="", tgt="", **kwargs):
    """
    Translates Amber Lipid21 residues (PA, PC, OL) to POPC using charmmlipid2amber.csv.
    """
    import os
    csv_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "charmmlipid2amber.csv")
    amber_to_charmm, _ = get_mapping(csv_path)
    
    is_gro = src.endswith(".gro")
    with open(src, "r") as f:
        lines = f.readlines()
        
    out_lines = []
    
    for i, line in enumerate(lines):
        if is_gro:
            if i < 2 or i == len(lines)-1:
                out_lines.append(line)
                continue
            resname = line[5:10].strip()
            atomname = line[10:15].strip()
            if resname in ["PA", "PC", "OL"]:
                key = (resname, atomname)
                if key in amber_to_charmm:
                    c_res, c_atom = amber_to_charmm[key]
                    new_resname = c_res.rjust(5)
                    new_atomname = c_atom.rjust(5)
                    # Renumber residues uniquely in GRO format (columns 1-5, 1-indexed based on 134 atoms per lipid)
                    residue_number = ((i - 2) // 134) + 1
                    res_num_str = str(residue_number).rjust(5)
                    out_lines.append(res_num_str + new_resname + new_atomname + line[15:])
                else:
                    out_lines.append(line)
            else:
                out_lines.append(line)
        else:
            if not line.startswith("ATOM"):
                out_lines.append(line)
                continue
            resname = line[17:21].strip()
            atomname = line[12:16].strip()
            if resname in ["PA", "PC", "OL"]:
                key = (resname, atomname)
                if key in amber_to_charmm:
                    c_res, c_atom = amber_to_charmm[key]
                    new_resname = c_res.ljust(4)
                    new_atomname = c_atom.center(4)
                    out_lines.append(line[:12] + new_atomname + " " + new_resname + line[21:])
                else:
                    out_lines.append(line)
            else:
                out_lines.append(line)
                
    with open(tgt if tgt else src, "w") as f:
        f.writelines(out_lines)
    return True

def charmm2amber(src="", tgt="", **kwargs):
    """
    Reads coordinate file. Converts POPC back to PA, PC, OL using mapping.
    """
    import os
    csv_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "charmmlipid2amber.csv")
    _, charmm_to_amber = get_mapping(csv_path)
        
    is_gro = src.endswith(".gro")
    with open(src, "r") as f:
        lines = f.readlines()
        
    out_lines = []
    
    for i, line in enumerate(lines):
        if is_gro:
            if i < 2 or i == len(lines)-1:
                out_lines.append(line)
                continue
            resname = line[5:10].strip()
            atomname = line[10:15].strip()
            if resname in ["POPC", "POP"]:
                key = ("POPC", atomname)
                if key in charmm_to_amber:
                    a_res, a_atom = charmm_to_amber[key]
                    new_resname = a_res.rjust(5)
                    new_atomname = a_atom.rjust(5)
                    out_lines.append(line[:5] + new_resname + new_atomname + line[15:])
                else:
                    out_lines.append(line)
            else:
                out_lines.append(line)
        else:
            if not line.startswith("ATOM"):
                out_lines.append(line)
                continue
            resname = line[17:21].strip()
            atomname = line[12:16].strip()
            if resname in ["POPC", "POP"]:
                key = ("POPC", atomname)
                if key in charmm_to_amber:
                    a_res, a_atom = charmm_to_amber[key]
                    new_resname = a_res.ljust(4)
                    new_atomname = a_atom.center(4)
                    out_lines.append(line[:12] + new_atomname + " " + new_resname + line[21:])
                else:
                    out_lines.append(line)
            else:
                out_lines.append(line)
                
    with open(tgt if tgt else src, "w") as f:
        f.writelines(out_lines)
    return True
