import functools
import logging
import os
import shutil
import subprocess

import broker
import groerrors
import protein
import recipes
import settings
import utils

logger = logging.getLogger('pymemdyn.gromacs')


class Gromacs(object):
    def __init__(self, *args, **kwargs):
        self.logger = logging.getLogger('pymemdyn.gromacs.Gromacs')
        self.broker = kwargs.get("broker") or broker.Printing()
        self._membrane_complex = None
        self.wrapper = Wrapper()
        logging.basicConfig(filename='GROMACS.log',
                            format='%(asctime)s %(message)s',
                            datefmt='%m/%d/%Y %I:%M:%S',
                            level=logging.DEBUG)
        if "membrane_complex" in kwargs.keys():
            self.set_membrane_complex(kwargs["membrane_complex"])
            self.tpr = \
                self.membrane_complex.complex.proteins.pdb.replace(".pdb", 
                                                                   ".tpr")
        

    def set_membrane_complex(self, value):
        """
        set_membrane_complex: Sets the membrane object
        """
        self._membrane_complex = value

    def get_membrane_complex(self):
        return self._membrane_complex
    membrane_complex = property(get_membrane_complex, set_membrane_complex)

    def count_lipids(self, **kwargs):
        """
        count_lipids: Counts the lipids in source and writes a target with N4 tags
        """
        src = open(kwargs["src"], "r")
        tgt = open(kwargs["tgt"], "w")

        half = self.membrane_complex.complex.prot_z / 2
        self.membrane_complex.membrane.lipids_up = 0
        self.membrane_complex.membrane.lipids_down = 0
        self.membrane_complex.membrane.n_wats = 0

        if hasattr(self.membrane_complex.complex, "waters"):
            # Careful, some waters belong to crystal, not solvent
            self.membrane_complex.membrane.n_wats -= \
            self.membrane_complex.complex.waters.number

        for line in src:
            if len(line.split()) > 2:
                # This N4 and POP are the lipid markers
                if line.split()[2] == "N4" and line.split()[3] == "POP":
                    tgt.write(line)
                    if (float(line.split()[7]) >= half):
                        self.membrane_complex.membrane.lipids_up += 1
                    elif (float(line.split()[7]) < half):
                        self.membrane_complex.membrane.lipids_down += 1
                elif line.split()[2] == "OW":  # Water marker
                    self.membrane_complex.membrane.n_wats += 1

        src.close()
        tgt.close()

        return True

    def get_charge(self, **kwargs):
        """
        get_charge: Gets the total charge of a system using gromacs grompp command
        """
        out, err = self.wrapper.run_command({"gromacs": "grompp",
                                             "options": kwargs})

        # Now we are looking for this line:
        # System has non-zero total charge: 6.000002e+00
        charge = 0
        for line in err.decode("utf-8").split("\n"):
            if "total charge" in line:
                # In gromacs the charge is not displayed in scientific notation.
                # so this will result in giving a charge of 5, for a charge of 5.99999
                print(line)
                charge = abs(int(round(float(line.split()[-1]))))
                break

        self.membrane_complex.complex.negative_charge = 0
        self.membrane_complex.complex.positive_charge = 0
        if charge > 0:
            self.membrane_complex.complex.positive_charge = charge
            self.membrane_complex.complex.negative_charge = 0
        elif charge < 0:
            self.membrane_complex.complex.negative_charge = charge
            self.membrane_complex.complex.positive_charge = 0
        else:
            self.membrane_complex.complex.negative_charge = 0
            self.membrane_complex.complex.positive_charge = 0

        return True

    def get_ndx_groups(self, **kwargs):
        """
        get_ndx_groups: Run make_ndx and set the total number of groups found
        """
        out, err = self.wrapper.run_command({"gromacs": "make_ndx",
                                             "options": kwargs,
                                             "input": "q\n"})

        for line in out.decode("utf-8").split("\n"):
            if "atoms" in line:
                self.n_groups = int(line.split()[0])

        return True

    def get_ndx_sol(self, **kwargs):
        """
        get_ndx_sol: Run make_ndx and set the last number id for SOL found
        """
        out, err = self.wrapper.run_command({"gromacs": "make_ndx",
                                             "options": kwargs,
                                             "input": "q\n"})

        for line in out.decode("utf-8").split("\n"):
            # print(line)
            logger.debug(line)
            if "SOL" in line:
                self.n_sol = int(line.split()[0])

        return True

    def make_ndx(self, **kwargs):
        """
        make_ndx: Wraps the make_ndx command tweaking the input to reflect the
        characteristics of the complex
        """
        if not (self.get_ndx_groups(**kwargs)): return False
        n_group = self.n_groups

        if not (self.get_ndx_sol(**kwargs)): return False
        n_sol = self.n_sol

        # Create the solution with no crystal water, crossing fingers.
        input = "r SOL \n"
        input += "name {0} SOL\n".format(n_sol)
        input += "del {0}\n".format(n_sol)

        # Create the "wation" group (always present)
        n_group += 1
        input += " r SOL | r HOH | r CHL | r Cl* | r SOD | r Na* \n"
        input += "name {0} wation\n".format(n_group)

        # Create the "protlig" group
        n_group += 1
        input += " \"Protein\" | r LIG | r ALO \n"
        input += "name {0} protlig\n".format(n_group)

        # Create the "membr" group
        n_group += 1
        input += " r POP* | r CHO | r LIP\n"
        input += "name {0} membr\n".format(n_group)

        # This makes a separate group for each chain (if more than one)
        if type(self.membrane_complex.complex.proteins) == protein.Oligomer:
            for chain in self.membrane_complex.complex.proteins.chains:
                n_group += 1
                input += "a {0}-{1}\n".format(
                    self.membrane_complex.complex.proteins.points[chain][0],
                    # start point of chain X
                    self.membrane_complex.complex.proteins.points[chain][1])  
                    # end point of chain X
                input += "name {0} Protein_chain_{1}\n".format(n_group, chain)

        if hasattr(self.membrane_complex.membrane, "ions"):
            # TODO This makes the group ions 
            # n_group += 1
            # input += "1 || r LIG\nname {0} protlig\n".format(n_group)
            pass
        input += "q\n"

        # Now the wrapper itself
        out, err = self.wrapper.run_command({"gromacs": "make_ndx",
                                             "options": kwargs,
                                             "input": input})

        logging.debug("make_ndx command")
        logging.debug(err.decode().strip('\n'))
        logging.debug(out.decode().strip('\n'))

        self.logger.debug("make_ndx command")
        # self.logger.debug(err.decode().strip('\n'))
        # self.logger.debug(out.decode().strip('\n'))

        return True

    def make_topol_lipids(self, **kwargs):
        """
        make_topol_lipids: Add lipid positions to topol.top
        """
        topol = open("topol.top", "a")
        topol.write("; Number of POPC molecules with higher z-coord value:\n")
        topol.write("POPC " + str(self.membrane_complex.membrane.lipids_up))
        topol.write("\n; Number of POPC molecules with lower z-coord value:\n")
        topol.write("POPC " + str(self.membrane_complex.membrane.lipids_down))
        topol.write("\n; Total number of SOL molecules:\n")
        
        if hasattr(self.membrane_complex.complex, "HOH"):
            topol.write("SOL " + str(self.membrane_complex.membrane.n_wats - self.membrane_complex.complex.HOH._n_waters) + "\n")
        else:
            topol.write("SOL " + str(self.membrane_complex.membrane.n_wats) + "\n")

        topol.close()

        return True

    def manual_log(self, command, output):
        """
        manual_log: Redirect the output to file in command["options"]["log"]
        Some commands can't be logged via flag, so one has to catch and
        redirect stdout and stderr
        """
        log = open(command["options"]["log"], "w")
        log.writelines(str(output.decode('utf8').strip('\n')))
        # print(str(output.decode('utf8').strip('\n')))
        log.close()

        return True

    def relax(self, **kwargs):
        """
        relax: Relax a protein
        """
        if not os.path.isdir(kwargs["tgt_dir"]): os.makedirs(kwargs["tgt_dir"])
        posres = kwargs.get("posres", [])


        if type(self.membrane_complex.complex.proteins) == protein.Monomer:
            posres.append("posre.itp")
        elif type(self.membrane_complex.complex.proteins) == protein.Oligomer:
            for chain in self.membrane_complex.complex.proteins.chains:
                posres.extend(["posre_Protein_chain_"+chain+".itp"])

        # for var, value in vars(kwargs["membrane_complex"]).items():            
        #     if isinstance(value, protein.Ligand) or isinstance(value, protein.CrystalWaters) or isinstance(value, protein.Ions):
        #         posres.append(f"posre_{var}.itp")

        # if hasattr(self.membrane_complex.complex, "waters") and \
        #         self.membrane_complex.complex.waters:
        #     posres.append("posre_hoh.itp")
        # if hasattr(self.membrane_complex.complex, "ions") and \
        #         self.membrane_complex.complex.ions:
        #     posres.append("posre_ion.itp")

        for posre in posres:
            new_posre = open(os.path.join(kwargs["tgt_dir"], posre), "w")

            for line in open(os.path.join(kwargs["src_dir"], posre), "r"):
                if line.split()[-3:] == ["1000", "1000", "1000"]:
                    new_posre.write(" ".join(line.split()[:2] + \
                                             [str(kwargs["const"])] * 3))
                    new_posre.write("\n")
                else:
                    new_posre.write(line)
            new_posre.close()

        new_mdp = open(os.path.join(kwargs["tgt_dir"], f"eq{str(kwargs['const'])}.mdp"), "w")
        src_mdp = open(os.path.join(kwargs["src_dir"], kwargs["mdp"]), "r")

        for line in src_mdp:
            if (line.startswith("gen_vel")):
                new_mdp.write("gen_vel             = no\n")
            else:
                new_mdp.write(line)
        src_mdp.close()
        new_mdp.close()

        # utils.make_topol(target_dir=kwargs["tgt_dir"],
        #                  working_dir=os.getcwd(),
        #                  complex=self.membrane_complex.complex)

    def restrain_ca(self, **kwargs):
        """
        restrain_ca: Restrain only protein CA by index
        """
        # Get C-alpha atoms from index
        index_file = kwargs.get("index")
        index = open(index_file, "r")

        ca_ids = set()
        section = False
        for line in index:
            if line.startswith("[ C-alpha ]"):
                section = True
            elif line.startswith("[") and section:
                break
            elif section:
                ids = line.split()
                for id in ids:
                    ca_ids.add(int(id.strip()))
        
        self.logger.debug(f'Extracted C-alpha IDs: {sorted(ca_ids)}')

        index.close()

        # keep only CA atom restraints
        src_files = kwargs["src_files"]
        
        self.logger.debug(f'scr_files: {src_files}')
        for src in src_files:
            if os.path.exists(src):
                self.logger.debug(f'processing CA-restraints: {src}')
                with open(src, 'r') as file:
                    lines = file.readlines()

                ca_restraints = ['[ position_restraints ]\n']
                for line in lines:
                    if line:
                        parts = line.split()
                        if parts:
                            if parts[0].isdigit(): 
                                final_nr = int(parts[0])
                                if int(parts[0]) in ca_ids:
                                    ca_restraints.append(line)
            
                self.logger.debug(f'ca_restr.: {ca_restraints}')

                ca_corr = {x - final_nr for x in ca_ids}
                ca_ids = {x for x in ca_corr if x > 0}

                self.logger.debug(f'correction factor: {final_nr}')
                self.logger.debug(f'corrected CA-IDs: {sorted(ca_ids)}')
               
                with open(src, 'w') as file:
                    file.writelines(ca_restraints)
            else:
                self.logger.debug(f'{src} does not exist')

        return True

    def run_recipe(self, debugFast=False):
        """
        run_recipe: Run recipe for the complex
        """
        if not hasattr(self, "recipe"):
            self.select_recipe(debugFast=debugFast)

        self.repo_dir = self.wrapper.repo_dir

        for n, command_name in enumerate(self.recipe.steps):
            command = self.recipe.recipe[command_name]

            if command_name in self.recipe.breaks.keys():
                command["options"] = self.set_options(command["options"],
                                     self.recipe.breaks[command_name])

            # NOW RUN IT !
            self.broker.dispatch("{0} Step ({1}/{2}): {3}.".format(
                self.recipe.__class__.__name__,
                n + 1, len(self.recipe.steps),
                command_name))
            self.logger.info("{0} Step ({1}/{2}): {3}.".format(
                self.recipe.__class__.__name__,
                n + 1, len(self.recipe.steps),
                command_name))
            if ("gromacs") in command:
                # Either run a Gromacs pure command...
                if hasattr(self, "queue"): command["queue"] = self.queue
                out, err = self.wrapper.run_command(command)
                logging.debug(" ".join(self.wrapper.generate_command(command)))
                logging.debug(err.decode().strip('\n'))
                logging.debug(out.decode().strip('\n'))
                self.logger.debug(" ".join(self.wrapper.generate_command(command)))
                # self.logger.debug(err.decode().strip('\n'))
                self.logger.debug(out.decode().strip('\n'))
                ## Next line tests Gromacs output checking for known errors
                groerrors.GromacsMessages(gro_err=err,
                                          command=command["gromacs"])

                # Some commands are unable to log via flag, so we catch and
                # redirect stdout and stderr
                if command["gromacs"] in ["energy"]:
                    self.manual_log(command, out)
                    
            else:
                # ...or run a local function
                logging.debug(command)
                self.logger.debug(command)
                try:
                    f = getattr(self, command["command"])
                except AttributeError:
                    # Fallback to the utils module
                    f = getattr(utils, command["command"])

                if ("options") in command:
                    f(**command["options"])
                else:
                    f()
                logging.debug("FUNCTION: " + str(f.__doc__).strip())
                self.logger.debug("FUNCTION: " + str(f.__doc__).strip())

        return True

    def select_recipe(self, stage="", debugFast=False, full_relax=True):
        """
        select_recipe: Select the appropriate recipe for the complex
        """
        recipe = ""
        stage = stage or "Init"

        if self.membrane_complex:
            if not any(isinstance(var, protein.Ligand) for var in vars(self.membrane_complex.complex).values()):
                recipe += "Basic"
            else:
                recipe += "Ligand"

        recipe += stage  # This kwarg carries the proper recipe:
        # Init, Minimization, Equilibration...

        if hasattr(recipes, recipe):
            self.recipe = getattr(recipes, recipe)(debugFast=debugFast, membrane_complex=self.membrane_complex.complex, full_relax=full_relax)
        elif hasattr(recipes, "Basic" + stage):
            # Fall back to Basic recipe if no specific where found
            self.recipe = getattr(recipes, "Basic" + stage)(debugFast=debugFast, membrane_complex=self.membrane_complex.complex, full_relax=full_relax)

        return True

    def set_box_sizes(self):
        """
        set_box_sizes: Set length values for different boxes
        """
        self.membrane_complex.complex.set_nanom()
        self.membrane_complex.trans_box_size = \
            [str(self.membrane_complex.complex.gmx_prot_xy),
             str(self.membrane_complex.complex.gmx_prot_xy),
             str(self.membrane_complex.membrane.gmx_bilayer_z)]
        self.membrane_complex.bilayer_box_size = \
            [str(self.membrane_complex.membrane.gmx_bilayer_x),
             str(self.membrane_complex.membrane.gmx_bilayer_y),
             str(self.membrane_complex.membrane.gmx_bilayer_z)]
        self.membrane_complex.embeded_box_size = \
            [str(self.membrane_complex.membrane.gmx_bilayer_x),
             str(self.membrane_complex.membrane.gmx_bilayer_y),
             str(self.membrane_complex.complex.gmx_emb_z)]
        self.membrane_complex.protein_box_size = \
            [str(self.membrane_complex.complex.gmx_prot_xy),
             str(self.membrane_complex.complex.gmx_prot_xy),
             str(self.membrane_complex.complex.gmx_prot_z)]

        return True

    def set_chains(self, **kwargs):
        """
        set_chains: Set the REAL points of a dimer after protonation
        """
        src = kwargs.get("src")

        if type(self.membrane_complex.complex.proteins) == protein.Oligomer:
            points = self.membrane_complex.complex.proteins.points
            with open(src, "r") as pdb_fp:
                for line in pdb_fp:
                    if len(line) > 21 and line.startswith(("ATOM", "HETATM")):
                        atom_serial = int(line[6:11])
                        chain_id = line[21]
                        if chain_id != ' ':
                            if points[chain_id]:
                                points[chain_id] = [
                                    min(atom_serial, points[chain_id][0]),
                                    max(atom_serial, points[chain_id][1])]
                            else:
                                points[chain_id] = [atom_serial, atom_serial]
                        

            self.membrane_complex.complex.proteins.points = points

        return True

    def set_grompp(self, **kwargs):
        """
        set_grompp: Copy template files to working dir
        """
        for repo_src in kwargs.keys():
            shutil.copy(os.path.join(
                self.repo_dir, kwargs[repo_src]),
                repo_src)

        return True

    def clean_itp(self, **kwargs):
        """
        clean_itp: Cut a itp file to be usable later as with restraints
        """
        src_files = kwargs["src_files"]
        
        for src in src_files:
            if os.path.exists(src):
                with open(src, 'r') as file:
                    lines = file.readlines()

                tgt_lines = []
                toggle_keep = True
                for line in lines:
                    if line.startswith("#ifdef POSRES"):
                        toggle_keep = False
                    
                    if toggle_keep == True:
                        tgt_lines.append(line)
                    
                    if line.startswith("#endif"):
                        toggle_keep = True

                with open(src, 'w') as file:
                    file.writelines(tgt_lines)

        return True

    def set_itp(self, **kwargs):
        """
        set_itp: Cut a top file to be usable later as itp
        """
        src = open(kwargs["src"], "r")
        tgt = open(kwargs["tgt"], "w")

        get_name = False

        for line in src:
            if line.startswith("#include"):
                pass
            elif line.startswith("; Include Position restraint file"):
                break
            else:
                tgt.write(line)

            if get_name and not line.startswith(";"):
                self.membrane_complex.complex.proteins.name = line.split()[0]
                get_name = False

            if line.startswith("[ moleculetype ]"):
                get_name = True

        tgt.close()
        src.close()

        return True

    def set_options(self, options, breaks):
        """
        set_options: Set break options from recipe
        """
        for option, value in breaks.items():
            # This is a hack to get the attribute recursively,
            # feeding getattr with dot-splitted string thanks to reduce
            # Here we charge some commands with options calculated
            new_option = functools.reduce(getattr,
                                value.split("."),
                                self)
            options[option] = new_option

        return options

    def set_popc(self, tgt=""):
        """
        set_popc: Create a pdb file only with the lipid bilayer (POP), no waters.
        Set some measures on the fly (height of the bilayer)
        """
        tgt = open(tgt, "w")
        pops = []

        for line in open(self.membrane_complex.membrane.pdb, "r"):
            if len(line.split()) > 7:
                if line.split()[3] == "POP":
                    tgt.write(line)
                    pops.append(float(line.split()[7]))

        tgt.close()
        self.membrane_complex.membrane.bilayer_z = max(pops) - min(pops)

        return True

    def set_protein_height(self, **kwargs):
        """
        set_protein_height: Get the z-axis center from a pdb file for membrane or
        solvent alignment
        """
        z_tot = []
        z_mem = []
        
        for line in open(kwargs["src"], "r"):
            if line[:4] == "ATOM":
                z_tot.append(float(line[46:54]))
                
                if line[17:20] == "POP" and line[13:15] == "N4":
                    z_mem.append(float(line[46:54]))
        
        if "solvate2" in kwargs.keys():
            correction = (sum(z_tot)/len(z_tot)-sum(z_mem)/len(z_mem))/10
            if correction <= 0:
                self.membrane_complex.complex.center_z = correction
            else:
                self.membrane_complex.complex.center_z = 0
            # Correction allows protpopc to be aligned correctly with water.gro.
            # Editconf translates the system, but it seems to only misalign 
            # water.gro and protpopc.pdb when the protein center is intracellular.
            
            # For any misalignment issues, look here
                
        else:
            self.membrane_complex.complex.center_z = (sum(z_tot)/len(z_tot))/10
        
        return True

    def set_protein_size(self, **kwargs):
        """
        set_protein_size: Get the protein maximum base width from a pdb file
        """
        for line in open(kwargs["src"], "r"):
            if line.startswith("CRYST1"):
                if kwargs["dir"] == "xy":
                    self.membrane_complex.complex.prot_xy = \
                        max(float(line.split()[1]),  # xyprotein
                            float(line.split()[2]))
                elif kwargs["dir"] == "z":
                    self.membrane_complex.complex.prot_z = \
                        float(line.split()[3])  # hprotein
                    self.set_box_sizes()

                break

        return True

    def set_stage_init(self, **kwargs):
        """
        set_stage_init: Copy a set of files from source to target dir
        """
        if not os.path.isdir(kwargs["tgt_dir"]): os.mkdir(kwargs["tgt_dir"])

        if "src_files" in kwargs.keys():
            for src_file in kwargs["src_files"]:
                if (os.path.isfile(os.path.join(kwargs["src_dir"], src_file))):
                    shutil.copy(os.path.join(kwargs["src_dir"], src_file),
                                os.path.join(kwargs["tgt_dir"],
                                            os.path.split(src_file)[1]))

        if "repo_files" in kwargs.keys():
            for repo_file in kwargs["repo_files"]:
                shutil.copy(os.path.join(self.repo_dir, repo_file),
                            os.path.join(kwargs["tgt_dir"], repo_file))

        return True

    def set_steep(self, **kwargs):
        """
        set_steep: Copy the template steep.mdp to target dir
        """
        shutil.copy(os.path.join(self.repo_dir, "steep.mdp"),
                    "steep.mdp")
        
        return True

    def set_water(self, **kwargs):
        """
        set_water: Create a water layer for a box
        """
        start = (self.membrane_complex.membrane.bilayer_zw - \
                 self.membrane_complex.complex.prot_z) / 2
        end = start + self.membrane_complex.complex.prot_z

        src = open(self.membrane_complex.membrane.pdb, "r")
        tgt = open(kwargs["tgt"], "w")

        res = "NULL"
        for line in src:
            if len(line.split()) > 7:
                if ((line.split()[2] == "OW") and
                        ((float(line.split()[7]) > end) or
                             (float(line.split()[7]) < start))):
                    res = line.split()[4]
                if ((line.split()[4] != res) and
                        (line.split()[3] == "SOL")):
                    tgt.write(line)

        tgt.close()
        src.close()

        return True


class Wrapper(object):
    def __init__(self, *args, **kwargs):
        self.work_dir = os.getcwd()
        # The gromacs to be used
        self.gromacs_dir = settings.GROMACS_PATH
        # The directory where all the files live
        self.repo_dir = settings.TEMPLATES_DIR

    def _common_io(self, src, tgt):
        """
        _common_io: Autoexpand many Gromacs commands that use -f for input
        and -o for the output file
        """
        return ["-f", src, "-o", tgt]

    def generate_command(self, kwargs):
        """
        generate_command: Receive some variables in kwargs, generate
        the appropriate command to be run. Return a set in the form of
        a string "command -with flags"
        """
        try:
            mode = kwargs["gromacs"]
        except KeyError:
            raise

        if "src" in kwargs["options"].keys():
            src = self._setDir(kwargs["options"]["src"])
        if "tgt" in kwargs["options"].keys():
            tgt = self._setDir(kwargs["options"]["tgt"])
        options = kwargs["options"]

        command = [os.path.join(self.gromacs_dir, "gmx")]
        command.extend([mode])
        if "queue" in kwargs.keys():
            if hasattr(kwargs["queue"], mode):
                # If we got a queue enabled for this command, use it
                command = list(kwargs["queue"].command)  # Already a list
                kwargs["queue"].make_script(
                    workdir=kwargs["options"]["dir"],
                    options=self._mode_mdrun(options))

        # Standard -f input -o output
        if mode in ["pdb2gmx", "editconf", "grompp", "trjconv",
                    "make_ndx", "genrestr", "energy"]:
            command.extend(self._common_io(src, tgt))

            if (mode == "pdb2gmx"):  # PDB2GMX
                command.extend(self._mode_pdb2gmx(options))
            if (mode == "editconf"):  # EDITCONF
                command.extend(self._mode_editconf(options))
            if (mode == "grompp"):  # GROMPP
                command.extend(self._mode_grompp(options))
            if (mode == "trjconv"):  # TRJCONV
                command.extend(self._mode_trjconv(options))
            if (mode == "make_ndx"):  # MAKE_NDX
                pass
            if (mode == "genrestr"):  # GENRSTR
                command.extend(self._mode_genrest(options))
            if (mode == "energy"):  # ENERGY
                pass

        else:
            if (mode == "eneconv"):  # ENECONV
                command.extend(self._mode_eneconv(options))
            if (mode == "solvate"):  # SOLVATE
                command.extend(self._mode_solvate(options))
            if (mode == "genion"):  # GENION
                command.extend(self._mode_genion(options))
            if (mode == "rms"):  # RMS
                command.extend(self._mode_rms(options))
            if (mode == "rmsf"):  # RMSF
                command.extend(self._mode_rmsf(options))
            if (mode == "trjcat"):  # TRJCAT
                command.extend(self._mode_trjcat(options))
            if (mode == "mdrun"):  # MDRUN_SLURM
                pass

        return command

    def _mode_editconf(self, kwargs):
        """
        _mode_editconf: Wrap the editconf command options
        """
        command = []

        if "dist" in kwargs.keys():
            command.extend(["-d", str(kwargs["dist"])])
        if "box" in kwargs.keys():
            command.extend(["-box"])
            command.extend(kwargs["box"])
        if "angles" in kwargs.keys():
            command.extend(["-angles"])
            command.extend(kwargs["angles"])
            command.extend(["-bt", kwargs["bt"]])
        if "center" in kwargs.keys():
            command.extend(["-center"])
            command.extend(kwargs["center"])
        if "translate_z" in kwargs.keys():
            command.extend(["-translate", "0", "0"])
            command.extend([str(kwargs["translate_z"])])
        

        return command

    def _mode_eneconv(self, kwargs):
        """
        mode_eneconv: Wrap the eneconv command options
        """
        src_files = utils.make_cat(kwargs["dir1"],
                                   kwargs["dir2"],
                                   kwargs["name"])
        command = ["-f"]
        command.extend(src_files)
        command.extend(["-o", self._setDir(kwargs["tgt"]),
                        "-settime"])

        return command

    def _mode_rms(self, kwargs):
        """
        _mode_rms: Wrap the rms command options
        """
        return ["-s", self._setDir(kwargs["src"]),
                "-f", self._setDir(kwargs["src2"]),
                "-o", self._setDir(kwargs["tgt"])]

    def _mode_rmsf(self, kwargs):
        """
        _mode_rmsf: Wrap the rmsf command options and report per
        residue RMSF
        """
        return ["-s", self._setDir(kwargs["src"]),
                "-f", self._setDir(kwargs["src2"]),
                "-o", self._setDir(kwargs["tgt"]),
                "-res"]

    def _mode_genion(self, kwargs):
        """
        _mode_genion: Wrap the genion command options
        """
        if (kwargs["np"]) == 0 and (kwargs["nn"] == 0):
            # Genion refuses to run at all if no ion are provided, so provide
            # one of each
            kwargs["np"] = 1
            kwargs["nn"] = 1
        command = ["-s", kwargs["src"],
                   "-o", kwargs["tgt"],
                   "-p", self._setDir(kwargs["src2"]),
                   "-n", kwargs["index"],
                   "-np", str(kwargs["np"]),
                   "-nn", str(kwargs["nn"]),
                   "-pname", "NA+",
                   "-nname", "CL-"]

        return command

    def _mode_genrest(self, kwargs):
        """
        _mode_genrest: Wrap the genrest command options
        """
        command = ["-fc"] + kwargs["forces"]

        if "index" in kwargs.keys():
            command.extend(["-n", kwargs["index"]])

        return command

    def _mode_grompp(self, kwargs):
        """
        _mode_grompp: Wrap the grompp command options
        """
        command = ["-maxwarn", " 3",
                   "-c", self._setDir(kwargs["src2"]),
                   "-r", self._setDir(kwargs["src2"]),
                   "-p", self._setDir(kwargs["top"]),
                   "-po", self._setDir("mdout.mdp")]
        if "index" in kwargs.keys():
            command.extend(["-n", self._setDir(kwargs["index"])])

        if "tgt_top" in kwargs.keys():
            command.extend(["-pp", self._setDir(kwargs["tgt_top"])]) # DEBUGGING

        return command

    def _mode_mdrun(self, kwargs):
        """
        _mode_mdrun: Wrap the mdrun command options
        """
        command = ["-s", kwargs["src"],
                   "-o", kwargs["tgt"],
                   "-e", kwargs["energy"],
                   "-c", kwargs["conf"],
                   "-g", kwargs["log"]]

        if ("traj" in kwargs.keys()):
            command.extend(["-x", kwargs["traj"]])

        if ("cpi" in kwargs.keys()):
            command.extend(["-cpi", kwargs["cpi"]])

        return command

    def _mode_pdb2gmx(self, kwargs):
        """
        _mode_pdb2gmx: Wrap the pdb2gmx command options
        """
        return ["-p", self._setDir(kwargs["top"]),
                "-i", self._setDir("posre.itp"),
                "-ignh", "-ff", "oplsaa", "-water", "spc"]

    def _mode_solvate(self, kwargs):
        '''
        _mode_solvate: Wrap the solvate command options
        '''
        return ["-cp", self._setDir(kwargs["cp"]),
                "-cs", self._setDir(kwargs["cs"]),
                "-p", self._setDir(kwargs["top"]),
                "-o", self._setDir(kwargs["tgt"])]

    def _mode_trjcat(self, kwargs):
        """
        _mode_trjcat: Wrap the trjcat command options
        """
        src_files = utils.make_cat(kwargs["dir1"],
                                   kwargs["dir2"],
                                   kwargs["name"])
        command = ["-f"]
        command.extend(src_files)
        command.extend(["-o", self._setDir(kwargs["tgt"]),
                        "-settime"])

        return command

    def _mode_trjconv(self, kwargs):
        """
        _mode_trjconv: Wrap the trjconv command options
        """
        command = ["-s", self._setDir(kwargs["src2"]),
                   "-pbc", kwargs["pbc"]]
        if "ur" in kwargs.keys():
            command.extend(["-ur", kwargs["ur"]])
        if "skip" in kwargs.keys():
            command.extend(["-skip", kwargs["skip"]])
        if "trans" in kwargs.keys():
            command.extend(["-trans"])
            command.extend([str(x) for x in kwargs["trans"]])
        else:
            command.extend(["-center"])

        return command

    def run_command(self, kwargs):
        """
        run_command: Run a command that comes in kwargs in a subprocess, and
        return the output as (output, errors)
        """
        command = self.generate_command(kwargs)

        if ("input" in kwargs.keys()):
            p = subprocess.Popen(command,
                                 stdin=subprocess.PIPE,
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE)
            gro_out, gro_errs = p.communicate(kwargs["input"].encode())
        else:
            p = subprocess.Popen(command,
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE)

            gro_out, gro_errs = p.communicate()

        return gro_out, gro_errs

    def _setDir(self, filename):
        """
        _setDir: Expand a filename with the work dir to save code space
        """
        return os.path.join(self.work_dir, filename)
