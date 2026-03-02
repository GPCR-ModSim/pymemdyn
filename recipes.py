"""This module describes the commandline or python commands for all the 
phases of pymemdyn. It consists of: 

    - Init
    - Minimization
    - Equilibration
    - Relaxation
    - Collecting results
    
"""

import os

import protein



class BasicInit(object):
    def __init__(self, **kwargs):
        # First we make a list of ordered steps
        self.steps = ["pdb2gmx", "set_itp", "clean_itp", "concat", "set_protein_height", "editconf",
                      "set_protein_size", "editconf2", "set_protein_size2",
                      "set_popc", "editconf3", "editconf4", "make_topol",
                      "editconf5", "solvate", "set_protein_height2", "set_water", 
                      "editconf6", "editconf7", "editconf8", 
                      "solvate2", "count_lipids", "make_topol2",
                      "make_topol_lipids", "make_ffoplsaanb", "set_grompp",
                      "set_chains", "make_ndx", "grompp", "trjconv",
                      "get_charge", "genion", "grompp2", "trjconv2",
                      "grompp3", "trjconv3", "clean_pdb"]

        # And then we define each step
        self.recipe = \
            {"pdb2gmx": {"gromacs": "pdb2gmx",  # 1
                         "options": {"src": "",
                                     "tgt": "proteinopls.pdb",
                                     "top": "protein.top"}},             
             
             "set_itp": {"command": "set_itp",  # 2
                         "options": {"src": "protein.top",
                                     "tgt": "protein.itp"}},

             "clean_itp": {"command": "clean_itp",
                           "options": {"src_files": []}},

             "concat": {"command": "concat",  # 3
                        "options": {"src": "proteinopls.pdb",
                                    "tgt": ""}},

             "set_protein_height": {"command": "set_protein_height",  # 4
                                    "options": {"src": "proteinopls.pdb"}},

             "editconf": {"gromacs": "editconf",  # 5
                          "options": {"src": "proteinopls.pdb",
                                      "tgt": "proteinopls.pdb",
                                      "dist": ""}},

             "set_protein_size": {"command": "set_protein_size",  # 6
                                  "options": {"src": "proteinopls.pdb",
                                              "dir": "xy"}},

             "editconf2": {"gromacs": "editconf",  # 7
                           "options": {"src": "proteinopls.pdb",
                                       "tgt": "proteinopls.pdb",
                                       "dist": ""}},

             "set_protein_size2": {"command": "set_protein_size",  # 8
                                   "options": {"src": "proteinopls.pdb",
                                               "dir": "z"}},

             "set_popc": {"command": "set_popc",  # 9
                          "options": {"tgt": "popc.pdb"}},

             "editconf3": {"gromacs": "editconf",  # 10
                           "options": {"src": "proteinopls.pdb",
                                       "tgt": "proteinopls.pdb",
                                       "box": "",
                                       "angles": ["90", "90", "120"],
                                       "bt": "tric"}},

             "editconf4": {"gromacs": "editconf",  # 11
                           "options": {"src": "popc.pdb",
                                       "tgt": "popc.gro",
                                       "box": ""}},

             "make_topol": {"command": "make_topol",  # 12
                            "options": {}},

             "editconf5": {"gromacs": "editconf",  # 13
                           "options": {"src": "proteinopls.pdb",
                                       "tgt": "proteinopls.pdb",
                                       "translate_z": ""}},

             "solvate": {"gromacs": "solvate",  # 14
                        "options": {"cp": "proteinopls.pdb",
                                    "cs": "popc.gro",
                                    "tgt": "protpopc.pdb",
                                    "top": "topol.top"}},

             "set_protein_height2": {"command": "set_protein_height",  # 15
                                    "options": {"src": "protpopc.pdb",
                                                "solvate2": "True"}}, 
             
             "set_water": {"command": "set_water",  # 16
                           "options": {"tgt": "water.pdb"}},

             "editconf6": {"gromacs": "editconf",  # 17
                           "options": {"src": "water.pdb",
                                       "tgt": "water.gro",
                                       "box": ""}},

             "editconf7": {"gromacs": "editconf",  # 18
                           "options": {"src": "protpopc.pdb",
                                       "tgt": "protpopc.pdb",
                                       "box": "",
                                       "angles": ["90", "90", "120"],
                                       "bt": "tric"}},
             
             "editconf8": {"gromacs": "editconf",  # 19
                           "options": {"src": "protpopc.pdb",
                                       "tgt": "protpopc.pdb",
                                       "translate_z": ""}},

             "solvate2": {"gromacs": "solvate",  # 20
                         "options": {"cp": "protpopc.pdb",
                                     "cs": "water.gro",
                                     "tgt": "tmp.pdb",
                                     "top": "topol.top"}},

             "count_lipids": {"command": "count_lipids",  # 21
                              "options": {"src": "tmp.pdb",
                                          "tgt": "popc.pdb"}},

             "make_topol2": {"command": "make_topol",  # 22
                             "options": {}},

             "make_topol_lipids": {"command": "make_topol_lipids"},  # 23

             "make_ffoplsaanb": {"command": "make_ffoplsaanb",  # 24
                                 "options": {}},

             "set_grompp": {"command": "set_grompp",  # 25
                            "options": {"steep.mdp": "steep.mdp",
                                        "popc.itp": "popc.itp",
                                        "spc.itp": "spc.itp",
                                        "ions.itp": "ions.itp",
                                        "ffoplsaabon_mod.itp": "ffoplsaabon_mod.itp",
                                        "ffoplsaa_mod.itp": "ffoplsaa_mod.itp"}},

             "set_chains": {"command": "set_chains",  # 26
                            "options": {"src": "proteinopls.pdb"}},

             "make_ndx": {"command": "make_ndx",  # 27
                          "options": {"src": "tmp.pdb",
                                      "tgt": "index.ndx"}},

             "grompp": {"gromacs": "grompp",  # 28
                        "options": {"src": "steep.mdp",
                                    "src2": "tmp.pdb",
                                    "tgt": "topol.tpr",
                                    "top": "topol.top",
                                    "index": "index.ndx"}},

             "trjconv": {"gromacs": "trjconv",  # 29
                         "options": {"src": "tmp.pdb",
                                     "src2": "topol.tpr",
                                     "tgt": "tmp.pdb",
                                     "pbc": "mol",
                                     "index": "index.ndx"},
                         "input": "1\n0\n"},

             "get_charge": {"command": "get_charge",  # 30
                            "options": {"src": "steep.mdp",
                                        "src2": "tmp.pdb",
                                        "tgt": "topol.tpr",
                                        "top": "topol.top",
                                        "index": "index.ndx"}},

             "genion": {"gromacs": "genion",  # 31
                        "options": {"src": "topol.tpr",
                                    "tgt": "output.pdb",
                                    "src2": "topol.top",
                                    "index": "index.ndx",
                                    "np": "",
                                    "nn": ""},
                        "input": " SOL \n"},

             "grompp2": {"gromacs": "grompp",  # 32
                         "options": {"src": "steep.mdp",
                                     "src2": "output.pdb",
                                     "tgt": "topol.tpr",
                                     "top": "topol.top"}},

             "trjconv2": {"gromacs": "trjconv",  # 33
                          "options": {"src": "output.pdb",
                                      "src2": "topol.tpr",
                                      "tgt": "output.pdb",
                                      "trans": [],
                                      "pbc": "mol"},
                          "input": "0\n"},

             "grompp3": {"gromacs": "grompp",  # 34
                         "options": {"src": "steep.mdp",
                                     "src2": "output.pdb",
                                     "tgt": "topol.tpr",
                                     "top": "topol.top"}},

             "trjconv3": {"gromacs": "trjconv",  # 35
                          "options": {"src": "output.pdb",
                                      "src2": "topol.tpr",
                                      "tgt": "hexagon.pdb",
                                      "ur": "compact",
                                      "pbc": "mol"},
                          "input": "1\n0\n"},
             
             "clean_pdb": {"command": "clean_pdb", # 36
                           "options": {"src": "hexagon.pdb",
                                       "tgt": "hexagon.pdb"}}
             }

        self.breaks = \
            {"pdb2gmx": {"src": "membrane_complex.complex.proteins.pdb_his"},
             "concat": {"tgt": "membrane_complex.complex"},
             "editconf": {"dist": "membrane_complex.box_height"},
             "editconf2": {"dist": "membrane_complex.box_width"},
             "editconf3": {"box": "membrane_complex.trans_box_size"},
             "editconf4": {"box": "membrane_complex.bilayer_box_size"},
             "editconf5": {"translate_z": "membrane_complex.complex.center_z"},
             "editconf6": {"box": "membrane_complex.embeded_box_size"},
             "editconf7": {"box": "membrane_complex.protein_box_size"},
             "editconf8": {"translate_z": "membrane_complex.complex.center_z"},
             "genion": {"nn": "membrane_complex.complex.positive_charge",
                        "np": "membrane_complex.complex.negative_charge"},
             "make_topol": {"complex": "membrane_complex.complex"},
             "make_topol2": {"complex": "membrane_complex.complex"},
             "make_ffoplsaanb": {"complex": "membrane_complex.complex"},
             "trjconv2": {"trans": "membrane_complex.complex.trans"}
             }

        if kwargs["membrane_complex"].proteins.chains:
            for chain in kwargs["membrane_complex"].proteins.chains:
                self.recipe["clean_itp"]["options"]["src_files"].append(
                                                    f"protein_Protein_chain_{chain}.itp")

        if kwargs["debugFast"] or False:
            self.recipe["set_grompp"]["options"]["steep.mdp"] = "steepDEBUG.mdp"


# This recipe modifies the previous one taking ligands into account
class LigandInit(BasicInit):
    def __init__(self, **kwargs):
        super(LigandInit, self).__init__(**kwargs)

        if kwargs["membrane_complex"]:
            for var, value in vars(kwargs["membrane_complex"]).items():
                if isinstance(value, protein.Ligand):

                    self.steps.insert(9, f"genrestr_{var}")
                    self.recipe[f"genrestr_{var}"] = \
                        {"gromacs": "genrestr", 
                        "options": {"src": "",
                                    "tgt": f"posre_{var}.itp",
                                    "index": f"ligand_{var}_ha.ndx",
                                    "forces": ["1000", "1000", "1000"]},
                        "input": "3\n"}

                    self.steps.insert(9, f"make_ndx_{var}")
                    self.recipe[f"make_ndx_{var}"] = \
                        {"gromacs": "make_ndx",
                        "options": {"src": "",
                                    "tgt": f"ligand_{var}_ha.ndx",
                                    var: True},
                        "input": "! a H*\nq\n"}

                    self.breaks[f"make_ndx_{var}"] = \
                        {"src": f"membrane_complex.complex.{var}.pdb"}
                    self.breaks[f"genrestr_{var}"] = \
                        {"src": f"membrane_complex.complex.{var}.pdb"}


##########################################################################
#                 Minimization                                           #
##########################################################################

class BasicMinimization(object):
    def __init__(self, **kwargs):
        self.steps = ["set_stage_init", "mdrun"]
        self.recipe = {
            "set_stage_init": {"command": "set_stage_init",  # 1
                               "options": {"src_dir": "",
                                           "src_files": ["topol.tpr"],
                                           "tgt_dir": "Rmin",
                                           "repo_files": ["steep.mdp"]}},

            "mdrun": {"gromacs": "mdrun",  # 2
                      "options": {"dir": "Rmin",
                                  "src": "topol.tpr",
                                  "tgt": "traj.trj",
                                  "energy": "ener.edr",
                                  "conf": "confout.gro",
                                  "log": "md.log"}},
        }
        self.breaks = {}

        if kwargs["debugFast"] or False:
            self.recipe["set_stage_init"]["options"]["repo_files"] = \
                ["eqDEBUG.mdp"]


##########################################################################
#                 Equilibration                                          #
##########################################################################

class BasicEquilibration(object):
    def __init__(self, **kwargs):
        self.steps = ["clean_itp", "editconf", "make_ndx", "set_grompp", "set_stage_init", "grompp", 
                      "set_stage_init2", "mdrun"]
        self.recipe = {
            "clean_itp": {"command": "clean_itp",
                           "options": {"src_files": ["protein.itp"]}},
            
            "editconf": {"gromacs": "editconf",  # 1
                         "options": {"src": "Rmin/confout.gro",
                                     "tgt": "min.pdb"}},

            "make_ndx": {"command": "make_ndx",  # 2
                         "options": {"src": "min.pdb",
                                     "tgt": "index.ndx"}},

            "set_grompp": {"command": "set_grompp",  # 3
                           "options": {"eq.mdp": "eq.mdp"}},

            "set_stage_init": {"command": "set_stage_init",  # 4
                               "options": {"src_dir": "",
                                           "src_files": ["eq.mdp"],
                                           "tgt_dir": "eq"}},            
                                           
            "grompp": {"gromacs": "grompp",  # 5
                       "options": {"src": "eq.mdp",
                                   "src2": "min.pdb",
                                   "top": "topol.top",
                                   "tgt": "topol.tpr",
                                   "tgt_top": "eq/processed.top", # DEBUGGING
                                   "index": "index.ndx"}},

            "set_stage_init2": {"command": "set_stage_init",  # 6
                                "options": {"src_dir": "",
                                            "src_files": ["topol.tpr",
                                                          "posre.itp"],
                                            "tgt_dir": "eq"}},

            "mdrun": {"gromacs": "mdrun",  # 7
                      "options": {"dir": "eq",
                                  "src": "topol.tpr",
                                  "tgt": "traj.trr",
                                  "energy": "ener.edr",
                                  "conf": "confout1000.gro",
                                  "traj": "traj.xtc",
                                  "log": "md_eq1000.log"}},
        }
        self.breaks = {}

        if kwargs["membrane_complex"].proteins.chains:
            for chain in kwargs["membrane_complex"].proteins.chains:
                self.recipe["set_stage_init2"]["options"]["src_files"].append(
                                                    f"posre_Protein_chain_{chain}.itp")
        if kwargs["membrane_complex"]:    
            for var, value in vars(kwargs["membrane_complex"]).items():
                if isinstance(value, protein.Ligand) or isinstance(value, protein.CrystalWaters) or isinstance(value, protein.Ions):
                    self.recipe["set_stage_init2"]["options"]["src_files"].append(
                                                        f"posre_{var}.itp")

        if kwargs["debugFast"] or False:
            self.recipe["grompp"]["options"]["src"] = "Rmin/eqDEBUG.mdp"
            self.recipe["set_stage_init"]["options"]["src_files"] = \
                ["eqDEBUG.mdp"]
            
                
#class LigandEquilibration(BasicEquilibration):
#    def __init__(self, **kwargs):
#        super(LigandEquilibration, self).__init__(**kwargs)

        # if kwargs["membrane_complex"]:
        #     for var, value in vars(kwargs["membrane_complex"]).items():
        #         if isinstance(value, protein.Ligand):
        #             self.steps.insert(2, "genrestr")
        #             self.recipe["genrestr"] = \
        #                 {"gromacs": "genrestr",
        #                 "options": {"src": "Rmin/topol.tpr",
        #                             "tgt": "protein_ca200.itp",
        #                             "index": "index.ndx",
        #                             "forces": ["200", "200", "200"]},
        #                 "input": "3\n"}
                    

##########################################################################
#                 Relaxation                                             #
##########################################################################

class BasicRelax(object):
    def __init__(self, **kwargs):
        self.steps = []
        self.recipe = {}
        for const in range(800, 0, -200):
            self.steps.extend(["relax{0}".format(const),
                               "set_stage_init{0}".format(const),
                               "grompp{0}".format(const),
                               "mdrun{0}".format(const)])
            tgt_dir = "eq/{0}".format(const)
            src_dir = "eq"
            self.recipe["relax%d" % const] = \
                {"command": "relax",  # 1, 5, 9, 13
                 "options": {"const": const,
                             "src_dir": src_dir,
                             "tgt_dir": tgt_dir,
                             "posres": [],
                             "mdp": "eq.mdp"}}
            
            self.recipe["set_stage_init%d" % const] = \
                {"command": "set_stage_init",  # 2, 6, 10, 14
                 "options": {"src_dir": "",
                             "src_files": ["topol.top",
                                           "ffoplsaa_mod.itp",
                                           "ffoplsaanb_mod.itp",
                                           "ffoplsaabon_mod.itp",
                                           "protein.itp",
                                           "popc.itp",
                                           "ions.itp",
                                           "spc.itp"],
                             "tgt_dir": "eq/{0}".format(const)}}   

            self.recipe["grompp%d" % const] = \
                {"gromacs": "grompp",  # 3, 7, 11, 15
                 "options": {"src": os.path.join(tgt_dir, "eq{0}.mdp".format(const)),
                             "src2": os.path.join(src_dir, "confout{0}.gro".format(const+200)),
                             "top": os.path.join(tgt_dir, "topol.top"),
                             "tgt": os.path.join(tgt_dir, "topol.tpr"),
                             "tgt_top": os.path.join(tgt_dir, "processed.top"), # DEBUGGING
                             "index": "index.ndx"}}

            self.recipe["mdrun%d" % const] = \
                {"gromacs": "mdrun",  # 4, 8, 12, 16
                 "options": {"dir": tgt_dir,
                             "src": "topol.tpr",
                             "tgt": "traj.trr",
                             "energy": "ener.edr",
                             "conf": "../confout{0}.gro".format(const),
                             "traj": "traj.xtc",
                             "log": "md_eq{0}.log".format(const)}}        

            for chain in kwargs["membrane_complex"].proteins.chains:
                self.recipe[f"set_stage_init{const}"]["options"]["src_files"].append(f"protein_Protein_chain_{chain}.itp")

            for var, value in vars(kwargs["membrane_complex"]).items():            
                if isinstance(value, protein.Ligand) or isinstance(value, protein.CrystalWaters) or isinstance(value, protein.Ions):
                    self.recipe[f"relax{const}"]["options"]["posres"].append(f"posre_{var}.itp")
                    self.recipe[f"set_stage_init{const}"]["options"]["src_files"].append(f"{var}.itp")

        # TODO: add section for copying .itp for oligomers (posre is handled in relax)
       
        self.steps.append("trjconv_to_pdb_halfway") # 17

        self.recipe["trjconv_to_pdb_halfway"] = {
            "gromacs": "trjconv",
            "options": {
                "src": "eq/confout200.gro",
                "src2": "eq/200/topol.tpr",
                "tgt": "eq/confout200.pdb",
                "pbc": "mol",
                "ur": "compact"
            },
            "input": "1\n0\n"
        }

        self.breaks = {}

        if kwargs["debugFast"] or False:
            for i in [x for x in self.recipe.keys() if x.startswith("relax")]:
                self.recipe[i]["options"]["mdp"] = "eqDEBUG.mdp"

# class LigandRelax(BasicRelax):
#     def __init__(self, **kwargs):
#         super(LigandRelax, self).__init__(**kwargs)
#         if kwargs["membrane_complex"]:
#             for var, value in vars(kwargs["membrane_complex"]).items():
#                 if isinstance(value, protein.Ligand):
#                     self.recipe["relax800"]["options"]["posres"].append(f"posre_{var}.itp")


##########################################################################
#                 C-Alpha Atoms Relaxation                               #
##########################################################################

class BasicCARelax(object):
    def __init__(self, **kwargs):
        self.steps = ["set_stage_init", "set_stage_init2", "restrain_ca", "grompp", "mdrun"]
        self.recipe = {
            "set_stage_init": {"command": "set_stage_init",  # 1
                               "options": {"src_dir": "eq",
                                           "tgt_dir": "eqProd",
                                           "src_files": ["confout200.gro"],
                                           "repo_files": ["eqCA.mdp"]}},

            "set_stage_init2": {"command": "set_stage_init",  # 1
                                "options": {"src_dir": "eq/200",
                                            "tgt_dir": "eqProd",
                                            "src_files": ["topol.top",
                                                          "ffoplsaa_mod.itp",
                                                          "ffoplsaanb_mod.itp",
                                                          "ffoplsaabon_mod.itp",
                                                          "protein.itp",
                                                          "posre.itp",
                                                          "popc.itp",
                                                          "ions.itp",
                                                          "spc.itp"]}},

            "restrain_ca": {"command": "restrain_ca", 
                            "options": {"src_files": ["eqProd/posre.itp"],
                                        "index": "index.ndx"}},
            
            # "genrestr": {"gromacs": "genrestr",  # 2
            #              "options": {"src": "Rmin/topol.tpr",
            #                          "tgt": "posre.itp",
            #                          "index": "index.ndx",
            #                          "forces": ["200"] * 3},
            #              "input": "3\n"},

            "grompp": {"gromacs": "grompp",  # 3
                       "options": {"src": "eqProd/eqCA.mdp",
                                   "src2": "eqProd/confout200.gro",
                                   "top": "eqProd/topol.top",
                                   "tgt": "eqProd/topol.tpr",
                                   "tgt_top": "eqProd/processed.top", # DEBUGGING
                                   "index": "index.ndx"}},

            "mdrun": {"gromacs": "mdrun",  # 4
                      "options": {"dir": "eqProd",
                                  "src": "topol.tpr",
                                  "tgt": "traj.trr",
                                  "energy": "ener.edr",
                                  "conf": "confout.gro",
                                  "traj": "traj.xtc",
                                  "log": "md_eqProd.log"}},
        }

        self.breaks = {}

        for chain in kwargs["membrane_complex"].proteins.chains:
            self.recipe[f"set_stage_init2"]["options"]["src_files"].extend([f"protein_Protein_chain_{chain}.itp", f"posre_Protein_chain_{chain}.itp"])
            self.recipe[f"restrain_ca"]["options"]["src_files"].extend([f"eqProd/posre_Protein_chain_{chain}.itp"])

        for var, value in vars(kwargs["membrane_complex"]).items():            
            if isinstance(value, protein.Ligand) or isinstance(value, protein.CrystalWaters) or isinstance(value, protein.Ions):
                self.recipe[f"set_stage_init2"]["options"]["src_files"].extend([f"{var}.itp", f"posre_{var}.itp"])

        if kwargs["debugFast"] or False:
            self.recipe["set_stage_init"]["options"]["src_files"] = \
                ["confout.gro", "eqDEBUG.mdp"]
            self.recipe["grompp"]["options"]["src"] = "eqProd/eqDEBUG.mdp"

        if kwargs["full_relax"] != True:
                self.steps.remove("mdrun")


##########################################################################
#                 Ballesteros-Weinstein Restrained Relaxation            #
##########################################################################

class BasicBWRelax(object):
    def __init__(self, **kwargs):
        self.steps = ["getbw","set_stage_init", "grompp", "mdrun"]
        self.recipe = {
            "getbw": {"command": "getbw",  # 1
                      "options": {"src": "proteinopls.pdb"}},

            "set_stage_init": {"command": "set_stage_init",  # 2
                               "options": {"src_dir": "eq",
                                           "tgt_dir": "eqProd",
                                           "src_files": ["confout200.gro",
                                                         "../disre.itp"],
                                           "repo_files": ["dres.mdp"]}},

            "grompp": {"gromacs": "grompp",  # 4
                       "options": {"src": "eqProd/dres.mdp",
                                   "src2": "eqProd/confout200.gro",
                                   "top": "topol.top",
                                   "tgt": "eqProd/topol.tpr",
                                   "tgt_top": "eqProd/processed.top", # DEBUGGING
                                   "index": "index.ndx"}},

            "mdrun": {"gromacs": "mdrun",  # 4
                      "options": {"dir": "eqProd",
                                  "src": "topol.tpr",
                                  "tgt": "traj.trr",
                                  "energy": "ener.edr",
                                  "conf": "confout.gro",
                                  "traj": "traj.xtc",
                                  "log": "md_eqBW.log"}},
        }

        self.breaks = {}

        if kwargs["debugFast"] or False:
            self.recipe["set_stage_init"]["options"]["src_files"] = \
                ["confout.gro", "eqDEBUG.mdp"]
            self.recipe["grompp"]["options"]["src"] = "eqProd/eqDEBUG.mdp"

        if kwargs["full_relax"] != True:
                self.steps.remove("mdrun")

##########################################################################
#                 Collect All Results & Output                           #
##########################################################################
class BasicCollectResults(object):
    def __init__(self, **kwargs):
        """
        This recipe navigates through the output of *GROMACS*, generating and
        collecting various files and logs, and finally putting it all together
        in a tar file.

        Provides a list of *steps* as the order to execute the recipe.
        dict *recipe* as the variables to pass to functions.
        dict *breaks* as points where object calling can put their vars.
        """
        self.breaks = {}
        self.steps = ["trjcat", "trjconv", "trjconv_to_pdb", "rms1", "rms2",
                      "rms3", "rmsf", "tot_ener", "temp", "pressure",
                      "volume", "set_stage_init", "grompp",
                      "set_end", "set_end_2", "set_end_3", 
                      "set_end_4", "set_end_5", "set_end_6",
                      "set_end_7", "tar_it"]

        self.recipe = {"trjcat":
                           {"gromacs": "trjcat",  # 1
                            "options": {"dir1": "eq",
                                        "dir2": "eqProd",
                                        "name": "traj.xtc",
                                        "tgt": "traj_EQ.xtc"},
                            "input": "c\n" * 6},

                       "trjconv":
                           {"gromacs": "trjconv",  # 2
                                   "options": {"src": "traj_EQ.xtc",
                                               "src2": "topol.tpr",
                                               "tgt": "traj_pymol.xtc",
                                               "ur": "compact",
                                               "skip": "2",
                                               "pbc": "mol"},
                                   "input": "1\n0\n"},

                       "trjconv_to_pdb":
                           {"gromacs": "trjconv", 
                                   "options": {"src": "eqProd/confout.gro",
                                               "src2": "eqProd/topol.tpr",
                                               "tgt": "eqProd/confout.pdb",
                                               "ur": "compact",
                                               "pbc": "mol"},
                                   "input": "1\n0\n"},

                       "rms1":
                           {"gromacs": "rms",  # 4
                                 "options": {"src": "eq/topol.tpr",
                                             "src2": "traj_EQ.xtc",
                                             "tgt": "rmsd-all-atom-vs-start.xvg"},
                                 "input": "1\n1\n"},

                       "rms2":
                           {"gromacs": "rms",  # 5
                                 "options": {"src": "eq/topol.tpr",
                                             "src2": "traj_EQ.xtc",
                                             "tgt": "rmsd-backbone-vs-start.xvg"},
                                 "input": "4\n4\n"},

                       "rms3":
                           {"gromacs": "rms",  # 6
                                 "options": {"src": "eq/topol.tpr",
                                             "src2": "traj_EQ.xtc",
                                             "tgt": "rmsd-calpha-vs-start.xvg"},
                                 "input": "3\n3\n"},

                       "rmsf":
                           {"gromacs": "rmsf",  # 7
                                 "options": {"src": "eq/topol.tpr",
                                             "src2": "traj_EQ.xtc",
                                             "tgt": "rmsf-per-residue.xvg"},
                                 "input": "1\n"},

                        "set_stage_init": {"command": "set_stage_init",  # 2
                                           "options": {"tgt_dir": ".",
                                                       "repo_files": ["prod.mdp"]}},

                        "grompp": {"gromacs": "grompp",  # 3
                                "options": {"src": "prod.mdp",
                                            "src2": "eq/confout200.gro",
                                            "top": "topol.top",
                                            "tgt": "topol.tpr",
                                            "tgt_top": "prod.top", # DEBUGGING
                                            "index": "index.ndx"}},

                       "set_end":
                           {"command": "set_stage_init",  # 12
                                   "options": {"src_dir": "eqProd",
                                               "src_files": ["confout.gro",
                                                             "processed.top",
                                                             "dres.mdp",
                                                             "eqCA.mdp"],
                                               "repo_files": ["prod.mdp",
                                                              "README.md",
                                                              "load_gpcr.pml"],
                                               "tgt_dir": "finalOutput"}},

                       "set_end_2":
                           {"command": "set_stage_init",  # 14
                                     "options": {"src_dir": "",
                                                 "src_files": [
                                                     "log.log",
                                                     "hexagon.pdb",
                                                     "index.ndx", 
                                                     "traj_EQ.xtc",
                                                     "ener_EQ.edr",
                                                     "traj_pymol.xtc",
                                                     "prod.top"],
                                                 "tgt_dir": "finalOutput"}},

                       "set_end_3":
                           {"command": "set_stage_init",  # 15
                                     "options": {"src_dir": "",
                                                 "src_files": ["tot_ener.xvg",
                                                               "tot_ener.log",
                                                               "temp.xvg",
                                                               "temp.log",
                                                               "pressure.xvg",
                                                               "pressure.log",
                                                               "volume.xvg",
                                                               "volume.log",                                                     
                                                               "rmsd-all-atom-vs-start.xvg",
                                                               "rmsd-calpha-vs-start.xvg",
                                                               "rmsd-backbone-vs-start.xvg",
                                                               "rmsf-per-residue.xvg"],
                                                 "tgt_dir": "finalOutput/reports"}},

                       "set_end_4":
                           {"command": "set_stage_init",  # 16
                                     "options": {"src_dir": "eq",
                                                 "src_files": ["md_eq1000.log",
                                                               "confout1000.gro",
                                                               "confout800.gro",
                                                               "confout600.gro",
                                                               "confout400.gro",
                                                               "confout200.gro",],
                                                 "tgt_dir": "finalOutput/logs"}},

                       "set_end_5":
                           {"command": "set_stage_init",  # 17
                                     "options": {"src_dir": "eq",
                                                 "src_files": [
                                                     "{0}/md_eq{0}.log".format(
                                                         const)
                                                     for const in
                                                     range(800, 0, -200)],
                                                 "tgt_dir": "finalOutput/logs"}},

                       "set_end_6":
                           {"command": "set_stage_init",  # 18
                                     "options": {"src_dir": "eqProd",
                                                 "src_files": ["md_eqBW.log",
                                                               "md_eqCA.log"],
                                                 "tgt_dir": "finalOutput/logs"}},
                       "set_end_7":
                           {"command": "set_stage_init",
                                     "options": {"src_dir": "",
                                                 "src_files": ["eqProd/confout.pdb",
                                                               "eq/confout200.pdb"],
                                                 "tgt_dir": "finalOutput"}},

                       "tar_it":
                           {"command": "tar_out", # 19
                                  "options": {"src_dir": "finalOutput",
                                              "tgt": "MD_output.tgz"}}, 
                       }

        options = {"tot_ener": "Total-Energy\n",
                   "temp": "Temperature\n",
                   "pressure": "Pressure\n",
                   "volume": "Volume\n"} # 8 to 11

        for option, gro_key in options.items():
            self.recipe[option] = \
                {"gromacs": "energy",
                 "options": {"src": "ener_EQ.edr",
                             "tgt": "{0}.xvg".format(option),
                             "log": "{0}.log".format(option)},
                 "input": gro_key}

        if kwargs["full_relax"] != True:
            self.recipe["trjcat"]["options"]["dir2"] = ""
            self.steps.remove("set_stage_init")
            self.steps.remove("grompp")

class BasicCACollectResults(BasicCollectResults):
    def __init__(self, **kwargs):
        super(BasicCACollectResults, self).__init__(**kwargs)
        self.steps.insert(2, "eneconv")
        self.recipe["eneconv"] = \
                {"gromacs": "eneconv", 
                 "options": {"dir1": "eq",
                         "dir2": "eqProd",
                         "name": "ener.edr",
                         "tgt": "ener_EQ.edr"},
                 "input": "y\nc\nc\nc\nc\nc\nc\n"}


class BasicBWCollectResults(BasicCollectResults):
    def __init__(self, **kwargs):
        super(BasicBWCollectResults, self).__init__(**kwargs)
        self.steps.insert(14, "set_end_BW")
        self.recipe["set_end_BW"] = \
                {"command": "set_stage_init",
                          "options": {"src_dir": "",
                                      "src_files": ["tot_ener2.xvg",
                                                    "tot_ener2.log",
                                                    "temp2.xvg",
                                                    "temp2.log",
                                                    "pressure2.xvg",
                                                    "pressure2.log",
                                                    "volume2.xvg",
                                                    "volume2.log"],                                                     
                                      "tgt_dir": "finalOutput/reports"}}

        self.steps.insert(10, "volume2")
        self.steps.insert(10, "pressure2")
        self.steps.insert(10, "temp2")
        self.steps.insert(10, "tot_ener2")
            
        options2 = {"tot_ener2": "Total-Energy\n",
                        "temp2": "Temperature\n",
                        "pressure2": "Pressure\n",
                        "volume2": "Volume\n"}

        for option, gro_key in options2.items():
            self.recipe[option] = \
                {"gromacs": "energy",
                 "options": {"src": "eqProd/ener.edr",
                             "tgt": "{0}.xvg".format(option),
                             "log": "{0}.log".format(option)},
                 "input": gro_key}

        self.steps.insert(2, "eneconv")
        self.recipe["eneconv"] = \
                {"gromacs": "eneconv", 
                 "options": {"dir1": "eq",
                             "dir2": "",
                             "name": "ener.edr",
                             "tgt": "ener_EQ.edr"},
                 "input": "y\nc\nc\nc\nc\nc\nc\n"}
        
