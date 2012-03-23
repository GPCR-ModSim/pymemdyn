import os

class MonomerRecipe(object):
    def __init__(self, *args):
        self.recipe = \
        [{"gromacs": "pdb2gmx", #0
          "options": {"src": "",
                      "tgt": "proteinopls.pdb",
                      "top": "protein.top"}},
         {"command": "set_itp", #1
          "options": {"src": "protein.top",
                      "tgt": "protein.itp"}},
         {"gromacs": "editconf", #2
          "options": {"src": "proteinopls.pdb",
                      "tgt": "proteinopls.pdb",
                      "dist": ""}},
         {"command": "set_protein_size", #3
          "options": {"src": "proteinopls.pdb",
                      "dir": "xy"}},
         {"gromacs": "editconf", #4
          "options": {"src": "proteinopls.pdb",
                      "tgt": "proteinopls.pdb",
                      "dist": ""}},
         {"command": "set_protein_size", #5
          "options": {"src": "proteinopls.pdb",
                      "dir": "z"}},
         {"command": "set_popc", #6
          "options": {"tgt": "popc.pdb"}},
         {"gromacs": "editconf", #7
          "options": {"src": "proteinopls.pdb",
                      "tgt": "proteinopls.pdb",
                      "box": "",
                      "angles": ["90", "90", "120"],
                      "bt": "tric"}},
         {"gromacs": "editconf", #8
          "options": {"src": "popc.pdb",
                      "tgt": "popc.pdb",
                      "box": ""}},
         {"command": "make_topol",
          "options": {"protein" :1}}, #9
         {"gromacs": "editconf", #10
          "options": {"src": "proteinopls.pdb",
                      "tgt": "proteinopls.pdb",
                      "translate": ["0", "0", "0"]}},
         {"gromacs": "genbox", #11
          "options": {"cp": "proteinopls.pdb",
                      "cs": "popc.pdb",
                      "tgt": "protpopc.pdb",
                      "top": "topol.top"}},
          {"command": "set_water", #12
           "options": {"tgt": "water.pdb"}},
          {"gromacs": "editconf", #13
           "options": {"src": "water.pdb",
                       "tgt": "water.pdb",
                       "box": ""}},
          {"gromacs": "editconf", #14
           "options": {"src": "protpopc.pdb",
                       "tgt": "protpopc.pdb",
                       "box": "",
                       "angles": ["90", "90", "120"],
                       "bt": "tric"}},
          {"gromacs": "genbox", #15
           "options": {"cp": "protpopc.pdb",
                       "cs": "water.pdb",
                       "tgt": "tmp.pdb",
                       "top": "topol.top"}},
          {"command": "count_lipids", #16
           "options": {"src": "tmp.pdb",
                       "tgt": "popc.pdb"}},
          {"command": "make_topol",#17
           "options": {"protein": 1}},
          {"command": "make_topol_lipids"}, #18
          {"command": "set_grompp", #19
           "options": {"steep.mdp": "steep.mdp",
                       "popc.itp": "popc.itp",
                       "ffoplsaanb_mod.itp": "ffoplsaanb_mod.itp",
                       "ffoplsaabon_mod.itp": "ffoplsaabon_mod.itp",
                       "ffoplsaa_mod.itp": "ffoplsaa_mod.itp"}},
          {"gromacs": "grompp", #20
           "options": {"src": "steep.mdp",
                       "src2": "tmp.pdb",
                       "tgt": "topol.tpr",
                       "top": "topol.top"}},
          {"gromacs": "trjconv", #21
           "options": {"src": "tmp.pdb",
                       "src2": "topol.tpr",
                       "tgt": "tmp.pdb",
                       "pbc": "mol"},
           "input": "1\n0\n"},
          {"command": "get_charge",
           "options": {"src": "steep.mdp",
                       "src2": "tmp.pdb",
                       "tgt": "topol.tpr",
                       "top": "topol.top"}}, #22
          {"gromacs": "genion", #23
           "options": {"src": "topol.tpr",
                       "src2": "topol.top",
                       "tgt": "output.pdb",
                       "np": "",
                       "nn": ""},
           "input": "SOL\n"},
          {"gromacs": "grompp", #24
           "options": {"src": "steep.mdp",
                       "src2": "output.pdb",
                       "tgt": "topol.tpr",
                       "top": "topol.top"}},
          {"gromacs": "trjconv", #25
           "options": {"src": "output.pdb",
                       "src2": "topol.tpr",
                       "tgt": "output.pdb",
                       "trans": [],
                       "pbc": "mol"},
           "input": "0\n"},
          {"gromacs": "grompp", #26
           "options": {"src": "steep.mdp",
                       "src2": "output.pdb",
                       "tgt": "topol.tpr",
                       "top": "topol.top"}},
          {"gromacs": "trjconv", #27
           "options": {"src": "output.pdb",
                       "src2": "topol.tpr",
                       "tgt": "hexagon.pdb",
                       "ur": "compact",
                       "pbc": "mol"},
           "input": "1\n0\n"}
        ]

        self.breaks = {0: {"src": "membrane_complex.complex.monomer.pdb"},
                       2: {"dist": "membrane_complex.box_height"},
                       4: {"dist": "membrane_complex.box_width"},
                       7: {"box": "membrane_complex.trans_box_size"},
                       8: {"box": "membrane_complex.bilayer_box_size"},
                       13: {"box": "membrane_complex.embeded_box_size"},
                       14: {"box": "membrane_complex.protein_box_size"},
                       23: {"np": "membrane_complex.complex.positive_charge",
                            "nn": "membrane_complex.complex.negative_charge"},
                       25: {"trans": "membrane_complex.complex.trans"}
                      }

        if "debug" in args:
            self.recipe[19]["options"]["steep.mdp"] = "steepDEBUG.mdp"

# This recipe modifies the previous one taking a ligand into account
class MonomerLigandRecipe(MonomerRecipe):
    def __init__(self):
        super(MonomerLigandRecipe, self).__init__()
        self.recipe[19]["options"]["ffoplsaanb_mod.itp"] =\
            "ffoplsaanb_mod_lig.itp"
        self.recipe[19]["options"]["ffoplsaabon_mod.itp"] =\
            "ffoplsaabon_mod_lig.itp"
        self.recipe[19]["options"]["ffoplsaa_mod.itp"] = "ffoplsaa_mod_lig.itp"
        self.recipe[17] =\
            {"command": "make_topol",
             "options": {"ligand": ""}}
        self.recipe[9] =\
            {"command": "make_topol",
             "options": {"ligand": ""}}

        self.recipe.insert(9,
            {"gromacs": "genrestr",
             "options": {"src": "",
                         "tgt": "posre_lig.itp",
                         "index": "ligand_ha.ndx",
                         "forces": ["1000", "1000", "1000"]},
             "input": "2\n"})

        self.recipe.insert(9,
            {"gromacs": "make_ndx",
             "options": {"src": "",
                         "tgt": "ligand_ha.ndx",
                         "ligand": True},
             "input": "! a H*\nq\n"})

        self.recipe.insert(2,
            {"command": "concat",
             "options": {"src": "proteinopls.pdb",
                         "tgt": ""}})

        self.breaks = {0: {"src": "membrane_complex.complex.monomer.pdb"},
                       2: {"tgt": "membrane_complex.complex.ligand.pdb"},
                       3: {"dist": "membrane_complex.box_height"},
                       5: {"dist": "membrane_complex.box_width"},
                       8: {"box": "membrane_complex.trans_box_size"},
                       9: {"box": "membrane_complex.bilayer_box_size"},
                       10: {"src": "membrane_complex.complex.ligand.pdb"},
                       11: {"src": "membrane_complex.complex.ligand.pdb"},
                       12: {"ligand": "membrane_complex.complex.ligand.pdb",
                            "protein": "membrane_complex.complex.monomer.pdb"},
                       #13: {"box": "membrane_complex.embeded_box_size"},
                       14: {"box": "membrane_complex.protein_box_size"},
                       16: {"box": "membrane_complex.embeded_box_size"},
                       17: {"box": "membrane_complex.protein_box_size"},
                       20: {"ligand": "membrane_complex.complex.ligand.pdb",
                            "protein": "membrane_complex.complex.monomer.pdb"},
                       26: {"np": "membrane_complex.complex.positive_charge",
                            "nn": "membrane_complex.complex.negative_charge"},
                       28: {"trans": "membrane_complex.complex.trans"},
                      }

class BasicMinimization(object):
    def __init__(self, *args):
        self.recipe = \
        [{"command": "set_stage_init", #0
          "options": {"src_dir": "",
                      "src_files": ["topol.tpr"],
                      "tgt_dir": "Rmin",
                      "repo_files": ["eq.mdp"]}},
         {"gromacs": "mdrun", #1
          "options": {"dir": "Rmin",
                      "src": "topol.tpr",
                      "tgt": "traj.trj",
                      "energy": "ener.edr",
                      "conf": "confout.gro",
                      "log": "md.log"}},
        ]
        self.breaks = {}

        if "debug" in args:
            self.recipe[0]["options"]["repo_files"] = ["eqDEBUG.mdp"]

class LigandMinimization(BasicMinimization):
    def __init__(self):
        super(LigandMinimization, self).__init__()

class BasicEquilibration(object):
    def __init__(self, args):
        self.recipe = \
        [{"gromacs": "editconf", #0
          "options": {"src": "Rmin/confout.gro",
                      "tgt": "min.pdb"}},
         {"command": "make_ndx", #1
          "options": {"src": "min.pdb",
                      "tgt": "index.ndx"}},
         {"gromacs": "grompp", #2
          "options": {"src": "Rmin/eq.mdp",
                      "src2": "min.pdb",
                      "top": "topol.top",
                      "tgt": "topol.tpr",
                      "index":"index.ndx"}},
         {"command": "set_stage_init", #3
          "options": {"src_dir": "Rmin",
                      "src_files": ["topol.tpr", "eq.mdp"],
                      "tgt_dir": "eq"}},
         {"command": "set_stage_init", #4
          "options": {"src_dir": "",
                      "src_files": ["posre.itp"],
                      "tgt_dir": "eq"}},
         {"gromacs": "mdrun", #5
          "options": {"dir": "eq",
                      "src": "topol.tpr",
                      "tgt": "traj.trj",
                      "energy": "ener.edr",
                      "conf": "confout.gro",
                      "traj": "traj.xtc",
                      "log": "md_eq1000.log"}},
        ]
        self.breaks = {}

        if "debug" in args:
            self.recipe[2]["options"]["src"] = "Rmin/eqDEBUG.mdp"
            self.recipe[3]["options"]["src_files"] =\
                ["topol.tpr", "eqDEBUG.mdp"]

class LigandEquilibration(BasicEquilibration):
    def __init__(self):
        super(LigandEquilibration, self).__init__()
        self.recipe[4]["options"]["src_files"].append("posre_lig.itp")
        self.recipe.insert(2,
            {"gromacs": "genrestr",
             "options": {"src": "Rmin/topol.tpr",
                         "tgt": "protein_ca200.itp",
                         "index": "index.ndx",
                         "forces": ["200", "200", "200"]},
             "input": "3\n"})

class BasicRelax(object):
    def __init__(self, *args):
        self.recipe = []
        for const in range(800, 0, -200):
            tgt_dir = "eq/{0}".format(const)
            src_dir = "eq"
            self.recipe += [
            {"command": "relax", #0, 3, 6, 9
             "options": {"const": const,
                         "src_dir": src_dir,
                         "tgt_dir": tgt_dir,
                         "mdp": "eq.mdp",
                         "posres": ["posre.itp"]}},
            {"gromacs": "grompp", #1, 4, 7, 10
             "options": {"src": os.path.join(tgt_dir, "eq.mdp"),
                         "src2": os.path.join(src_dir, "confout.gro"),
                         #top": os.path.join(tgt_dir, "topol.top"),
                         "top": "topol.top",
                         "tgt": os.path.join(tgt_dir, "topol.tpr"),
                         "index": "index.ndx"}},
             #TODO ese conf de abaixo hai que copialo, non vale asi
            {"gromacs": "mdrun", #2, 5, 8, 11
             "options": {"dir": tgt_dir,
                         "src": "topol.tpr",
                         "tgt": "traj.trr",
                         "energy": "ener.edr",
                         "conf": "../confout.gro",
                         "traj": "traj.xtc",
                         "log": "md_eq{0}.log".format(const)}},
            ]
        self.breaks = {}

        if "debug" in args:
            for i in range(0, 12, 3): #All relax lines
                self.recipe[i]["options"]["mdp"] = "eqDEBUG.mdp"

class LigandRelax(BasicRelax):
    def __init__(self):
        super(LigandRelax, self).__init__()
        self.recipe[4]["options"]["posres"].append("posre_lig.itp")

class CAEquilibrate(object):
    def __init__(self, *args):
        self.recipe = \
            [{"command": "set_stage_init", #0
              "options": {"src_dir": "eq",
                          "tgt_dir": "eqCA",
                          "src_files": ["confout.gro", "eq.mdp"]}},
             {"gromacs": "genrestr", #1
              "options": {"src": "Rmin/topol.tpr",
                          "tgt": "posre.itp",
                          "index": "index.ndx",
                          "forces": ["200"] * 3},
              "input": "3\n"},
             {"gromacs": "grompp", #2
              "options": {"src": "eqCA/eq.mdp",
                          "src2": "eqCA/confout.gro",
                          "top": "topol.top",
                          "tgt": "eqCA/topol.tpr",
                          "index": "index.ndx"}},
             {"gromacs": "mdrun", #3
              "options": {"dir": "eqCA",
                          "src": "topol.tpr",
                          "tgt": "traj.trr",
                          "energy": "ener.edr",
                          "conf": "confout.gro",
                          "traj": "traj.xtc",
                          "log": "md_eqCA.log"}},
             {"gromacs": "trjcat", #4
              "options": {"dir1": "eq",
                          "dir2": "eqCA",
                          "name": "traj.xtc",
                          "tgt": "traj_EQ.xtc"},
              "input": "c\n" * 6},
             {"gromacs": "eneconv", #5
              "options": {"dir1": "eq",
                          "dir2": "eqCA",
                          "name": "ener.edr",
                          "tgt": "ener_EQ.edr"},
              "input": "c\n" * 6}
             ]

        self.breaks = {}

        if "debug" in args:
            self.recipe[0]["options"]["src_files"] =\
                ["confout.gro", "eqDEBUG.mdp"]
            self.recipe[2]["options"]["src"] = "eqCA/eqDEBUG.mdp"

