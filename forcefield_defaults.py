import os

ROOT_DIR = os.path.dirname(os.path.realpath(__file__))
TEMPLATES_ROOT = os.path.join(ROOT_DIR, "templates")

FORCEFIELD_DEFAULTS = {
    "oplsaa": {
        "proteinff": "oplsaa",
        "waterff": "spc",
        "lipidff": "berger",
        "ligandff": "ligpargen",
        "templates_dir": os.path.join(TEMPLATES_ROOT, "oplsaa"),
        "ff_nb": os.path.join(TEMPLATES_ROOT, "oplsaa/ffoplsaanb_base.itp"),
        "set_grompp_files": { "steep.mdp": "steep.mdp",
                                        "ff_files": ["ions.itp", "spc.itp", "ffoplsaanb_base.itp", "ffoplsaabon_mod.itp","ffoplsaa_mod.itp", "ffoplsaanb_lip.itp"]},
    },
    "oplsaam": {
        "proteinff": "oplsaam",
        "waterff": "tip3p_opls",
        "lipidff": "berger",
        "ligandff": "ligpargen",
        "templates_dir": os.path.join(TEMPLATES_ROOT, "oplsaam"),
        "ff_nb": os.path.join(TEMPLATES_ROOT, "oplsaam/ffoplsaamnb_base.itp"),
        "set_grompp_files": { "steep.mdp": "steep.mdp",
                             "ff_files": ["ions.itp", "tip3p_opls.itp", "ffoplsaamnb_base.itp", "ffoplsaambon_mod.itp", "ffoplsaam_mod.itp", "ffoplsaamnb_lip.itp"]}
    },
    "amber14sb": {
        "proteinff": "amber14sb",
        "waterff": "tip3p_amber",
        "lipidff": "amberlipid21",
        "ligandff": "openff-2.3.0.offxml",
        "templates_dir": os.path.join(TEMPLATES_ROOT, "amber14sb"),
        "ff_nb": os.path.join(TEMPLATES_ROOT, "amber14sb/ffnonbonded.itp"),
        "set_grompp_files": { "steep.mdp": "steep.mdp",
                             "ff_files": ["ions.itp", "tip3p_amber.itp", "ffbonded.itp", "forcefield.itp", "amberlipid21_atomtypes.itp"]}
    }
}

LIPIDFF_TEMPLATES = {
    "amberlipid21": {
        "templates_dir": os.path.join(TEMPLATES_ROOT, "membrane/amberlipid21"),
        "set_grompp_files": { "popc.itp": "popc.itp" }
    },
    "berger": {
        "templates_dir": os.path.join(TEMPLATES_ROOT, "membrane/berger"),
        "set_grompp_files": { "popc.itp": "popc.itp" }
    }
}

WATERFF_TEMPLATES = {
    # The chosen water file is copied as water.itp, which every topol.top includes.
    "spc": {
        "templates_dir": os.path.join(TEMPLATES_ROOT, "water"),
        "set_grompp_files": {"water.itp": "spc.itp"}
    },
    "tip3p_opls": {
        "templates_dir": os.path.join(TEMPLATES_ROOT, "water"),
        "set_grompp_files": {"water.itp": "tip3p_opls.itp"}
    },
    "tip3p_amber": {
        "templates_dir": os.path.join(TEMPLATES_ROOT, "water"),
        "set_grompp_files": {"water.itp": "tip3p_amber.itp"}
    }
}

def get_forcefield_defaults(proteinff):
    if proteinff is None:
        return None
    return FORCEFIELD_DEFAULTS.get(proteinff.lower())


def get_set_grompp_files(proteinff=None, lipidff=None, waterff=None):
    """
    Return a mapping of files that should be copied for set_grompp
    based on the selected protein, lipid, and water forcefields.
    """
    file_sources = {}

    def add_file(dest, src, templates_dir):
        file_sources[dest] = {"src": src, "dir": templates_dir}

    def normalize_set_grompp(entry, templates_dir):
        if not entry:
            return
        for key, value in entry.items():
            if key == "ff_files" and isinstance(value, list):
                for filename in value:
                    add_file(filename, filename, templates_dir)
            else:
                add_file(key, value, templates_dir)

    protein = get_forcefield_defaults(proteinff)
    if protein:
        normalize_set_grompp(protein.get("set_grompp_files", {}), protein["templates_dir"])

    if lipidff:
        lipid = LIPIDFF_TEMPLATES.get(lipidff.lower())
        if lipid:
            normalize_set_grompp(lipid.get("set_grompp_files", {}), lipid["templates_dir"])

    if waterff:
        water = WATERFF_TEMPLATES.get(waterff.lower())
        if water:
            normalize_set_grompp(water.get("set_grompp_files", {}), water["templates_dir"])

    return file_sources
