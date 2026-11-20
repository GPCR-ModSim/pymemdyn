import os

import forcefield_defaults
import gromacs


class Membrane(object):
    """
    Set the characteristics of the membrane in the complex.
    """
    def __init__(self, *args, **kwargs):
        lipidforcefield = kwargs.get("lipidforcefield") or "berger"
        lipid_template = forcefield_defaults.LIPIDFF_TEMPLATES.get(lipidforcefield.lower())
        if lipid_template:
            repo_dir = lipid_template["templates_dir"]
        else:
            repo_dir = kwargs.get("repo_dir") or gromacs.Wrapper().repo_dir

        # Bilayer PDB and sizes (Angstrom) selected by lipid forcefield.
        # x, y, zw: box of the bilayer with waters
        # z: bilayer without waters (max - min z coordinate of the lipids)
        bilayer_map = {
            "berger": {"pdb": "x4bilayer.pdb",
                       "x": 61.842, "y": 61.842, "zw": 154.053,
                       "z": 57.78},   # Membrane thickness: 39.71
            "amberlipid21": {"pdb": "lipid21bilayer.pdb",
                             "x": 64.329, "y": 64.329, "zw": 153.242,
                             "z": 56.36}   # Membrane thickness: 38.76
        }
        
        # Unknown lipid forcefield: fall back to the Berger bilayer
        bilayer = bilayer_map.get(lipidforcefield.lower(), bilayer_map["berger"])

        if "pdb" not in kwargs.keys():
            self.pdb = os.path.join(repo_dir, bilayer["pdb"])
        else:
            self.pdb = kwargs["pdb"]

        self.bilayer_x = bilayer["x"]    # xbilayer
        self.bilayer_y = bilayer["y"]    # ybilayer
        self.bilayer_zw = bilayer["zw"]  # hbilayerw i.e. with waters
        self.bilayer_z = bilayer["z"]    # zbilayer i.e. without waters
        self.set_nanom()

    def set_nanom(self):
        """
        Convert some measurements to nanometers to comply with GROMACS units.
        """
        nanometer = 10
        self.gmx_bilayer_x = self.bilayer_x / nanometer
        self.gmx_bilayer_y = self.bilayer_y / nanometer
        self.gmx_bilayer_z = self.bilayer_z / (nanometer - 1) #Loose the lipids
        self.gmx_bilayer_zw = self.bilayer_zw / nanometer

        return True