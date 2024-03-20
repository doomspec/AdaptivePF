from math import sin, cos, pi

"""
This file provides the function to produce geometry object for the chemistry packages to use.
H2, LiH, H2O, N2 are provided.
"""


def get_H2_geo(bond_len):
    atom_1 = 'H'
    atom_2 = 'H'
    coordinate_1 = (0.0, 0.0, 0.0)
    coordinate_2 = (bond_len, 0.0, 0.0)
    geometry = [(atom_1, coordinate_1), (atom_2, coordinate_2)]
    return geometry


def get_H6_geo(bond_len):
    atom_1 = 'H'
    atom_2 = 'H'
    atom_3 = 'H'
    atom_4 = 'H'
    atom_5 = 'H'
    atom_6 = 'H'
    coordinate_1 = (0.0, 0.0, 0.0)
    coordinate_2 = (bond_len, 0.0, 0.0)
    coordinate_3 = (bond_len * 2, 0.0, 0.0)
    coordinate_4 = (bond_len * 3, 0.0, 0.0)
    coordinate_5 = (bond_len * 4, 0.0, 0.0)
    coordinate_6 = (bond_len * 5, 0.0, 0.0)
    geometry = [(atom_1, coordinate_1), (atom_2, coordinate_2), (atom_3, coordinate_3),
                (atom_4, coordinate_4), (atom_5, coordinate_5), (atom_6, coordinate_6)]
    return geometry


def get_H4_geo(bond_len):
    atom_1 = 'H'
    atom_2 = 'H'
    atom_3 = 'H'
    atom_4 = 'H'
    coordinate_1 = (0.0, 0.0, 0.0)
    coordinate_2 = (bond_len, 0.0, 0.0)
    coordinate_3 = (bond_len * 2, 0.0, 0.0)
    coordinate_4 = (bond_len * 3, 0.0, 0.0)
    geometry = [(atom_1, coordinate_1), (atom_2, coordinate_2), (atom_3, coordinate_3),
                (atom_4, coordinate_4)]
    return geometry


def get_LiH_geo(bond_len):
    atom_1 = 'H'
    atom_2 = 'Li'
    coordinate_1 = (0.0, 0.0, 0.0)
    coordinate_2 = (bond_len, 0.0, 0.0)
    geometry = [(atom_1, coordinate_1), (atom_2, coordinate_2)]
    return geometry


def get_N2_geo(bond_len):
    atom_1 = 'N'
    atom_2 = 'N'
    coordinate_1 = (0.0, 0.0, 0.0)
    coordinate_2 = (bond_len, 0.0, 0.0)
    geometry = [(atom_1, coordinate_1), (atom_2, coordinate_2)]
    return geometry


def get_H2O_geo(bond_len):
    atom_1 = 'H'
    atom_2 = 'H'
    atom_3 = 'O'
    coordinate_1 = (0.0, 0.0, 0.0)
    coordinate_2 = (2 * cos(104.5 / 180 * pi) * bond_len, 0.0, 0.0)
    coordinate_3 = (1 * cos(104.5 / 180 * pi) * bond_len,
                    sin(104.5 / 180 * pi) * bond_len, 0.0)
    geometry = [(atom_1, coordinate_1), (atom_2, coordinate_2),
                (atom_3, coordinate_3)]
    return geometry


def get_H2S_geo(bond_len):
    atom_1 = 'H'
    atom_2 = 'H'
    atom_3 = 'S'
    coordinate_1 = (0.0, 0.0, 0.0)
    coordinate_2 = (2 * cos(104.5 / 180 * pi) * bond_len, 0.0, 0.0)
    coordinate_3 = (1 * cos(104.5 / 180 * pi) * bond_len,
                    sin(104.5 / 180 * pi) * bond_len, 0.0)
    geometry = [(atom_1, coordinate_1), (atom_2, coordinate_2),
                (atom_3, coordinate_3)]
    return geometry


def get_CH2_geo(bond_len):
    atom_1 = 'H'
    atom_2 = 'H'
    atom_3 = 'C'
    coordinate_1 = (0.0, 0.0, 0.0)
    coordinate_2 = (2 * cos(104.5 / 180 * pi) * bond_len, 0.0, 0.0)
    coordinate_3 = (1 * cos(104.5 / 180 * pi) * bond_len,
                    sin(104.5 / 180 * pi) * bond_len, 0.0)
    geometry = [(atom_1, coordinate_1), (atom_2, coordinate_2),
                (atom_3, coordinate_3)]
    return geometry


def get_C2H2_geo(bond_len):
    atom_1 = 'H'
    atom_2 = 'C'
    atom_3 = 'C'
    atom_4 = 'H'
    coordinate_1 = (0.0, 0.0, 0.0)
    coordinate_2 = (1.06, 0.0, 0.0)
    coordinate_3 = (bond_len + 1.06, 0.0, 0.0)
    coordinate_4 = (bond_len + 2.12, 0.0, 0.0)
    geometry = [(atom_1, coordinate_1), (atom_2, coordinate_2),
                (atom_3, coordinate_3), (atom_4, coordinate_4)]
    return geometry


equilibrium_geometry_dict = {"H2": 0.74, "H6": 1.5, "H4": 1.5, "LiH": 1.4, "H2O": 0.96, "N2": 1.1}
geometry_generator_dict = {"H2": get_H2_geo, "H6": get_H6_geo, "H4": get_H4_geo, "LiH": get_LiH_geo, "H2O": get_H2O_geo,
                           "N2": get_N2_geo}
