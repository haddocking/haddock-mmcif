from pathlib import Path
import subprocess
import shlex
import os
from haddock2mmcif.modules.utils import backmap

AA_DICTIONARY = {
    "CYS": "C",
    "ASP": "D",
    "SER": "S",
    "GLN": "Q",
    "LYS": "K",
    "ILE": "I",
    "PRO": "P",
    "THR": "T",
    "PHE": "F",
    "ASN": "N",
    "GLY": "G",
    "HIS": "H",
    "LEU": "L",
    "ARG": "R",
    "TRP": "W",
    "ALA": "A",
    "VAL": "V",
    "GLU": "E",
    "TYR": "Y",
    "MET": "M",
}


TOOL_PATH = Path(Path(__file__).parents[3], "tools")


class PDB:
    def __init__(self, pdb_f):
        self.pdb_file = Path(pdb_f)
        self.atom_list = []
        self.seq_dic = {}
        self.map_dic = {}
        self.interface_dic = {}

    def load(self):
        pdb_dic = {}
        seq_id = 0
        with open(self.pdb_file) as fh:
            for line in fh.readlines():
                if line.startswith("ATOM"):
                    atom_id = line[12:16].strip()
                    chain = line[21]
                    resname = line[17:20]
                    resnum = int(line[22:26])
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    element = line[77:79].strip()
                    if chain not in pdb_dic:
                        pdb_dic[chain] = {}
                        self.map_dic[chain] = {}
                        seq_id = 0

                    if resnum not in pdb_dic[chain]:
                        pdb_dic[chain][resnum] = resname
                        seq_id += 1  # this needs to be before the assign below
                        self.map_dic[chain][seq_id] = resnum
                        # self.map_dic[chain][resnum] = seq_id

                    atom_data = chain, seq_id, element, atom_id, x, y, z
                    self.atom_list.append(atom_data)

        for chain in pdb_dic:
            self.seq_dic[chain] = []
            for three_letter_res in pdb_dic[chain].values():
                try:
                    one_letter_res = AA_DICTIONARY[three_letter_res]
                except KeyError:
                    # this aminoacid is not in the dictionary,
                    #  overwrite for now
                    one_letter_res = "A"
                self.seq_dic[chain].append(one_letter_res)

    def get_interface(self, cutoff=5.0):
        cmd = f"{TOOL_PATH}/contact-chainID {self.pdb_file} {cutoff}"
        out = subprocess.check_output(shlex.split(cmd))  # nosec
        for line in out.decode("utf-8").split(os.linesep):
            data = line.split()
            if data:
                res_a, chain_a, _, res_b, chain_b, _, _ = data
                res_a = int(res_a)
                res_b = int(res_b)

                if chain_a not in self.interface_dic:
                    self.interface_dic[chain_a] = []

                if chain_b not in self.interface_dic:
                    self.interface_dic[chain_b] = []

                mapped_res_a = backmap(self.map_dic[chain_a], res_a)
                if mapped_res_a not in self.interface_dic[chain_a]:
                    self.interface_dic[chain_a].append(mapped_res_a)

                mapped_res_b = backmap(self.map_dic[chain_b], res_b)
                if mapped_res_b not in self.interface_dic[chain_b]:
                    self.interface_dic[chain_b].append(mapped_res_b)
