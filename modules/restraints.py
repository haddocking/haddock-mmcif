import re
import logging
restlog = logging.getLogger('log')


class UnambigRestraint():

    def __init__(self, tbl_file):
        self.tbl_file = tbl_file
        self.tbl_list = []

    def load(self):
        """Load the unambiguous restraints."""
        resid_regex = r"resid\s(\d*)"
        segid_regex = r"segid\s*(\w*)"
        distances_regex = r"(\d*.?\d*)\s(\d*.?\d*)\s(\d*.?\d*)$"
        with open(self.tbl_file, 'r') as fh:
            for line in fh.readlines():
                if line.startswith('assign'):
                    res_i, res_j = re.findall(resid_regex, line)
                    segid_i, segid_j = re.findall(segid_regex, line)
                    distances = re.findall(distances_regex, line)[0]
                    distance, lower_bound, upper_bound = distances
                    restraint = (int(res_i), segid_i, int(res_j), segid_j,
                                 float(distance), float(lower_bound),
                                 float(upper_bound))
                    self.tbl_list.append(restraint)


class AmbigRestraint():

    def __init__(self, tbl_file):
        self.tbl_file = tbl_file
        self.tbl_dic = {}

    def load(self):
        """Load the ambig restraints."""
        tbl_regex = r"resid\s*(\d*).*segid\s*(\w*)"
        with open(self.tbl_file, 'r') as fh:
            for line in fh.readlines():
                if line.startswith('assign'):
                    # this is a special line
                    # do something special
                    match = re.search(tbl_regex, line)
                    active_res = int(match.group(1))
                    active_segid = str(match.group(2))
                    active = (active_res, active_segid)
                    self.tbl_dic[active] = []
                else:
                    # can be a valid line or not
                    match = re.search(tbl_regex, line)
                    if match:
                        pas_res = int(match.group(1))
                        pas_segid = str(match.group(2))
                        passive = (pas_res, pas_segid)
                        self.tbl_dic[active].append(passive)
