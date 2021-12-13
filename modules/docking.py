import ihm.model
import logging
dockinglog = logging.getLogger('log')


class DockingModel(ihm.model.Model):
    """Subclass to save memory."""
    # ======================================================================
    # IMPORTANT #
    # To add the atoms, the class module needs a list containing the following:
    #  [(<ihm.AsymUnit object...108edf5b0>, 1, 'C', 'CA', 1.0, 2.0, 3.0), ...]
    # Which means that the AsymUnit object will be copied many times and this
    #  will eat a lot of memory.
    # To avoid this we subclass IHM Model class and override get_atoms function
    # ======================================================================
    def __init__(self, assymetric_dic, atom_list, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.asym_unit_map = assymetric_dic
        self.atom_list = atom_list

    def get_atoms(self):
        for asym, seq_id, type_symbol, atom_id, x, y, z in self.atom_list:
            yield ihm.model.Atom(asym_unit=self.asym_unit_map[asym],
                                 type_symbol=type_symbol, seq_id=seq_id,
                                 atom_id=atom_id, x=x, y=y, z=z)
