"""Interface that gives useful infos, stereochemistry and chemical groups from a compound's name in the clipboard."""

from .molinterface import find_smiles_patterns, find_rings, format_rings, open_url,find_chemical_groups, get_compound_info, display_image1_in_tkinter, display_image2_in_tkinter, highlight_chemical_groups, update_gui_labels, update_molecule_info, clipboard_monitor, toggle_image_display, main, smarts_patterns


__version__ = "0.0.1"
