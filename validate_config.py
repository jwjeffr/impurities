#!/usr/bin/python

import json
import sys


def main():

    config_file_name = sys.argv[1]

    with open(config_file_name, 'r') as file:
        config = json.load(file)
    
    # first, validate that subdictionaries with systems as dictionaries only have keys from system

    for key in ["Structure", "Dominant Order Parameter", "Atoms", "System Colors", "System Line Styles", "Type Maps"]:
        assert list(config[key].keys()) == config["Systems"]

    # then, make sure that atom lists are consistent for each system

    for system in config["Systems"]:

        assert config["Atoms"][system] == list(config["Type Maps"][system].values())
    
    # then, check that all atom types have a defined color

    all_types = [config["Atoms"][system] for system in config["Systems"]]
    
    # flatten list
    all_types = [item for sublist in all_types for item in sublist]
    all_types = list(set(all_types))

    for type_ in all_types:
        _ = config["Atom Colors"][type_]

    # lastly, check that all atom properties have the same keys
    assert config["Atom Colors"].keys() == config["Atom Radii"].keys() == config["Atom Abbreviations"].keys()

    print(f'{config_file_name} is a valid config')

if __name__ == '__main__':

    main()
