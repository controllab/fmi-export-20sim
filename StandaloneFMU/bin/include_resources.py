# -*- coding: utf-8 -*-

"""
This script is called during the export of a 20-sim model to an FMU. It copies
data files used by a 20-sim model to the resources directory of the generated FMU.

The script expects two input arguments:
    1. The full path to the ModelConfiguration.xml file
    2. The full path to the resources directory
"""

import logging
import os
from shutil import copy
from sys import argv
import xml.etree.ElementTree as etree

DATA_FILE_EXTENSIONS = ('.txt', '.dat', '.csv')

class ScriptExit( Exception ): pass

def setup_logging(log_filename='output'):
    """
    Configure logging to write messages to console and file at the DEBUG level.
    """
    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
        datefmt='%Y-%m-%d %H:%M',
        filename='{}.log'.format(log_filename)
    )
    # Create a handler for dispatching log messages to the console
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.DEBUG)
    # Create a simpler format for the console
    formatter = logging.Formatter('%(name)-12s: %(levelname)-8s %(message)s')
    console_handler.setFormatter(formatter)
    # Register the console handler to the logging
    logging.getLogger('').addHandler(console_handler)

def check_arguments(argv):
    if len(argv) != 3:
        # First argument is always the script's name
        expected_args_count = len(argv) - 1
        logging.error('Incorrect number of arguments, expected two, received {}'.format(expected_args_count))
        return

if __name__ == "__main__":
    try:
        # TODO Choose appropriate location for this script's log file.
        setup_logging(argv[1])

        # Check arguments given to this script are as expected
        check_arguments(argv)

        model_conf = argv[1]
        res_dir = argv[2]

        if not os.path.exists(model_conf):
            logging.error('Model configuration file not found: {}'.format(model_conf))
            raise ScriptExit
        if not os.path.isdir(res_dir):
            logging.error('Resources directory not found: {}'.format(res_dir))
            raise ScriptExit

        # Get the location on disk of the 20-sim model from the ModelConfiguration.xml
        tree = etree.parse(model_conf)
        root = tree.getroot()
        originalModelNode = root.findall('.//originalModel')
    
        if len(originalModelNode) != 1:
            logging.error('Invalid content, unable to find model path information in {}'.format(model_conf))
            raise ScriptExit
        
        model_name = originalModelNode[0].text
        model_location = os.path.dirname(model_name)

        # Find data files specified in the ModelConfiguration.xml as model variables
        data_files = list()
        for model_var in root.iter('modelVariable'):
            type = model_var.find('type').text
            if type == 'string':
                value = model_var.find('value').text
                if value.endswith(DATA_FILE_EXTENSIONS):
                    # A data file without a path is assumed to be in the same dir as the
                    # model, so prepend the model's location to it
                    if not os.path.dirname(value):
                        value = ''.join((model_location, os.sep, value))
            
                    data_files.append(value)

        # Copy found data files to resources directory
        for df in data_files:
            new_file = copy(df, res_dir)
            if not new_file:
                logging.error('Unable to copy to data file {} to the resources directory {}'.format(df, res_dir))
    except ScriptExit:
        pass
    except Exception as err:
        logging.error('An unexpected error occurred: {}\n'.format(err))

