#!/usr/bin/env python
#



import logging
import os 

def setup(config):
    '''
    Builds directories in config file 
    Download the appropriate supplement data.
    Only needs to be done once.
    '''

    log = logging.getLogger('biccn')
    config = config
    # directories

    resourcedir = os.path.expanduser(config.get('setup', 'resourcedir'))
    outputdir = os.path.expanduser(config.get('setup', 'outputdir'))
    figuredir =os.path.expanduser(config.get('setup', 'figuredir'))

    dirs_to_make = [metadir, 
                    f'{cachedir}/biccn', 
                    tempdir, 
                    resourcedir, 
                    outputdir, 
                    figuredir]

    for direc in dirs_to_make:
        try:
            os.makedirs(direc)
        except FileExistsError:
            pass