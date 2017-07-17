# coding: utf-8

from __future__ import unicode_literals, division

import subprocess
import shutil
import os

from custodian.custodian import Job

__author__ = "Maja-Olivia Lenz"
__version__ = "0.1"
__maintainer__ = "Maja-Olivia Lenz"
__email__ = "lenz@fhi-berlin.mpg.de"
__status__ = "Development"
__date__ = "13/07/17"

class FHIaimsJob(Job):
    """
    A basic FHI-aims Job.
    """
    
    def __init__(self,aims_cmd, control_in='control.in', geometry_in='geometry.in',
                 output_file='aims.out', backup=True):
        """
        Initialized a basic FHIaims job.
        """
        self.aims_cmd = aims_cmd
        self.control_in = control_in
        self.geometry_in = geometry_in
        self.output_file = output_file
        self.backup = backup

    def setup(self):
        """
        Allows for pre-processing. Here: performs backup if needed.
        """
        if self.backup:
            shutil.copy(self.control_in, '{}.orig'.format(self.control_in))
            shutil.copy(self.geometry_in, '{}.orig'.format(self.geometry_in))

    def run(self):
        """
        Performs FHI aims job run.
        """
        if self._check_success():
            return None

        with open(self.output_file,'w') as fout:
            p = subprocess.Popen(self.aims_cmd, stdout=fout)
        return p

    def postprocess(self):
        """
        Postprocessing like renaming files, zip files etc. could be done
        """
        pass

    def _check_success(self):
        if os.path.exists(self.output_file):
            with open(self.output_file,'r') as f:
                for line in f:
                    if 'Have a nice day' in line:
                        return True
        return False
