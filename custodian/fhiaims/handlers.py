#from __future__ import unicode_literals, division

from custodian.custodian import ErrorHandler
from custodian.utils import backup

from ase.io import read
from ase.io.aims import write_aims
from ase.calculators.aims import Aims
from ase.calculators.calculator import Parameters

import numpy as np
from collections import Counter

__author__ = "Maja-Olivia Lenz"
__version__ = "0.1"
__maintainer__ = "Maja-Olivia Lenz"
__email__ = "lenz@fhi-berlin.mpg.de"
__status__ = "Development"
__date__ = "13/07/17"

class FHIaimsErrorHandler(ErrorHandler):
    """
    Error handler for FHIaims jobs.
    Does not use ansible/interpreter/Modder but fixes in place.
    This error handler is monitoring the running job for errors.
    """
    is_monitor = False

    error_msgs = {
        "trusted_descent" : "trusted_descent",
        "check_cpu_consistency" : "check_cpu_consistency_matrix", 
        "check_for_close_encounters" : "check_for_close_encounters",
        "scf_solver": "scf_solver: SCF cycle not converged.",
        "ill-conditioned": "Attention: Ill-conditioned overlap matrix detected",
        "too_many_neighbors": "get_neighbor_lists: Too many neighbors"
    }

    def __init__(self,output_file='aims.out',control_in='control.in',
                 geometry_in='geometry.in', param_file='parameters.ase'):
        """
        Initializes with output file name where errors are read from.
        """
        self.output_file = output_file
        self.control_in = control_in
        self.geometry_in = geometry_in
        self.param_file = param_file

        self.errors = set()
        self.error_count = Counter()

    def check(self):
        """
        Checks output file for errors. Returns Boolean.
        """
        self.errors = set()
        with open(self.output_file,'r') as f:
            for line in f:
                for err,msg in FHIaimsErrorHandler.error_msgs.items():
                    if msg in line:
                        self.errors.add(err)
        return len(self.errors) > 0

    def correct(self):
        """
        Corrects detected errors. Returns dict.
        """
        backup([self.output_file, self.geometry_in, self.control_in, self.param_file])
        actions = []
       
        for e in self.errors: 
            if "trusted_descent" == e:
                act = self._fix_trusted_descent(self.error_count[e])
                actions.append(act)
                self.error_count[e] += 1
            elif "check_cpu_consistency" == e:
                if self.error_count[e] == 0:
                    self._set_control('check_cpu_consistency', '.false.')
                    actions.append('check_cpu_consistency false')
                else:
                    actions.append(None)
                self.error_count[e] += 1
            elif "check_for_close_encounters" == e:
                # duplicate parameters, atoms too close?
                act = self._fix_close_encounters(self.error_count[e])
                actions.append(act)
                self.error_count[e] += 1
            elif "scf_solver" == e:
                act = self._fix_scf_not_converged(self.error_count[e])
                actions.append(act)
                self.error_count[e] += 1
            elif "ill-conditioned" == e:
                act = self._fix_ill_conditioned(self.error_count[e])
                actions.append(act)
                self.error_count[e] += 1 
            elif "too_many_neighbors" == e:
                act = self._fix_close_encounters(self.error_count[e])
                actions.append(act)
                self.error_count[e] += 1
            else:
                # Unimplemented error
                actions.append(None)
 
        return {"errors": list(self.errors), "actions": actions}

    # ======= internal functions for error fixing ===================

    def _fix_scf_not_converged(self, error_count):
        if error_count == 0:
            key = 'sc_iter_limit'
            modifier = 2.5
            action = 'Increase max. limit of scf iterations'
            self._modify_control(key,modifier)
        elif error_count == 1:
            keys = ('mixer', 'charge_mix_param', 'occupation_type')
            vals = ('pulay', 0.05, 'gaussian 0.05')
            action = 'Correction for metals 1'
            self._set_control(keys,vals)            
        elif error_count in range(2,4):
            mix_param = 10**(-error_count-1)
            keys = ('mixer', 'charge_mix_param')
            vals = ('linear', mix_param)
            action = 'Correction for metals 2'
            self._set_control(keys,vals)            
        elif error_count == 4:
            #  Standard strategy for metals failed, try Slab
            keys = ('mixer','charge_mix_param','occupation_type','preconditioner')
            vals = ('pulay', 0.2, 'gaussian 0.01', 'kerker off')
            action = 'Slab'
            self._set_control(keys,vals)            
        else:
            return None

        return action

    def _fix_ill_conditioned(self, error_count):
        factor = 1.5
        if error_count == 0:
            # increase whole cell if material too dense
            atoms = read(self.geometry_in, format='aims')
            cell = atoms.get_cell() * factor
            atoms.set_cell(cell)
            # write new geometry file
            write_aims(self.geometry_in,atoms)
            return 'increased cell'
        elif error_count == 1:
            # restore original cell size
            atoms = read(self.geometry_in, format='aims')
            cell = atoms.get_cell() * 1.0 / factor
            atoms.set_cell(cell)
            write_aims(self.geometry_in,atoms)
            # try instead:
            keys = ('basis_threshold','override_illconditioning')
            vals = (1e-6, '.true.')
            self._set_control(keys, vals)
            return 'override ill-conditioning'

    def _fix_trusted_descent(self, error_count):
        if error_count >= 4:
            return None
        length = 0.015-error_count*0.005    # default 0.025
        self._set_control('harmonic_length_scale',length)
        return 'increased harmonic_length_scale'

    def _fix_close_encounters(self, error_count):
        if error_count >= 4:
            return None
        atoms = read(self.geometry_in,format='aims')

        # find smallest distance between two atoms and indices
        dists = atoms.get_all_distances()
        dists = np.triu(dists)  # only upper triag. because of symmetry
        smallest = np.min( dists[np.nonzero(dists)] ) # smallest non zero dist.
        indices = zip(*np.where( dists == smallest ))

        # increase smallest distance by factor, keep center of bond (fix=0.5)
        factor = 1.5
        for ind in indices:
            atoms.set_distance(ind[0], ind[1], smallest*factor, fix=0.5)

        # write new geometry file
        write_aims(self.geometry_in,atoms)
        return 'increased smallest distance between atoms'
           
    def _set_control(self, keys, values):
        # set a value in control in file        

        # read parameters.ase
        parameters = Parameters.read(self.param_file)

        # MODIFY 
        if isinstance(keys,(tuple,list)):
            for key, value in zip(keys, values):
                parameters[key] = value
        elif isinstance(keys,str):
            parameters[keys] = values
        else:
            print 'ERROR in _set_control: key_list must be string, tuple or list'
            return None
        
        # set calculator object with new parameters
        calc = Aims(xc='',species_dir='')
        calc.parameters = parameters
        # read atoms object to write new input
        atoms = read(self.geometry_in, format='aims') 
        # write new parameters file:
        calc.parameters.write(self.param_file)
        # write modified control.in
        Aims.write_control(calc,atoms,self.control_in)
        Aims.write_species(calc,atoms,self.control_in)

    def _modify_control(self, key, modifier, mode='multiply')
        # by default multiplies the key's value by modifier

        # read parameters.ase
        parameters = Parameters.read(self.param_file)

        if mode == 'multiply':
            parameter[key] = parameter[key] * modifier

        # set calculator object with new parameters
        calc = Aims(xc='',species_dir='')
        calc.parameters = parameters
        # read atoms object to write new input
        atoms = read(self.geometry_in, format='aims')
        # write new parameters file:
        calc.parameters.write(self.param_file)
        # write modified control.in
        Aims.write_control(calc,atoms,self.control_in)
        Aims.write_species(calc,atoms,self.control_in)
