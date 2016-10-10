#!/usr/bin/env python

#  This file is part of MAUS: http://micewww.pp.rl.ac.uk:8080/projects/maus
#
#  MAUS is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  MAUS is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with MAUS.  If not, see <http://www.gnu.org/licenses/>.

import os
import sys
import math
import json
import ROOT
import numpy.random
import copy

import xboa.common
from xboa.bunch import Bunch
from xboa.hit import Hit

import maus_cpp.simulation
import maus_cpp.globals
import maus_cpp.field

class Tracking(object): # pylint: disable=R0902
    """Plotting class handles generation of particles and tracking"""
    def __init__(self, analysis_config, maus_config):
        """
        Initialise the tracking class
        - configuration should be a valid datacard set in string representation
        """
        self.analysis_config = analysis_config
        self.maus_config = json.loads(maus_config)
        self.lattice_file = self.analysis_config.geometry
        pid = abs(self.analysis_config.fit_pid)
        self.m_mu = xboa.common.pdg_pid_to_mass[pid]
        self.z_start = self.analysis_config.z_tku-1e-9
        self.pz = None
        #self.bz_us = 0.e-3
        #self.beta = 333.
        self.verbose_level = 5
        self.modules = maus_cpp.mice_module.MiceModule(self.lattice_file)

    def scale_mice_module(self, mod_name, scale_factor):
        if not self._recursive_set_scale_factor(mod_name, self.modules, scale_factor):
            raise KeyError("Failed to find module "+str(mod_name))

    def reset_modules(self):
        self.modules = maus_cpp.mice_module.MiceModule(self.lattice_file)

    def misalign_mice_module(self, dx, dxp, dy, dyp, dz, mod_name):
        self._recursive_misalign_mice_module(dx, dxp, dy, dyp, dz, mod_name, self.modules)

    def _recursive_misalign_mice_module(self, dx, dxp, dy, dyp, dz, mod_name, mice_mod):
        if mice_mod.get_name() == mod_name:
            position = mice_mod.get_property("Position", "hep3vector")
            position["x"] += dx
            position["y"] += dy
            position["z"] += dz
            mice_mod.set_property("Position", "hep3vector", position)
            rotation = mice_mod.get_property("Rotation", "hep3vector")
            rotation["x"] += dxp
            rotation["y"] += dyp
            mice_mod.set_property("Rotation", "hep3vector", rotation)
        children = [self._recursive_misalign_mice_module(dx, dxp, dy, dyp, dz, mod_name, mod) for mod in mice_mod.get_children()]
        mice_mod.set_children(children)
        return mice_mod

    def get_scale_factor(self, mod_name):
        return self._recursive_get_scale_factor(mod_name, self.modules)

    def _recursive_get_scale_factor(self, mod_name, mice_mod):
        if mice_mod.get_name() == mod_name:
            try:
                return mice_mod.get_property("ScaleFactor", "double")
            except KeyError:
                return 1.
        else:
            for mod in mice_mod.get_children():
                scale = self._recursive_get_scale_factor(mod_name, mod)
                if scale != None:
                    return scale
            return None                  

    def _recursive_set_scale_factor(self, mod_name, mice_mod, scale_factor):
        if mice_mod.get_name() == mod_name:
            pol = mice_mod.get_property("ScaleFactor", "double")
            pol = pol/abs(pol)
            mice_mod.set_property("ScaleFactor", "double", scale_factor*pol)
            return True
        children = mice_mod.get_children()
        for mod in children:
            if self._recursive_set_scale_factor(mod_name, mod, scale_factor):
                mice_mod.set_children(children)
                return True
        return False


    def do_tracking(self, delta_x, delta_y, delta_px, delta_py, print_fields = False):
        """Do the tracking"""
        print "Tracking setup"
        if maus_cpp.globals.has_instance():
            maus_cpp.globals.death()
        self._fiddle_configuration()
        maus_cpp.globals.set_monte_carlo_mice_modules(self.modules)
        hit_list = self.make_hit_list([delta_x, delta_y, delta_px, delta_py])
        print "Tracking", len(hit_list), "particles..."
        pos = self.track_particle(hit_list)
        print "Tracked particles"
        if print_fields:
            print maus_cpp.field.str(True)
        return pos

    def make_hit_list(self, deltas):
        defaults = {'pz':self.pz, 'mass':self.m_mu, 'charge':1., 'pid':-13, 'z':self.z_start}
        hit_list = [Hit.new_from_dict(defaults, 'energy')]
        delta_keys = ['x', 'y', 'px', 'py']
        for parity in [-1., +1.]:
            for i in range(4):
                defaults[delta_keys[i]] = parity*deltas[i]
                hit_list.append(Hit.new_from_dict(defaults, 'energy'))
                defaults[delta_keys[i]] = 0.
        return hit_list

    def track_particle(self, hit_list):
        """Track a particle with initial position x,y"""
        primaries = [self._get_primary(hit) for hit in hit_list]
        event_list = json.loads(maus_cpp.simulation.track_particles(json.dumps(primaries)))
        pos = [[hit["position"] for hit in event["virtual_hits"]] for event in event_list]
        for primary in primaries:
            primary['primary']['momentum']['x'] *= 1.
            primary['primary']['momentum']['y'] *= 1.
            primary['primary']['momentum']['z'] *= 1.
        event_list = json.loads(maus_cpp.simulation.track_particles(json.dumps(primaries)))
        event_list = [event["virtual_hits"] for event in event_list]
        return event_list

    def _fiddle_configuration(self):
        """Switch off physics processes and force keep steps to True"""
        self.maus_config["physics_processes"] = "none"
        self.maus_config["keep_steps"] = True
        self.maus_config["verbose_level"] = self.verbose_level
        maus_cpp.globals.birth(json.dumps(self.maus_config))

    def _get_primary(self, hit):
        """Build a primary with initial position (x, y, z_start) and large
        momentum in the z-direction"""
        return {
          'primary':{
            'position':{'x':hit['x'], 'y':hit['y'], 'z':hit['z']},
            'momentum':{'x':hit['px'], 'y':hit['py'], 'z':hit['pz']},
            'particle_id':hit['pid'],
            'energy':hit['energy'],
            'time':hit['t'],
            'random_seed':0
          }
        }


