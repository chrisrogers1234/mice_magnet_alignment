#!/usr/bin/env python

#  This file is part of MAUS: http://micewww.pp.rl.ac.uk/projects/maus
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

#pylint: disable = W0613, E1101, R0913, R0914

"""
Example to demonstrate use of the optics routines for matching.

In this example, MAUS will build a set of solenoid field maps for MICE Step IV 
then find a "matched" solution for the field maps, given coil current limits. We
use the ROOT TMinuit routines for this.

"Matched" solution means that the beta function is flat in the upstream solenoid
and symmetric about the centre of the focus coil (flat at z=0 mm).

Two classes are defined:
  - Lattice handles the Lattice things (controls magnet currents, wraps optics
    routines)
  - LatticeOptimiser handles the optimisation
"""

import json
import os
import math

import numpy
import ROOT

import xboa.Common

import maus_cpp.field
import maus_cpp.globals
import maus_cpp.simulation
import maus_cpp.polynomial_map
import maus_cpp.covariance_matrix
import Configuration

class Lattice:
    """
    Magnetic lattice and functions to generate and plot transfer matrices on
    said lattice
    """

    def __init__(self, momentum, z_start, coil_currents, lattice_file):
        """
        Initialises some data

        - momentum: beam momentum for the transfer matrix
        - z_start: start z position of the beam
        - coil_currents: python dictionary that maps the MiceModule name of the
                         coil to the scale factor that should be applied
        """
        self.mom = momentum
        self.z_in = z_start
        self.mass = xboa.Common.pdg_pid_to_mass[13]
        self.optics_model = None
        self.first = True
        self.coil_currents = coil_currents
        self.lattice_file = lattice_file
        self.initial_beta = None
        self.initial_alpha = 0.
        self.tm_list = []
        self.z_list = []
        self.setup_lattice(coil_currents)

    def setup_lattice(self, coil_currents):
        """
        Setup the lattice for optics work

        - coil_currents: python dictionary that maps the MiceModule name
          of the coil to the scale factor that should be applied

        Sets up data cards and loads MAUS with those data cards, then builds the
        geometry and field model, then builds the transfer matrices
        """
        config_str = Configuration.Configuration().\
                                          getConfigJSON(command_line_args=False)
        config_json = json.loads(config_str)
        config_json["simulation_reference_particle"] = {
          "random_seed": 0,
          "energy":(self.mass**2+self.mom**2)**0.5,
          "particle_id":-13,
          "time": 0.0,
          "position":{"x":0.0, "y":0.0, "z":self.z_in},
          "momentum":{"x":0.0, "y":0.0, "z":1.0}
        }
        config_json["simulation_geometry_filename"] = self.lattice_file
        config_json["physics_processes"] = "none"
        config_json["verbose_level"] = 1
        config_str = json.dumps(config_json)
        if maus_cpp.globals.has_instance():
            maus_cpp.globals.death()
        maus_cpp.globals.birth(config_str)
        # now set up the mice modules
        self.set_magnet_scale_factors(coil_currents)

    def set_magnet_scale_factors(self, coil_currents, mice_modules = None):
        """
        Set magnet scale factors

        - coil_currents: python dictionary that maps the MiceModule name of the
                         coil to the scale factor that should be applied
        """
        root_module = False
        if mice_modules == None:
            root_module = True
            if self.first:
                for key, value in coil_currents.iteritems():
                    print key, value
            self.coil_currents = coil_currents
            # first get the mice modules
            mice_modules = maus_cpp.mice_module.MiceModule(self.lattice_file)
            # get the children of the mice_modules
            # note this makes a copy (list) of the existing children on mice_modules
        children = mice_modules.get_children()
        # find the coil module and reset the scale factor property
        # note we enforce polarity; so if the scale factor was negative it
        # should stay negative
        for child in children:
            for coil_name in coil_currents.keys():
                new_scale_factor = coil_currents[coil_name]
                # force Read mode after first iteration (optimisation)
                if not self.first:
                    child.set_property("FieldMapMode", "string", "Read")
                if child.get_name().find(coil_name) > -1:
                    polarity = child.get_property("ScaleFactor", "double")
                    polarity /= abs(polarity)+1e-9
                    coil_scale = abs(new_scale_factor)*polarity
                    child.set_property("ScaleFactor", "double", coil_scale)
                    #print "Setting", coil_name, "ScaleFactor to", coil_scale
            # recurse
            self.set_magnet_scale_factors(coil_currents, child)
        # set the mice_modules children to the new values
        mice_modules.set_children(children)
        # update the monte carlo with the new mice modules
        if root_module:
            maus_cpp.globals.set_monte_carlo_mice_modules(mice_modules)
        # reinitialise the optics model (transfer matrices)
        self.first = False

    def calculate_transfer_matrices(self):
        self.tm_list = []
        self.z_list = []
        tracks = []
        dim = 4
        deltas = [1.]*dim
        tm_dict_of_lists = {} # map station to hits on that station
        z_dict = {}
        for track_number, index in enumerate(self.get_index_2(dim)):
            values = [deltas[i]*(index[i]-1) for i in range(dim)]
            energy = (values[1]**2+values[3]**2+self.mom**2+self.mass**2)**0.5
            tracks.append({'primary':
                    {
                        'position':{'x':values[0], 'y':values[2], 'z':self.z_in},
                        'momentum':{'x':values[1], 'y':values[3], 'z':self.mom},
                        'particle_id':-13,
                        'energy':energy,
                        'time':0.,
                        'random_seed':0
                    }
                }
            )
        tracks = maus_cpp.simulation.track_particles(json.dumps(tracks))
        tracks = json.loads(tracks)
        for a_track in tracks:
            for vhit in a_track['virtual_hits']:
                station = vhit['station_id']
                if station not in tm_dict_of_lists:
                    tm_dict_of_lists[station] = []
                tm_dict_of_lists[station].append([
                        vhit['position']['x'], vhit['momentum']['x'],
                        vhit['position']['y'], vhit['momentum']['y'],
                ])
                z_dict[station] = vhit['position']['z']
        station_list = sorted(tm_dict_of_lists.keys())
        hits0 = tm_dict_of_lists[station_list[0]]
        for a_station in station_list:
            hits = tm_dict_of_lists[a_station]
            if len(hits) != len(tracks):
                continue
            tm = maus_cpp.polynomial_map.PolynomialMap.least_squares_fit(hits0, hits, 1)
            self.tm_list.append(tm)
            self.z_list.append(z_dict[a_station])

    def get_index_2(self, dim):
        index = [1]*dim
        yield index
        for i in range(4):
            index = [1]*dim
            index[i] = 0
            yield index
        for i in range(4):
            index = [1]*dim
            index[i] = 2
            yield index

    def plot_file(self, name):
        """
        Get the plot file name
        - name: describes what was plotted
        - image_format: format specifier e.g. png, eps, root
        Returns a string with the file name 
        """
        return name+'.png'
#+'_p='+str(round(self.mom, 2))+\
        '_FC='+str(round(self.coil_currents["FocusCoil"], 2))+\
        '_M1U='+str(round(self.coil_currents["MatchCoil1U"], 2))+\
        '_M2U='+str(round(self.coil_currents["MatchCoil2U"], 2))+\
        '_M1D='+str(round(self.coil_currents["MatchCoil1D"], 2))+\
        '_M2D='+str(round(self.coil_currents["MatchCoil2D"], 2))+'.png'

    @staticmethod
    def plot_field(canvas=None, line_color=1, graph_name='bz on axis'):
        """
        Plot the field on-axis
        - canvas: if None, generate a new canvas. Else plot on the specified
                  canvas
        - line_color: specify a line_color for the graph
        - graph_name: specify a graph name - will show up in e.g. canvas legend
        """
        z_list = range(-5000, 5000, 50)
        bz_list = []
        for z_coord in z_list:
            bz_field = maus_cpp.field.get_field_value(0, 0, z_coord, 0)[2]*1e3
            bz_list.append(bz_field)
        hist, graph = xboa.Common.make_root_graph(graph_name,
                                                  z_list, "z [mm]",
                                                  bz_list, "B_{z} [T]")
        if canvas == None:
            canvas = xboa.Common.make_root_canvas("bz")
            hist.Draw()
        canvas.cd()
        graph.Draw()
        graph.SetLineColor(line_color)
        canvas.Update()
        return canvas, graph

    def get_beta(self, covariance_matrix):
        """
        Return the transverse beta function, given by (Var(x)+Var(y))*p/(2m)
        - covariance_matrix: use this matrix to calculate the 

        Note that for all calculations transverse emittance is 1 
        """
        return (covariance_matrix[0][0]+\
                covariance_matrix[2][2])*self.mom/2./self.mass

    def get_alpha(self, covariance_matrix):
        """
        Return the transverse alpha function
      
        Alpha is given by
            (Cov(x, px)+Cov(y, py))*p/(2m)

        Note that for all calculations transverse emittance is 1 
        """
        return -(covariance_matrix[0][1]+\
                 covariance_matrix[2][3])/2./self.mass

    def get_phi(self, z_position):
        delta_z_list = [abs(z - z_position) for z in self.z_list]
        closest = delta_z_list.index(min(delta_z_list))
        tm = self.tm_list[closest].get_coefficients_as_matrix()
        tm_2d = [[tm[i][j] for i in range(0, 2)] for j in range(1, 3)]
        cosphi = (tm_2d[0][0]+tm_2d[1][1])/2.
        tm_2d[0][0] -= cosphi
        tm_2d[1][1] -= cosphi
        det_2d = numpy.linalg.det(numpy.array(tm_2d))
        if det_2d < 0:
            return -5.
        sinphi = det_2d**0.5
        phi =  math.atan2(sinphi, cosphi)
        return phi

    def plot_beta(self, canvas=None, line_color=1, graph_name='beta on axis'):
        """
        Plot the beta function as a function of z-position
        - canvas: if None, generate a new canvas. Else plot on the specified
                  canvas
        - line_color: specify a line_color for the graph
        - graph_name: specify a graph name - will show up in e.g. canvas legend
        """
        self.calculate_transfer_matrices()
        z_list = range(int(self.z_in), int(abs(self.z_in)), 100)
        ellipse_list = [self.transport_covariance_matrix(z) for z in z_list]
        beta_list = [self.get_beta(ellipse) for ellipse in ellipse_list]
        hist, graph = xboa.Common.make_root_graph(graph_name,
                                              z_list, "z [mm]",
                                              beta_list, "#beta_{#perp}   [mm]",
                                              ymin=0.)
        hist.SetTitle("")
        if canvas == None:
            canvas = xboa.Common.make_root_canvas("beta")
            hist.Draw()
        canvas.cd()
        hist.GetYaxis().SetTitleOffset(1.5)
        graph.Draw("SAME")
        graph.SetLineColor(line_color)
        canvas.Update()
        return canvas, graph

    def plot_rms(self, canvas=None, line_color=1, graph_name='beta on axis'):
        """
        Plot the beta function as a function of z-position
        - canvas: if None, generate a new canvas. Else plot on the specified
                  canvas
        - line_color: specify a line_color for the graph
        - graph_name: specify a graph name - will show up in e.g. canvas legend
        """
        self.calculate_transfer_matrices()
        z_list = range(int(self.z_in), int(abs(self.z_in)), 100)
        ellipse_list = [self.transport_covariance_matrix(z) for z in z_list]
        beta_list = [self.get_beta(ellipse) for ellipse in ellipse_list]
        hist, graph = xboa.Common.make_root_graph(graph_name,
                                              z_list, "z [mm]",
                                              beta_list, "#beta_{#perp}   [mm]",
                                              ymin=0.)
        hist.SetTitle("")
        if canvas == None:
            canvas = xboa.Common.make_root_canvas("beta")
            hist.Draw()
        canvas.cd()
        hist.GetYaxis().SetTitleOffset(1.5)
        graph.Draw("SAME")
        graph.SetLineColor(line_color)
        canvas.Update()
        return canvas, graph


    def plot_beam_ellipse(self, z_target, axis, canvas=None, line_color=1):
        """
        Plot the beta function as a function of z-position
        - canvas: if None, generate a new canvas. Else plot on the specified
                  canvas
        - line_color: specify a line_color for the graph
        - graph_name: specify a graph name - will show up in e.g. canvas legend
        """
        self.calculate_transfer_matrices()
        ellipse = self.transport_covariance_matrix(z_target)
        tfunc = xboa.common.make_root_ellipse_function([0., 0.], ellipse, [15], xmin=-250, xmax=250, ymin=-100, ymax=100)
        ellipse = [[ellipse[i][j] for j in range(2*axis, 2*axis+1)] for i in range(2*axis, 2*axis+1)]
        print ellipse
        tfunc.SetLineColor(line_color)
        if canvas == None:
            canvas = xboa.Common.make_root_canvas("ellipse")
            tfunc.Draw()
        else:
            canvas.cd()
            tfunc.Draw("SAME")
        canvas.Update()
        return canvas, tfunc

    def plot_phi(self, canvas=None, line_color=1, graph_name='phase advance'):
        """
        Plot the beta function as a function of z-position
        - canvas: if None, generate a new canvas. Else plot on the specified
                  canvas
        - line_color: specify a line_color for the graph
        - graph_name: specify a graph name - will show up in e.g. canvas legend
        """
        self.calculate_transfer_matrices()
        phi_list = [self.get_phi(z) for z in self.z_list]
        hist, graph = xboa.Common.make_root_graph(graph_name,
                                              self.z_list, "z [mm]",
                                              phi_list, "#phi [rad]")
        if canvas == None:
            canvas = xboa.Common.make_root_canvas("phi")
            hist.Draw()
        canvas.cd()
        hist.GetYaxis().SetTitleOffset(1.5)
        graph.Draw()
        graph.SetLineColor(line_color)
        canvas.Update()
        return canvas, graph

    def transport_measured_covariance_matrix(self, z_position, measured_covariance_matrix):
        numpy_ellipse_in = numpy.array(measured_covariance_matrix)
        delta_z_list = [abs(z - z_position) for z in self.z_list]
        closest = delta_z_list.index(min(delta_z_list))
        tm = self.tm_list[closest]
        numpy_tm = numpy.array(tm.get_coefficients_as_matrix())[0:4, 1:5]
        #numpy_ellipse_in = numpy.array([[ellipse_start.get_element(i, j) for i in range(3, 7)] for j in range(3, 7)])
        # V_out = M^T V_in M
        numpy_ellipse_out = numpy.dot(numpy_ellipse_in, numpy.transpose(numpy_tm))
        numpy_ellipse_out = numpy.dot(numpy_tm, numpy_ellipse_out)
        return numpy_ellipse_out


    def transport_covariance_matrix(self, z_position):
        """
        Transport the covariance matrix to nearest virtual plane to z_position
        """
        initial_bz = maus_cpp.field.get_field_value(0, 0, self.z_in, 0)[2]
        # beta function for constant bz solenoid and dbeta/dz = 0.
        if self.initial_beta == None:
            self.initial_beta = self.mom/initial_bz*2./xboa.Common.constants['c_light']
        delta_z_list = [abs(z - z_position) for z in self.z_list]
        closest = delta_z_list.index(min(delta_z_list))
        tm = self.tm_list[closest]
        numpy_tm = numpy.array(tm.get_coefficients_as_matrix())[0:4, 1:5]
        ellipse_start = maus_cpp.covariance_matrix.create_from_penn_parameters(
             mass = self.mass,
             momentum = self.mom,
             emittance_t = 1.,
             beta_t = self.initial_beta,
             alpha_t = self.initial_alpha,
             emittance_l = 1.,
             beta_l = 1.,
             bz = initial_bz,
        )
        numpy_ellipse_in = numpy.array([[ellipse_start.get_element(i, j) for i in range(3, 7)] for j in range(3, 7)])
        # V_out = M^T V_in M
        numpy_ellipse_out = numpy.dot(numpy_ellipse_in, numpy.transpose(numpy_tm))
        numpy_ellipse_out = numpy.dot(numpy_tm, numpy_ellipse_out)
        return numpy_ellipse_out.tolist()

    def print_tm(self):
        for z, tm in zip(self.z_list, self.tm_list):
            print "TM at z", z
            tm_el = tm.get_coefficients_as_matrix()
            for row in tm_el:
                for el in row:
                    print str(round(el, 3)).rjust(15),
                print

    def get_tm(self, z_pos):
        delta_z_list = [abs(z-z_pos) for z in self.z_list]
        z_tm_list = zip(delta_z_list, self.tm_list)
        delta_z, tm = min(z_tm_list)
        return tm

def tm_hack(lattice):
    import decoupled_transfer_matrix
    #lattice.print_tm()
    tm_00 = lattice.get_tm(-1900.)
    tm_01 = lattice.get_tm(0.)
    tm_02 = lattice.get_tm(1900.)
    lattice.z_in = 0.
    lattice.calculate_transfer_matrices()
    tm_11 = lattice.get_tm(0.)
    tm_12 = lattice.get_tm(1900.)

    tm_00 = numpy.array(tm_00.get_coefficients_as_matrix())[0:4, 1:5]
    tm_11 = numpy.array(tm_11.get_coefficients_as_matrix())[0:4, 1:5]

    tm_01 = numpy.array(tm_01.get_coefficients_as_matrix())[0:4, 1:5]
    tm_02 = numpy.array(tm_02.get_coefficients_as_matrix())[0:4, 1:5]
    tm_12 = numpy.array(tm_12.get_coefficients_as_matrix())[0:4, 1:5]

    print "tm_{00}"
    print tm_00
    print "tm_{11}"
    print tm_11
    print
    print "tm_{01}"
    print tm_01
    print "tm_{02}"
    print tm_02
    print "tm_{12}"
    print tm_12
    print "tm_{12} tm_{01}"
    print numpy.dot(tm_12, tm_01)
    print "tm_{12}^{-1}"
    print numpy.linalg.inv(tm_12)

    numpy.set_printoptions(linewidth=200)
    decoupled_transfer_matrix.DecoupledTransferMatrix.det_tolerance = 0.01
    tm_01_decoupled = decoupled_transfer_matrix.DecoupledTransferMatrix(tm_01)
    tm_12_decoupled = decoupled_transfer_matrix.DecoupledTransferMatrix(tm_12)

    print "Decoupled tm_01"
    print numpy.real(tm_01_decoupled.t)
    print "Decoupled tm_12"
    print numpy.real(tm_12_decoupled.t)
    print "Decoupled tm_12^{-1}"
    print numpy.real(numpy.linalg.inv(tm_12_decoupled.t))


def main():
    """
    Plot beta and bz before and after matching the beamline
    """
    lattice_file = os.path.expandvars(
                "Step_IV_Coils.dat"
              )
    lattice_optimiser = LatticeOptimiser(200., 50., 0., 0., lattice_file)
    lattice = lattice_optimiser.lattice
    lattice.initial_beta = 400.
    lattice.initial_alpha = +0.5


    lattice_optimiser.setup_minuit(100.)
    canvas, tfunc = lattice.plot_beam_ellipse(-2900., 0, line_color=4)
    tfunc 
    canvas, tfunc = lattice.plot_beam_ellipse(+2900., 0, line_color=8, canvas = canvas)

    #lattice_optimiser.optimise(100)
    #bz_canvas, bz_graph = lattice.plot_field(
    #                                      line_color=4,
    #                                      graph_name='')
    #bz_canvas.Print("bz.png")


    tm_hack(lattice)
    beta_canvas = None
    graph_list = []
#    for beta, color in [(None, 1), (500., 2), (2000., 8)]:
#        lattice.initial_beta = beta
    beta_canvas, beta_graph = lattice.plot_beta(canvas=beta_canvas,
                      line_color=1) # also calculates initial beta
    beta_graph.SetName('#beta = '+str(lattice.initial_beta)+' mm')
    graph_list.append(beta_graph)
#    make_root_legend(beta_canvas, graph_list, 0.12, 0.5, 0.32, 0.895)
    beta_canvas.Print("beta_vs_z_beta.png")
    raw_input()

    lattice.initial_beta = None
    mom_str = str(round(lattice.mom, 2))
    beta_canvas, beta_graph = lattice.plot_beta(
                                line_color=1,
                                graph_name='pz = '+mom_str+' MeV/c')
    graph_list = [beta_graph]
    for pz, color in [(140., 2), (240., 6), (300., 8)]:
        lattice.mom = pz
        beta_canvas, beta_graph = lattice.plot_beta(canvas=beta_canvas,
                          line_color=color,
                          graph_name='pz = '+str(pz)+' MeV/c')
        graph_list.append(beta_graph)
    make_root_legend(beta_canvas, graph_list, 0.12, 0.5, 0.32, 0.895)
    lattice.mom = 200.
    beta_canvas.Print("beta_vs_z_pz.png")


    lattice.calculate_transfer_matrices()
    ellipse_start = lattice.transport_covariance_matrix(-2600.)
    beta_start = lattice.get_beta(ellipse_start)
    ellipse_end = lattice.transport_covariance_matrix(2600.)
    beta_end = lattice_optimiser.lattice.get_beta(ellipse_end)
    phi_start = lattice_optimiser.lattice.get_phi(-2600.)
    phi_end = lattice_optimiser.lattice.get_phi(2600.)
    print maus_cpp.field.str(True)
    print "beta at start is ", beta_start, "mm"
    print "beta at end is  ", beta_end, "mm"
    print "phi at start is ", phi_start, "rad"
    print "phi at end is  ", phi_end, "rad"

    print "Finished - press <Enter> to end"
    raw_input()

# NOTES
# Check the coil geometry
# Use the Parzen formula to do periodic lattice function
# We probably need to run in a sort of doublet mode so we need to get 2pi phase advance out of the magnets...

if __name__ == "__main__":
    main()

