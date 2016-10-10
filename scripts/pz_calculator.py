import math
import xboa.common

class PzCalculator(object):
    def __init__(self, config):
        self.config = config
        self.c_light = xboa.common.constants["c_light"]
        self.mass = xboa.common.pdg_pid_to_mass[13]

    def get_separation_z(self, event, z):
        tku = event["data_in"]
        tkd = event["data_out"]
        rx = (tku[0] + tku[1]*(z-self.config.z_tku)) - \
             (tkd[0] + tkd[1]*(z-self.config.z_tkd))
        ry = (tku[2] + tku[3]*(z-self.config.z_tku)) - \
             (tkd[2] + tkd[3]*(z-self.config.z_tkd))
        r = (rx**2 + ry**2)**0.5
        return r

    def get_path_length(self, event):
        tku = event["data_in"]
        tkd = event["data_out"]
        ptu = (tku[1]**2 + tku[3]**2)**0.5
        zu = self.config.z_fc-self.config.z_tof1
        path_length_u = zu/math.cos(math.atan(ptu))

        ptd = (tkd[1]**2 + tkd[3]**2)**0.5
        zd = self.config.z_tof2 - self.config.z_fc
        path_length_d = zd/math.cos(math.atan(ptd))

        path_length = path_length_u + path_length_d
        return path_length

    def get_p_tot(self, event):
        if event["data_in"] == None or event["data_out"] == None or event["tof12"] == None:
            return None, None, None
        path_length = self.get_path_length(event)
        beta = path_length/(event["tof12"])/self.c_light
        if beta >= 1 or beta < 0:
            return None, None, None
        gamma = 1/(1-beta**2)**0.5
        p_tot = self.mass*beta*gamma
        miss_distance = self.get_separation_z(event, self.config.z_fc)
        return p_tot, path_length, miss_distance

