import numpy as np
from astropy.cosmology import FlatLambdaCDM
from astropy import constants as const
from astropy import units as u
import pandas as pd


class Component():
    def __init__(self, x, y, maj, min, pos, flux, date, mjd, year, delta_x_est=0, delta_y_est=0,
                 component_number=-1, is_core=False, redshift=0, scale=60 * 60 * 10 ** 3):
        self.x = x
        self.y = y
        self.mjd = mjd
        self.maj = maj
        self.min = min
        self.pos = pos
        self.flux = flux
        self.date = date
        self.mjd = mjd
        self.year = year
        self.component_number = component_number
        self.is_core = is_core
        self.delta_x_est = self.x
        self.delta_y_est = self.y
        self.distance_to_core = np.sqrt(self.delta_x_est ** 2 + self.delta_y_est ** 2)
        self.redshift = redshift
        self.tb = 5.44e9 * self.flux * (1 + self.redshift) / self.maj / self.min
        self.scale = scale

    def set_distance_to_core(self, core_x, core_y):
        self.delta_x_est = self.x - core_x
        self.delta_y_est = self.y - core_y
        self.distance_to_core = np.sqrt(self.delta_x_est ** 2 + self.delta_y_est ** 2)

    def assign_component_number(self, number):
        self.component_number = number


class ComponentCollection():
    def __init__(self, components=[], name=""):
        self.components = components
        if len(components) > 0:
            self.redshift = components[0].redshift
            self.scale = components[0].scale
        else:
            self.redshift = 0
            self.scale = 1
        self.name = name
        self.mjds = []
        self.year = []
        self.dist = []
        self.dist_err = []
        self.time = []
        self.xs = []
        self.ys = []

        for comp in components:
            self.year.append(comp.year)
            self.dist.append(comp.distance_to_core * self.scale)
            self.dist_err.append(comp.maj / 2 * self.scale)
            self.xs.append(comp.delta_x_est)
            self.ys.append(comp.delta_y_est)

    def length(self):
        return len(self.components)

    def get_speed(self):
        if self.length() > 2:
            cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

            def reduced_chi2(fit, x, y, yerr, N, n):
                return 1. / (N - n) * np.sum(((y - fit) / yerr) ** 2.)

            t_mid = (np.min(self.year) + np.max(self.year)) / 2.
            time = np.array(self.year) - t_mid

            linear_fit, cov_matrix = np.polyfit(time, self.dist, 1, cov='scaled')  # ,w=1./dist_err)

            speed = linear_fit[0]
            speed_err = np.sqrt(cov_matrix[0, 0])
            y0 = linear_fit[1] - t_mid*speed
            y0_err = np.sqrt(cov_matrix[1, 1])
            beta_app = speed * (np.pi / (180 * self.scale * u.yr)) * (
                    cosmo.luminosity_distance(self.redshift) / (const.c.to('pc/yr') * (1 + self.redshift)))
            beta_app_err = speed_err * (np.pi / (180 * self.scale * u.yr)) * (
                    cosmo.luminosity_distance(self.redshift) / (const.c.to('pc/yr') * (1 + self.redshift)))
            d_crit = np.sqrt(1 + beta_app ** 2)
            d_crit_err = (1 + beta_app) ** (-0.5) * beta_app * beta_app_err
            dist_0_est = linear_fit[1] - speed * t_mid
            t_0 = - linear_fit[1] / speed + t_mid
            sum_x = time / np.array(self.dist_err) ** 2
            sum_x2 = time ** 2 / np.array(self.dist_err) ** 2
            sum_err = 1. / np.array(self.dist_err) ** 2
            Delta = np.sum(sum_err) * np.sum(sum_x2) - (np.sum(sum_x)) ** 2
            t_0_err = np.sqrt((cov_matrix[1, 1] / speed ** 2) + (linear_fit[1] ** 2 * cov_matrix[0, 0] / speed ** 4) +
                              2 * linear_fit[1] / speed ** 3 * np.sum(sum_x) / Delta)
            red_chi_sqr = reduced_chi2(linear_fit[0] * time + linear_fit[1], time, self.dist, self.dist_err, len(time),
                                       len(linear_fit))

        else:
            speed = 0
            speed_err = 0
            y0 = 0
            y0_err = 0
            beta_app = 0
            beta_app_err = 0
            d_crit = 0
            d_crit_err = 0
            dist_0_est = 0
            t_0 = 0
            t_0_err = 0
            red_chi_sqr = 0

        return {"speed": speed, "speed_err": speed_err, "y0": y0, "y0_err": y0_err,
                "beta_app": beta_app, "beta_app_err": beta_app_err, "d_crit": d_crit, "d_Crit_err": d_crit_err,
                "dist_0_est": dist_0_est, "t_0": t_0, "t_0_err": t_0_err, "red_chi_sqr": red_chi_sqr}
