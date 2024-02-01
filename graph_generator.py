import matplotlib as mpl
import matplotlib.pyplot as plt
from kivy.metrics import dp
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.colors as colors
from astropy.io import fits
from astropy.modeling import models, fitting
import os
from mpl_toolkits.axes_grid1 import make_axes_locatable
from kinematics import Component
from astropy.time import Time

#optimized draw on Agg backend
mpl.rcParams['path.simplify'] = True
mpl.rcParams['path.simplify_threshold'] = 1.0
mpl.rcParams['agg.path.chunksize'] = 1000

#define some matplotlib figure parameters
mpl.rcParams['font.family'] = 'Quicksand'
#mpl.rcParams['axes.spines.top'] = False
#mpl.rcParams['axes.spines.right'] = False
mpl.rcParams['axes.linewidth'] = 1.0

font_size_axis_title=dp(13)
font_size_axis_tick=dp(12)

class KinematicPlot(object):
    def __init__(self):

        super().__init__()
        self.fig, self.ax = plt.subplots(1, 1)
        self.fig.subplots_adjust(left=0.13,top=0.96,right=0.93,bottom=0.2)


    def plot_kinematics(self,component_collection,color):
        if component_collection.length()>0:
            self.ax.scatter(component_collection.year,component_collection.dist,c=color,marker=".")
        self.ax.set_xlabel('Time [year]', fontsize=font_size_axis_title)
        self.ax.set_ylabel('Distance from Core [mas]', fontsize=font_size_axis_title)

    def plot_fluxs(self,component_collection,color):
        if component_collection.length() > 0:
            self.ax.plot(component_collection.year, component_collection.fluxs, c=color, label=component_collection.name,marker=".")
        self.ax.set_xlabel('Time [year]', fontsize=font_size_axis_title)
        self.ax.set_ylabel('Flux Density [Jy]', fontsize=font_size_axis_title)

    def plot_tbs(self,component_collection,color):
        if component_collection.length() > 0:
            self.ax.plot(component_collection.year, component_collection.tbs, c=color, label=component_collection.name,marker=".")
        self.ax.set_xlabel('Time [year]', fontsize=font_size_axis_title)
        self.ax.set_ylabel('Brightness Temperature [K]', fontsize=font_size_axis_title)
        self.ax.set_yscale("log")



    def set_limits(self,x,y):
        self.ax.set_xlim(x)
        self.ax.set_ylim(y)
    def plot_linear_fit(self,x_min,x_max,slope,y0,color,label=""):
        def y(x):
            return slope*x+y0
        self.ax.plot([x_min,x_max],[y(x_min),y(x_max)],color,label=label)




class FitsImage(object):
    """class that generate Matplotlib graph."""

    def __init__(self,clean_image_file,model_image_file=""):
        """Create empty structure plot. 
        
        """       
        super().__init__()

        self.clean_image_file=clean_image_file
        self.model_image_file=model_image_file
        self.fig, self.ax = plt.subplots(1, 1)

        self.components=[]

        #component default color
        self.component_color = "black"

        # Define some variables (create gui for it later)
        fit_noise = True  # if True, the noise value and rms deviation will be fitted as described in the PhD-thesis of Moritz Böck (https://www.physik.uni-wuerzburg.de/fileadmin/11030400/Dissertation_Boeck.pdf); if False, the noise frome difmap will be used

        # Contour plot
        contour = True  # if True, a contour plot will be done
        contour_color = 'grey'  # input: array of color-strings; if None, the contour-colormap (contour_cmap) will be used
        contour_cmap = None  # matplotlib colormap string
        contour_alpha = 1  # transparency
        contour_width = 0.5  # contour linewidth



        # Image colormap
        self.im_colormap = False  # if True, a image colormap will be done
        im_color = 'inferno'  # string for matplotlib colormap

        # Plot mode
        individual_plots = True  # if True, all epochs will be plotted to individual pdf files
        combined_plot = False  # if True, all epochs will be fitted to one combined pdf file

        # Overplot Gaussian-components
        if model_image_file=="":
            overplot_gauss = False  # if True all Gaussian components are plotted
        else:
            overplot_gauss = True

        gauss_linewidth = 0.5  # linewidth of the ellipse
        gauss_color = 'black'  # color of the ellipse

        # Overplot clean-components
        overplot_clean = False  # if True all clean components are plotted
        clean_alpha = 1  # float for sympol transparency
        clean_linewidth = 0.5  # clean linewidth of the symbol

        # Read clean files in
        hdu_list = fits.open(clean_image_file)

        # Set name
        self.name = hdu_list[0].header["OBJECT"]

        self.freq = float(hdu_list[0].header["CRVAL3"]) #frequency in Hertz


        # Unit selection and adjustment
        degpp=abs(hdu_list[0].header["CDELT1"])#degree per pixel

        if degpp>0.01:
            unit = 'deg'
            scale = 1.
        elif degpp>6.94e-6:
            unit = 'arcmin'
            scale = 60.
        elif degpp>1.157e-7:
            scale = 60. * 60.
            unit = 'arcsec'
        else:
            scale = 60. * 60. * 1000.
            unit = 'mas'

        # Convert Pixel into unit
        X = np.linspace(0, hdu_list[0].header["NAXIS1"], hdu_list[0].header["NAXIS1"],
                        endpoint=False)  # NAXIS1: number of pixels at R.A.-axis
        for j in range(len(X)):
            X[j] = (X[j] - hdu_list[0].header["CRPIX1"]) * hdu_list[0].header[
                "CDELT1"] * scale  # CRPIX1: reference pixel, CDELT1: deg/pixel
        X[int(hdu_list[0].header["CRPIX1"])] = 0.0

        Y = np.linspace(0, hdu_list[0].header["NAXIS2"], hdu_list[0].header["NAXIS2"],
                        endpoint=False)  # NAXIS2: number of pixels at Dec.-axis
        for j in range(len(Y)):
            Y[j] = (Y[j] - hdu_list[0].header["CRPIX2"]) * hdu_list[0].header[
                "CDELT2"] * scale  # CRPIX2: reference pixel, CDELT2: deg/pixel
        Y[int(hdu_list[0].header["CRPIX2"])] = 0.0

        extent = np.max(X), np.min(X), np.min(Y), np.max(Y)

        # Get image data
        image_data = hdu_list[0].data
        Z = image_data[0, 0, :, :]

        # Fit noise level (according to https://www.physik.uni-wuerzburg.de/fileadmin/11030400/Dissertation_Boeck.pdf)
        if fit_noise == True:
            Z1 = Z.flatten()
            bin_heights, bin_borders = np.histogram(Z1 - np.min(Z1) + 10 ** (-5), bins="auto")
            bin_widths = np.diff(bin_borders)
            bin_centers = bin_borders[:-1] + bin_widths / 2.
            bin_heights_err = np.where(bin_heights != 0, np.sqrt(bin_heights), 1)

            t_init = models.Gaussian1D(np.max(bin_heights), np.median(Z1 - np.min(Z1) + 10 ** (-5)), 0.001)
            fit_t = fitting.LevMarLSQFitter()
            t = fit_t(t_init, bin_centers, bin_heights, weights=1. / bin_heights_err)

            # Set contourlevels to mean value + 3 * rms_noise * 2 ** x
            levs1 = t.mean.value + np.min(Z1) - 10 ** (-5) + 3 * t.stddev.value * np.logspace(0, 100, 100,
                                                                                              endpoint=False, base=2)
            levs = t.mean.value + np.min(Z1) - 10 ** (-5) - 3 * t.stddev.value * np.logspace(0, 100, 100,
                                                                                             endpoint=False, base=2)
            levs = np.flip(levs)
            levs = np.concatenate((levs, levs1))

        else:
            levs1 = 3 * hdu_list[0].header["NOISE"] * np.logspace(0, 100, 100, endpoint=False, base=2)
            levs = np.flip(-levs1)
            levs = np.concatenate((levs, levs1))



        # Image colormap
        if self.im_colormap == True:
            col = self.ax.imshow(Z, cmap=im_color, norm=colors.SymLogNorm(linthresh=levs1[0], linscale=0.5, vmin=levs[99],
                                                                     vmax=0.5 * np.max(Z), base=10.), extent=extent,
                            origin='lower')
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="5%", pad=0.05)
            cbar = fig.colorbar(col, use_gridspec=True, cax=cax)
            cbar.set_label('Flux Density [Jy]')

        # Contour plot
        if contour == True:
            if self.im_colormap == True:
                self.ax.contour(X, Y, Z, linewidths=contour_width, levels=levs, colors='grey', alpha=contour_alpha,
                           cmap=contour_cmap)
            else:
                self.ax.contour(X, Y, Z, linewidths=contour_width, levels=levs, colors=contour_color, alpha=contour_alpha,
                           cmap=contour_cmap)

        # Set beam parameters
        beam_maj = hdu_list[0].header["BMAJ"] * scale
        beam_min = hdu_list[0].header["BMIN"] * scale
        beam_pa = hdu_list[0].header["BPA"]

        #plot limits
        ra_min = np.min(X)
        ra_max = np.max(X)
        dec_min = np.min(Y)
        dec_max = np.max(Y)

        # Set beam ellipse, sourcename and observation date positions
        size_x = np.absolute(ra_max) + np.absolute(ra_min)
        size_y = np.absolute(dec_max) + np.absolute(dec_min)
        if size_x > size_y:
            ell_x = ra_max - beam_maj
            ell_y = dec_min + beam_maj
            name_x = ra_max - size_x * 0.05
            name_y = dec_max - size_x * 0.05
            date_x = ra_min + size_x * 0.05
            date_y = name_y
        else:
            ell_x = ra_max - beam_maj
            ell_y = dec_min + beam_maj
            name_x = ra_max - size_y * 0.05
            name_y = dec_max - size_y * 0.05
            date_x = ra_min + size_y * 0.05
            date_y = name_y

        # Plot beam
        beam = Ellipse([ell_x, ell_y], beam_maj, beam_min, -beam_pa + 90, fc='grey')
        self.ax.add_artist(beam)

        #plot_date
        date=self.get_date(hdu_list)

        self.ax.set_title(date, fontsize=font_size_axis_title)

        """
        if self.im_colormap == True:
            self.ax.text(date_x, date_y, date, color='grey', ha='right', va='top')
        else:
            self.ax.text(date_x, date_y, date, color='black', ha='right', va='top')
        """

        """
        # Plot name
        if self.im_colormap == True:
            self.ax.text(name_x, name_y, self.name, color='grey', ha='left', va='top')
        else:
            self.ax.text(name_x, name_y, self.name, color='black', ha='left', va='top')
        """


        # Read modelfit files in
        if (overplot_gauss == True) or (overplot_clean == True):
            model_df=self.getComponentInfo(model_image_file)

            #sort in gauss and clean components
            model_gauss_df = model_df[model_df["Major_axis"] > 0.].reset_index()
            model_clean_df = model_df[model_df["Major_axis"] == 0.].reset_index()

            # Overplot clean components
            if overplot_clean == True:
                c_x = model_clean_df["Delta_x"]
                c_y = model_clean_df["Delta_y"]
                c_flux = model_clean_df["Flux"]

                for j in range(len(c_x)):
                    if c_flux[j] < 0.:
                        self.ax.plot(c_x[j] * scale, c_y[j] * scale, marker='+', color='red', alpha=clean_alpha,
                                linewidth=0.2, zorder=2)
                    else:
                        self.ax.plot(c_x[j] * scale, c_y[j] * scale, marker='+', color='green', alpha=clean_alpha,
                                linewidth=0.2, zorder=2)

            # Overplot Gaussian components
            if overplot_gauss == True:

                g_x = model_gauss_df["Delta_x"]
                g_y = model_gauss_df["Delta_y"]
                g_maj = model_gauss_df["Major_axis"]
                g_min = model_gauss_df["Minor_axis"]
                g_pos = model_gauss_df["PA"]
                g_flux = model_gauss_df["Flux"]
                g_date = model_gauss_df["Date"]
                g_mjd = model_gauss_df["mjd"]
                g_year = model_gauss_df["Year"]

                for j in range(len(g_x)):
                    #plot component
                    component_plot = self.plotComponent(g_x[j],g_y[j],g_maj[j],g_min[j],g_pos[j],scale)
                    component=Component(g_x[j],g_y[j],g_maj[j],g_min[j],g_pos[j],g_flux[j],g_date[j],g_mjd[j],g_year[j],scale=scale,freq=self.freq)
                    self.components.append([component_plot,component])

        hdu_list.close()

        self.xmin,self.xmax = ra_min, ra_max
        self.ymin,self.ymax = dec_min, dec_max
        
        self.fig.subplots_adjust(left=0.13,top=0.96,right=0.93,bottom=0.2)

        # Plot look tuning
        self.ax.set_aspect('equal', adjustable='box', anchor='C')
        self.ax.set_xlim(ra_min, ra_max)
        self.ax.set_ylim(dec_min, dec_max)
        self.ax.invert_xaxis()
        self.ax.set_xlabel('Relative R.A. [' + unit + ']',fontsize=font_size_axis_title)
        self.ax.set_ylabel('Relative DEC. [' + unit + ']',fontsize=font_size_axis_title)
        self.fig.tight_layout()

    def plotComponent(self,x,y,maj,min,pos,scale):

        # Plotting ellipses
        comp = Ellipse([x * scale, y * scale], maj * scale, min * scale, -pos + 90,
                       fill=False, zorder=2, color=self.component_color, lw=0.5)
        ellipse=self.ax.add_artist(comp)

        # Plotting axes of the ellipses
        maj1_x = x - np.sin(-np.pi / 180 * pos) * maj * 0.5
        maj1_y = y + np.cos(-np.pi / 180 * pos) * maj * 0.5
        maj2_x = x + np.sin(-np.pi / 180 * pos) * maj * 0.5
        maj2_y = y - np.cos(-np.pi / 180 * pos) * maj * 0.5

        min1_x = x - np.sin(-np.pi / 180 * (pos + 90)) * min * 0.5
        min1_y = y + np.cos(-np.pi / 180 * (pos + 90)) * min * 0.5
        min2_x = x + np.sin(-np.pi / 180 * (pos + 90)) * min * 0.5
        min2_y = y - np.cos(-np.pi / 180 * (pos + 90)) * min * 0.5

        line1=self.ax.plot([maj1_x * scale, maj2_x * scale], [maj1_y * scale, maj2_y * scale], color=self.component_color, lw=0.5)
        line2=self.ax.plot([min1_x * scale, min2_x * scale], [min1_y * scale, min2_y * scale], color=self.component_color, lw=0.5)


        return [ellipse,line1,line2]
    def get_date(self,hdu_list):

        # Plot date
        time = hdu_list[0].header["DATE-OBS"]
        time = time.split("T")[0]
        time = time.split("/")
        if len(time) == 1:
            date = time[0]
        elif len(time) == 3:
            if len(time[0]) < 2:
                day = "0" + time[0]
            else:
                day = time[0]
            if len(time[1]) < 2:
                month = "0" + time[1]
            else:
                month = time[1]
            if len(time[2]) == 2:
                if 45 < int(time[2]) < 100:
                    year = "19" + time[2]
                elif int(time[2]) < 46:
                    year = "20" + time[2]
            elif len(time[2]) == 4:
                year = time[2]
            date = year + "-" + month + "-" + day
        return date

    #gets components from .fits file
    def getComponentInfo(self,filename):
        data_df = pd.DataFrame()
        hdu_list = fits.open(filename)
        comp_data = hdu_list[1].data
        comp_data1 = np.zeros((len(comp_data), len(comp_data[0])))
        date = np.array([])
        year = np.array([])
        mjd = np.array([])
        for j in range(len(comp_data)):
            comp_data1[j, :] = comp_data[j]
            date1=self.get_date(hdu_list)
            date = np.append(date, date1)
            t = Time(date1)
            year = np.append(year, t.jyear)
            mjd = np.append(mjd, t.mjd)
        comp_data1_df = pd.DataFrame(data=comp_data1,
                                     columns=["Flux", "Delta_x", "Delta_y", "Major_axis", "Minor_axis", "PA",
                                              "Typ_obj"])
        comp_data1_df["Date"] = date
        comp_data1_df["Year"] = year
        comp_data1_df["mjd"] = mjd
        comp_data1_df.sort_values(by=["Delta_x", "Delta_y"], ascending=False, inplace=True)
        if data_df.empty:
            data_df = comp_data1_df
        else:
            data_df = pd.concat([data_df, comp_data1_df], axis=0, ignore_index=True)
        return data_df

    def change_plot_lim(self,x_min,x_max,y_min,y_max):
        self.ax.set_xlim(x_min, x_max)
        self.ax.set_ylim(y_min, y_max)


                
