import matplotlib as mpl
import matplotlib.pyplot as plt
from kivy.metrics import dp
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from matplotlib.collections import LineCollection
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


class ImageData(object):
    def __init__(self,
                 fits_file,
                 model="",
                 lin_pol=[],
                 evpa=[],
                 pol_from_stokes=True,
                 stokes_q="",
                 stokes_u=""):

        self.file_path = fits_file

        # Read clean files in
        hdu_list=fits.open(fits_file)
        self.hdu_list = hdu_list

        #read stokes data from input files if defined
        if stokes_q != "":
            stokes_q = fits.open(fits_file)[0].data[0, 0, :, :]
        else:
            stokes_q=[]

        if stokes_u != "":
            stokes_u = fits.open(fits_file)[0].data[0, 0, :, :]
        else:
            stokes_u=[]

        # Set name
        self.name = hdu_list[0].header["OBJECT"]

        self.freq = float(hdu_list[0].header["CRVAL3"])  # frequency in Hertz

        # Unit selection and adjustment
        self.degpp = abs(hdu_list[0].header["CDELT1"])  # degree per pixel

        if self.degpp > 0.01:
            self.unit = 'deg'
            self.scale = 1.
        elif self.degpp > 6.94e-6:
            self.unit = 'arcmin'
            self.scale = 60.
        elif self.degpp > 1.157e-7:
            self.scale = 60. * 60.
            self.unit = 'arcsec'
        else:
            self.scale = 60. * 60. * 1000.
            self.unit = 'mas'

        # Convert Pixel into unit
        self.X = np.linspace(0, hdu_list[0].header["NAXIS1"], hdu_list[0].header["NAXIS1"],
                        endpoint=False)  # NAXIS1: number of pixels at R.A.-axis
        for j in range(len(self.X)):
            self.X[j] = (self.X[j] - hdu_list[0].header["CRPIX1"]) * hdu_list[0].header[
                "CDELT1"] * self.scale  # CRPIX1: reference pixel, CDELT1: deg/pixel
        self.X[int(hdu_list[0].header["CRPIX1"])] = 0.0

        self.Y = np.linspace(0, hdu_list[0].header["NAXIS2"], hdu_list[0].header["NAXIS2"],
                        endpoint=False)  # NAXIS2: number of pixels at Dec.-axis
        for j in range(len(self.Y)):
            self.Y[j] = (self.Y[j] - hdu_list[0].header["CRPIX2"]) * hdu_list[0].header[
                "CDELT2"] * self.scale  # CRPIX2: reference pixel, CDELT2: deg/pixel
        self.Y[int(hdu_list[0].header["CRPIX2"])] = 0.0

        self.extent = np.max(self.X), np.min(self.X), np.min(self.Y), np.max(self.Y)

        self.image_data = hdu_list[0].data
        self.Z = self.image_data[0, 0, :, :]

        #read in polarization input

        # check if FITS file contains more than just Stokes I
        only_stokes_i = False
        if hdu_list[0].data.shape[0] == 1:
            only_stokes_i = True
        if (np.shape(self.Z) == np.shape(stokes_q) and np.shape(self.Z) == np.shape(stokes_u) and
                        np.shape(stokes_q) == np.shape(stokes_u)):
            only_stokes_i = True #in this case override the polarization data with the data that was input to Q and U

        if only_stokes_i:
            pols=1
            # Check if linpol/evpa/stokes_i have same dimensions!
            dim_wrong = True
            if pol_from_stokes:
                if (np.shape(self.Z) == np.shape(stokes_q) and np.shape(self.Z) == np.shape(stokes_u) and
                        np.shape(stokes_q) == np.shape(stokes_u)):
                    dim_wrong = False
                    self.stokes_q=stokes_q
                    self.stokes_u=stokes_u
                else:
                    self.lin_pol = np.zeros(np.shape(self.Z))
                    self.evpa = np.zeros(np.shape(self.Z))
            else:
                if (np.shape(self.Z) == np.shape(lin_pol) and np.shape(self.Z) == np.shape(evpa) and
                        np.shape(lin_pol) == np.shape(evpa)):
                    dim_wrong = False
                    self.lin_pol=lin_pol
                    self.evpa=evpa
                else:
                    self.lin_pol=np.zeros(np.shape(self.Z))
                    self.evpa=np.zeros(np.shape(self.Z))
        else:
            pols=3
            dim_wrong=False
            self.stokes_q=hdu_list[0].data[1,0,:,:]
            self.stokes_u=hdu_list[0].data[2,0,:,:]

        if pol_from_stokes and not dim_wrong:
            self.lin_pol = np.sqrt(self.stokes_q ** 2 + self.stokes_u ** 2)
            self.evpa = 0.5 * np.arctan2(self.stokes_u, self.stokes_q)


        # Set beam parameters
        try:
            self.beam_maj = hdu_list[0].header["BMAJ"] * self.scale
            self.beam_min = hdu_list[0].header["BMIN"] * self.scale
            self.beam_pa = hdu_list[0].header["BPA"]
        except:
            self.beam_maj = 0
            self.beam_min = 0
            self.beam_pa = 0

        self.date = get_date(fits_file)

        if model!="":
            #TODO basic checks if file is valid
            self.model=getComponentInfo(model)
        else:
            self.model=None

        hdu_list.close()


class FitsImage(object):
    """class that generate Matplotlib graph."""

    def __init__(self,
                 clean_image_file, #path to a .fits file (cleaned/final image)
                 model_image_file="", #path to a model file (.fits (or .mod to be implemented))
                 stokes_i_sigma_cut=3, #sigma_cut for stokes_i_contours
                 plot_mode="stokes_i", #possible modes "stokes_i", "lin_pol", "frac_pol"
                 im_colormap=True, #Choose whether to do colormap or not
                 contour=True, #Choose whether to do contour plot or not
                 contour_color = 'grey',  # input: array of color-strings; if None, the contour-colormap (contour_cmap) will be used
                 contour_cmap = None,  # matplotlib colormap string
                 contour_alpha = 1,  # transparency
                 contour_width = 0.5,  # contour linewidth
                 im_color='inferno', # string for matplotlib colormap
                 plot_beam=True, #choose whether to plot beam or not
                 ###HERE STARTS POLARIZATION INPUT
                 lin_pol=[],  # 2d list/array of lin pol data
                 evpa=[],  # 2d list/array of evpa data
                 pol_from_stokes=True,
                 # choose whether to use lin_pol & evpa input (False) OR stokes_q and stokes_u input (TRUE)
                 stokes_q="",  # filepath to stokes q fits image
                 stokes_u="", # filepath to stokes u fits image
                 plot_evpa=True, #decide whether to plot EVPA or not
                 evpa_len=6,  # choose length of EVPA in pixels
                 lin_pol_sigma_cut=3,  # choose lowest sigma contour for Lin Pol plot
                 evpa_distance=5,  # choose distance of EVPA vectors to draw in pixels
                 rotate_evpa=0,  # rotate EVPAs by a given angle in degrees (North through East)
                 evpa_color="white", #set EVPA color for plot
                 rcparams={}  # option to modify matplotlib look
                 ):

        super().__init__()

        #read image
        self.clean_image_file=clean_image_file
        self.clean_image = ImageData(clean_image_file,
                                     model=model_image_file,
                                     lin_pol=lin_pol,
                                     evpa=evpa,
                                     pol_from_stokes=pol_from_stokes,
                                     stokes_q=stokes_q,
                                     stokes_u=stokes_u)

        #set parameters
        self.name = self.clean_image.name
        self.freq = self.clean_image.freq
        image_data = self.clean_image.image_data
        X = self.clean_image.X
        Y = self.clean_image.Y
        Z = self.clean_image.Z
        unit = self.clean_image.unit
        scale = self.clean_image.scale
        degpp = self.clean_image.degpp
        extent = self.clean_image.extent
        date=self.clean_image.date
        # Set beam parameters
        beam_maj = self.clean_image.beam_maj
        beam_min = self.clean_image.beam_min
        beam_pa = self.clean_image.beam_pa
        self.evpa_color=evpa_color

        #plot limits
        ra_max,ra_min,dec_min,dec_max=extent

        self.model_image_file=model_image_file
        self.fig, self.ax = plt.subplots(1, 1)

        self.components=[]

        #component default color
        self.component_color = "black"

        fit_noise = True  # if True, the noise value and rms deviation will be fitted as described in the PhD-thesis of Moritz Böck (https://www.physik.uni-wuerzburg.de/fileadmin/11030400/Dissertation_Boeck.pdf); if False, the noise frome difmap will be used

        # Image colormap
        self.im_colormap = im_colormap  # if True, a image colormap will be done

        # Overplot Gaussian-components
        if model_image_file=="":
            overplot_gauss = False  # if True all Gaussian components are plotted
        else:
            overplot_gauss = True

        # Overplot clean-components
        overplot_clean = False  # if True all clean components are plotted
        clean_alpha = 1  # float for sympol transparency

        #get sigma levs
        levs, levs1 = get_sigma_levs(Z,stokes_i_sigma_cut)

        # Image colormap
        if self.im_colormap == True and plot_mode=="stokes_i":
            self.plotColormap(Z,im_color,levs,levs1,extent)
            contour_color="white"


        if (plot_mode=="lin_pol" or plot_mode=="frac_pol") and np.sum(self.clean_image.lin_pol)!=0:

            levs_linpol, levs1_linpol = get_sigma_levs(self.clean_image.lin_pol, lin_pol_sigma_cut)

            if plot_mode=="lin_pol":
                self.plotColormap(self.clean_image.lin_pol,im_color,levs_linpol,levs1_linpol,extent,
                                  label="Linear Polarized Intensity [Jy/beam]")
            if plot_mode=="frac_pol":
                plot_lin_pol = np.array(self.clean_image.lin_pol)
                plot_frac_pol = plot_lin_pol / np.array(self.clean_image.Z)
                plot_frac_pol = np.ma.masked_where((plot_lin_pol < levs1_linpol[0]) | (self.clean_image.Z<levs1[0]),
                                                  plot_frac_pol)  # Check if this is actually the right thing to do

                self.plotColormap(plot_frac_pol,im_color,np.zeros(100),[0.01],extent,
                                  label="Fractional Linear Polarization")
                self.evpa_color="black"
                contour_color="grey"

            if plot_evpa:
                self.plotEvpa(self.clean_image.evpa, rotate_evpa, evpa_len, evpa_distance, levs1_linpol, levs1)

        # Contour plot
        if contour == True:
            self.ax.contour(X, Y, Z, linewidths=contour_width, levels=levs, colors=contour_color,
                            alpha=contour_alpha,
                            cmap=contour_cmap)

        # Set beam ellipse, sourcename and observation date positions
        size_x = np.absolute(ra_max) + np.absolute(ra_min)
        size_y = np.absolute(dec_max) + np.absolute(dec_min)
        if size_x > size_y:
            ell_x = ra_max - beam_maj
            ell_y = dec_min + beam_maj
        else:
            ell_x = ra_max - beam_maj
            ell_y = dec_min + beam_maj

        if plot_beam:
            # Plot beam
            beam = Ellipse([ell_x, ell_y], beam_maj, beam_min, -beam_pa + 90, fc='grey')
            self.ax.add_artist(beam)

        self.ax.set_title(date, fontsize=font_size_axis_title)

        # Read modelfit files in
        if (overplot_gauss == True) or (overplot_clean == True):
            model_df = getComponentInfo(model_image_file)

            # sort in gauss and clean components
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
                    # plot component
                    component_plot = self.plotComponent(g_x[j], g_y[j], g_maj[j], g_min[j], g_pos[j], scale)
                    component = Component(g_x[j], g_y[j], g_maj[j], g_min[j], g_pos[j], g_flux[j], g_date[j],
                                          g_mjd[j], g_year[j], scale=scale, freq=self.freq)
                    self.components.append([component_plot, component])

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

    def plotColormap(self,
                     Z, #2d data array to plot
                     im_color, #colormap to use
                     levs, #sigma levs output
                     levs1, #sigma levs output
                     extent, #plot lims x_min,x_max,y_min,y_max
                     label="Flux Density [Jy]" #label for colorbar
                     ):
        col = self.ax.imshow(Z, cmap=im_color, norm=colors.SymLogNorm(linthresh=levs1[0], linscale=0.5, vmin=levs[99],
                                                                      vmax=0.5 * np.max(Z), base=10.), extent=extent,
                             origin='lower')
        divider = make_axes_locatable(self.ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cbar = self.fig.colorbar(col, use_gridspec=True, cax=cax)
        cbar.set_label(label)

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


    def change_plot_lim(self,x_min,x_max,y_min,y_max):
        self.ax.set_xlim(x_min, x_max)
        self.ax.set_ylim(y_min, y_max)

    def plotEvpa(self,evpa,rotate_evpa,evpa_len,evpa_distance,levs1_linpol,levs1_i):

        evpa_len=evpa_len*self.clean_image.degpp*self.clean_image.scale

        stokes_i=self.clean_image.Z
        # plot EVPA
        evpa = evpa + rotate_evpa / 180 * np.pi

        # create mask where to plot EVPA (only where stokes i and lin pol have plotted contours)
        mask = np.zeros(np.shape(stokes_i), dtype=bool)
        mask[:] = (self.clean_image.lin_pol > levs1_linpol[0]) * (stokes_i > levs1_i[0])
        XLoc, YLoc = np.where(mask)

        y_evpa = evpa_len * np.cos(evpa[mask])
        x_evpa = evpa_len * np.sin(evpa[mask])

        SelPix = range(0, len(stokes_i), evpa_distance)

        lines = []
        for i in range(0, len(XLoc)):
            if XLoc[i] in SelPix and YLoc[i] in SelPix:
                Xpos = -float(self.clean_image.X[XLoc[i]]) #TODO check where the minus comes from!
                Ypos = -float(self.clean_image.Y[YLoc[i]]) #TODO check where the minus comes from!
                X0 = float(Xpos - x_evpa[i] / 2.)
                X1 = float(Xpos + x_evpa[i] / 2.)
                Y0 = float(Ypos - y_evpa[i] / 2.)
                Y1 = float(Ypos + y_evpa[i] / 2.)
                lines.append(((Y0, X0), (Y1, X1)))
        lines = tuple(lines)

        ##TODO something is still wrong here with the sign and orientation of the EVPA, please double check!


        # plot the evpas
        evpa_lines = LineCollection(lines, colors=self.evpa_color, linewidths=2)
        self.ax.add_collection(evpa_lines)


# takes a an image (2d) array as input and calculates the sigma levels for plotting, sigma_contour_limit denotes the sigma level of the lowest contour
def get_sigma_levs(image,  # 2d array/list
                   sigma_contour_limit=3  # choose the lowest sigma contour to plot
                   ):
    Z1 = image.flatten()
    bin_heights, bin_borders = np.histogram(Z1 - np.min(Z1) + 10 ** (-5), bins="auto")
    bin_widths = np.diff(bin_borders)
    bin_centers = bin_borders[:-1] + bin_widths / 2.
    bin_heights_err = np.where(bin_heights != 0, np.sqrt(bin_heights), 1)

    t_init = models.Gaussian1D(np.max(bin_heights), np.median(Z1 - np.min(Z1) + 10 ** (-5)), 0.001)
    fit_t = fitting.LevMarLSQFitter()
    t = fit_t(t_init, bin_centers, bin_heights, weights=1. / bin_heights_err)
    noise = t.stddev.value

    # Set contourlevels to mean value + 3 * rms_noise * 2 ** x
    levs1 = t.mean.value + np.min(Z1) - 10 ** (-5) + sigma_contour_limit * t.stddev.value * np.logspace(0, 100, 100,
                                                                                                        endpoint=False,
                                                                                                        base=2)
    levs = t.mean.value + np.min(Z1) - 10 ** (-5) - sigma_contour_limit * t.stddev.value * np.logspace(0, 100, 100,
                                                                                                       endpoint=False,
                                                                                                       base=2)
    levs = np.flip(levs)
    levs = np.concatenate((levs, levs1))

    return levs, levs1

#gets components from .fits file
def getComponentInfo(filename):

    #TODO also include reading .mod files

    data_df = pd.DataFrame()
    hdu_list = fits.open(filename)
    comp_data = hdu_list[1].data
    comp_data1 = np.zeros((len(comp_data), len(comp_data[0])))
    date = np.array([])
    year = np.array([])
    mjd = np.array([])
    for j in range(len(comp_data)):
        comp_data1[j, :] = comp_data[j]
        date1=get_date(filename)
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

def get_date(filename):

    hdu_list=fits.open(filename)
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
                
