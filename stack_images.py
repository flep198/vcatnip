#this script can stack VLBI .fits images in stokes parameters and in linear polarization and EVPA (weighted and unweighted)
#it can also do some beam folding operations required prior to stacking

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from astropy.modeling import models, fitting
from matplotlib.collections import LineCollection  
import pexpect
from pexpect import replwrap  
import os
from numpy import linalg

def stack_images(image_array, #input images to be stacked
        weighted=False, #choose whether to use the weighted option
        weights_array=[] #put in weights to use for stacking
        ):
    
    #make sure everything is numpy
    for ind,image in enumerate(image_array):
        image_array[ind]=np.array(image)

    for ind,weights in enumerate(weights_array):
        weights_array[ind]=np.array(weights)
    
    #get some general properties
    dim=image_array[0].shape
    stacked_image=np.empty(shape=dim)
    n_images=len(image_array)
    
    #check if all input images have the same dimension
    wrong_size=False
    for image in image_array:
        if image.shape != dim:
            wrong_size=True
    
    #check that weights_array has same dimension as the image_array
    wrong_weights=False
    if weighted:
        for weights in weights_array:
            if weights.shape != dim:
                wrong_weights = True
        if len(weights_array) != n_images:
            wrong_weights = True 
    else:
        #setup a weights_array with equal weights if weighting is not used
        weights_array=np.empty((n_images,dim[0],dim[1]))
        for i in range(n_images):
            weights_array[i]=np.ones(dim)


    if wrong_size:
        raise Exception("Error! Your input arrays do not all have the same dimension!")
    elif wrong_weights:
        raise Exception("Error! You have selected a weighted stack but your weights are not the same dimension as your images!")
    else:
        for i in range(len(image_array[0])):
            for j in range(len(image_array[0][0])):
                pixel_value=0
                weight_sum=0
                for k in range(n_images):
                    pixel_value+=weights_array[k][i][j]*image_array[k][i][j]
                    weight_sum+=weights_array[k][i][j]
                pixel_value=pixel_value/weight_sum
                stacked_image[i,j]=pixel_value
        
        return stacked_image

#aligns multiple 2d arrays on the brightest pixel centered on a given pixel position
#if align_by is set, the image array will be shifted according to the data in align_by

def align_images(image_array, #input images to be aligned
        center=[-1,-1], #position where to put the brightest pixel, default or negative values put the brightest pixel to the center
        align_by=[] #can input another set of image (equal dimensions as image_array) to align by.
        ):

    #check if there was input to the align_by field
    do_align_by = False
    for ind,image in enumerate(align_by):
        align_by[ind]=np.array(image)
    align_by=np.array(align_by)

    if len(align_by)>0:
        do_align_by = True
    
    #make sure everything is numpy
    for ind,image in enumerate(image_array):
        image_array[ind]=np.array(image)
    image_array=np.array(image_array)
    center=np.array(center)

    if not do_align_by:
        align_by=image_array

    #get some general properties
    dim=image_array[0].shape
    n_images=len(image_array)
    
    #set center pixel
    if np.any(center < 0):
        #if negative values are given for center point (or default) we use the center pixel
        center=np.array([int(round(dim[0]/2-1)),int(round(dim[1]/2)-1)])
    
    #check if all input images have the same dimension
    wrong_size=False
    for image in image_array:
        if image.shape != dim:
            wrong_size=True

    #check if align_by images have the same size as the image array
    wrong_size_align_by = False
    if align_by.shape != image_array.shape:
        wrong_size_align_by=True
 
    if wrong_size:
        raise Exception("Error! The input images have different sizes!")
    elif wrong_size_align_by:
        raise Exception("Error! The size of the align_by images differs from the image size!")
    else:
        for ind,image in enumerate(image_array):
            
            #find indeces of brightest pixel
            brightest_pixel=np.unravel_index(np.argmax(align_by[ind], axis=None), dim)
            #find shift vector
            shift=center-brightest_pixel

            #do the shift
            image_array[ind]=np.roll(image,(shift[0],shift[1]),axis=(0,1))
        
        return image_array


#this function takes file paths to fits files (e.g. from DIFMAP or CASA exportfits) and creates a stacked fits file from them
#by stacking I,Q,U and all IFs independent from each other. With the "align" option, the images will be centered at the brightest pixel in total intensity (aligned)
#it can be either passed fits_files including all polarizations or

def stack_fits(fits_files, #a list of filepaths to fits files (either full polarization or just stokes I)
        stokes_q_fits=[], #a list of filepaths to fits files containing stokes Q
        stokes_u_fits=[], #a list of filepaths to fits files conataining stokes U
        export_fits=False, #choose whether to write an output fits file (stacked)
        output_file="stacked.fits", #choose file name for output file
        overwrite=True, #choose whether to overwrite an already existing image or not
        align=True #choose whether to align the images on the brightest pixel in Stokes I before stacking (all pols will be aligned according to the brightest pixel in Stokes I)
        ):
    
    #check if there is more than one fits file
    wrong_len=False
    if len(fits_files)<1:
        wrong_len=True
    
    #check if fits files are in STOKES format
    not_stokes=True
    for fits_file in fits_files:
        file=fits.open(fits_file)
        for i in range(10):
            try:
                if "STOKES" in file[0].header["CTYPE"+str(i)]:
                    not_stokes=False
            except:
                pass
    
    #check if FITS file contains more than just Stokes I
    only_stokes_i=False
    for fits_file in fits_files:
        file=fits.open(fits_file)
        if file[0].data.shape[0]==1:
            only_stokes_i=True

 
    if wrong_len:
        raise Exception("Error! Please put in more than one fits file, otherwise stacking makes no sense!")
    elif not_stokes:
        raise Exception("Error! Your fits-files are not in STOKES format. This is currently not implemented!")
    else:
        
        if only_stokes_i and (len(stokes_q_fits)!=len(fits_files) or len(stokes_u_fits)!=len(fits_files)):
            print("Warning! Only Stokes I input given!")
            print("-> will produce only Stokes I stacked image")
            pols=1 #do only stokes I in this case
        else:
            pols=3 #do all three polarization Stokes I, Q, U

        dim=fits.open(fits_files[0])[0].data.shape
        output_stacked=np.empty((pols,dim[1],dim[2],dim[3]))

        for pol in range(pols): #iterate over polarizations
            for spw in range(dim[1]): #iterate over spws/IFs
                data_to_stack=np.empty((len(fits_files),dim[2],dim[3]))
                if only_stokes_i:
                    if pol==0:
                        for ind,fits_file in enumerate(fits_files):
                            data_to_stack[ind]=fits.open(fits_file)[0].data[0][spw]
                    if pol==1:
                        for ind,fits_file in enumerate(stokes_q_fits):
                            data_to_stack[ind]=fits.open(fits_file)[0].data[0][spw]
                    if pol==2:
                        for ind,fits_file in enumerate(stokes_u_fits):
                            data_to_stack[ind]=fits.open(fits_file)[0].data[0][spw]
                else:
                    for ind,fits_file in enumerate(fits_files):
                        data_to_stack[ind]=fits.open(fits_file)[0].data[pol][spw]
                
                #store stokes i data for later use
                if pol==0:
                    i_data_to_stack=data_to_stack

                #align polarizations according to stokes i peak
                if align:
                    output_stacked[pol][spw]=stack_images(align_images(data_to_stack,align_by=i_data_to_stack))
                else:
                    output_stacked[pol][spw]=stack_images(data_to_stack)

        if export_fits:
            file=fits.open(fits_files[0])
            file[0].data=output_stacked
            file.writeto(output_file,overwrite=overwrite)
        
        return output_stacked

#this funciton takes file paths to fits files (e.g. from DIFMAP) and creates a stacked fits file from them
#by first calculation linear polarization P and EVPA and stacking P and EVPA and NOT Q,U. If weighted is set to true,
#the EVPA stack is weighted with the linear Polarization. The "align" option centers all polarizations on the Stokes I peak before stacking

def stack_pol_fits(fits_files, #list of file paths to fits files with full polarization or Stokes I only data
        stokes_q_fits=[], #list of file paths to fits files with Stokes Q data
        stokes_u_fits=[], #list of file paths to fits files with Stokes U data
        weighted=False, #choose whether to weight the EVPA stacking by the level of linear polarization
        align=True #choose whether to align the images on the brightest pixel in Stokes I before stacking (all pols will be aligned according to the brightes pixel in Stokes I)
        ):

    #check if there is more than one fits file
    wrong_len=False
    if len(fits_files)<1:
        wrong_len=True
    
    #check if fits files are in STOKES format
    not_stokes=True
    for fits_file in fits_files:
        file=fits.open(fits_file)
        for i in range(10):
            try:
                if "STOKES" in file[0].header["CTYPE"+str(i)]:
                    not_stokes=False
            except:
                pass

    #check if FITS file contains more than just Stokes I
    only_stokes_i=False
    for fits_file in fits_files:
        file=fits.open(fits_file)
        if file[0].data.shape[0]==1:
            only_stokes_i=True

    if wrong_len:
        raise Exception("Error! Please put in more than one fits file, otherwise stacking makes no sense!")
    elif not_stokes:
        raise Exception("Error! Your fits-files are not in STOKES format. This is currently not implemented!")
    elif only_stokes_i and (len(stokes_q_fits)!=len(fits_files) or len(stokes_u_fits)!=len(fits_files)):
        raise Exception("Error! At least one of your input files contains only Stokes I but no polarization information and you did not specify a matching Stokes Q and U path!")
    else:
    
        dim=fits.open(fits_files[0])[0].data.shape
        output_stacked=np.empty((3,dim[1],dim[2],dim[3]))

        #read in all fits files
        image_storage=np.empty((len(fits_files),3,dim[1],dim[2],dim[3]))
        for ind,fits_file in enumerate(fits_files):
           
            if only_stokes_i: #if the fits files only contain one polarization, it reads from the different input files
                image_storage[ind]=[fits.open(fits_file)[0].data[0],fits.open(stokes_q_fits[ind])[0].data[0],fits.open(stokes_u_fits[ind])[0].data[0]]
            else: #if the fits file contains multiple polarization it will read from the single fits file
                image_storage[ind]=fits.open(fits_file)[0].data[0:3]     

        #calculate linear polarization and EVPAs from Stokes Q and U for all input images and IFs/spw
        pol_images=np.empty((dim[1],3,len(fits_files),dim[2],dim[3])) #-> here we will only save stokes i,lin.pol and evpa
        for ind,image in enumerate(image_storage):
            for spw in range(dim[1]):
                stokes_i=image[0][spw]
                stokes_q=image[1][spw]
                stokes_u=image[2][spw]
                pol_images[spw][0][ind] = stokes_i #Stokes I
                pol_images[spw][1][ind] = np.sqrt(stokes_q**2+stokes_u**2) #linear polarization
                pol_images[spw][2][ind] = 0.5*np.arctan2(stokes_u,stokes_q) #EVPA
    
        #do the actual stacking
        for spw in range(dim[1]):
            
            if align:
                output_stacked[0][spw]=stack_images(align_images(pol_images[spw][0])) #stack stokes_i
                output_stacked[1][spw]=stack_images(align_images(pol_images[spw][1],align_by=pol_images[spw][0])) #stack lin.pol.
                if weighted:
                    output_stacked[2][spw]=stack_images(align_images(pol_images[spw][2],align_by=pol_images[spw][0]),weighted=True,weights_array=pol_images[spw][1])    
                else:
                    output_stacked[2][spw]=stack_images(align_images(pol_images[spw][2],align_by=pol_images[spw][0])) #stack EVPA without weights

            else:
                output_stacked[0][spw]=stack_images(pol_images[spw][0]) #stack stokes_i
                output_stacked[1][spw]=stack_images(pol_images[spw][1]) #stack lin.pol.
                if weighted:
                    output_stacked[2][spw]=stack_images(pol_images[spw][2],weighted=True,weights_array=pol_images[spw][1])    
                else:
                    output_stacked[2][spw]=stack_images(pol_images[spw][2]) #stack EVPA without weights

        return output_stacked
        
#takes a an image (2d) array as input and calculates the sigma levels for plotting, sigma_contour_limit denotes the sigma level of the lowest contour
def get_sigma_levs(image, #2d array/list
        sigma_contour_limit=3 #choose the lowest sigma contour to plot
        ):
    
    Z1 = image.flatten()
    bin_heights, bin_borders = np.histogram(Z1 - np.min(Z1) + 10 ** (-5), bins = "auto")
    bin_widths = np.diff(bin_borders)
    bin_centers = bin_borders[:-1] + bin_widths / 2.
    bin_heights_err = np.where(bin_heights != 0, np.sqrt(bin_heights), 1)
    
    t_init = models.Gaussian1D(np.max(bin_heights), np.median(Z1 - np.min(Z1) + 10 ** (-5)), 0.001)
    fit_t = fitting.LevMarLSQFitter()
    t = fit_t(t_init, bin_centers, bin_heights, weights = 1. / bin_heights_err)
    noise=t.stddev.value

    #Set contourlevels to mean value + 3 * rms_noise * 2 ** x
    levs1 = t.mean.value + np.min(Z1) - 10 ** (-5) + sigma_contour_limit * t.stddev.value * np.logspace(0, 100, 100, endpoint = False, base = 2)
    levs = t.mean.value + np.min(Z1) - 10 ** (-5) - sigma_contour_limit * t.stddev.value * np.logspace(0, 100, 100, endpoint = False, base = 2)
    levs = np.flip(levs)
    levs = np.concatenate((levs, levs1))
    
    return levs,levs1[0]

#function to plot polarization image either with stokes q and u as input or with lin. pol. and evpa
def plot_pol_image(stokes_i, #2d list/array of stokes i data
        lin_pol=[], #2d list/array of lin pol data
        evpa=[], #2d list/array of evpa data
        pol_from_stokes=False, #choose whether to use lin_pol & evpa input (False) OR stokes_q and stokes_u input (TRUE)
        stokes_q=[], #2d list/array of stokes_q data
        stokes_u=[], #2d list/array of stokes_u data
        evpa_len=4, #choose length of EVPA in pixels
        stokes_i_sigma_cut=3, #choose lowest sigma contour for Stokes I plot
        lin_pol_sigma_cut=3, #choose lowest sigma contour for Lin Pol plot
        evpa_distance=10, #choose distance of EVPA vectors to draw
        rotate_evpa=0, #rotate EVPAs by a given angle in degrees (North through East) 
        cmap="inferno", #choose matplotlib colormap for plot
        rcparams={} #option to modify matplotlib look
        ):
    
    #update plot look
    plt.rcParams.update(rcparams)

    #Check if linpol/evpa/stokes_i have same dimensions!
    dim_wrong=True
    if pol_from_stokes:
        if np.shape(stokes_i) == np.shape(stokes_q) and np.shape(stokes_i) == np.shape(stokes_u):
            dim_wrong=False
    else:
        if np.shape(stokes_i) == np.shape(lin_pol) and np.shape(stokes_i) == np.shape(evpa):
            dim_wrong=False


    X = np.linspace(0,len(stokes_i),len(stokes_i))
    Y = np.linspace(0,len(stokes_i[0]),len(stokes_i[0]))

    #plot total intensity (stokes i)

    plot_levs_i, lowest_contour_i = get_sigma_levs(stokes_i,sigma_contour_limit=stokes_i_sigma_cut)
    fig, ax = plt.subplots()

    ax.imshow(stokes_i, norm = colors.SymLogNorm(linthresh = lowest_contour_i, 
        linscale = 0.5, 
        vmin = 1 * lowest_contour_i, 
        vmax = 0.5 * np.max(stokes_i),base=10), 
        origin="lower",
        cmap=cmap)
    ax.contour(X,Y,stokes_i,levels=plot_levs_i,colors="white")

    plt.show()

    if dim_wrong:
        raise Exception("Warning! Polarization data has wrong dimensions, could not produce polarization plot.")
    else:
        #calculate EVPA and lin pol from stokes q and u
        if pol_from_stokes:
            lin_pol=np.sqrt(stokes_q**2+stokes_u**2)
            evpa=0.5*np.arctan2(stokes_u,stokes_q)
   
        lin_pol[np.isnan(lin_pol)]=0.0

        #plot linear polarization

        plot_levs_lin_pol, lowest_contour_lin_pol = get_sigma_levs(lin_pol,sigma_contour_limit=lin_pol_sigma_cut) #TRY USING stokes i levels and noise here
        fig, ax = plt.subplots()

        ax.imshow(lin_pol, norm = colors.SymLogNorm(linthresh = lowest_contour_lin_pol, 
            linscale = 0.5, 
            vmin = 1 * lowest_contour_lin_pol, 
            vmax = 0.5 * np.max(lin_pol),base=10), 
            origin="lower",
            cmap=cmap)
        ax.contour(X,Y,lin_pol,levels=plot_levs_lin_pol,colors="white")

        #plot EVPA
        evpa=evpa+rotate_evpa/180*np.pi

        #create mask where to plot EVPA (only where stokes i and lin pol have plotted contours)
        mask = np.zeros(np.shape(stokes_i),dtype=bool)
        mask[:] = (lin_pol > lowest_contour_lin_pol)*(stokes_i > lowest_contour_i)
        XLoc, YLoc = np.where(mask)

        y_evpa=evpa_len*np.cos(evpa[mask])
        x_evpa=evpa_len*np.sin(evpa[mask])

        SelPix = range(0,len(stokes_i),evpa_distance)
    
        lines=[]
        for i in range(0,len(XLoc)):
            if XLoc[i] in SelPix and YLoc[i] in SelPix:
                Xpos = float(XLoc[i])
                Ypos = float(YLoc[i])
                X0 = float(Xpos-x_evpa[i]/2.)
                X1 = float(Xpos+x_evpa[i]/2.)
                Y0 = float(Ypos-y_evpa[i]/2.)
                Y1 = float(Ypos+y_evpa[i]/2.)
                lines.append(((Y0,X0),(Y1,X1)))    
        lines = tuple(lines)

        #plot the evpas
        evpa_lines = LineCollection(lines,colors='r',linewidths=2)
        ax.add_collection(evpa_lines)

        plt.show()

#this method is intended to determine the (smallest) common beam of multiple fits-images

def get_common_beam(fits_files,
        ppe=100, #sample points used per input ellipse
        plot_beams=False, #makes a simple plot of all the beams for double checking
        tolerance=0.0001 #adjust the numeric tolarance for determining the common beam
        ):
       
    fig = plt.figure()
    ax = fig.add_subplot()

    sample_points=np.empty(shape=(ppe*len(fits_files),2))
        
    for ind,file_path in enumerate(fits_files):
        file=fits.open(file_path)
        bmaj = float(file[0].header["BMAJ"]*60*60*1000)
        bmin = float(file[0].header["BMIN"]*60*60*1000)
        posa = float(file[0].header["BPA"])/180*np.pi
               
        #sample ellipse points
        
        ellipse_angles=np.linspace(0,2*np.pi,ppe)
        X=-bmin/2*np.sin(ellipse_angles)
        Y=bmaj/2*np.cos(ellipse_angles)
        
        #rotate them according to position angle
        X_rot=-X*np.cos(posa)-Y*np.sin(posa)
        Y_rot=X*np.sin(posa)+Y*np.cos(posa)

        for i in range(ppe):
            sample_points[ind*ppe+i]=np.array([X_rot[i],Y_rot[i]])
        if plot_beams:
            plt.plot(X_rot,Y_rot,c="k")

    #find minimum ellipse
    (center,radii,rotation)=getMinVolEllipse(sample_points,tolerance=tolerance)  

       
    #find out bmaj, bmin and posa
    bmaj_ind=np.argmax(radii)
    
    if bmaj_ind==0:
        bmaj=2*radii[0]
        bmin=2*radii[1]
        posa=-np.arcsin(rotation[1][0])/np.pi*180-90
    else:
        bmaj=2*radii[1]
        bmin=2*radii[0]
        posa=-np.arcsin(rotation[1][0])/np.pi*180

    #make posa from -90 to +90
    if posa>90:
        posa=posa-180
    elif posa<-90:
        posa=posa+180

    # plot ellipsoid
    if plot_beams:
        from matplotlib import patches
        ellipse=patches.Ellipse(center,bmin,bmaj,angle=posa,fill=False,zorder=2,linewidth=2,color="r")
        ax.add_patch(ellipse)

        ax.axis("equal")
        plt.show()
     

    return bmaj,bmin,posa


#gets common ellipse from point selection (needed for smallest common beam calculation), adapted from https://github.com/minillinim/ellipsoid/blob/master/ellipsoid.py 
def getMinVolEllipse(P=None, tolerance=0.1):
        """ Find the minimum volume ellipsoid which holds all the points

        Based on work by Nima Moshtagh
        http://www.mathworks.com/matlabcentral/fileexchange/9542
        and also by looking at:
        http://cctbx.sourceforge.net/current/python/scitbx.math.minimum_covering_ellipsoid.html
        Which is based on the first reference anyway!

        Here, P is a numpy array of 2 dimensional points like this:
        P = [[x,y], <-- one point per line
             [x,y],
             [x,y]]

        Returns:
        (center, radii, rotation)

        """

        (N, d) = np.shape(P)
        d = float(d)

        # Q will be our working array
        Q = np.vstack([np.copy(P.T), np.ones(N)])
        QT = Q.T

        # initializations
        err = 1.0 + tolerance
        u = (1.0 / N) * np.ones(N)

        # Khachiyan Algorithm
        while err > tolerance:
            V = np.dot(Q, np.dot(np.diag(u), QT))
            M = np.diag(np.dot(QT , np.dot(linalg.inv(V), Q)))    # M the diagonal vector of an NxN matrix
            j = np.argmax(M)
            maximum = M[j]
            step_size = (maximum - d - 1.0) / ((d + 1.0) * (maximum - 1.0))
            new_u = (1.0 - step_size) * u
            new_u[j] += step_size
            err = np.linalg.norm(new_u - u)
            u = new_u

        # center of the ellipse
        center = np.dot(P.T, u)

        # the A matrix for the ellipse
        A = linalg.inv(
                       np.dot(P.T, np.dot(np.diag(u), P)) -
                       np.array([[a * b for b in center] for a in center])
                       ) / d
        # Get the values we'd like to return
        U, s, rotation = linalg.svd(A)
        radii = 1.0/np.sqrt(s)

        return (center, radii, rotation)

#this method folds images with a defined beam (either custom input or common beam), needs .uvf,.mod and .fits files from the original images as inputs
#requires the definition of "difmap_path", since it runs difmap in the background, so you need to point difmap_path to where your difmap executable sits

def fold_with_beam(fits_files, #array of file paths to fits images input
        difmap_path, #path to where the difmap executable is located
        bmaj=-1, #beam major axis to fold with (in mas)
        bmin=-1, #beam minor axis to fold with (in mas)
        posa=-1, #beam position angle to fold with (in deg)
        channel="i", #polarization channel to use for folding (default stokes "i"), possible options "q","u","v"
        output_dir="output", #Name and path to the output directory
        n_pixel=2048, #number of pixels in output image
        pixel_size=-1, #pixel size in mas (default uses 1/10 of bmin)
        use_common_beam=False, #choose whether to fold with the common beam of the input arrays (TRUE) or whether to use the input bmaj,bmin,posa beam (FALSE, default)  
        mod_files=[], #optional input array of file paths to .mod files for fits_files
        uvf_files=[] #optional input array of file paths to .uvf files for fits_files
        ):

        #check if custom beam is correctly defined when using it
        if not use_common_beam and (bmaj<0 or bmin<0 or posa<0):
            raise Exception("Please define a sensible custom beam using 'bmaj', 'bmin' and 'posa' kwargs or select choose use_common_beam=True to use the common beam.")
        elif use_common_beam:
            bmaj,bmin,posa=get_common_beam(fits_files) 
        
        #set pixel_size if not manually set:
        if pixel_size<0:
            pixel_size=bmin/10

        #check if there was input to mod_files
        if len(fits_files)!=len(mod_files):
            mod_files=[]
            print("No or insufficient number of mod files defined. Will try to guess their names from .fits file names")
            for fits_file in fits_files:
                mod_files=np.append(mod_files,fits_file[:-5]+".mod")
        
        #check if there was input to uvf_files
        if len(fits_files)!=len(uvf_files):
            uvf_files=[]
            print("No or insufficient number of uvf files defined. Will try to guess their names from .fits file names")
            for fits_file in fits_files:
                uvf_files=np.append(uvf_files,fits_file[:-5]+".uvf")

        #create output directory
        if not os.path.isdir(output_dir):
            os.makedirs(output_dir,exist_ok=True)
        
        #add difmap to PATH
        if difmap_path != None and not difmap_path in os.environ['PATH']:
            os.environ['PATH'] = os.environ['PATH'] + ':{0}'.format(difmap_path) 
            
        # Initialize difmap call
        child = pexpect.spawn('difmap', encoding='utf-8', echo=False)
        child.expect_exact("0>",None, 2)

        def send_difmap_command(command,prompt="0>"):
            child.sendline(command)
            child.expect_exact(prompt, None, 2)

        for ind, fits_file in enumerate(fits_files):
            send_difmap_command("obs " + uvf_files[ind])
            send_difmap_command("uvw 0,-1") #use natural weighting
            send_difmap_command("select " + channel)
            send_difmap_command("rmod " + mod_files[ind])
            send_difmap_command("maps " + str(n_pixel) + "," + str(pixel_size))
            send_difmap_command("restore " + str(bmaj) + "," + str(bmin) + "," + str(posa))
            send_difmap_command("save " + output_dir + "/" + '.'.join(fits_file.split("/")[-1].split(".")[0:-1])+"_convolved")
        
        os.system("rm -rf difmap.log*")
        
        print("Convolving complete!")
    


if __name__ == '__main__':
    
    #Some basic testing

    import glob
    #stack1=stack_fits(["0506+056.u.2023_07_01.icn.fits","0506+056.u.2023_07_01.icn.fits"],stokes_q_fits=["0506+056.u.2023_07_01.icn.fits","0506+056.u.2023_07_01.icn.fits"],stokes_u_fits=["0506+056.u.2023_07_01.icn.fits","0506+056.u.2023_07_01.icn.fits"])
    #plot_pol_image(stack1[0][0],lin_pol=stack1[1][0],evpa=stack1[2][0],stokes_i_sigma_cut=5)#),pol_from_stokes=True,stokes_q=stack1[1][0],stokes_u=stack1[2][0])
    #stack2=stack_fits(glob.glob("*.fits"))
    #plot_pol_image(stack2[0][0])#,lin_pol=stack2[1][0],evpa=stack2[2][0])

    print(fold_with_beam(glob.glob("*.fits"),difmap_path="/opt/uvf_difmap_2.5p/",use_common_beam=True))


    #stack1=stack_fits(["BR224AKu_imaged/polarization/0506_final_pol.fits","BR224BKu_imaged/polarization/0506_final_pol.fits"])
    #plot_pol_image(stack1[0][0],pol_from_stokes=True,stokes_q=stack1[1][0],stokes_u=stack1[2][0])

    #stack2=stack_pol_fits(["BR224AKu_imaged/polarization/0506_final_pol.fits","BR224BKu_imaged/polarization/0506_final_pol.fits"],weighted=True)
    #plot_pol_image(stack2[0][0],lin_pol=stack2[1][0],evpa=stack2[2][0])

    #file=fits.open("3C279_final_pol.fits")[0].data
    #plot_pol_image(file[0][0],pol_from_stokes=True,stokes_q=file[1][0],stokes_u=file[2][0])
    #print(align_images([[[1,4,3,4],[1,2,3,3],[1,2,3,3],[1,2,3,3]],[[1,2,3,3],[1,2,3,4],[1,2,3,3],[5,2,3,3]]]))
