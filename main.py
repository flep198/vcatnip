from kivy.app import App
from kivy.uix.popup import Popup
from kivy.uix.togglebutton import ToggleButton
from kivy.uix.tabbedpanel import TabbedPanel
from kivy.uix.boxlayout import BoxLayout
from kivy.uix.label import Label
from kivy.properties import ObjectProperty, StringProperty, ListProperty
from kivy.uix.button import Button
from kivy.uix.checkbox import CheckBox
from kivy.uix.anchorlayout import AnchorLayout
from graph_widget import MatplotFigure
from kivy.utils import platform
from kinematics import ComponentCollection
import numpy as np
from graph_generator import FitsImage, ImageData, KinematicPlot, get_date, getComponentInfo
from astroquery.ipac.ned import Ned
from astropy.io import fits
import os
import pandas as pd
import glob
from stack_images import stack_fits, stack_pol_fits, get_common_beam, fold_with_beam
import subprocess


#avoid conflict between mouse provider and touch (very important with touch device)
#no need for android platform
if platform != 'android':
    from kivy.config import Config
    Config.set('input', 'mouse', 'mouse,disable_on_activity')

class FileChoosePopup(Popup):
    load = ObjectProperty()

class FileSavePopup(Popup):
    load = ObjectProperty()

    def load_filepath_to_view(self, selection):
        if selection!=[]:
            self.ids.filename_save.text = str(selection[0])

class FileImportPopup(Popup):
    load = ObjectProperty()

class ModelFits(TabbedPanel):
    file_info = StringProperty("No file chosen")
    the_popup = ObjectProperty(None)
    save_dialog = ObjectProperty(None)
    modelfit_filepaths = ListProperty([])
    clean_filepaths = ListProperty([])
    stokes_q_filepaths = ListProperty([])
    stokes_u_filepaths = ListProperty([])
    uvf_filepaths = ListProperty([])
    casa_clean_model_filepaths = ListProperty([])
    plots = []
    components = []
    active_component_ind=None
    component_colors=["green","red","blue","yellow","orange","purple"]
    figure_widgets=[]
    core_component_ind=None
    name=""
    component_collections = []

    #### STACKING variables
    stacking_single_plots = []
    stacking_single_plot_buttons = []
    stacking_single_plot_checkboxes = []
    final_stack_image = ""

    #### PLOTTING variables
    plotting_single_plots = []
    plotting_single_plots_data = []
    plotting_single_plot_buttons = []
    plotting_single_plot_checkboxes = []

    #### GENERAL FUNCTIONS ACROSS TABS

    def show_popup(self,title,text,button_text):
        popup = Popup(title=title, auto_dismiss=False, size_hint=(0.35,0.35))
        box_layout=BoxLayout(orientation="vertical")
        box_layout.add_widget(Label(text=text,size_hint_y=0.8))
        close_button=Button(text=button_text,size_hint_y=0.2)
        close_button.bind(on_press=popup.dismiss)
        box_layout.add_widget(close_button)
        popup.add_widget(box_layout)
        popup.open()

    def sort_fits_by_date(self,fits_files):
        fits_files = np.array(fits_files)
        if len(fits_files)>0:
            date=[]

            for filepath in fits_files:
                plot_data=ImageData(filepath)
                plot=FitsImage(plot_data)
                date=np.append(date,get_date(filepath))

            args = date.argsort()
            fits_files = fits_files[args]
        return fits_files.tolist()

    def sort_uvf_by_date(self,uvf_files):
        uvf_files = np.array(uvf_files)
        if len(uvf_files) > 0:
            date = []

            for filepath in uvf_files:
                date = np.append(date, fits.open(filepath)[0].header["DATE-OBS"])
            args = date.argsort()
            uvf_files = uvf_files[args]
        return uvf_files.tolist()
    #### END OF GENERAL FUNCTIONS USED ACROSS TABS

    #### START OF IMPORT FUNCTIONS

    def load_uvf(self,selection):
        self.plots = []
        # sort by date
        self.uvf_filepaths = self.sort_uvf_by_date(selection)
        try:
            self.the_popup.dismiss()
        except:
            pass
        self.file_info = str(len(self.uvf_filepaths)) + " Files selected"
        self.ids.get_file_uvf.text = self.file_info

    def load_modelfit(self, selection):
        self.plots=[]
        # sort by date
        self.modelfit_filepaths = self.sort_fits_by_date(selection)
        try:
            self.the_popup.dismiss()
        except:
            pass
        self.file_info=str(len(self.modelfit_filepaths))+" Files selected"
        self.ids.get_file_modelfit.text = self.file_info

    def load_clean(self, selection):
        self.plots=[]
        self.clean_filepaths = self.sort_fits_by_date(selection)
        try:
            self.the_popup.dismiss()
        except:
            pass
        self.file_info=str(len(self.clean_filepaths))+" Files selected"
        self.ids.get_file_clean.text = self.file_info

    def load_stokes_q(self, selection):
        self.plots=[]
        self.stokes_q_filepaths = self.sort_fits_by_date(selection)
        try:
            self.the_popup.dismiss()
        except:
            pass
        self.file_info=str(len(self.stokes_q_filepaths))+" Files selected"
        self.ids.get_file_stokes_q.text = self.file_info

    def load_stokes_u(self, selection):
        self.plots=[]
        self.stokes_u_filepaths = self.sort_fits_by_date(selection)
        try:
            self.the_popup.dismiss()
        except:
            pass
        self.file_info=str(len(self.stokes_u_filepaths))+" Files selected"
        self.ids.get_file_stokes_u.text = self.file_info

    def load_casa_clean_model(self,selection):
        self.plots=[]
        self.casa_clean_model_filepaths = self.sort_fits_by_date(selection)
        try:
            self.the_popup.dismiss()
        except:
            pass
        self.file_info = str(len(self.casa_clean_model_filepaths)) + " Files selected"
        self.ids.get_file_casa_clean_model.text = self.file_info


    #### END OF IMPORT FUNCTIONS

    #### START OF KINEMATIC FUNCTIONS

    def open_popup(self,type):
        if type == "modelfit":
            self.the_popup = FileChoosePopup(load=self.load_modelfit)
        elif type == "uvf":
            self.the_popup = FileChoosePopup(load=self.load_uvf)
        elif type == "clean":
            self.the_popup = FileChoosePopup(load=self.load_clean)
        elif type == "stokes_q":
            self.the_popup = FileChoosePopup(load=self.load_stokes_q)
        elif type == "stokes_u":
            self.the_popup = FileChoosePopup(load=self.load_stokes_u)
        elif type == "casa_clean_model":
            self.the_popup = FileChoosePopup(load=self.load_casa_clean_model)
        self.the_popup.open()

    def open_save_dialog(self):
        self.save_dialog = FileSavePopup()
        self.save_dialog.open()

    def open_import_dialog(self):
        self.import_dialog = FileImportPopup()
        self.import_dialog.open()

    def current_color(self,ind):
        if ind is not None:
            ind=int(ind % len(self.component_colors))
            return self.component_colors[ind]
        else:
            return None

    def create_kinematic_plots(self):
        #sort files by date
        date=[]
        fits_images=[]

        if len(self.clean_filepaths)!=len(self.modelfit_filepaths) and len(self.clean_filepaths)>0:
            self.show_popup("Warning","Please use an equal number of modelfit and clean images","Continue")
            self.modelfit_filepaths=[]
            self.clean_filepaths=[]
        elif len(self.clean_filepaths) == 0 and len(self.modelfit_filepaths) == 0:
            self.show_popup("Warning", "No data selected", "Continue")
        elif len(self.clean_filepaths) == len(self.modelfit_filepaths):
            #create plots for view page
            for ind,filepath in enumerate(self.modelfit_filepaths):
                plot_data=ImageData(self.clean_filepaths[ind],model=filepath)
                plot=FitsImage(plot_data,overplot_gauss=True)
                fits_images=np.append(fits_images,plot)
            self.show_popup("Information", "File loading completed. Have fun doing kinematics!", "Continue")
        elif len(self.clean_filepaths) == 0 and len(self.modelfit_filepaths)>0:
            # create plots for view page
            for filepath in self.modelfit_filepaths:
                plot_data = ImageData(filepath,model=filepath)
                plot = FitsImage(plot_data,overplot_gauss=True)
                fits_images = np.append(fits_images, plot)
            self.show_popup("Warning",
                            "No clean images imported, using only modelfit images.\n Have fun doing kinematics!",
                            "Continue")
        elif len(self.clean_filepaths) > 0 and len(self.modelfit_filepaths) ==0:
            self.show_popup("Warning","No modelfits imported","Continue")



        for plot in fits_images:
            self.plots.append(plot)
            self.name=plot.name
            new_figure=MatplotFigure(size_hint_x=None,
            width=100)
            new_figure.figure=plot.fig
            new_figure.touch_mode="pan"
            self.ids.figures.add_widget(new_figure)
            self.figure_widgets.append(new_figure)

    def change_plot_lims(self):
        #try since the text might be no numbers

        try:
            ra_min = float(self.ids.input_ra_min.text)
            ra_max = float(self.ids.input_ra_max.text)
            dec_min = float(self.ids.input_dec_min.text)
            dec_max = float(self.ids.input_dec_max.text)

            for plot in self.plots:
                plot.change_plot_lim(ra_max,ra_min,dec_min,dec_max)
                plot.ax.figure.canvas.draw_idle()
                plot.ax.figure.canvas.flush_events()
        except:
            pass



    def add_component(self):
        count= len(self.components)
        while "Component " + str(count) in self.components:
            count+=1
        self.components.append("Component " + str(count))

        button=ToggleButton(
            text=self.components[-1],
            size_hint_y=None,
            height=50,
            group="components",
        )
        button.bind(on_release=self.set_active_component)
        self.ids.component_list.add_widget(button)

        #create list entry in kinematic table

        box_id="kinematic"+str(count)
        box=BoxLayout(orientation="horizontal",
                      size_hint_y=None,
                      height=50)

        #Name Field
        label_id=box_id+"_id"
        component_id_label=Label(text=str(count))
        box.add_widget(component_id_label)

        #Speed Field
        component_id_label=Label(text="")
        box.add_widget(component_id_label)

        #Beta_app Field
        component_id_label = Label(text="")
        box.add_widget(component_id_label)

        #Doppler Field
        component_id_label = Label(text="")
        box.add_widget(component_id_label)

        self.ids.kinematic_list.add_widget(box)


    def remove_component(self):
        #find active button
        ind_to_remove=self.active_component_ind

        #reset core component if core component is removed
        if ind_to_remove == self.core_component_ind:

            #reset core component ind to None
            self.core_component_ind = None

            #reset information in component elements
            for plot in self.plots:
                for comp in plot.components:
                    comp[1].is_core=True
                    comp[1].delta_x_est=comp[1].x
                    comp[1].delta_y_est=comp[1].y
                    comp[1].set_distance_to_core(0,0)


        self.active_component_ind=None
        button_to_remove = next( (t for t in ToggleButton.get_widgets('components') if t.state=='down'), None)
        #remove button
        if button_to_remove:
            self.ids.component_list.remove_widget(button_to_remove)
            for plot in self.plots:
                for comp in plot.components:
                    if comp[1].component_number == ind_to_remove:
                        comp[1].component_number=-1 #reset component assignment
                        comp[0][0].set_color(plot.component_color)
                        comp[0][1][0].set_color(plot.component_color)
                        comp[0][2][0].set_color(plot.component_color)
                        plot.ax.figure.canvas.draw_idle()
                        plot.ax.figure.canvas.flush_events()

        #remove component from kinematic table
        for box in self.ids.kinematic_list.children:
            for row in box.children:
                if row.text == str(ind_to_remove):
                    self.ids.kinematic_list.remove_widget(box)

        self.update_kinematic_plot()

    def load_filepath_to_view(self, selection):
        if selection!=[]:
            self.ids.filename_save.text = str(selection[0])

    def choose_component(self,figure_widget,x,y):

        #check what plot is currently displayed
        for plot in self.plots:
            if figure_widget.figure == plot.fig:
                current_plot=plot
        fig_comps=current_plot.components
        final_dist=0

        for ind,comp in enumerate(fig_comps):
            x_ellipse,y_ellipse=comp[0][0].get_center()
            #calculate distance between click and component center
            dist=(x-x_ellipse)**2+(y-y_ellipse)**2
            if dist<final_dist or final_dist==0:
                final_ind=ind
                final_comp=comp
                final_dist=dist

        current_color=self.current_color(self.active_component_ind)

        if self.active_component_ind is not None:
            #component is already selected -> unselect it
            if final_comp[0][1][0].get_color() == current_color:
                final_comp[0][0].set_color(current_plot.component_color)
                final_comp[0][1][0].set_color(current_plot.component_color)
                final_comp[0][2][0].set_color(current_plot.component_color)
                #reset component assignment
                final_comp[1].assign_component_number(-1)
            #check if component is still unassigend -> assign it
            elif final_comp[0][1][0].get_color() == current_plot.component_color:
                final_comp[0][0].set_color(current_color)
                final_comp[0][1][0].set_color(current_color)
                final_comp[0][2][0].set_color(current_color)
                #assign component to number for kinematics
                final_comp[1].assign_component_number(self.active_component_ind)
            #unselect all previously selected components of this color
            for comp in fig_comps:
                if comp != final_comp and comp[0][1][0].get_color() == current_color:
                    comp[0][0].set_color(current_plot.component_color)
                    comp[0][1][0].set_color(current_plot.component_color)
                    comp[0][2][0].set_color(current_plot.component_color)
                    comp[1].assign_component_number(-1)

        #update kinematic plot
        self.update_kinematic_plot()

    def update_core_dist(self):

        for plot in self.plots:
            core_set=False
            for comp in plot.components:
                if comp[1].component_number == self.core_component_ind:
                    comp[1].is_core = True
                    core_x = comp[1].x
                    core_y = comp[1].y
                    core_set = True
                else:
                    comp[1].is_core = False
            if core_set:
                for comp in plot.components:
                    comp[1].set_distance_to_core(core_x, core_y)

    def set_core_component(self):

        self.core_component_ind = self.active_component_ind

        #remove "Core" from all buttons
        for button in ToggleButton.get_widgets('components'):
            if " (Core)" in button.text:
                button.text=button.text.replace(" (Core)","")

        #Add "Core" to button text
        active_button = next((t for t in ToggleButton.get_widgets('components') if t.state == 'down'), None)
        if active_button:
            active_button.text=active_button.text+" (Core)"

        #update kinematics
        self.update_kinematic_plot()

    def set_active_component(self,button=""):
        # find active button
        active_button = next((t for t in ToggleButton.get_widgets('components') if t.state == 'down'), None)
        if active_button:
            #set active component
            for i,comp in enumerate(self.components):
                if comp in active_button.text:
                    self.active_component_ind=i
        else:
            self.active_component_ind = None


    def update_kinematic_plot(self,do_export=False,export_path=""):

        #currently creates a new plot everytime => maybe rewrite to just update existing plot?
        self.kplot=KinematicPlot()
        self.ids.kinematic_plot.figure=self.kplot.fig
        self.component_collections=[]
        self.update_core_dist()

        t_max=0
        if len(self.plots)>0:
            t_min=self.plots[0].components[0][1].year
            d_max=0
            tb_max=0
            tb_min=self.plots[0].components[0][1].tb
            flux_max=0
            for i in range(len(self.components)):
                collection=[]
                for plot in self.plots:
                    for comp in plot.components:
                        if comp[1].component_number == i:
                            collection.append(comp[1])
                            # determine plot limits
                            if comp[1].year > t_max:
                                t_max = comp[1].year
                            if comp[1].year < t_min:
                                t_min = comp[1].year
                            if comp[1].distance_to_core*comp[1].scale > d_max:
                                d_max = comp[1].distance_to_core*comp[1].scale
                            if comp[1].tb > tb_max:
                                tb_max = comp[1].tb
                            if comp[1].tb < tb_min:
                                tb_min = comp[1].tb
                            if comp[1].flux > flux_max:
                                flux_max = comp[1].flux

                comp_collection=ComponentCollection(collection,name=self.components[i])
                self.component_collections.append(comp_collection)

                # find out which plot is currently selected:
                active_button = next((t for t in ToggleButton.get_widgets('kinematic_select') if t.state == 'down'), None)

                if active_button.text == "Flux Density":
                    self.kplot.plot_fluxs(comp_collection,self.current_color(i))
                    self.kplot.set_limits([t_min - 0.1 * (t_max - t_min), t_max + 0.1 * (t_max - t_min)],
                                          [0, 1.2 * flux_max])
                elif active_button.text == "Brightness Temperature":
                    self.kplot.plot_tbs(comp_collection,self.current_color(i))
                    self.kplot.set_limits([t_min - 0.1 * (t_max - t_min), t_max + 0.1 * (t_max - t_min)],
                                          [tb_min, 10 * tb_max])
                else:
                    self.kplot.plot_kinematics(comp_collection,self.current_color(i))
                    self.kplot.set_limits([t_min - 0.1 * (t_max - t_min), t_max + 0.1 * (t_max - t_min)],
                                          [0, 1.2 * d_max])

                fit_data = comp_collection.get_speed()
                if len(collection)>2 and active_button.text == "Kinematic":
                    self.kplot.plot_linear_fit(t_min-0.1*(t_max-t_min),t_max+0.1*(t_max-t_min),
                                            fit_data["speed"],
                                            fit_data["y0"],self.current_color(i),
                                            label=self.components[i])

                #write data to table
                for row in self.ids.kinematic_list.children:
                    for label in row.children:
                        if label.text == str(i):
                            final_row=row
                try:
                    labels=final_row.children
                    labels[2].text="{:.2f}".format(fit_data["speed"])+" +/- "+"{:.2f}".format(fit_data["speed_err"])
                    labels[1].text="{:.2f}".format(fit_data["beta_app"])+" +/- "+"{:.2f}".format(fit_data["beta_app_err"])
                    labels[0].text="{:.2f}".format(fit_data["d_crit"])+" +/- "+"{:.2f}".format(fit_data["d_crit_err"])
                except:
                    pass

            self.kplot.ax.legend()
            self.kplot.ax.set_title(self.name)
            self.kplot.fig.tight_layout()
            self.kplot.ax.figure.canvas.draw_idle()
            self.kplot.ax.figure.canvas.flush_events()

            #export the plot as pdf and png
            if do_export:
                self.kplot.fig.savefig(export_path+".pdf")
                self.kplot.fig.savefig(export_path+".png",dpi=300)

    def get_redshift(self):

        #get redshift from NED, otherwise set to zero
        try:
            self.redshift = np.average(Ned.get_table(self.name, table="redshifts")["Published Redshift"])
        except:
            self.redshift = 0.00

        self.ids.redshift.text = "{:.5f}".format(self.redshift)
        self.update_redshift()

    def update_redshift(self):

        try:
            self.redshift = float(self.ids.redshift.text)
        except:
            pass
        #set redshift of components
        for plot in self.plots:
            for comp in plot.components:
                comp[1].redshift=self.redshift
                comp[1].tb=1.22e12/15**2 * comp[1].flux * (1 + self.redshift) / comp[1].maj / comp[1].min /(comp[1].scale)**2

        if len(self.components)>0:
            self.update_kinematic_plot()

    #used to save/export the kinematic results, writes a directory which includes two .csv files for the component info
    #and a sub-folder /fits with the fits files, and pdf and png files for the plots
    def save_kinematics(self,save_text,selection):
        save_path=str(selection[0])

        #set default file name if name was not specified
        if save_text==str(selection[0]):
            save_path=save_path+"/Kinematic_"+"VCAT_export.vcat"
        else:
            save_path=save_text

        #create save directories
        os.makedirs(save_path,exist_ok=True)
        os.makedirs(save_path+"/modelfit_files",exist_ok=True)
        os.makedirs(save_path+"/clean_fits",exist_ok=True)

        #copy fits files to common directory
        for file in self.modelfit_filepaths:
            os.system("cp " + file + " " + save_path + "/modelfit_filess/")

        for file in self.clean_filepaths:
            os.system("cp " + file + " " + save_path + "/clean_fits/")

        #export Table as seen on the screen with kinematic results
        export_infos=[]
        for coll in self.component_collections:
            export_infos.append(coll.get_speed())
        df = pd.DataFrame(export_infos)
        df.to_csv(save_path + '/kinematic_fit.csv',index=False)

        #export Table with component identification
        export_infos=[]
        for plot in self.plots:
            for comp in plot.components:
                export_infos.append(comp[1].get_info())
        df = pd.DataFrame(export_infos)
        df.to_csv(save_path + '/component_info.csv', index=False)

        #export plots as pdf and png

        #get selector button
        kinematic_button = self.ids.kinematic_choose
        flux_button = self.ids.flux_choose
        tb_button = self.ids.tb_choose

        k_state_initial=kinematic_button.state
        f_state_initial=flux_button.state
        t_state_initial=tb_button.state

        kinematic_button.state = 'down'
        flux_button.state = 'normal'
        tb_button.state = 'normal'
        self.update_kinematic_plot(do_export=True,export_path=save_path+"/kinematic_plot")

        kinematic_button.state = 'normal'
        flux_button.state = 'down'
        tb_button.state = 'normal'
        self.update_kinematic_plot(do_export=True, export_path=save_path + "/flux_density_plot")

        kinematic_button.state = 'normal'
        flux_button.state = 'normal'
        tb_button.state = 'down'
        self.update_kinematic_plot(do_export=True, export_path=save_path + "/tb_plot")

        #reset buttons to initial state
        kinematic_button.state = k_state_initial
        flux_button.state = f_state_initial
        tb_button.state = t_state_initial
        self.update_kinematic_plot()

        self.show_popup("Export Info","Export successful to \n" + save_path,"Continue")

    #used to reimport data that was exported with the function above. Needs a directory path
    def import_kinematics(self,directory):

        #TODO reset everything before importing stuff

        #check if it is a valid vcat kinematics folder
        if (os.path.isfile(directory[0]+"/component_info.csv") and os.path.isfile(directory[0]+"/kinematic_fit.csv") and
                os.path.exists(directory[0]+"/modelfit_files") and os.path.exists(directory[0]+"/clean_fits")):

            #import the fits files
            selection=glob.glob(directory[0]+"/modelfit_files/*")
            self.load_modelfit(selection)
            selection = glob.glob(directory[0] + "/clean_fits/*")
            self.load_modelfit(selection)
            self.create_kinematic_plots()

            #add components
            comp_info=pd.read_csv(directory[0]+"/component_info.csv")
            ncomps=len(np.unique(comp_info[comp_info["component_number"]>=0]["component_number"]))
            for i in range(ncomps):
                self.add_component()

            #set core
            core_ind=comp_info[comp_info["is_core"]==True]["component_number"].values[0]
            for t in ToggleButton.get_widgets('components'):
                if "Component " +str(core_ind) in t.text:
                    t.state="down"
                    self.set_active_component()
                    self.set_core_component()
                    t.state='normal'
            self.set_active_component()


            #now identify them with each other
            for plot in self.plots:
                for comp in plot.components:
                    #match plot component with component in list
                    filter_df=comp_info[(round(comp_info["x"],15)==round(comp[1].x,15))]
                    ind=filter_df[round(comp_info["y"],15)==round(comp[1].y,15)]["component_number"].values[0]

                    if ind>=0: #only due this if component was assigned
                        #activate correct button
                        for t in ToggleButton.get_widgets('components'):
                            if "Component " + str(ind) in t.text:
                                t.state='down'
                                self.set_active_component()
                                t.state='normal'

                        current_color = self.current_color(self.active_component_ind)
                        comp[0][0].set_color(current_color)
                        comp[0][1][0].set_color(current_color)
                        comp[0][2][0].set_color(current_color)
                        # assign component to number for kinematics
                        comp[1].assign_component_number(self.active_component_ind)

                        self.set_active_component()

            self.update_kinematic_plot()

            #set redshift
            self.ids.redshift.text="{:.5f}".format(comp_info["redshift"].values[0])
            self.update_kinematic_plot()

            message="Imported Data successfully."

        else:
            #throw error message
            message="Please select a valid .vcat kinematic folder!"

        #create popup for import info
        self.show_popup("Importing Kinematics",message,"Continue")


    #### END OF KINEMATIC FUNCTIONS

    #### START OF STACKING FUNCTIONS

    def load_stacking_files(self):

        #TODO some basic background checks on the files to see if they are valid and if they exist (write some popups)
        #TODO also check for polarizations. If there is only Stokes I input, grey out the polarization stacking options

        for ind,file in enumerate(self.clean_filepaths):
            if len(self.stokes_u_filepaths)>ind and len(self.stokes_q_filepaths)>ind:
                if len(self.modelfit_filepaths)>ind:
                    model=self.modelfit_filepaths[ind]
                else:
                    model=""
                plot_data=ImageData(file,model=model,stokes_u=self.stokes_u_filepaths[ind],stokes_q=self.stokes_q_filepaths[ind])
                image=FitsImage(plot_data,plot_mode="frac_pol",plot_evpa=True)
            else:
                #try to load model from clean .fits file
                if len(self.modelfit_filepaths)>ind:
                    model=self.modelfit_filepaths[ind]
                else:
                    model=""
                plot_data=ImageData(file,model=model)
                image=FitsImage(plot_data,plot_mode="frac_pol",plot_evpa=True)

            #check if polarization information was given:
            if np.sum(image.clean_image.stokes_q)==0 or np.sum(image.clean_image.stokes_u)==0:
                self.ids.weighted_check.disabled=True
                self.ids.stack_pol_check.disabled=True
                self.ids.stack_stokes_check.disabled=True
                self.ids.stacked_image_linpol.disabled=True
                self.ids.stacked_image_fracpol.disabled=True

            self.stacking_single_plots.append(image)
            button = ToggleButton(
                text=str(get_date(file)),
                size_hint_y=None,
                size_hint_x=0.9,
                height=50,
                group="stacking_single_plots",
            )
            button.bind(on_release=self.change_stacking_single_plot)

            #add checkbox
            checkbox=CheckBox(active=True,size_hint_x=0.1)

            self.stacking_single_plot_checkboxes.append(checkbox)
            self.stacking_single_plot_buttons.append(button)
            self.ids.stacking_single_list.add_widget(button)
            self.ids.stacking_single_list.add_widget(checkbox)

            #select first button by default
            try:
                self.stacking_single_plot_buttons[0].state="down"
                self.change_stacking_single_plot(self.stacking_single_plot_buttons[0])
            except:
                pass
        if self.ids.weighted_check.disabled:
            self.show_popup("Warning", "Only doing Stokes I stacking! \n No valid polarization data detected",
                            "Continue")
        else:
            self.show_popup("Success","Loaded data including polarization!","Continue")


    #set image input as currently displayed plot
    def change_stacking_single_plot(self,button=""):
        ind = np.where(np.array(self.stacking_single_plot_buttons)==button)[0][0]
        self.ids.stacking_single_plot.figure = self.stacking_single_plots[ind].fig


    def stack_images(self):

        #find files to stack
        files_to_stack=[]
        files_to_stack_q=[]
        files_to_stack_u=[]
        files_to_stack_uvf=[]
        files_to_stack_models=[]
        files_to_stack_casa_clean_models=[]

        do_beam_restore = self.ids.restore_beam_check.active #TODO modify this for text beam input and "get custom beam"-button
        models_loaded = False
        fold_polarization_beams=True

        for ind,box in enumerate(self.stacking_single_plot_checkboxes):
            if box.active:
                files_to_stack.append(self.clean_filepaths[ind])
                #try finding polarization
                try:
                    files_to_stack_q.append(self.stokes_q_filepaths[ind])
                    files_to_stack_u.append(self.stokes_u_filepaths[ind])
                except:
                    if len(fits.open(files_to_stack[ind])[0].data)==1:
                        fold_polarization_beams=False
                        self.show_popup("Warning","No sufficient polarization info found! \n Only staking Stokes I.",
                                    "Continue")
                #try importing models and uvf files for beam folding
                try:
                    files_to_stack_uvf.append(self.uvf_filepaths[ind])
                except:
                    if do_beam_restore:
                        self.show_popup("Warning","No uvf file loaded! \n Cannot restore to common beam",
                                    "Continue")
                        do_beam_restore = False
                try:
                    files_to_stack_models.append(self.modelfit_filepaths[ind])
                    models_loaded=True
                except:
                    pass
                #load casa clean models
                try:
                    files_to_stack_casa_clean_models.append(self.casa_clean_model_filepaths[ind])
                except:
                    pass

        #check if we need to restore the beam sizes first before proceeding
        if do_beam_restore:
            difmap_path = self.ids.difmap_path.text ##TODO fix this, so difmap executable is found automatically!
            mod_file_paths=[]
            mod_file_paths_q=[]
            mod_file_paths_u=[]

            for i in range(len(files_to_stack)):
                if models_loaded:
                    image=ImageData(files_to_stack[i],model=files_to_stack_models[i])
                else:
                    image=ImageData(files_to_stack[i])
                mod_file_paths.append("tmp/mod_files/"+image.date+".mod")

            if fold_polarization_beams:
                try:
                    #DIFMAP style
                    for file in files_to_stack_q:
                        image_q=ImageData(file,model_save_dir="tmp/mod_files_q/")
                        mod_file_paths_q.append("tmp/mod_files_q/"+image_q.date+".mod")
                    if len(files_to_stack_q)==0:
                        raise Exception()
                except:

                    # TRY to import CASA clean model
                    try:
                        for file in files_to_stack_casa_clean_models:
                            image_model=ImageData(file,model_save_dir="tmp/",is_casa_model=True)
                            mod_file_paths_q.append("tmp/mod_files_q/"+image_model.date+".mod")
                    except:
                        self.show_popup("Error","Stokes Q .fits file does not contain clean model!","Continue")
                        fold_polarization_beams=False

                try:
                    #DIFMAP style
                    for file in files_to_stack_u:
                        image_u=ImageData(file,model_save_dir="tmp/mod_files_u/")
                        mod_file_paths_u.append("tmp/mod_files_u/"+image_u.date+".mod")
                    if len(files_to_stack_u)==0:
                        raise Exception()
                except:
                    try:
                        # TRY to import CASA clean model
                        for file in files_to_stack_casa_clean_models:
                            image_model=ImageData(file,model_save_dir="tmp/",is_casa_model=True)
                            mod_file_paths_u.append("tmp/mod_files_u/"+image_model.date+".mod")
                    except:
                        self.show_popup("Error","Stokes U .fits file does not contain clean model!","Continue")
                        fold_polarization_beams=False

            degpp=image.degpp*image.scale
            npix=len(image.X)

            #Restore Stokes I
            try:
                restore_maj=float(self.ids.beam_maj.text)
                restore_min=float(self.ids.beam_min.text)
                restore_pa=float(self.ids.beam_pos.text)
            except:
                self.show_popup("Error","Please set beam to restore with!","Continue")

            fold_with_beam(files_to_stack, difmap_path=difmap_path, bmaj=restore_maj, bmin=restore_min, posa=restore_pa,
                           output_dir="tmp/restored_fits", n_pixel=npix*4,
                           pixel_size=degpp/2, mod_files=mod_file_paths,
                           uvf_files=files_to_stack_uvf)

            if fold_polarization_beams:
                #Restores Stokes Q
                fold_with_beam(files_to_stack, #input stokes I files here here to use the common beam from Stokes I!
                               difmap_path=difmap_path, bmaj=restore_maj, bmin=restore_min, posa=restore_pa,
                               channel="q", output_dir="tmp/restored_fits_q",
                               n_pixel=npix * 4,
                               pixel_size=degpp / 2, mod_files=mod_file_paths_q,
                               uvf_files=files_to_stack_uvf)
                for i in range(len(files_to_stack)):
                    files_to_stack_q.append("tmp/restored_fits_q/"+".".join(
                        files_to_stack[i].split("/")[-1].split(".")[0:-1]) + "_convolved.fits")

                #Restore Stokes U
                fold_with_beam(files_to_stack, #input stokes I files here here to use the common beam from Stokes I!
                               difmap_path=difmap_path, bmaj=restore_maj, bmin=restore_min, posa=restore_pa,
                               channel="u", output_dir="tmp/restored_fits_u",
                               n_pixel=npix * 4,
                               pixel_size=degpp / 2, use_common_beam=True, mod_files=mod_file_paths_u,
                               uvf_files=files_to_stack_uvf)
                for i in range(len(files_to_stack)):
                    files_to_stack_u.append("tmp/restored_fits_u/"+".".join(
                        files_to_stack[i].split("/")[-1].split(".")[0:-1]) + "_convolved.fits")

            #change file names:
            for i in range(len(files_to_stack)):
                files_to_stack[i]="tmp/restored_fits/"+'.'.join(files_to_stack[i].split("/")[-1].split(".")[0:-1])+"_convolved.fits"



        align=self.ids.do_alignment_check.active
        weighted=self.ids.weighted_check.active
        if self.ids.stack_stokes_check.active:
            output_stacked = stack_fits(files_to_stack,align=align,stokes_q_fits=files_to_stack_q,stokes_u_fits=files_to_stack_u)
            if len(output_stacked) == 1:
                stack_image = ImageData(files_to_stack[0],pol_from_stokes=True)
                stack_image.Z = output_stacked[0][0]
            elif len(output_stacked) > 1: #check for polarization
                stack_image = ImageData(files_to_stack[0], stokes_i=output_stacked[0][0],
                                        pol_from_stokes=True,stokes_q=output_stacked[1][0],stokes_u=output_stacked[2][0])
                stack_image.Z = output_stacked[0][0]
                stack_image.stokes_q = output_stacked[1][0]
                stack_image.stokes_u = output_stacked[2][0]
                stack_image.lin_pol = np.sqrt(output_stacked[1][0]**2+output_stacked[2][0]**2)
                stack_image.evpa=0.5*np.arctan2(output_stacked[2][0],output_stacked[1][0])

        elif self.ids.stack_pol_check.active:
            output_stacked=stack_pol_fits(files_to_stack,weighted=weighted,align=align,stokes_u_fits=files_to_stack_u,stokes_q_fits=files_to_stack_q)
            #create new image data with header info from first fits file
            stack_image = ImageData(files_to_stack[0],pol_from_stokes=False)
            stack_image.Z = output_stacked[0][0]
            if len(output_stacked) > 1:  # check for polarization
                stack_image.lin_pol = output_stacked[1][0]
                stack_image.evpa = output_stacked[2][0]
        #create plot
        self.final_stack_image=stack_image
        self.ids.stacked_image_i.state="down"
        StackPlot = FitsImage(stack_image,title="Stacked Image")
        self.ids.stacked_image.figure = StackPlot.fig

    def change_final_stack_plot(self,mode,button):
        if self.final_stack_image!="":
            if mode in ["lin_pol","frac_pol"]:
                plot_evpa=True
            else:
                plot_evpa=False
            StackPlot = FitsImage(self.final_stack_image,plot_mode=mode,plot_evpa=plot_evpa,title="Stacked Image")
            self.ids.stacked_image.figure = StackPlot.fig
            button.state="down"

    def fill_in_common_beam(self):
        files_to_stack=[]
        for ind,box in enumerate(self.stacking_single_plot_checkboxes):
            if box.active:
                files_to_stack.append(self.clean_filepaths[ind])
        if len(files_to_stack)>0:
            maj,min,pos=get_common_beam(files_to_stack)
            self.ids.beam_maj.text="{:.4f}".format(maj)
            self.ids.beam_min.text = "{:.4f}".format(min)
            self.ids.beam_pos.text = "{:.4f}".format(pos)
        else:
            self.show_popup("Error","Load stacking files first!","Continue")

    ### BEGIN PLOT VIEW

    def load_plotting_files(self):
    ### TODO disable pol buttons if there is no polarization in input
        for ind,clean_file in enumerate(self.clean_filepaths):

            clean_path=self.clean_filepaths[ind]

            try:
                model_path=self.modelfit_filepaths[ind]
            except:
                model_path=""
            try:
                stokes_u_path=self.stokes_u_filepaths[ind]
            except:
                stokes_u_path= ""
            try:
                stokes_q_path = self.stokes_q_filepaths[ind]
            except:
                stokes_q_path = ""

            try:# import CASA STYLE models
                ImageData(self.casa_clean_model_filepaths[ind],model_save_dir="tmp/",is_casa_model=True)
            except:
                pass
            plot_data=ImageData(clean_path,model=model_path,stokes_q=stokes_q_path,stokes_u=stokes_u_path)
            self.plotting_single_plots_data.append(plot_data)
            plot=FitsImage(plot_data)
            self.plotting_single_plots.append(plot)

            button = ToggleButton(
                text=str(plot_data.date),
                size_hint_y=None,
                size_hint_x=0.9,
                height=50,
                group="plotting_single_plots",
            )
            button.bind(on_release=self.change_plotting_single_plot)

            # add checkbox
            checkbox = CheckBox(active=True, size_hint_x=0.1)

            self.plotting_single_plot_checkboxes.append(checkbox)
            self.plotting_single_plot_buttons.append(button)
            self.ids.plotting_single_list.add_widget(button)
            self.ids.plotting_single_list.add_widget(checkbox)

    def change_plotting_single_plot(self,button=""):
        ind = np.where(np.array(self.plotting_single_plot_buttons)==button)[0][0]
        self.ids.plotting_single_plot.figure = self.plotting_single_plots[ind].fig

        #Update Source information on right side
        image=self.plotting_single_plots[ind].clean_image
        self.ids.source_name.text = str(image.name)
        self.ids.frequency.text = "{:.2f}".format(image.freq*1e-9) + " GHz"
        self.ids.obs_date.text = str(image.date)
        self.ids.integrated_flux_image.text = "{:.2f}".format(image.integrated_flux_image*1000) +  " mJy"
        self.ids.integrated_flux_clean.text = "{:.2f}".format(image.integrated_flux_clean*1000) +  " mJy"
        self.ids.peak_flux.text = "{:.2f}".format(np.max(image.Z*1000)) + " mJy/beam"
        self.ids.integrated_pol_flux_image.text = "{:.2f}".format(image.integrated_pol_flux_image * 1000) + " mJy"
        self.ids.integrated_pol_flux_clean.text = "{:.2f}".format(image.integrated_pol_flux_clean * 1000) + " mJy"
        self.ids.pol_peak_flux.text = "{:.2f}".format(np.max(image.lin_pol * 1000)) + " mJy/beam"
        self.ids.stokes_i_noise.text = "{:.2f}".format(image.noise*1000) + " mJy/beam"
        self.ids.pol_noise.text = "{:.2f}".format(image.pol_noise*1000) + " mJy/beam"
        self.ids.evpa.text = "{:.2f}".format(image.evpa_average/np.pi*180) + "°"
        self.ids.fractional_polarization.text = (
                "{:.2f}".format(image.integrated_pol_flux_clean/image.integrated_flux_clean*100)+"%")

    #TODO IMPLEMENT THIS!
    def replot(self):

        #get parameters from input fields
        active_button = next((t for t in ToggleButton.get_widgets('plot_mode_select') if t.state == 'down'), None)
        try:
            if active_button.text == "Lin. Pol.":
                plot_mode = "lin_pol"
            elif active_button.text == "Frac. Pol.":
                plot_mode = "frac_pol"
            else:
                plot_mode = "stokes_i"
        except:
            plot_mode = "stokes_i"
            self.ids.plot_mode_i.state = "down"
        try:
            xlim = [float(self.ids.ra_min.text), float(self.ids.ra_max.text)]
        except:
            xlim = []
        try:
            ylim = [float(self.ids.dec_min.text), float(self.ids.dec_max.text)]
        except:
            ylim = []

        for ind, data in enumerate(self.plotting_single_plots_data):
            self.plotting_single_plots[ind] = FitsImage(data,
                                                        plot_mode=plot_mode,
                                                        stokes_i_sigma_cut=float(self.ids.stokes_i_sigma_cut.text),
                                                        im_colormap = self.ids.im_colormap.active,
                                                        contour = self.ids.contour.active,
                                                        contour_color = self.ids.contour_color.text,
                                                        contour_cmap = self.ids.contour_cmap.text,
                                                        contour_alpha = float(self.ids.contour_alpha.text),
                                                        contour_width = float(self.ids.contour_width.text),
                                                        im_color = self.ids.im_color.text,
                                                        plot_beam = self.ids.plot_beam.active,
                                                        overplot_gauss = self.ids.overplot_gauss.active,
                                                        overplot_clean = self.ids.overplot_clean.active,
                                                        component_color = self.ids.component_color.text,
                                                        xlim = xlim,
                                                        ylim = ylim,
                                                        plot_evpa = self.ids.plot_evpa.active,
                                                        evpa_width = float(self.ids.evpa_width.text),
                                                        evpa_len = float(self.ids.evpa_len.text),
                                                        lin_pol_sigma_cut = float(self.ids.lin_pol_sigma_cut.text),
                                                        evpa_distance = float(self.ids.evpa_distance.text),
                                                        rotate_evpa = float(self.ids.rotate_evpa.text),
                                                        evpa_color = self.ids.evpa_color.text,
                                                        title = self.ids.title.text,
                                                        rcparams = self.ids.rcparams.text)

        #find active button
        active_button = next((t for t in ToggleButton.get_widgets('plotting_single_plots') if t.state == 'down'), None)
        self.change_plotting_single_plot(button=active_button)



class VCAT(App):

    def build(self):
        self.screen=ModelFits()
        return self.screen

    #### START OF KINEMATIC FUNCTIONS

    def set_touch_mode(self,touch_mode):
        if touch_mode == "zoombox":
            for figure in self.screen.figure_widgets:
                figure.touch_mode = "pan"
                self.screen.ids.figure_scroll.do_scroll_x = False
                #self.screen.ids.figure_scroll.update_canvas()
        if touch_mode == "pointselect":
            for figure in self.screen.figure_widgets:
                figure.touch_mode = "pointselect"
                self.screen.ids.figure_scroll.do_scroll_x = True

    def home_plot(self):
        for figure in self.screen.figure_widgets:
            figure.home()

    def add_component(self):
        self.screen.add_component()

    def remove_component(self):
        self.screen.remove_component()

    def choose_component(self,figure_widget,x,y):
        self.screen.choose_component(figure_widget,x,y)

    def set_core_component(self):
        self.screen.set_core_component()

    def get_redshift(self):
        self.screen.get_redshift()

    def update_redshift(self):
        self.screen.update_redshift()

    def update_kinematic_plot(self):
        self.screen.update_kinematic_plot()

    def save_kinematics(self,save_text,selection):
        self.screen.save_kinematics(save_text,selection)

    def import_kinematics(self,directory):
        self.screen.import_kinematics(directory)

    #### END OF KINEMATIC FUNCTIONS

    #### START OF STACKING FUNCTIONS

if __name__ == "__main__":
    VCAT().run()