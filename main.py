from kivy.app import App
from kivy.uix.popup import Popup
from kivy.uix.togglebutton import ToggleButton
from kivy.uix.tabbedpanel import TabbedPanel
from kivy.uix.boxlayout import BoxLayout
from kivy.uix.label import Label
from kivy.properties import ObjectProperty, StringProperty, ListProperty
from kivy.uix.button import Button
from kivy.uix.checkbox import CheckBox
from graph_widget import MatplotFigure
from kivy.utils import platform
import numpy as np
from astroquery.ipac.ned import Ned
from astropy.io import fits
from astropy.time import Time
import os
import pandas as pd
import glob
from matplotlib.patches import Ellipse
from vcat.stacking_helpers import stack_fits, stack_pol_fits, fold_with_beam, modelfit_ehtim_pol, modelfit_ehtim_full_pol, modelfit_difmap
from vcat.kinematics import ComponentCollection, Component
from vcat import FitsImage, ImageData, KinematicPlot, ImageCube
from vcat.helpers import get_date, get_common_beam, write_mod_file_from_components, get_residual_map, Jy2JyPerBeam, JyPerBeam2Jy
from mojave_db_access import upload_csv_to_MOJAVE, download_kinematic_from_MOJAVE, query_models, check_password


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

class ExportMoviePopup(Popup):
    load = ObjectProperty()

    def load_filepath_to_view(self, selection):
        if selection!=[]:
            self.ids.export_movie_save.text = str(selection[0])

class ModelfitExportPopup(Popup):
    load = ObjectProperty()

    def load_filepath_to_view(self, selection):
        if selection!=[]:
            self.ids.export_modelfit_save.text = str(selection[0])

class PlotExportPopup(Popup):
    load = ObjectProperty()

    def load_filepath_to_view(self, selection):
        if selection!=[]:
            self.ids.export_plot_save.text = str(selection[0])

class MOJAVEExportPopup(Popup):
    load = ObjectProperty()

class MOJAVEImportPopup(Popup):
    load = ObjectProperty()

class MOJAVEPasswordPopup(Popup):
    load = ObjectProperty()

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
    component_colors=["#ef476f","#ffd166","#06d6a0","#118ab2","#073b4c"]
    figure_widgets=[]
    core_component_ind=None
    name=""
    component_collections = []
    mojave_password=""

    #### MODELFIT variables
    modelfit_plots = []
    modelfit_plot_buttons = []
    current_modelfit_residual_plot = []
    new_modelfit_component_coords = []
    current_new_ellipse_plot = []

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
    noise_method=""

    #### GENERAL FUNCTIONS ACROSS TABS

    def set_mojave_password(self,password):
        #TODO implement a check if the password is correct, if not, open MOJAVEPasswordPopup again.
        if check_password(password):
            self.mojave_password=password
        else:
            self.mojave_password_dialog = MOJAVEPasswordPopup()
            self.mojave_password_dialog.open()

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

        try:
            if len(fits_files)>0:
                date=[]

                for filepath in fits_files:
                    date=np.append(date,get_date(filepath))

                args = date.argsort()
                fits_files = fits_files[args]
        except:
            pass

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

        #try attaching the new files to previous ones
        combined = selection + self.uvf_filepaths
        # sort by date
        self.uvf_filepaths = self.sort_uvf_by_date(combined)

        try:
            self.the_popup.dismiss()
        except:
            pass
        self.file_info = str(len(self.uvf_filepaths)) + " Files selected"
        self.ids.get_file_uvf.text = self.file_info

    def load_modelfit(self, selection):

        # try attaching the new files to previous ones
        combined = selection + self.modelfit_filepaths

        # sort by date
        try:
            self.modelfit_filepaths = self.sort_fits_by_date(combined)
        except:
            pass
        try:
            self.the_popup.dismiss()
        except:
            pass
        self.file_info=str(len(self.modelfit_filepaths))+" Files selected"
        self.ids.get_file_modelfit.text = self.file_info

    def load_clean(self, selection):

        # try attaching the new files to previous ones
        combined = selection + self.clean_filepaths

        self.clean_filepaths = self.sort_fits_by_date(combined)
        try:
            self.the_popup.dismiss()
        except:
            pass
        self.file_info=str(len(self.clean_filepaths))+" Files selected"
        self.ids.get_file_clean.text = self.file_info

    def load_stokes_q(self, selection):

        # try attaching the new files to previous ones
        combined = selection + self.stokes_q_filepaths

        self.stokes_q_filepaths = self.sort_fits_by_date(combined)
        try:
            self.the_popup.dismiss()
        except:
            pass
        self.file_info=str(len(self.stokes_q_filepaths))+" Files selected"
        self.ids.get_file_stokes_q.text = self.file_info

    def load_stokes_u(self, selection):

        # try attaching the new files to previous ones
        combined = selection + self.stokes_u_filepaths

        self.stokes_u_filepaths = self.sort_fits_by_date(combined)
        try:
            self.the_popup.dismiss()
        except:
            pass
        self.file_info=str(len(self.stokes_u_filepaths))+" Files selected"
        self.ids.get_file_stokes_u.text = self.file_info

    def load_casa_clean_model(self,selection):

        # try attaching the new files to previous ones
        combined = selection + self.casa_clean_model_filepaths

        self.casa_clean_model_filepaths = self.sort_fits_by_date(combined)
        try:
            self.the_popup.dismiss()
        except:
            pass
        self.file_info = str(len(self.casa_clean_model_filepaths)) + " Files selected"
        self.ids.get_file_casa_clean_model.text = self.file_info

    def reset_input_files(self):

        #reset file arrays
        self.casa_clean_model_filepaths=[]
        self.stokes_u_filepaths=[]
        self.stokes_q_filepaths=[]
        self.clean_filepaths=[]
        self.modelfit_filepaths=[]
        self.uvf_filepaths=[]

        #reset text in field
        self.ids.get_file_casa_clean_model.text=""
        self.ids.get_file_stokes_u.text=""
        self.ids.get_file_stokes_q.text=""
        self.ids.get_file_clean.text=""
        self.ids.get_file_modelfit.text=""
        self.ids.get_file_uvf.text=""

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

    def open_movie_dialog(self):
        self.movie_dialog = ExportMoviePopup()
        self.movie_dialog.open()

    def open_import_dialog(self):
        self.import_dialog = FileImportPopup()
        self.import_dialog.open()

    def open_mojave_import_dialog(self):
        self.the_popup = MOJAVEImportPopup()
        self.the_popup.open()

        if self.mojave_password=="":
            self.mojave_password_dialog = MOJAVEPasswordPopup()
            self.mojave_password_dialog.open()

    def search_mojave_models(self,sourcename):

        models=query_models(sourcename,self.mojave_password)
        if len(models)>0:
            output_text="The following models are available in the MOJAVE DB:"
            for mod in models:
                output_text+="\n"+mod
        else:
            output_text="No models found in the MOJAVE DB"
        self.show_popup("Available MOJAVE models",output_text,"Continue")


    def current_color(self,ind):
        if ind is not None:
            ind=int(ind % len(self.component_colors))
            return self.component_colors[ind]
        else:
            return None

    def create_kinematic_plots(self):

        fits_images=[]
        sort_inds=[]
        #fits images to plot
        modelfit_files_to_plot=[]
        clean_files_to_plot=[]
        uvf_files_to_plot=[]

        #get already used modelfit files:
        already_used_fits=[]
        for plot in self.plots:
            already_used_fits.append(plot.model_image_file)
        for ind, file_path in enumerate(self.modelfit_filepaths):
            if file_path in already_used_fits:
                #in this case the file is already plotted
                pass
            else:
                #in this case we need to add the plot for the file
                modelfit_files_to_plot.append(file_path)
                sort_inds.append(ind)

                #try also attaching the clean files if they exist
                try:
                    clean_files_to_plot.append(self.clean_filepaths[ind])
                except:
                    pass

                # try also attaching the uvf files if they exist
                try:
                    uvf_files_to_plot.append(self.uvf_filepaths[ind])
                except:
                    pass

        if len(clean_files_to_plot)!=len(modelfit_files_to_plot) and len(clean_files_to_plot)>0:
            self.show_popup("Warning","Please use an equal number of modelfit and clean images","Continue")

        elif len(clean_files_to_plot) == 0 and len(modelfit_files_to_plot) == 0:
            self.show_popup("Warning", "No data selected", "Continue")
        #check if clean maps AND modelfits were provided
        elif len(clean_files_to_plot) == len(modelfit_files_to_plot):
            warn_uvf=False
            #create plots for view page
            for ind,filepath in enumerate(modelfit_files_to_plot):
                if not len(uvf_files_to_plot) == len(modelfit_files_to_plot):
                    plot_data=ImageData(clean_files_to_plot[ind],model=filepath,noise_method=self.noise_method,is_ehtim_model=self.ids.is_ehtim_model.active)
                    warn_uvf=True
                else:
                    plot_data = ImageData(clean_files_to_plot[ind], model=filepath, uvf_file=uvf_files_to_plot[ind],
                                          difmap_path=self.ids.difmap_path.text,noise_method=self.noise_method,is_ehtim_model=self.ids.is_ehtim_model.active)
                plot=FitsImage(plot_data,plot_model=True,adjust_comp_size_to_res_lim=True)
                fits_images=np.append(fits_images,plot)
            self.show_popup("Information", "File loading completed. Have fun doing kinematics!", "Continue")
            if warn_uvf:
                self.show_popup("Warning",
                                "No .uvf files loaded, will not be able to calculate good upper limits for TB!",
                                "Continue")

        #otherwise check if no clean maps were provided and only modelfit maps
        elif len(clean_files_to_plot) == 0 and len(modelfit_files_to_plot)>0:
            # create plots for view page
            warn_uvf=False
            for ind,filepath in enumerate(modelfit_files_to_plot):
                if not len(uvf_files_to_plot) == len(modelfit_files_to_plot):
                    plot_data = ImageData(filepath,model=filepath,noise_method=self.noise_method,is_ehtim_model=self.ids.is_ehtim_model.active)
                    warn_uvf=True
                else:
                    plot_data = ImageData(filepath, model=filepath, uvf_file=uvf_files_to_plot[ind],
                                          difmap_path=self.ids.difmap_path.text,noise_method=self.noise_method,is_ehtim_model=self.ids.is_ehtim_model.active)
                plot = FitsImage(plot_data,plot_model=True,adjust_comp_size_to_res_lim=True)
                fits_images = np.append(fits_images, plot)
            self.show_popup("Warning",
                            "No clean images imported, using only modelfit images.\n Have fun doing kinematics!",
                            "Continue")
            if warn_uvf:
                self.show_popup("Warning",
                                "No .uvf files loaded, will not be able to calculate upper limits for TB!",
                                "Continue")
        elif len(clean_files_to_plot) > 0 and len(modelfit_files_to_plot) ==0:
            self.show_popup("Warning","No modelfits imported","Continue")

        #remove all existing figures to later add them again together with the new ones
        for figure in self.figure_widgets:
            self.ids.figures.remove_widget(figure)

        frequencies_old = []
        #find frequencies that are already plotted
        for plot in self.plots:
            frequency="{:.0f}".format(plot.freq/1e9)+" GHz"
            if frequency not in frequencies_old:
                frequencies_old.append(frequency)

        frequencies = []
        for ind,plot in enumerate(fits_images):
            frequency="{:.0f}".format(plot.freq/1e9)+" GHz"
            if not frequency in frequencies and not frequency in frequencies_old:
                frequencies.append(frequency)
            self.plots.insert(sort_inds[ind],plot)
            self.name=plot.name
            new_figure = MatplotFigure(size_hint_x=None,
                                       width=100)
            new_figure.figure = plot.fig
            new_figure.touch_mode = "pan"
            self.figure_widgets.insert(sort_inds[ind],new_figure)

        #add widgets again together with new ones
        for new_figure in self.figure_widgets:
            self.ids.figures.add_widget(new_figure)

        #create frequency toggle switches if multiple frequencies were imported
        for ind,frequency in enumerate(frequencies):

            #add buttons for the final kinematic plot to switch between frequencies
            toggle_state="down" if ind==0 else "normal"
            button = ToggleButton(
                text=frequency,
                size_hint_y=1,
                group="plot_frequency_select",
                state=toggle_state
            )
            button.bind(on_release=self.update_kinematic_plot)
            self.ids.plot_frequency_select_buttons.add_widget(button)

            #add buttons for the scroll view to hide/show specific frequencies
            button = ToggleButton(
                text=frequency,
                size_hint_y=1,
                state="down"
            )
            button.bind(on_release=self.update_scroll_plots)
            self.ids.scroll_frequency_select_buttons.add_widget(button)

    def refresh_plot(self):

        try:
            ind=int(self.ids.refresh_plot_ind.text)

            #create new ImageData object
            try:
                plot_data=ImageData(self.clean_filepaths[ind], model=self.modelfit_filepaths[ind], uvf_file=self.uvf_filepaths[ind],
                          difmap_path=self.ids.difmap_path.text,noise_method=self.noise_method,is_ehtim_model=self.ids.is_ehtim_model.active)
            except:
                try:
                    plot_data=ImageData(self.modelfit_filepaths[ind],model=self.modelfit_filepaths[ind],
                                        uvf_file=self.uvf_filepaths[ind],noise_method=self.noise_method,is_ehtim_model=self.ids.is_ehtim_model.active)
                except:
                    try:
                        plot_data=ImageData(self.clean_filepaths[ind],model=self.modelfit_filepaths[ind],
                                            noise_method=self.noise_method,is_ehtim_model=self.ids.is_ehtim_model.active)
                    except:
                        plot_data=ImageData(self.modelfit_filepaths[ind],model=self.modelfit_filepaths[ind],
                                            noise_method=self.noise_method,is_ehtim_model=self.ids.is_ehtim_model.active)

            plot = FitsImage(plot_data, plot_model=True, adjust_comp_size_to_res_lim=True)
            #create Figure
            new_figure = MatplotFigure(size_hint_x=None,
                                       width=100)
            new_figure.figure = plot.fig
            new_figure.touch_mode = "pan"

            self.plots[ind]=plot

            # remove all plots from scroll view
            for widget in self.figure_widgets:
                self.ids.figures.remove_widget(widget)

            #add new widget
            self.figure_widgets[ind]= new_figure

            #plot all of them again
            for widget in self.figure_widgets:
                self.ids.figures.add_widget(widget)

            self.update_kinematic_plot()
        except:
            self.show_popup("Warning", "Please use an existing image number!", "Continue")

    #takes care of showing/hiding plots depending on the frequency selected.
    def update_scroll_plots(self,button_dummy):
        #remove all plots from scroll view
        for widget in self.figure_widgets:
            self.ids.figures.remove_widget(widget)

        #add back only the ones that are selected
        for ind,plot in enumerate(self.plots):
            for button in self.ids.scroll_frequency_select_buttons.children:
                if "{:.0f}".format(plot.freq/1e9)+" GHz" == button.text:
                    if button.state=="down":
                        self.ids.figures.add_widget(self.figure_widgets[ind])
                    else:
                        self.ids.figures.remove_widget(self.figure_widgets[ind])

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



    def add_component(self,number=-1):
        if number == -1:
            count = len(self.components)
            while "Component " + str(count) in self.components:
                count += 1
        else:
            count=number

        self.components.append("Component " + str(count))

        component_box=BoxLayout(
            size_hint_y=None,
            height=50
        )

        button=ToggleButton(
            text=self.components[-1],
            size_hint_y=1,
            size_hint_x=0.9,
            group="components",
        )
        button.bind(on_release=self.set_active_component)

        checkbox=CheckBox(active=False,size_hint_x=0.1)
        checkbox.bind(on_release=self.update_kinematic_plot)

        component_box.add_widget(button)
        component_box.add_widget(checkbox)
        self.ids.component_list.add_widget(component_box)


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
            #remove component button
            self.ids.component_list.remove_widget(button_to_remove.parent)

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
                for obj in button.parent.children:
                    obj.disabled=False

        #Add "Core" to button text
        active_button = next((t for t in ToggleButton.get_widgets('components') if t.state == 'down'), None)
        if active_button:
            active_button.text=active_button.text+" (Core)"
            for obj in active_button.parent.children:
                if obj != active_button:
                    obj.active = False
                    obj.disabled=True


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

        # find out which plot is currently selected:
        active_button = next((t for t in ToggleButton.get_widgets('kinematic_select') if t.state == 'down'), None)

        plots_to_fit=self.plots

        #check for chi-square plot
        #TODO Make it separate for every frequency!!
        if active_button != None and active_button.text == "Chi-Square":
            self.kplot.plot_chi_square(self.uvf_filepaths, self.modelfit_filepaths, self.ids.difmap_path.text)


        #calculate Core Shift based on modelfits TODO update this with new capability in VCAT
        if active_button != None and active_button.text == "Core Shift" and len(plots_to_fit)>0:
            t_max=0
            t_min=plots_to_fit[0].components[0][1].year

            #first find available frequencies
            frequencies=[]
            for plot in plots_to_fit:
                plot_freq="{:.0f}".format(plot.freq/1e9)
                if plot_freq not in frequencies:
                    frequencies.append(plot_freq)

            #get mulfitfreq component lists for core shift and do kinematic fit
            final_mf_fits = []
            for i in range(len(self.components)):
                collections=[]
                for frequency_to_plot in frequencies:
                    collection = []
                    for plot in plots_to_fit:
                        if "{:.0f}".format(plot.freq/1e9) == frequency_to_plot:
                            for comp in plot.components:
                                if comp[1].component_number == i:
                                    # determine plot limits
                                    if comp[1].year > t_max:
                                        t_max = comp[1].year
                                    if comp[1].year < t_min:
                                        t_min = comp[1].year
                                    collection.append(comp[1])
                    fit_data=ComponentCollection(collection,name=self.components[i]).get_speed()
                    collections.append(fit_data)
                final_mf_fits.append(collections)

            values = []

            #calculate core shift for selected components and plot it
            for i in range(len(final_mf_fits)):

                #check the checkboxes which components to use:
                use=False
                for child in self.ids.component_list.children:
                    use1=False
                    use2=False
                    for child2 in child.children:
                        if not use:
                            if isinstance(child2,CheckBox):
                                use1 = child2.active
                            if isinstance(child2,ToggleButton):
                                if "Component "+str(i) in child2.text:
                                    use2 = True
                            use = use1 and use2

                #cast list to float
                frequencies=[float(frequency) for frequency in frequencies]


                for j in range(len(frequencies)):
                    if use and j!=np.argmax(frequencies):
                        component_fits=final_mf_fits[i]
                        highest_freq_fits=component_fits[np.argmax(frequencies)]
                        current_freq_fits=component_fits[j]

                        #only do it if the fit worked
                        if not ((highest_freq_fits["speed"]==0 and highest_freq_fits["y0"]==0)
                                or (current_freq_fits["speed"]==0 and current_freq_fits["y0"]==0)):

                            #calculate core shift evolution
                            new_values=(current_freq_fits["speed"] - highest_freq_fits["speed"]) * np.array([t_max,t_min]) + (
                                        current_freq_fits["y0"] - highest_freq_fits["y0"])
                            #get extreme values to calculate plot lims
                            values=np.append(values,new_values)

                            #plot it
                            self.kplot.plot_linear_fit(t_min - 0.1 * (t_max - t_min), t_max + 0.1 * (t_max - t_min),
                                                       current_freq_fits["speed"]-highest_freq_fits["speed"],
                                                       current_freq_fits["y0"]-highest_freq_fits["y0"], self.current_color(i),
                                                       label=self.components[i]+", "+str(np.max(frequencies))+"-"+str(frequencies[j])+" GHz")
            #set plot lims
            if len(values)>0:
                plot_min_y = np.min(values)
                plot_max_y = np.max(values)
                self.kplot.set_limits([t_min - 0.1 * (t_max - t_min), t_max + 0.1 * (t_max - t_min)],
                    [plot_min_y-0.2,plot_max_y+0.2])


        # filter out which frequency to plot
        frequency_button = next((t for t in ToggleButton.get_widgets('plot_frequency_select') if t.state == 'down'),
                                None)
        if (len(ToggleButton.get_widgets('plot_frequency_select')) > 0):
            if (frequency_button == None):
                frequency_button = ToggleButton.get_widgets('plot_frequency_select')[0]
                frequency_button.state = "down"
            frequency_to_plot = frequency_button.text
            plots_to_fit = []
            for plot in self.plots:
                if "{:.0f}".format(plot.freq / 1e9) + " GHz" == frequency_to_plot:
                    plots_to_fit.append(plot)
        else:
            plots_to_fit = self.plots

        t_max=0
        if len(plots_to_fit)>0:
            t_min=plots_to_fit[0].components[0][1].year
            d_max=0
            tb_max=0
            tb_min=plots_to_fit[0].components[0][1].tb
            flux_max=0
            for i in range(len(self.components)):
                collection=[]
                for plot in plots_to_fit:
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
                            if comp[1].tb > tb_max and comp[1].tb != float('inf'):
                                tb_max = comp[1].tb
                            if comp[1].tb < tb_min:
                                tb_min = comp[1].tb
                            if comp[1].flux > flux_max:
                                flux_max = comp[1].flux

                comp_collection=ComponentCollection(collection,name=self.components[i])
                self.component_collections.append(comp_collection)

                if active_button != None and active_button.text == "Flux":
                    self.kplot.plot_fluxs(comp_collection,self.current_color(i))
                    self.kplot.set_limits([t_min - 0.1 * (t_max - t_min), t_max + 0.1 * (t_max - t_min)],
                                          [0, 1.2 * flux_max])
                elif active_button != None and active_button.text == "TB":
                    self.kplot.plot_tbs(comp_collection,self.current_color(i))
                    self.kplot.set_limits([t_min - 0.1 * (t_max - t_min), t_max + 0.1 * (t_max - t_min)],
                                          [tb_min, 10 * tb_max])
                elif active_button != None and active_button.text == "Kinematic":
                    self.kplot.plot_kinematics(comp_collection,self.current_color(i))
                    self.kplot.set_limits([t_min - 0.1 * (t_max - t_min), t_max + 0.1 * (t_max - t_min)],
                                          [0, 1.2 * d_max])
                elif active_button != None and active_button.text == "PA":
                    self.kplot.plot_pas(comp_collection, self.current_color(i))
                    self.kplot.set_limits([t_min - 0.1 * (t_max - t_min), t_max + 0.1 * (t_max - t_min)],
                                          [0, 1.2 * d_max])

                if (active_button != None) and (len(collection)>2) and (active_button.text == "Kinematic"):
                    fit_data = comp_collection.get_speed()[0]
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
                        labels[2].text="{:.2f}".format(fit_data["speed"])+" +/- "+"{:.2f}".format(abs(fit_data["speed_err"]))
                        labels[1].text="{:.2f}".format(fit_data["beta_app"])+" +/- "+"{:.2f}".format(abs(fit_data["beta_app_err"]))
                        labels[0].text="{:.2f}".format(fit_data["d_crit"])+" +/- "+"{:.2f}".format(abs(fit_data["d_crit_err"]))
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
        os.makedirs(save_path+"/uvf_files",exist_ok=True)

        #copy fits files to common directory
        for file in self.modelfit_filepaths:
            os.system("cp " + file + " " + save_path + "/modelfit_files/")

        for file in self.clean_filepaths:
            os.system("cp " + file + " " + save_path + "/clean_fits/")

        for file in self.uvf_filepaths:
            os.system("cp " + file + " " + save_path + "/uvf_files/")

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
        core_shift_choose = self.ids.core_shift_choose

        k_state_initial=kinematic_button.state
        f_state_initial=flux_button.state
        t_state_initial=tb_button.state
        core_shift_state_initial=core_shift_choose.state

        for frequency_btn in ToggleButton.get_widgets('plot_frequency_select'):

            frequency_btn.state="down"
            for frequency_btn2 in ToggleButton.get_widgets('plot_frequency_select'):
                if frequency_btn!=frequency_btn2:
                    frequency_btn2.state="normal"

            freq=frequency_btn.text.replace(" ","")
            kinematic_button.state = 'down'
            flux_button.state = 'normal'
            tb_button.state = 'normal'
            core_shift_choose.state = 'normal'
            self.update_kinematic_plot(do_export=True,export_path=save_path+"/kinematic_plot_"+freq)

            kinematic_button.state = 'normal'
            flux_button.state = 'down'
            tb_button.state = 'normal'
            core_shift_choose.state = 'normal'
            self.update_kinematic_plot(do_export=True, export_path=save_path + "/flux_density_plot_"+freq)

            kinematic_button.state = 'normal'
            flux_button.state = 'normal'
            tb_button.state = 'down'
            core_shift_choose.state = 'normal'
            self.update_kinematic_plot(do_export=True, export_path=save_path + "/tb_plot_"+freq)

            frequency_btn.state="normal"

        kinematic_button.state = 'normal'
        flux_button.state = 'normal'
        tb_button.state = 'normal'
        core_shift_choose.state = 'down'
        self.update_kinematic_plot(do_export=True, export_path=save_path + "/core_shift")

        #reset buttons to initial state
        kinematic_button.state = k_state_initial
        flux_button.state = f_state_initial
        tb_button.state = t_state_initial
        core_shift_choose.state = core_shift_state_initial
        self.update_kinematic_plot()

        self.show_popup("Export Info","Export successful to \n" + save_path,"Continue")

        #Optionally export the fits to MOJAVE database
        self.component_info_csv=save_path + '/component_info.csv'
        self.show_mojave_popup()

    #used to export movie from kinematics
    def save_movie(self,popup,save_text,selection):

        save_path = str(selection[0])

        # set default file name if name was not specified
        if save_text == str(selection[0]):
            save_path = save_path + "/movie.mp4"
        else:
            save_path = save_text

        #create ImageCube
        images=[]
        for plot in self.plots:
            image_data=plot.clean_image
            comps=[]
            for comp in plot.components:
                comps.append(comp[1])
            image_data.components=comps
            images.append(image_data)

        im_cube=ImageCube(images)

        #extract plot settings from popup
        if popup.ids.movie_plot_mode=="Stokes I":
            plot_mode="stokes_i"
        elif popup.ids.movie_plot_mode=="Lin. Pol.":
            plot_mode="lin_pol"
        elif popup.ids.movie_plot_mode=="Frac. Pol.":
            plot_mode="frac_pol"
        else:
            plot_mode="stokes_i"
        #TODO potentially add spix....
        duration = float(popup.ids.movie_duration.text)
        fps = float(popup.ids.movie_fps.text)

        interval=int(1000/fps)
        n_frames=int(fps*duration)

        args={
            "n_frames": n_frames,
            "interval": interval,
            "fps": fps,
            "freq": popup.ids.movie_freq.text,
            "plot_mode": plot_mode,
            "plot_components": popup.ids.movie_overplot_gauss.active,
            "im_colormap": popup.ids.movie_im_colormap.active,
            "im_color": popup.ids.movie_im_color.text,
            "contour": popup.ids.movie_contour.active,
            "stokes_i_sigma_cut": float(popup.ids.movie_stokes_i_sigma_cut.text),
            "contour_color": popup.ids.movie_contour_color.text,
            "contour_cmap": popup.ids.movie_contour_cmap.text,
            "contour_alpha": float(popup.ids.movie_contour_alpha.text),
            "contour_width": float(popup.ids.movie_contour_width.text),
            "component_color": popup.ids.movie_component_color.text,
            "rcparams": popup.ids.movie_rcparams.text,
            "background_color": popup.ids.movie_background_color.text,
            "plot_beam": popup.ids.movie_plot_beam.active,
            "xlim": [float(popup.ids.movie_ra_max.text),float(popup.ids.movie_ra_min.text)],
            "ylim": [float(popup.ids.movie_dec_min.text),float(popup.ids.movie_dec_max.text)],
            "plot_evpa": popup.ids.movie_plot_evpa.active,
            "rotate_evpa": float(popup.ids.movie_rotate_evpa.text),
            "evpa_len": float(popup.ids.movie_evpa_len.text),
            "evpa_distance": float(popup.ids.movie_evpa_distance.text),
            "evpa_color": popup.ids.movie_evpa_color.text,
            "evpa_width": float(popup.ids.movie_evpa_width.text),
            "lin_pol_sigma_cut": float(popup.ids.movie_lin_pol_sigma_cut.text),
            "save":save_path
              }

        #Check if we need to convolve with the common beam
        if popup.ids.movie_convolve.active:
            im_cube=im_cube.restore()

        #create movie
        im_cube.movie(**args)

        popup.dismiss()


    def show_mojave_popup(self):
        self.the_popup = MOJAVEExportPopup()
        self.the_popup.open()

        if self.mojave_password=="":
            self.mojave_password_dialog = MOJAVEPasswordPopup()
            self.mojave_password_dialog.open()

    def export_data_to_mojave(self,observer,password,source):
        upload_csv_to_MOJAVE(self.component_info_csv, observer, password, source)

    def import_data_from_mojave(self,source,observer,band,password):
        download_kinematic_from_MOJAVE(source,band,observer,password,self.ids.difmap_path.text,"/tmp/tmp_data_vcat")
        self.import_kinematics(["/tmp/tmp_data_vcat",""])

    #used to reimport data that was exported with the function above. Needs a directory path
    def import_kinematics(self,directory,check=True):

        #check if it is a valid vcat kinematics folder
        if (os.path.isfile(directory[0]+"/component_info.csv") and os.path.isfile(directory[0]+"/kinematic_fit.csv") and
                os.path.exists(directory[0]+"/modelfit_files") and os.path.exists(directory[0]+"/clean_fits")):

            #import the fits files
            selection=glob.glob(directory[0]+"/modelfit_files/*")
            self.load_modelfit(selection)
            selection = glob.glob(directory[0] + "/clean_fits/*")
            self.load_clean(selection)
            try:
                selection = glob.glob(directory[0]+"/uvf_files/*")
                self.load_uvf(selection)
            except:
                pass
            self.create_kinematic_plots()

            #add components
            comp_info=pd.read_csv(directory[0]+"/component_info.csv")
            ncomps=len(np.unique(comp_info[comp_info["component_number"]>=0]["component_number"]))
            for i in np.unique(comp_info[comp_info["component_number"]>=0]["component_number"]):
                self.add_component(number=i)

            #set core
            core_ind=comp_info[comp_info["is_core"]==True]["component_number"].values[0]
            for t in ToggleButton.get_widgets('components'):
                if "Component " +str(core_ind) in t.text:
                    t.state="down"
                    self.set_active_component()
                    self.set_core_component()
                    t.state='normal'
            self.set_active_component()

            frequencies=[]
            #now identify them with each other
            for plot in self.plots:
                for comp in plot.components:
                    print(comp[1].get_info())
                    try:
                        #match plot component with component in list
                        filter_df=comp_info[(round(comp_info["x"],15)==round(comp[1].x,15))]
                        ind=filter_df[round(comp_info["y"],15)==round(comp[1].y,15)]["component_number"].values[0]
                        print(ind)
                    except:
                        ind=-1
                    if ind>=0: #only do this if component was assigned
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

    def reset_kinematics(self):

        for i in range(np.max([len(self.components),len(self.ids.plot_frequency_select_buttons.children)])):
            #delete components
            for box in self.ids.kinematic_list.children:
                self.ids.kinematic_list.remove_widget(box)

            for box in self.ids.component_list.children:
                self.ids.component_list.remove_widget(box)

            #delete frequency buttons
            for box in self.ids.plot_frequency_select_buttons.children:
                self.ids.plot_frequency_select_buttons.remove_widget(box)

            for box in self.ids.scroll_frequency_select_buttons.children:
                self.ids.scroll_frequency_select_buttons.remove_widget(box)

        self.active_component_ind = None
        self.components = []

        #remove plots
        self.plots=[]

        #reset redshift
        self.ids.redshift.text = "0.0"

        #delete plots
        for figure in self.figure_widgets:
            self.ids.figures.remove_widget(figure)

        self.update_kinematic_plot()

    #### END OF KINEMATIC FUNCTIONS

    #### START OF MODELFIT FUNCTIONS

    def load_modelfit_files(self):

        # TODO some basic background checks on the files to see if they are valid and if they exist (write some popups)
        # TODO also check for polarizations. If there is only Stokes I input, grey out the polarization stacking options
        modelfit_pol_loaded=False

        if len(self.uvf_filepaths) != len(self.clean_filepaths):
            self.show_popup("Error","Please load an equal number of clean files and .uvf files!","Continue")
            return 0

        for ind, file in enumerate(self.clean_filepaths):
            uvf_file = self.uvf_filepaths[ind]
            if len(self.modelfit_filepaths) > ind:
                model = self.modelfit_filepaths[ind]
                plot_model = True
                plot_comp_ids = True
            else:
                model = ""
                plot_model = False
                plot_comp_ids = False
            if len(self.stokes_u_filepaths) > ind and len(self.stokes_q_filepaths) > ind:
                plot_data = ImageData(file, model=model, uvf_file=uvf_file, stokes_u=self.stokes_u_filepaths[ind],
                                      stokes_q=self.stokes_q_filepaths[ind], noise_method=self.noise_method,fit_comp_polarization=True,is_ehtim_model=self.ids.is_ehtim_model.active)
                image = FitsImage(plot_data, plot_mode="lin_pol", plot_model=plot_model, plot_comp_ids=plot_comp_ids,
                                  plot_evpa=True, evpa_color="black", contour_color="grey")

            else:
                plot_data = ImageData(file, model=model, uvf_file=uvf_file, noise_method=self.noise_method,fit_comp_polarization=True,is_ehtim_model=self.ids.is_ehtim_model.active)
                image = FitsImage(plot_data, plot_mode="lin_pol", plot_model=plot_model, plot_comp_ids=plot_comp_ids,
                                  plot_evpa=True, evpa_color="black", contour_color="grey")


            # check if polarization information was given:
            if np.sum(image.clean_image.stokes_q) == 0 or np.sum(image.clean_image.stokes_u) == 0:
                #TODO disable polarization options
                modelfit_pol_loaded=False
                pass
            else:
                modelfit_pol_loaded=True

            self.modelfit_plots.append(image)

            button = ToggleButton(
                text=str(get_date(file)),
                size_hint_y=None,
                size_hint_x=1,
                height=50,
                group="modelfit_plots",
            )
            button.bind(on_release=self.change_modelfit_plot)

            self.modelfit_plot_buttons.append(button)
            self.ids.modelfit_list.add_widget(button)

            # select first button by default
            try:
                self.modelfit_plot_buttons[0].state = "down"
                self.change_modelfit_plot(self.modelfit_plot_buttons[0])
            except:
                pass
        if not modelfit_pol_loaded:
            self.show_popup("Warning", "Data loaded successfully, \n but no valid polarization data detected",
                            "Continue")
        else:
            self.show_popup("Success", "Loaded data including polarization!", "Continue")

    def change_modelfit_plot(self,button=""):
        try:
            ind = np.where(np.array(self.modelfit_plot_buttons) == button)[0][0]
        except:
            for i1, but in enumerate(self.ids.modelfit_list.children):
                i1 = len(self.ids.modelfit_list.children) - 1 - i1
                if but.state == "down":
                    ind=i1
            try:
                ind
            except:
                return None

        self.ids.modelfit_image.figure = self.modelfit_plots[ind].fig

        self.ids.modelfit_component_list.clear_widgets()
        self.modelfit_component_buttons=[]

        #read out components attached to the ImageData object
        comps=self.modelfit_plots[ind].clean_image.components

        #add them to the component list
        for comp in comps:
            button = ToggleButton(
                text="Component " + str(comp.component_number),
                size_hint_y=None,
                size_hint_x=1,
                height=50,
                group="modelfit_components",
            )
            button.bind(on_release=self.change_modelfit_component_select)

            self.modelfit_component_buttons.append(button)
            self.ids.modelfit_component_list.add_widget(button)

        #create temporary working directory
        os.system("rm -rf " + "/tmp/tmp_data_vcat/")
        os.makedirs("/tmp/tmp_data_vcat", exist_ok=True)

        if self.ids.modelfit_method.text=="Stokes I/DIFMAP":
            #write .mod files
            write_mod_file_from_components(comps,"i","/tmp/tmp_data_vcat/i_model.mod")

            #calculate new residual plot
            get_residual_map(self.modelfit_plots[ind].clean_image.uvf_file,"/tmp/tmp_data_vcat/i_model.mod","/tmp/tmp_data_vcat/i_model.mod",
                             difmap_path=self.ids.difmap_path.text, save_location="/tmp/tmp_data_vcat/i_res.fits",
                             npix=len(self.modelfit_plots[ind].clean_image.X)*2,
                             pxsize=self.modelfit_plots[ind].clean_image.degpp*self.modelfit_plots[ind].clean_image.scale)

            #create residual plot
            res_image = ImageData("/tmp/tmp_data_vcat/i_res.fits")
            res_image.residual_map = res_image.Z
            res_plot = res_image.plot(plot_mode="residual",contour=False,title="Residual Map",component_color="green")

            #get current modelfit plot (right)
            modelfit_image = ImageData(self.modelfit_plots[ind].clean_image.fits_file).plot()

            for comp in comps:
                res_plot.plotComponent(comp.x,comp.y,comp.maj,comp.min,comp.pos,comp.scale,fillcolor="",id=comp.component_number)
                modelfit_image.plotComponent(comp.x,comp.y,comp.maj,comp.min,comp.pos,comp.scale,fillcolor="",id=comp.component_number)

        elif self.ids.modelfit_method.text=="Lin. Pol./ehtim":
            #write .mod files
            write_mod_file_from_components(comps, "q", "/tmp/tmp_data_vcat/q_model.mod")
            write_mod_file_from_components(comps, "u", "/tmp/tmp_data_vcat/u_model.mod")
            #calculate residual maps
            get_residual_map(self.modelfit_plots[ind].clean_image.uvf_file, "/tmp/tmp_data_vcat/q_model.mod","/tmp/tmp_data_vcat/q_model.mod",difmap_path=self.ids.difmap_path.text,
                             channel="q", save_location="/tmp/tmp_data_vcat/q_res.fits", npix=len(self.modelfit_plots[ind].clean_image.X)*2,
                             pxsize=self.modelfit_plots[ind].clean_image.degpp*self.modelfit_plots[ind].clean_image.scale)
            get_residual_map(self.modelfit_plots[ind].clean_image.uvf_file, "/tmp/tmp_data_vcat/u_model.mod","/tmp/tmp_data_vcat/u_model.mod",difmap_path=self.ids.difmap_path.text,
                             channel="u", save_location="/tmp/tmp_data_vcat/u_res.fits", npix=len(self.modelfit_plots[ind].clean_image.X)*2,
                             pxsize=self.modelfit_plots[ind].clean_image.degpp*self.modelfit_plots[ind].clean_image.scale)

            #create residual image
            res_image = ImageData(self.modelfit_plots[ind].clean_image.fits_file,
                                  stokes_q="/tmp/tmp_data_vcat/q_res.fits", stokes_u="/tmp/tmp_data_vcat/u_res.fits")

            res_plot = res_image.plot(plot_mode="lin_pol", plot_evpa=True,title="Residual Map")
            # get current modelfit plot (right)
            modelfit_image = ImageData(self.modelfit_plots[ind].clean_image.fits_file,
                                       stokes_q=self.modelfit_plots[ind].clean_image.stokes_q_path,
                                       stokes_u=self.modelfit_plots[ind].clean_image.stokes_u_path).plot(plot_mode="lin_pol",plot_evpa=True)

            for comp in comps:
                res_plot.plotComponent(comp.x,comp.y,comp.maj,comp.min,comp.pos,comp.scale,fillcolor="",id=comp.component_number,evpa=comp.evpa)
                modelfit_image.plotComponent(comp.x, comp.y, comp.maj, comp.min, comp.pos, comp.scale, fillcolor="",
                                       id=comp.component_number, evpa=comp.evpa)

        elif self.ids.modelfit_method.text=="Full Pol./ehtim":
            #write .mod files
            write_mod_file_from_components(comps, "i", "/tmp/tmp_data_vcat/i_model.mod")
            write_mod_file_from_components(comps, "q", "/tmp/tmp_data_vcat/q_model.mod")
            write_mod_file_from_components(comps, "u", "/tmp/tmp_data_vcat/u_model.mod")
            #calculate residual maps
            get_residual_map(self.modelfit_plots[ind].clean_image.uvf_file, "/tmp/tmp_data_vcat/i_model.mod","/tmp/tmp_data_vcat/i_model.mod",
                             difmap_path=self.ids.difmap_path.text,
                             channel="i", save_location="/tmp/tmp_data_vcat/i_res.fits",
                             npix=len(self.modelfit_plots[ind].clean_image.X) * 2,
                             pxsize=self.modelfit_plots[ind].clean_image.degpp * self.modelfit_plots[
                                 ind].clean_image.scale)
            get_residual_map(self.modelfit_plots[ind].clean_image.uvf_file, "/tmp/tmp_data_vcat/q_model.mod", "/tmp/tmp_data_vcat/q_model.mod",
                             difmap_path=self.ids.difmap_path.text,
                             channel="q", save_location="/tmp/tmp_data_vcat/q_res.fits", npix=len(self.modelfit_plots[ind].clean_image.X)*2,
                             pxsize=self.modelfit_plots[ind].clean_image.degpp*self.modelfit_plots[ind].clean_image.scale)
            get_residual_map(self.modelfit_plots[ind].clean_image.uvf_file, "/tmp/tmp_data_vcat/u_model.mod","/tmp/tmp_data_vcat/u_model.mod",
                             difmap_path=self.ids.difmap_path.text,
                             channel="u", save_location="/tmp/tmp_data_vcat/u_res.fits", npix=len(self.modelfit_plots[ind].clean_image.X)*2,
                             pxsize=self.modelfit_plots[ind].clean_image.degpp*self.modelfit_plots[ind].clean_image.scale)

            #create residual image
            res_image = ImageData("/tmp/tmp_data_vcat/i_res.fits",
                                  stokes_q="/tmp/tmp_data_vcat/q_res.fits", stokes_u="/tmp/tmp_data_vcat/u_res.fits")

            res_plot = res_image.plot(plot_mode="lin_pol", plot_evpa=True,title="Residual Map")
            # get current modelfit plot (right)
            modelfit_image = ImageData(self.modelfit_plots[ind].clean_image.fits_file,
                                       stokes_q=self.modelfit_plots[ind].clean_image.stokes_q_path,
                                       stokes_u=self.modelfit_plots[ind].clean_image.stokes_u_path).plot(plot_mode="lin_pol",plot_evpa=True)

            for comp in comps:
                res_plot.plotComponent(comp.x,comp.y,comp.maj,comp.min,comp.pos,comp.scale,fillcolor="",id=comp.component_number,evpa=comp.evpa)
                modelfit_image.plotComponent(comp.x, comp.y, comp.maj, comp.min, comp.pos, comp.scale, fillcolor="",
                                       id=comp.component_number, evpa=comp.evpa)
        #assign image
        self.ids.modelfit_residual_image.figure = res_plot.fig
        self.ids.modelfit_image.figure = modelfit_image.fig
        self.current_modelfit_residual_plot = res_plot


    def change_modelfit_component_select(self,button=""):
        for i1, but in enumerate(self.ids.modelfit_list.children):
            i1 = len(self.ids.modelfit_list.children) - 1 - i1
            if but.state=="down":
                for ind,button in enumerate(self.ids.modelfit_component_list.children):
                    ind=len(self.ids.modelfit_component_list.children)-1-ind
                    if button.state == "down":
                        #find current component
                        current_comp=self.modelfit_plots[i1].clean_image.components[ind]

                        #set parameters
                        self.ids.component_flux.text=f"{current_comp.flux*1e3:.2f}"
                        self.ids.component_pos.text=f"{current_comp.pos:.2f}"
                        self.ids.component_maj.text=f"{current_comp.maj*current_comp.scale:.2f}"
                        self.ids.component_min.text=f"{current_comp.min*current_comp.scale:.2f}"
                        self.ids.component_x.text=f"{current_comp.x*current_comp.scale:.2f}"
                        self.ids.component_y.text=f"{current_comp.y*current_comp.scale:.2f}"
                        self.ids.component_lin_pol.text=f"{current_comp.lin_pol*1e3:.2f}"
                        self.ids.component_evpa.text=f"{current_comp.evpa:.2f}"

    #Function to add new modelfit component
    def add_modelfit_component(self):
        for i1, but in enumerate(self.ids.modelfit_list.children):
            i1 = len(self.ids.modelfit_list.children) - 1 - i1
            if but.state == "down":
                image_data=self.modelfit_plots[i1].clean_image
                self.modelfit_plots[i1].clean_image.components.append(
                    Component(float(self.ids.component_x.text)/image_data.scale,
                              float(self.ids.component_y.text)/image_data.scale,
                              float(self.ids.component_maj.text)/image_data.scale,
                              float(self.ids.component_min.text)/image_data.scale,
                              float(self.ids.component_pos.text),
                              float(self.ids.component_flux.text)*1e-3,
                              image_data.date,image_data.mjd,Time(image_data.date).decimalyear,
                              lin_pol=float(self.ids.component_lin_pol.text) * 1e-3,
                              evpa=float(self.ids.component_evpa.text),
                              scale=image_data.scale,component_number=len(image_data.components)))

                #reload plots
                self.change_modelfit_plot(but)

    #Function to remove modelfit component
    def remove_modelfit_component(self,button=""):
        for i1, but in enumerate(self.ids.modelfit_list.children):
            i1 = len(self.ids.modelfit_list.children) - 1 - i1
            if but.state == "down":
                for ind,button in enumerate(self.ids.modelfit_component_list.children):
                    if button.state == "down":
                        ind = len(self.ids.modelfit_component_list.children) - 1 - ind
                        self.modelfit_plots[i1].clean_image.components.pop(ind)

                        # reload plots
                        self.change_modelfit_plot(but)

    def update_modelfit_component(self,button=""):
        for i1, but in enumerate(self.ids.modelfit_list.children):
            i1=len(self.ids.modelfit_list.children)-1-i1
            if but.state == "down":
                for ind,button in enumerate(self.ids.modelfit_component_list.children):
                    if button.state == "down":
                        ind = len(self.ids.modelfit_component_list.children) - 1 - ind
                        comp=self.modelfit_plots[i1].clean_image.components[ind]
                        try:
                            self.modelfit_plots[i1].clean_image.components[ind].flux=float(self.ids.component_flux.text)*1e-3
                            self.modelfit_plots[i1].clean_image.components[ind].x = float(self.ids.component_x.text)/comp.scale
                            self.modelfit_plots[i1].clean_image.components[ind].y = float(self.ids.component_y.text)/comp.scale
                            self.modelfit_plots[i1].clean_image.components[ind].pos = float(self.ids.component_pos.text)
                            self.modelfit_plots[i1].clean_image.components[ind].maj = float(self.ids.component_maj.text)/comp.scale
                            self.modelfit_plots[i1].clean_image.components[ind].min = float(self.ids.component_min.text)/comp.scale
                            self.modelfit_plots[i1].clean_image.components[ind].lin_pol = float(self.ids.component_lin_pol.text)*1e-3
                            self.modelfit_plots[i1].clean_image.components[ind].evpa = float(self.ids.component_evpa.text)
                        except:
                            pass
                        # reload plots
                        self.change_modelfit_plot(but)

    #TODO load modelfit from other image as start model
    def load_modelfit_components(self):
        pass

    def do_modelfit(self,button=""):
        #find currently selected image
        ind=-1
        final_but=0
        for i1, but in enumerate(self.ids.modelfit_list.children):
            i1 = len(self.ids.modelfit_list.children) - 1 - i1
            if but.state == "down":
                ind=i1
                final_but=but

        if ind==-1:
            self.show_popup("Error","Please select a file from the file list first!", "Continue")
            return None

        #get current components
        comps=self.modelfit_plots[ind].clean_image.components
        if len(comps)==0:
            self.show_popup("Error", "Please add components first!", "Continue")

        if self.ids.modelfit_method.text=="Stokes I/DIFMAP":
            write_mod_file_from_components(self.modelfit_plots[ind].clean_image.components, "i", export="/tmp/mod.mod", adv=True)
            self.modelfit_plots[ind].clean_image.components=modelfit_difmap(self.modelfit_plots[ind].clean_image.uvf_file,
                                                                            "/tmp/mod.mod",int(self.ids.modelfit_count.text),
                                                                            difmap_path=self.ids.difmap_path.text,
                                                                            components=self.modelfit_plots[ind].clean_image.components)
            self.modelfit_plots[ind].clean_image.fit_comp_polarization()
        elif self.ids.modelfit_method.text=="Lin. Pol./ehtim":
            if self.ids.modelfit_minimizer.text=="Dynesty Dynamic":
                minimizer="dynesty_dynamic"
            else:
                minimizer="scipy.optimize.minimize"
            self.modelfit_plots[ind].clean_image.components=modelfit_ehtim_pol(self.modelfit_plots[ind].clean_image.uvf_file,
                                                                           comps,int(self.ids.modelfit_count.text),
                                                                           npix=self.ids.modelfit_npix.text,
                                                                           fov=self.ids.modelfit_fov.text,
                                                                            minimizer=minimizer,
                                                                           nwalker=int(self.ids.nwalker_count.text),
                                                                            max_size=float(self.ids.max_model_size.text),
                                                                            max_flux=float(self.ids.max_model_flux.text),
                                                                            max_dist=float(self.ids.max_model_dist.text),
                                                                            circ_gauss=self.ids.fit_circ.active)

        elif self.ids.modelfit_method.text=="Full Pol./ehtim":
            if self.ids.modelfit_minimizer.text=="Dynesty Dynamic":
                minimizer="dynesty_dynamic"
            else:
                minimizer="scipy.optimize.minimize"
            self.modelfit_plots[ind].clean_image.components=modelfit_ehtim_full_pol(self.modelfit_plots[ind].clean_image.uvf_file,
                                                                           comps,int(self.ids.modelfit_count.text),
                                                                           npix=self.ids.modelfit_npix.text,
                                                                           fov=self.ids.modelfit_fov.text,
                                                                            minimizer=minimizer,
                                                                           nwalker=int(self.ids.nwalker_count.text),
                                                                            max_size=float(self.ids.max_model_size.text),
                                                                            max_flux=float(self.ids.max_model_flux.text),
                                                                            max_dist=float(self.ids.max_model_dist.text),
                                                                            circ_gauss=self.ids.fit_circ.active)


        self.change_modelfit_plot(final_but)

    def open_export_modelfit_popup(self,button=""):
        if len(self.modelfit_plots) > 0:
            self.export_modelfit_dialog = ModelfitExportPopup()
            self.export_modelfit_dialog.open()

    def set_residual_touch_mode(self,button="",mode="zoom"):
        if mode=="zoom":
            self.ids.modelfit_residual_image.touch_mode="cursor"
        elif mode=="add_component":
            self.ids.modelfit_residual_image.touch_mode="ellipseselect"
            self.new_modelfit_component_coords=[]


    def plot_new_modelfit_component(self,xdata,ydata):
        self.new_modelfit_component_coords.append([xdata,ydata])

        #find currently active plot
        ind=-1
        for i1, but in enumerate(self.ids.modelfit_list.children):
            i1 = len(self.ids.modelfit_list.children) - 1 - i1
            if but.state == "down":
                ind = i1
        if ind==-1:
            return None

        ax=self.current_modelfit_residual_plot.ax

        n_points=len(self.new_modelfit_component_coords)
        try:
            self.current_new_ellipse_plot.remove()
        except:
            pass
        if n_points==1:
            #plot point
            self.current_new_ellipse_plot  = ax.scatter(self.new_modelfit_component_coords[0][0],
                                                        self.new_modelfit_component_coords[0][1])
            ax.figure.canvas.draw_idle()
            ax.figure.canvas.flush_events()
        elif n_points==2:
            #plot line
            self.current_new_ellipse_plot,  = ax.plot([self.new_modelfit_component_coords[0][0],
                                                       self.new_modelfit_component_coords[1][0]],
                                                       [self.new_modelfit_component_coords[0][1],
                                                        self.new_modelfit_component_coords[1][1]])
            ax.figure.canvas.draw_idle()
            ax.figure.canvas.flush_events()
        elif n_points==3:
            # add component as ellipse and empty new_modelfit_component_coords
            #center of the ellipse
            center_x, center_y = self.new_modelfit_component_coords[0][0], self.new_modelfit_component_coords[0][1]

            #calculate one axis of the ellipse:
            maj_x, maj_y = (self.new_modelfit_component_coords[1][0] - center_x ,
                            self.new_modelfit_component_coords[1][1] - center_y)

            #calculate the second axis of the ellipse
            min_x, min_y = (self.new_modelfit_component_coords[2][0] - center_x ,
                            self.new_modelfit_component_coords[2][1] - center_y)

            maj=2*np.sqrt(maj_x**2+maj_y**2)
            min=2*np.sqrt(min_x**2+min_y**2)

            if maj<min:
                maj_old=maj
                maj=min
                min=maj_old

            #calculate position angle
            PA=np.arctan2(maj_x,maj_y)/np.pi*180

            if PA<-90:
                PA+=180
            elif PA>90:
                PA-=180

            #plot ellipse
            self.current_new_ellipse_plot=Ellipse((center_x,center_y),width=maj,height=min,angle=-(PA+90))
            ax.add_patch(self.current_new_ellipse_plot)

            #modify the text boxes
            self.ids.component_x.text=f"{center_x:.2f}"
            self.ids.component_y.text=f"{center_y:.2f}"
            self.ids.component_maj.text = f"{maj:.2f}"
            self.ids.component_min.text = f"{min:.2f}"
            self.ids.component_pos.text = f"{PA:.2f}"

            #extract flux and linpol/EVPA from image within ellipse
            # Get the indices of the ellipse in image space
            self.current_modelfit_residual_plot.clean_image.masking(mask_type="reset")
            self.current_modelfit_residual_plot.clean_image.masking(mask_type="ellipse",
                                                                    args={'e_args':[maj,min,PA],
                                                                          'e_xoffset': center_x,'e_yoffset': center_y})

            mask=np.array(self.current_modelfit_residual_plot.clean_image.mask,dtype=bool)

            if self.ids.modelfit_method.text == "Stokes I/DIFMAP":
                flux=np.sum(self.current_modelfit_residual_plot.clean_image.residual_map[mask])
            else:
                flux = np.sum(self.current_modelfit_residual_plot.clean_image.Z[mask])
            lin_pol=np.sum(self.current_modelfit_residual_plot.clean_image.lin_pol[mask])

            try:
                stokes_q = np.sum(self.current_modelfit_residual_plot.clean_image.stokes_q[mask])
                stokes_u = np.sum(self.current_modelfit_residual_plot.clean_image.stokes_u[mask])
                evpa=1/2*np.arctan2(stokes_u,stokes_q)
            except:
                evpa = 0

            bmaj=self.modelfit_plots[i1].clean_image.beam_maj
            bmin=self.modelfit_plots[i1].clean_image.beam_min
            pxincr=self.current_modelfit_residual_plot.clean_image.degpp*self.current_modelfit_residual_plot.clean_image.scale


            self.ids.component_flux.text = f"{JyPerBeam2Jy(flux,bmaj,bmin,pxincr)*1e3:.2f}"
            self.ids.component_lin_pol.text = f"{JyPerBeam2Jy(lin_pol,bmaj,bmin,pxincr) * 1e3:.2f}"
            self.ids.component_evpa.text = f"{evpa/np.pi*180:.2f}"


            self.new_modelfit_component_coords=[]

    #### END OF MODELFIT FUNCTIONS

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
                plot_data=ImageData(file,model=model,stokes_u=self.stokes_u_filepaths[ind],
                                    stokes_q=self.stokes_q_filepaths[ind],noise_method=self.noise_method,
                                    is_ehtim_model=self.ids.is_ehtim_model.active)
                image=FitsImage(plot_data,plot_mode="frac_pol",plot_evpa=True,evpa_color="black",contour_color="grey")
            else:
                #try to load model from clean .fits file
                if len(self.modelfit_filepaths)>ind:
                    model=self.modelfit_filepaths[ind]
                else:
                    model=""
                plot_data=ImageData(file,model=model,noise_method=self.noise_method,is_ehtim_model=self.ids.is_ehtim_model.active)
                image=FitsImage(plot_data,plot_mode="frac_pol",plot_evpa=True,evpa_color="black",contour_color="grey")


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
        if np.sum(image.clean_image.stokes_q)==0 or np.sum(image.clean_image.stokes_u)==0:
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

        #create ImageData and ImageCube
        im_datas=[]
        for i in range(len(files_to_stack)):
            try:
                fits_q=files_to_stack_q[i]
            except:
                fits_q=""
            try:
                fits_u=files_to_stack_u[i]
            except:
                fits_u=""
            try:
                uvf_file=files_to_stack_uvf[i]
            except:
                uvf_file=""
            try:
                model=files_to_stack_models[i]
            except:
                model=""
            try:
                model=files_to_stack_casa_clean_models[i]
                is_casa_model=True
            except:
                is_casa_model=False
                try:
                    model = files_to_stack_models[i]
                except:
                    model = ""


            im_datas.append(ImageData(fits_file=files_to_stack[i],
                                      stokes_q=fits_q,
                                      stokes_u=fits_u,
                                      uvf_file=uvf_file,
                                      model=model,
                                      is_casa_model=is_casa_model,
                                      noise_method=self.noise_method,
                                      is_ehtim_model=self.ids.is_ehtim_model.active
                                      ))
        im_cube=ImageCube(image_data_list=im_datas)

        #check if we need to restore the beam sizes first before proceeding
        if do_beam_restore:
            #Restore Stokes I
            try:
                restore_maj=float(self.ids.beam_maj.text)
                restore_min=float(self.ids.beam_min.text)
                restore_pa=float(self.ids.beam_pos.text)
            except:
                self.show_popup("Error","Please set beam to restore with!","Continue")

            im_cube=im_cube.restore(beam_maj=restore_maj,beam_min=restore_min,beam_posa=restore_pa)



        align=self.ids.do_alignment_check.active
        if align:
            im_cube=im_cube.center()

        stack_image=im_cube.stack(mode="all").images.flatten()[0]

        #create plot
        self.final_stack_image=stack_image
        self.ids.stacked_image_i.state="down"
        StackPlot = FitsImage(stack_image,title="Stacked Image")
        self.ids.stacked_image.figure = StackPlot.fig

    def change_final_stack_plot(self,mode,button):
        if self.final_stack_image!="":
            if mode in ["lin_pol","frac_pol"]:
                plot_evpa=True
                if mode =="frac_pol":
                    evpa_color="black"
                    contour_color="grey"
                else:
                    evpa_color="white"
                    contour_color="grey"
            else:
                plot_evpa=False
                evpa_color="white"
                contour_color="grey"
            StackPlot = FitsImage(self.final_stack_image,plot_mode=mode,plot_evpa=plot_evpa,evpa_color=evpa_color,
                                  contour_color=contour_color,title="Stacked Image")
            self.ids.stacked_image.figure = StackPlot.fig
            button.state="down"

    def fill_in_common_beam(self):
        files_to_stack=[]
        for ind,box in enumerate(self.stacking_single_plot_checkboxes):
            if box.active:
                files_to_stack.append(self.clean_filepaths[ind])
        if len(files_to_stack)>0:
            im_cube=ImageCube().import_files(files_to_stack)
            maj,min,pos=im_cube.get_common_beam()
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
                ImageData(self.casa_clean_model_filepaths[ind],model_save_dir="tmp/",is_casa_model=True,noise_method=self.noise_method,is_ehtim_model=self.ids.is_ehtim_model.active)
            except:
                pass
            plot_data=ImageData(clean_path,model=model_path,stokes_q=stokes_q_path,stokes_u=stokes_u_path,noise_method=self.noise_method,is_ehtim_model=self.ids.is_ehtim_model.active)
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
        self.ids.beam_text.text = "{:.2f}".format(image.beam_maj) + " mas x " + "{:.2f}".format(image.beam_min) + " mas, " + "{:.2f}".format(image.beam_pa)+"°"
        if image.integrated_flux_clean!=0:
            self.ids.fractional_polarization.text = (
                    "{:.2f}".format(image.integrated_pol_flux_clean/image.integrated_flux_clean*100)+"%")
        else:
            self.ids.fractional_polarization.text = "-"

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
                                                        plot_model = self.ids.overplot_gauss.active,
                                                        plot_clean = self.ids.overplot_clean.active,
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
                                                        background_color=self.ids.background_color.text,
                                                        rcparams = self.ids.rcparams.text,
                                                        )

        #find active button
        active_button = next((t for t in ToggleButton.get_widgets('plotting_single_plots') if t.state == 'down'), None)
        self.change_plotting_single_plot(button=active_button)

    def open_export_plot_dialog(self):
        if len(self.plotting_single_plots)>0:
            self.export_plot_dialog = PlotExportPopup()
            self.export_plot_dialog.open()

    def export_plotting_plot(self,export_name):
        # find active button
        active_button = next((t for t in ToggleButton.get_widgets('plotting_single_plots') if t.state == 'down'), None)
        ind = np.where(np.array(self.plotting_single_plot_buttons) == active_button)[0][0]
        active_plot = self.plotting_single_plots[ind]
        active_plot.export(export_name)

    def export_modelfit(self,export_name):
        active_button = next((t for t in ToggleButton.get_widgets('modelfit_plots') if t.state == 'down'), None)
        ind = np.where(np.array(self.modelfit_plot_buttons) == active_button)[0][0]
        active_plot = self.modelfit_plots[ind]

        # get current components
        comps = active_plot.clean_image.components
        if len(comps) == 0:
            self.show_popup("Error", "Please add components first!", "Continue")

        if self.ids.modelfit_method.text == "Stokes I/DIFMAP":
            write_mod_file_from_components(comps, "i", export=export_name,adv=False)
        elif self.ids.modelfit_method.text == "Lin. Pol./ehtim":
            if self.ids.modelfit_minimizer.text == "Dynesty Dynamic":
                minimizer = "dynesty_dynamic"
            else:
                minimizer = "scipy.optimize.minimize"
            #we run ehtim with niter=0 and an export name to export the current model
            self.modelfit_plots[ind].clean_image.components = modelfit_ehtim_pol(
                self.modelfit_plots[ind].clean_image.uvf_file,
                comps, 0,
                npix=self.ids.modelfit_npix.text,
                fov=self.ids.modelfit_fov.text,
                minimizer=minimizer,
                nwalker=int(self.ids.nwalker_count.text),
                max_size=float(self.ids.max_model_size.text),
                max_flux=float(self.ids.max_model_flux.text),
                max_dist=float(self.ids.max_model_dist.text),
                circ_gauss=self.ids.fit_circ.active,
                export_model=export_name,
                skip_fit=True)
        elif self.ids.modelfit_method.text == "Full Pol./ehtim":
            if self.ids.modelfit_minimizer.text == "Dynesty Dynamic":
                minimizer = "dynesty_dynamic"
            else:
                minimizer = "scipy.optimize.minimize"
            #we run ehtim with niter=0 and an export name to export the current model
            self.modelfit_plots[ind].clean_image.components = modelfit_ehtim_full_pol(
                self.modelfit_plots[ind].clean_image.uvf_file,
                comps, 0,
                npix=self.ids.modelfit_npix.text,
                fov=self.ids.modelfit_fov.text,
                minimizer=minimizer,
                nwalker=int(self.ids.nwalker_count.text),
                max_size=float(self.ids.max_model_size.text),
                max_flux=float(self.ids.max_model_flux.text),
                max_dist=float(self.ids.max_model_dist.text),
                circ_gauss=self.ids.fit_circ.active,
                export_model=export_name,
                skip_fit=True)

class VCAT(App):

    def build(self):
        self.screen=ModelFits()
        self.screen.noise_method=self.screen.ids.noise_method.text
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

    def save_movie(self,popup,save_text,selection):
        self.screen.save_movie(popup,save_text,selection)

    def import_kinematics(self,directory):
        self.screen.import_kinematics(directory)

    def export_plotting_plot(self,name,selection):
        self.screen.export_plotting_plot(name)

    def export_modelfit(self,name,selection):
        self.screen.export_modelfit(name)

    def export_data_to_mojave(self,observer,password,source):
        self.screen.export_data_to_mojave(observer,password,source)

    def import_data_from_mojave(self,source,observer,band,password):
        self.screen.import_data_from_mojave(source,observer,band,password)

    def set_mojave_password(self,password):
        self.screen.set_mojave_password(password)

    #### END OF KINEMATIC FUNCTIONS

    #### START OF MODELFIT FUNCTIONS

    def plot_new_modelfit_component(self,xdata,ydata):
        self.screen.plot_new_modelfit_component(xdata,ydata)

if __name__ == "__main__":
    VCAT().run()
