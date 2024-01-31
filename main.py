from kivy.app import App
from kivy.uix.popup import Popup
from kivy.uix.togglebutton import ToggleButton
from kivy.uix.tabbedpanel import TabbedPanel
from kivy.uix.boxlayout import BoxLayout
from kivy.uix.label import Label
from kivy.properties import ObjectProperty, StringProperty, ListProperty
from kivy.uix.button import Button
from graph_widget import MatplotFigure
from kivy.utils import platform
from kinematics import ComponentCollection
import numpy as np
from graph_generator import FitsImage, KinematicPlot
from astroquery.ipac.ned import Ned
from astropy.io import fits
import os
import pandas as pd
import glob


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
    filepaths = ListProperty([])
    plots = []
    components = []
    active_component_ind=None
    component_colors=["green","red","blue","yellow","orange","purple"]
    figure_widgets=[]
    core_component_ind=None
    name=""
    component_collections = []

    def current_color(self,ind):
        if ind is not None:
            ind=int(ind % len(self.component_colors))
            return self.component_colors[ind]
        else:
            return None

    def open_popup(self):
        self.the_popup = FileChoosePopup(load=self.load)
        self.the_popup.open()

    def open_save_dialog(self):
        self.save_dialog = FileSavePopup()
        self.save_dialog.open()

    def open_import_dialog(self):
        self.import_dialog = FileImportPopup()
        self.import_dialog.open()

    def load(self, selection):
        self.plots=[]
        self.filepaths = selection
        try:
            self.the_popup.dismiss()
        except:
            pass
        self.file_info=str(len(self.filepaths))+" Files selected"
        self.ids.get_file.text = self.file_info

        #sort files by date
        date=[]
        fits_images=[]

        #create plots for view page
        for filepath in self.filepaths:
            plot=FitsImage(filepath,filepath)
            date=np.append(date,plot.get_date(fits.open(filepath)))
            fits_images=np.append(fits_images,plot)

        #sort them according to obsdate
        args=date.argsort()
        fits_images=fits_images[args]

        for plot in fits_images:
            self.plots.append(plot)
            self.name=plot.name
            new_figure=MatplotFigure(size_hint_x=None,
            width=100)
            new_figure.figure=plot.fig
            new_figure.touch_mode="pan"
            self.ids.figures.add_widget(new_figure)
            self.figure_widgets.append(new_figure)

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

    def set_active_component(self):
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

    def show_popup(self,title,text,button_text):
        popup = Popup(title=title, auto_dismiss=False, size_hint=(0.35,0.35))
        box_layout=BoxLayout(orientation="vertical")
        box_layout.add_widget(Label(text=text,size_hint_y=0.8))
        close_button=Button(text=button_text,size_hint_y=0.2)
        close_button.bind(on_press=popup.dismiss)
        box_layout.add_widget(close_button)
        popup.add_widget(box_layout)
        popup.open()

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
        os.makedirs(save_path+"/fits",exist_ok=True)

        #copy fits files to common directory
        for file in self.filepaths:
            os.system("cp " + file + " " + save_path + "/fits/")

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
                os.path.exists(directory[0]+"/fits")):

            #import the fits files
            selection=glob.glob(directory[0]+"/fits/*")
            self.load(selection)

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
                    t.state='normal'
            self.set_core_component()
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


class VCAT(App):

    def build(self):
        self.screen=ModelFits()
        return self.screen

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

if __name__ == "__main__":
    VCAT().run()