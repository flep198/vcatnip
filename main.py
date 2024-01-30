from kivy.app import App
from kivy.uix.popup import Popup
from kivy.uix.togglebutton import ToggleButton
from kivy.uix.tabbedpanel import TabbedPanel
from kivy.uix.boxlayout import BoxLayout
from kivy.uix.label import Label
from kivy.uix.scrollview import ScrollView
from kivy.properties import ObjectProperty, StringProperty, ListProperty
from graph_widget import MatplotFigure
from kivy.utils import platform
from itertools import cycle
from kinematics import ComponentCollection, Component
import numpy as np
from kivy.lang import Builder
from graph_generator import FitsImage, KinematicPlot
from astroquery.ipac.ned import Ned
from astropy.io import fits

#avoid conflict between mouse provider and touch (very important with touch device)
#no need for android platform
if platform != 'android':
    from kivy.config import Config
    Config.set('input', 'mouse', 'mouse,disable_on_activity')

class FileChoosePopup(Popup):
    load = ObjectProperty()

class FileSavePopup(Popup):
    load = ObjectProperty

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

    def load(self, selection):
        self.plots=[]
        self.filepaths = selection
        self.the_popup.dismiss()
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

    def set_active_component(self,button):
        # find active button
        active_button = next((t for t in ToggleButton.get_widgets('components') if t.state == 'down'), None)
        if active_button:
            #set active component
            for i,comp in enumerate(self.components):
                if comp in active_button.text:
                    self.active_component_ind=i
        else:
            self.active_component_ind = None


    def update_kinematic_plot(self):

        #currently creates a new plot everytime => maybe rewrite to just update existing plot?
        self.kplot=KinematicPlot()
        self.ids.kinematic_plot.figure=self.kplot.fig

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
                    labels[2].text="{:.2f}".format(float(fit_data["speed"]))+" +/- "+"{:.2f}".format(float(fit_data["speed_err"]))
                    labels[1].text="{:.2f}".format(float(fit_data["beta_app"]))+" +/- "+"{:.2f}".format(float(fit_data["beta_app_err"]))
                    labels[0].text="{:.2f}".format(float(fit_data["d_crit"]))+" +/- "+"{:.2f}".format(float(fit_data["d_crit_err"]))
                except:
                    pass

            self.kplot.ax.legend()
            self.kplot.ax.figure.canvas.draw_idle()
            self.kplot.ax.figure.canvas.flush_events()

    def get_redshift(self):

        #get redshift from NED, otherwise set to zero
        try:
            self.redshift = np.average(Ned.get_table(self.name, table="redshifts")["Published Redshift"])
        except:
            self.redshift = 0.00

        self.ids.redshift.text = str(self.redshift)
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


if __name__ == "__main__":
    VCAT().run()