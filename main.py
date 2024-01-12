from kivy.app import App
from kivy.uix.popup import Popup
from kivy.uix.togglebutton import ToggleButton
from kivy.uix.tabbedpanel import TabbedPanel
from kivy.properties import ObjectProperty, StringProperty, ListProperty
from graph_widget import MatplotFigure
from kivy.utils import platform
from itertools import cycle
from kinematics import Component
import numpy as np
from kivy.lang import Builder
from graph_generator import FitsImage, KinematicPlot

#avoid conflict between mouse provider and touch (very important with touch device)
#no need for android platform
if platform != 'android':
    from kivy.config import Config
    Config.set('input', 'mouse', 'mouse,disable_on_activity')




class FileChoosePopup(Popup):
    load = ObjectProperty()

class ModelFits(TabbedPanel):
    file_info = StringProperty("No file chosen")
    the_popup = ObjectProperty(None)
    filepaths = ListProperty([])
    plots = []
    components = []
    active_component_ind=None
    component_colors=["green","red","blue","yellow","orange","purple"]
    figure_widgets=[]

    def current_color(self,ind):
        if ind is not None:
            ind=int(ind % len(self.component_colors))
            return self.component_colors[ind]
        else:
            return None

    def open_popup(self):
        self.the_popup = FileChoosePopup(load=self.load)
        self.the_popup.open()
    def load(self, selection):
        self.plots=[]
        self.filepaths = selection
        self.the_popup.dismiss()
        self.file_info=str(len(self.filepaths))+" Files selected"
        self.ids.get_file.text = self.file_info


        #create plots for view page
        for filepath in self.filepaths:
            plot=FitsImage(filepath,filepath)
            self.plots.append(plot)
            new_figure=MatplotFigure()
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


    def remove_component(self):
        #find active button
        button_to_remove = next( (t for t in ToggleButton.get_widgets('components') if t.state=='down'), None)
        #remove button
        if button_to_remove:
            #if button_to_remove.text in self.components:
            #    self.components.remove(button_to_remove.text)
            self.ids.component_list.remove_widget(button_to_remove)

        #need to change plots as well and delete component association.

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

        #TODO: Before creating the kinematics plot, run the entire kinematic analysis first! Then plot the results!
        #update kinematic plot
        self.update_kinematic_plot()

    def set_active_component(self,button):
        # find active button
        active_button = next((t for t in ToggleButton.get_widgets('components') if t.state == 'down'), None)
        if active_button:
            #set active component
            for i,comp in enumerate(self.components):
                if comp == active_button.text:
                    self.active_component_ind=i
        else:
            self.active_component_ind = None


    def update_kinematic_plot(self):
        #currently creates a new plot everytime => maybe rewrite to just update existing plot?
        self.kplot=KinematicPlot()
        self.ids.kinematic_plot.figure=self.kplot.fig

        #plot components
        for plot in self.plots:
            for comp in plot.components:
                color=comp[0][1][0].get_color()
                self.kplot.plot_component(comp[1],color)




class VCAT(App):

    def build(self):
        self.screen=ModelFits()
        return self.screen

    def set_touch_mode(self,touch_mode):
        for figure in self.screen.figure_widgets:
            figure.touch_mode = touch_mode

    def home_plot(self):
        for figure in self.screen.figure_widgets:
            figure.home()

    def add_component(self):
        self.screen.add_component()

    def remove_component(self):
        self.screen.remove_component()

    def choose_component(self,figure_widget,x,y):
        self.screen.choose_component(figure_widget,x,y)


if __name__ == "__main__":
    VCAT().run()