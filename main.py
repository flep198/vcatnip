from kivy.app import App
from kivy.uix.popup import Popup
from kivy.uix.togglebutton import ToggleButton
from kivy.uix.tabbedpanel import TabbedPanel
from kivy.properties import ObjectProperty, StringProperty, ListProperty
from graph_widget import MatplotFigure
from kivy.utils import platform


#avoid conflict between mouse provider and touch (very important with touch device)
#no need for android platform
if platform != 'android':
    from kivy.config import Config
    Config.set('input', 'mouse', 'mouse,disable_on_activity')

from kivy.lang import Builder
from graph_generator import FitsImage


class FileChoosePopup(Popup):
    load = ObjectProperty()

class ModelFits(TabbedPanel):
    file_info = StringProperty("No file chosen")
    the_popup = ObjectProperty(None)
    filepaths = ListProperty([])
    plots = []
    components = []
    kinematics = []

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
            self.plots.append(FitsImage(filepath,filepath))

        if len(self.plots)==1:
            self.ids.figure_center.figure = self.plots[0].fig
            self.ids.figure_center.touch_mode = "pointselect"
        elif len(self.plots)==2:
            figure_center=self.ids.figure_center
            self.ids.figure_left.figure = self.plots[0].fig
            self.ids.figure_left.touch_mode = "pointselect"
            self.ids.figure_right.figure = self.plots[1].fig
            self.ids.figure_right.touch_mode = "pointselect"
        else:
            self.ids.figure_left.figure = self.plots[0].fig
            self.ids.figure_left.touch_mode = "pointselect"
            self.ids.figure_center.figure = self.plots[1].fig
            self.ids.figure_center.touch_mode = "pointselect"
            self.ids.figure_right.figure = self.plots[2].fig
            self.ids.figure_right.touch_mode = "pointselect"

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
        self.ids.component_list.add_widget(button)

    def remove_component(self):
        #find active button
        button_to_remove = next( (t for t in ToggleButton.get_widgets('components') if t.state=='down'), None)
        #remove item
        if button_to_remove:
            if button_to_remove.text in self.components:
                self.components.remove(button_to_remove.text)
            self.ids.component_list.remove_widget(button_to_remove)

    def choose_component(self,x,y):
        #check what plot is currently displayed
        for plot in self.plots:
            if self.ids.figure_center.figure == plot.fig:
                current_plot=plot
        fig_comps=current_plot.components
        final_dist=0
        for ind,comp in enumerate(fig_comps):
            x_ellipse,y_ellipse=comp[0].get_center()
            #calculate distance between click and component center
            dist=(x-x_ellipse)**2+(y-y_ellipse)**2
            if dist<final_dist or final_dist==0:
                final_ind=ind
                final_comp=comp
                final_dist=dist

        if final_comp[1][0].get_color() == "red":
            final_comp[0].set_color("black")
            final_comp[1][0].set_color("black")
            final_comp[2][0].set_color("black")
        else:
            final_comp[0].set_color("red")
            final_comp[1][0].set_color("red")
            final_comp[2][0].set_color("red")

        for comp in fig_comps:
            if comp != final_comp and comp[1][0].get_color() == "red":
                comp[0].set_color("black")
                comp[1][0].set_color("black")
                comp[2][0].set_color("black")

        #current_plot.fig.canvas.draw()
        #current_plot.fig.canvas.flush_events()





class VCAT(App):

    def build(self):
        self.screen=ModelFits()
        return self.screen

    def set_touch_mode(self,touch_mode):
        self.screen.ids.figure_center.touch_mode = touch_mode

    def home_plot(self):
        self.screen.ids.figure_center.home()

    def add_component(self):
        self.screen.add_component()

    def remove_component(self):
        self.screen.remove_component()

    def choose_component(self,x,y):
        self.screen.choose_component(selfx,y)


if __name__ == "__main__":
    VCAT().run()