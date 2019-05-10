import sys, os
from pyforms.basewidget import BaseWidget
from pyforms.controls import ControlFile, ControlSlider, ControlButton, ControlText, ControlDir, ControlCombo
from AnyQt.QtWidgets import QFileDialog


# Fix a bug of ControlDir in some version of pyforms
class ControlDirFix(ControlDir):
    def click(self):
        value = QFileDialog.getExistingDirectory(self.parent, self._label)

        if value and len(value) > 0:
            self.value = value


class HelpWindow(BaseWidget):
    def __init__(self):
        BaseWidget.__init__(self)

        # Text Help
        # close button


class ProgressWindows(BaseWidget):
    def __init__(self):
        BaseWidget.__init__(self)
        # stop button
        # Text progress



class TREEasyWindow(BaseWidget):

    def __init__(self):
        super(TREEasyWindow, self).__init__('TREEasy')

        # Forms
        self._directory = ControlDirFix()
        self._directory.label = 'Gene/Locus File Directory'
        self._species_file = ControlFile('Species Name File')
        self._gene_file = ControlFile("Gene Name File")
        self._type = ControlCombo("Gene/Locus Type")

        self._root = ControlText("Root Taxon(s)")

        self._type.add_item('CDS')
        self._type.add_item('nonCDS')

        self._bootstrap = ControlSlider("Bootstrap Value")
        self._bootstrap.value = 50
        self._bootstrap.min = 0
        self._bootstrap.max = 100
        self._cross = ControlText("Cross Value")
        self._cross.value = "3"
        self._network_number = ControlSlider("Maximum Network Number")
        self._network_number.value = 3
        self._network_number.min = 0
        self._network_number.max = 5
        self._threads_number = ControlSlider("Threads Number")
        self._threads_number.value = 8
        self._threads_number.min = 1
        self._threads_number.max = 32

        self._start = ControlButton('Start')
        self._start.value = self.__process_start

        self._formset = [
            '_directory',
            '_species_file',
            '_gene_file',
            '_root',
            ('_cross', '_network_number'),
            '_bootstrap',
            ('_type', '_threads_number', '_start')
        ]

        self.mainmenu = [
            {
                'File': [
                    {'Example': self.__example}
                ]
            },
            {
                'Help': [
                    {'Help': self.__help_event},
                    {'About': self.__about_event}
                ]
            }
        ]

    def __process_start(self):
        seq_file = self._directory.value
        typename = self._type.value
        roottaxon = self._root.value
        species_namefile = self._species_file.value
        gene_namefile = self._gene_file.value
        boot_value = self._bootstrap.value
        cross_value = self._cross.value
        thread_number = self._threads_number.value
        network_number = self._network_number.value
        # disable start button
        self._start.enabled = False
        os.system('python TREEasy.py -d %s -c %s -s %s -g %s -b %s -r %s -n %s -k %s -t %s' %
                  (seq_file, typename, species_namefile, gene_namefile, boot_value, roottaxon, network_number,
                   cross_value, thread_number))
        self._start.enabled = True
        pass

    def __help_event(self):
        pass

    def __about_event(self):
        pass

    def __example(self):
        pass


if __name__ == '__main__':
    from pyforms import start_app

    start_app(TREEasyWindow, geometry=[100, 100, 600, 480])
