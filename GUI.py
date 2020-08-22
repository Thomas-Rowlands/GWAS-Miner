import os

from PyQt5 import uic
from PyQt5.QtWidgets import QApplication, QPushButton, QFileDialog, QWidget
import Phenotype_Finder
from rtgo import ReadyThready

class MainForm:
    def __init__(self):
        self.Form, self.Window = uic.loadUiType("main_form.ui")
        self.app = QApplication([])
        self.window = self.Window()
        self.form = self.Form()
        self.form.setupUi(self.window)
        # Set widget attributes/visibility
        self.setup_interface()
        # Assign event handlers
        self.add_handlers()

    def setup_interface(self):
        self.form.status_lbl.setHidden(True)
        self.form.running_progress_bar.setHidden(True)

    def add_handlers(self):
        self.form.study_directory_btn.clicked.connect(self.get_file_directory)
        self.form.run_nlp_btn.clicked.connect(self.begin_processing)

    def get_file_directory(self):
        self.form.study_directory_input.setText(QFileDialog.getExistingDirectory(caption='Select Folder'))
        self.form.study_file_listview.clear()
        for file in os.listdir(self.form.study_directory_input.text()):
            if file.endswith(".json"):
                self.form.study_file_listview.addItem(file)

    def begin_processing(self):
        self.form.running_progress_bar.setHidden(False)
        Phenotype_Finder.process_studies(self.form.study_directory_input.text())

    def update_progress_text(self, text):
        self.form.status_lbl.setText(text)
        self.form.status_lbl.setHidden(False)

    def open(self):
        self.window.show()
        self.app.exec_()
