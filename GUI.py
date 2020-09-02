import os
import re
import sys
import traceback
from io import BytesIO

from PyQt5 import uic
from PyQt5.QtCore import QRunnable, pyqtSlot, QObject, pyqtSignal, QThreadPool, Qt
from PyQt5.QtGui import QPixmap
from PyQt5.QtSvg import QGraphicsSvgItem, QSvgWidget
from PyQt5.QtWidgets import QApplication, QFileDialog, QPushButton, QGraphicsScene, QLabel, QFrame

import Phenotype_Finder

from svglib.svglib import svg2rlg
from reportlab.graphics import renderPM


class MainForm:
    def __init__(self):
        self.Form, self.Window = uic.loadUiType("res/gwas_miner.ui")
        self.app = QApplication([])
        self.app.setAttribute(Qt.AA_EnableHighDpiScaling, True)
        self.window = self.Window()
        self.form = self.Form()
        self.form.setupUi(self.window)
        # Set widget attributes/visibility
        self.__setup_interface()
        # Assign event handlers
        self.__add_handlers()
        self.threadpool = QThreadPool()
        self.is_cancelled = False
        self.is_running = False
        self.worker = None
        self.run_worker(Phenotype_Finder.load_nlp_object, None, self.study_processing_finished_callback, True)
        self.dependency_svgs = []
        self.dependency_index = 0

    def __setup_interface(self):
        """
        Initialise user interface elements
        """
        self.form.stackedWidget.setGeometry(0, 0, self.window.width(), self.window.height())
        self.navigate_to_page(0)
        self.form.status_lbl.setHidden(True)
        self.form.loading_svg.setHidden(True)
        splash_loading_scene = QGraphicsScene()
        main_loading_scene = QGraphicsScene()
        splash_loading_svg = QGraphicsSvgItem("res/loading.svg")
        splash_loading_svg.setScale(0.4)
        main_loading_svg = QGraphicsSvgItem("res/loading.svg")
        main_loading_svg.setScale(0.4)
        splash_loading_scene.addItem(splash_loading_svg)
        main_loading_scene.addItem(main_loading_svg)
        self.form.test_loading_graphicsview.setScene(splash_loading_scene)
        self.form.loading_svg.setScene(main_loading_scene)

    def __add_handlers(self):
        """
        Assign event handlers for interactive elements.
        """
        self.form.study_directory_btn.clicked.connect(self.get_file_directory)
        self.form.run_nlp_btn.clicked.connect(self.run_nlp_btn_handler)
        self.form.visualise_back_btn.clicked.connect(self.visualise_back_btn_handler)
        self.form.analyse_study_btn.clicked.connect(self.analyse_study_btn_handler)
        self.form.dependency_next_btn.clicked.connect(self.get_next_dependency)
        self.form.dependency_previous_btn.clicked.connect(self.get_previous_dependency)

    def get_next_dependency(self):
        if self.dependency_index < len(self.dependency_svgs):
            self.dependency_index += 1
            self.render_dependency_svg(self.dependency_index)

    def get_previous_dependency(self):
        if self.dependency_index > 0:
            self.dependency_index -= 1
            self.render_dependency_svg(self.dependency_index)

    def navigate_to_page(self, page_number):
        """
        Sets the current index of the stack widget
        """
        self.form.stackedWidget.setCurrentIndex(page_number)

    def visualise_back_btn_handler(self):
        """
        Event handler for visualise back button.
        Navigates to the main page.
        """
        self.navigate_to_page(1)

    def render_study_visualisation(self, file, disable_controls=False):
        """
        Begin rendering study visualisations using Displacy on a worker thread.
        """
        if disable_controls:
            self.toggle_controls(False)
        study = Phenotype_Finder.prepare_study(self.form.study_directory_input.text(), file)
        self.form.loading_svg.setHidden(False)
        self.worker = Worker(Phenotype_Finder.get_study_visualisations, [study])
        self.worker.signals.progress_update.connect(self.set_progress_text)
        self.worker.signals.error.connect(self.set_progress_text)
        self.worker.signals.finished.connect(self.visualisation_finished_callback)
        self.threadpool.start(self.worker)

    def analyse_study_btn_handler(self):
        """
        Event handler for analyse study button.
        Begins the study visualisation rendering process.
        """
        if self.form.study_file_listwidget.selectedItems():
            self.render_study_visualisation(self.form.study_file_listwidget.selectedItems()[0].text())

    def get_file_directory(self):
        """
        Opens a folder selection dialog and adds each file within the selected folder
        to the study list widget.
        """
        path = QFileDialog.getExistingDirectory(caption='Select Study Directory')
        self.form.study_directory_input.setText(path)
        self.form.study_file_listwidget.clear()
        if self.validate_directory(path):
            for file in os.listdir(self.form.study_directory_input.text()):
                if file.endswith("_maintext.json"):
                    self.form.study_file_listwidget.addItem(file)

    def validate_directory(self, path):
        """
        Check the input directory exists.
        """
        try:
            valid = os.listdir(path)
            if not path or not valid:
                self.set_progress_text("Please select a valid directory.")
                return False
        except FileNotFoundError as f:
            self.set_progress_text("Invalid directory.")
            return False
        return True

    def toggle_controls(self, state):
        """
        Sets the enabled state of every push button in the application.
        """
        for widget in self.form.__dict__.items():
            if isinstance(widget[1], QPushButton):
                widget[1].setEnabled(state)

    def run_worker(self, func, args, finished_callback, disable_controls=False):
        """
        Executes a provided function on a worker thread.
        """
        if disable_controls:
            self.toggle_controls(False)
            self.navigate_to_page(0)
        else:
            self.form.loading_svg.setHidden(False)
        self.worker = Worker(func, args)
        self.worker.signals.progress_update.connect(self.set_progress_text)
        self.worker.signals.error.connect(self.set_progress_text)
        self.worker.signals.finished.connect(finished_callback)
        self.threadpool.start(self.worker)

    def run_nlp_btn_handler(self):
        """
        Begin processing the files identified within the study list widget.
        """
        if self.is_running:
            self.is_running = False
            Phenotype_Finder.is_cancelled = True
            self.form.run_nlp_btn.setEnabled(False)
            self.set_progress_text("Cancelling...")
            self.form.run_nlp_btn.setText("Please \nWait...")
        else:
            if not self.validate_directory(self.form.study_directory_input.text()):
                return
            self.is_running = True
            Phenotype_Finder.is_cancelled = False
            self.form.run_nlp_btn.setText("Stop \nProcessing")
            self.run_worker(Phenotype_Finder.process_studies, (self.form.study_directory_input.text(), None),
                            self.study_processing_finished_callback, True)

    def convert_svg_textpath(self, svg):
        skip = True
        coords = []
        result = svg
        #TODO Reduce the height of the SVG and scale all Y values to compensate
        # Get the starting coordinate of the relationship link
        for match in re.findall(r'(?<=d="M)([0-9,. ]+)', result):
            skip = not skip
            if skip:
                continue
            coords.append(match.split(" "))
        for coord in coords:
            line_start_coord = coord[1].split(",")
            line_stop_coord = coord[2].split(",")
            x = (float(line_start_coord[0]) + float(line_stop_coord[0])) / 2
            y = (float(line_start_coord[1]) + float(line_stop_coord[1])) / 2
            result = re.sub(r"(<textPath[^>]+.)", F"<tspan class=\"displacy-label\" fill=\"currentColor\" x=\"{x - 10}\" y=\"{y - 8}\">", result, count=1)
        result = result.replace("</textPath>", "</tspan>")
        # result = re.sub(r"(height=\"[0-9\.]+\")", "height=\"300\"", result, count=1)
        # result = re.sub(r"(height: [0-9\.]+px)", "height: 300px", result, count=1)
        return result


    def convert_svg_to_png(self, svg):
        return_val = None
        rlg = svg2rlg(BytesIO(bytes(self.convert_svg_textpath(svg), encoding="utf-8")))
        return_val = renderPM.drawToString(rlg, fmt="PNG")
        return return_val

    def render_dependency_svg(self, index):
        qp = QPixmap()
        qp.loadFromData(self.convert_svg_to_png(self.dependency_svgs[index]))
        self.form.dependency_image_label.setPixmap(qp)
        self.update_dependency_index_text()
        self.form.dependency_scrollarea.horizontalScrollBar().setValue(self.form.dependency_scrollarea.horizontalScrollBar().maximum() / 2)
        self.form.dependency_scrollarea.verticalScrollBar().setValue(self.form.dependency_scrollarea.verticalScrollBar().maximum())

    def update_dependency_index_text(self):
        self.form.dependency_index_label.setText(F"{self.dependency_index + 1} of {len(self.dependency_svgs) + 1}")

    def reformat_html(self, html):
        result = html.replace("line-height: 2.5;", "line-height: 1.15;").replace("-apple-system", "calibri")
        result = result.replace("padding: 0.45em 0.6em;", "padding: 30px;")
        result = result.replace("border-radius: 0.35em;", "border-radius: 5px;")
        return result

    def visualisation_finished_callback(self, response):
        entity_html = self.reformat_html(response.html[0])
        self.form.entity_visualisation_browser.setHtml(entity_html)
        self.dependency_svgs.clear()
        self.dependency_index = 0
        self.dependency_svgs = response.html[1]
        self.render_dependency_svg(0)
        self.navigate_to_page(2)
        self.form.loading_svg.hide()
        self.form.status_lbl.setText("Ready")

    def study_processing_finished_callback(self, status):
        if self.is_cancelled:
            self.form.status_lbl.setText("Cancelled.")
        else:
            self.form.status_lbl.setText("Complete.")
        self.form.loading_svg.setHidden(True)
        self.toggle_controls(True)
        self.form.run_nlp_btn.setText("Begin \nProcessing")
        self.navigate_to_page(1)

    def set_progress_text(self, progress):
        self.form.status_lbl.setText(progress)
        self.form.status_lbl.setHidden(False)

    def update_results_files(self, result):
        if result.status:
            self.form.result_file_listwidget.addItem(result.text)
            self.form.result_tab_widget.setTabText(0, F"Succeeded ({self.form.result_file_listwidget.count()})")
        else:
            self.form.results_failed_listwidget.addItem(result.text)
            self.form.result_tab_widget.setTabText(1, F"Failed ({self.form.results_failed_listwidget.count()})")

    def open(self):
        self.window.show()
        self.app.exec_()


class WorkerSignals(QObject):
    """
    Defines the signals available from a running worker thread.

    """
    finished = pyqtSignal(object)
    error = pyqtSignal(str)
    progress_update = pyqtSignal(str)  # Message


class QtFinishedResponse:
    def __init__(self, status, text, html=None):
        self.status = status  # Bool
        self.text = text  # String
        self.html = html


class Worker(QRunnable):
    """
    Worker thread

    Inherits from QRunnable to handler worker thread setup, signals and wrap-up.

    :param callback: The function callback to run on this worker thread. Supplied args and
                     kwargs will be passed through to the runner.
    :type callback: function
    :param args: Arguments to pass to the callback function
    :param kwargs: Keywords to pass to the callback function

    """

    def __init__(self, fn, args):
        super(Worker, self).__init__()
        # Store constructor arguments (re-used for processing)
        self.fn = fn
        self.args = args
        self.signals = WorkerSignals()
        self.progress_callback = self.signals.progress_update
        self.finished_callback = self.signals.finished

    @pyqtSlot()
    def run(self):
        """
        Initialise the runner function with passed args, kwargs.
        """

        # Retrieve args/kwargs here; and fire processing using them
        try:
            if not self.args:
                self.fn(self.progress_callback, self.finished_callback)
            else:
                self.fn(*self.args, self.progress_callback, self.finished_callback)
        except:
            traceback.print_exc()
            exctype, value = sys.exc_info()[:2]
            self.signals.error.emit(traceback.format_exc())
