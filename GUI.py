import os
import sys
import traceback

from PyQt5 import uic
from PyQt5.QtCore import QRunnable, pyqtSlot, QObject, pyqtSignal, QThreadPool, Qt
from PyQt5.QtWidgets import QApplication, QPushButton, QFileDialog, QWidget
import Phenotype_Finder


class MainForm:
    def __init__(self):
        self.Form, self.Window = uic.loadUiType("main_form.ui")
        self.app = QApplication([])
        self.app.setAttribute(Qt.AA_EnableHighDpiScaling, True)
        self.window = self.Window()
        self.form = self.Form()
        self.form.setupUi(self.window)
        # Set widget attributes/visibility
        self.setup_interface()
        # Assign event handlers
        self.add_handlers()
        self.threadpool = QThreadPool()

    def setup_interface(self):
        self.form.status_lbl.setHidden(True)
        self.form.running_progress_bar.setHidden(True)

    def add_handlers(self):
        self.form.study_directory_btn.clicked.connect(self.get_file_directory)
        self.form.run_nlp_btn.clicked.connect(self.begin_processing)

    def get_file_directory(self):
        self.form.study_directory_input.setText(QFileDialog.getExistingDirectory(caption='Select Study Directory'))
        self.form.study_file_listwidget.clear()
        for file in os.listdir(self.form.study_directory_input.text()):
            if file.endswith(".json"):
                self.form.study_file_listwidget.addItem(file)

    def begin_processing(self):
        worker = Worker(Phenotype_Finder.process_studies, (self.form.study_directory_input.text(), None))
        worker.signals.progress_update.connect(self.set_progress_text)
        worker.signals.error.connect(self.set_progress_text)
        worker.signals.finished.connect(self.update_results_files)
        self.form.running_progress_bar.setHidden(False)
        self.threadpool.start(worker)

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
    result = pyqtSignal(str)
    progress_update = pyqtSignal(str)  # Message


class QtFinishedResponse:
    def __init__(self, status, text):
        self.status = status
        self.text = text


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
            result = self.fn(*self.args, self.progress_callback, self.finished_callback)
            self.signals.result.emit(result)
        except:
            traceback.print_exc()
            exctype, value = sys.exc_info()[:2]
            self.signals.error.emit(traceback.format_exc())
