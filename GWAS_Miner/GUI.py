import os
import sys
from io import BytesIO

from PyQt5 import uic
from PyQt5.QtCore import QRunnable, pyqtSlot, QObject, pyqtSignal, QThreadPool, Qt, QSize
from PyQt5.QtGui import QPixmap, QIcon
from PyQt5.QtSvg import QGraphicsSvgItem
from PyQt5.QtWidgets import QApplication, QFileDialog, QPushButton, QGraphicsScene, QTableWidgetItem, QListWidgetItem, \
    QBoxLayout, QLabel, QHeaderView, QCheckBox, QHBoxLayout, QWidget
from reportlab.graphics import renderPM
from svglib.svglib import svg2rlg

from GWAS_Miner import GWASMiner
from functools import partial


class MainForm:
    def __init__(self):
        self.Form, self.Window = uic.loadUiType("GWAS_Miner/res/gwas_miner.ui")
        QApplication.setAttribute(Qt.AA_EnableHighDpiScaling, True)
        self.app = QApplication([])
        self.window = self.Window()
        self.form = self.Form()
        self.form.setupUi(self.window)
        self.center()
        # Set widget attributes/visibility
        self.__setup_interface()
        # Assign event handlers
        self.__add_handlers()
        self.threadpool = QThreadPool()
        self.is_cancelled = False
        self.is_running = False
        self.worker = None
        self.run_worker(GWASMiner.load_nlp_object, None, self.initial_loading_finished_callback, True)
        self.dependency_svgs = []
        self.dependency_index = 0
        self.previous_page = 0

    def __load_style(self, theme):
        with open(F'GWAS_Miner/res/{theme}') as file:
            self.window.setStyleSheet(file.read())

    def __setup_interface(self):
        """
        Initialise user interface elements
        """
        theme = "light_theme.qss"
        if GWASMiner.config.get(section="preferences", option="theme").lower() == "dark":
            theme = "dark_theme.qss"
        self.__load_style(theme)
        self.form.stackedWidget.setGeometry(0, 0, self.window.width(), self.window.height())
        self.navigate_to_page(0)
        self.form.status_lbl.setHidden(True)
        self.set_loading_visible(True)
        splash_loading_scene = QGraphicsScene()
        main_loading_scene = QGraphicsScene()
        splash_loading_svg = QGraphicsSvgItem("GWAS_Miner/res/loading.svg")
        splash_loading_svg.setScale(0.4)
        main_loading_svg = QGraphicsSvgItem("GWAS_Miner/res/loading.svg")
        main_loading_svg.setScale(0.4)
        splash_loading_scene.addItem(splash_loading_svg)
        main_loading_scene.addItem(main_loading_svg)
        self.form.test_loading_graphicsview.setScene(splash_loading_scene)
        self.form.loading_svg.setScene(main_loading_scene)

    def center(self):
        frame_gm = self.window.frameGeometry()
        screen = QApplication.desktop().screenNumber(
            QApplication.desktop().cursor().pos())
        center_point = QApplication.desktop().screenGeometry(screen).center()
        frame_gm.moveCenter(center_point)
        self.window.move(frame_gm.topLeft())

    def __add_handlers(self):
        """
        Assign event handlers for interactive elements.
        """
        self.form.study_directory_btn.clicked.connect(self.get_file_directory)
        self.form.run_nlp_btn.clicked.connect(self.run_nlp_btn_handler)
        self.form.visualise_back_btn.clicked.connect(self.visualise_back_btn_handler)
        self.form.dependency_next_btn.clicked.connect(self.get_next_dependency)
        self.form.dependency_previous_btn.clicked.connect(self.get_previous_dependency)
        self.form.update_ontology_cache_action.triggered.connect(self.update_ontology_cache_handler)
        self.form.quit_action.triggered.connect(self.quit_action_handler)
        self.form.settings_action.triggered.connect(self.settings_action_handler)
        self.form.settings_save_btn.clicked.connect(self.settings_save_handler)
        self.form.settings_cancel_btn.clicked.connect(self.settings_cancel_handler)
        self.form.result_view_back_btn.clicked.connect(self.result_back_btn_click_handler)
        self.form.file_select_all_checkbox.stateChanged.connect(self.set_study_checkbox_states)

    def set_study_checkbox_states(self):
        state = self.form.file_select_all_checkbox.isChecked()
        table = self.form.study_file_tablewidget
        for i in range(table.rowCount()):
            table.cellWidget(i, 2).findChild(QCheckBox).setChecked(state)

    def result_back_btn_click_handler(self):
        self.navigate_to_page(1)
        self.form.result_viewer_textbrowser.clear()

    def view_result(self, file):
        with open(F"output/{file}_result.json", "r", encoding="utf-8") as f_in:
            self.form.result_viewer_textbrowser.setText(f_in.read())
        self.navigate_to_page(4)

    def settings_cancel_handler(self):
        self.navigate_to_page(self.previous_page)

    def settings_save_handler(self):
        previous_theme = GWASMiner.config.get("preferences", "theme")
        new_theme = self.form.theme_combobox.currentText()
        GWASMiner.config.set("preferences", "theme", new_theme)
        if previous_theme != new_theme:
            self.__load_style(F"{new_theme.lower()}_theme.qss")
        GWASMiner.save_config()
        self.navigate_to_page(self.previous_page)

    def settings_action_handler(self):
        if GWASMiner.config.get("preferences", "theme").lower() == "light":
            self.form.theme_combobox.setCurrentText("Light")
        else:
            self.form.theme_combobox.setCurrentText("Dark")
        self.navigate_to_page(3)

    @staticmethod
    def quit_action_handler():
        sys.exit("Quitting...")

    def set_splash_loading_text(self, text):
        self.form.loading_label.setText(text)

    def update_ontology_cache_handler(self):
        from GWAS_Miner import Ontology
        self.set_splash_loading_text("Updating ontology data...")
        self.run_worker(Ontology.update_ontology_cache, None, self.ontology_updated_callback, disable_controls=True)

    def ontology_updated_callback(self, response=None):
        """
        Return to the Main frame.
        """
        self.navigate_to_page(1)

    def get_next_dependency(self):
        """
        Display the next dependency visualisation.
        """
        if self.dependency_index < len(self.dependency_svgs):
            self.dependency_index += 1
            self.render_dependency_svg(self.dependency_index)

    def get_previous_dependency(self):
        """
        Display the previous dependency visualisation.
        """
        if self.dependency_index > 0:
            self.dependency_index -= 1
            self.render_dependency_svg(self.dependency_index)

    def navigate_to_page(self, page_number):
        """
        Sets the current index of the stack widget
        0 = Loading, 1 = main, 2 = study visualisation, 3 = settings
        """
        self.previous_page = self.form.stackedWidget.currentIndex()
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
        study = GWASMiner.prepare_study(self.form.study_directory_input.text(), file)
        self.set_loading_visible(True)
        self.worker = Worker(GWASMiner.get_study_visualisations, [study])
        self.worker.signals.progress_update.connect(self.set_progress_text)
        self.worker.signals.error.connect(self.set_progress_text)
        self.worker.signals.finished.connect(self.visualisation_finished_callback)
        self.threadpool.start(self.worker)

    def set_loading_visible(self, state):
        if not state:
            if not self.is_running:
                self.form.loading_svg.setHidden(True)
        else:
            self.form.loading_svg.setHidden(False)

    def get_file_directory(self):
        """
        Opens a folder selection dialog and adds each file within the selected folder
        to the study list widget.
        """
        path = QFileDialog.getExistingDirectory(caption='Select Study Directory')
        if not path:
            return
        self.form.study_directory_input.setText(path)
        self.form.study_file_tablewidget.setRowCount(0)
        if self.validate_directory(path):
            self.form.file_select_all_checkbox.setChecked(True)
            for file in os.listdir(self.form.study_directory_input.text()):
                if file.endswith("_maintext.json"):
                    # Add analysis button for the study.
                    analyse_btn = QPushButton()
                    analyse_icon = QIcon("GWAS_Miner/res/icons/statistics.png")
                    analyse_btn.setIcon(analyse_icon)
                    analyse_btn.clicked.connect(partial(self.render_study_visualisation, file))  # partial is essential
                    analyse_widget = QWidget()
                    analyse_layout = QHBoxLayout(analyse_widget)
                    analyse_layout.setAlignment(Qt.AlignCenter)
                    analyse_layout.addWidget(analyse_btn)
                    analyse_layout.setContentsMargins(0, 0, 0, 0)
                    # Add study label text.
                    study = QTableWidgetItem()
                    study.setTextAlignment(Qt.AlignCenter)
                    study.setData(Qt.DisplayRole, file)
                    # Add checkbox for study record.
                    check_box_widget = QWidget()
                    check_box_layout = QHBoxLayout(check_box_widget)
                    check_box = QCheckBox()
                    check_box.setChecked(True)
                    check_box_layout.addWidget(check_box)
                    check_box_layout.setAlignment(Qt.AlignCenter)
                    check_box_layout.setContentsMargins(0, 0, 0, 0)

                    row_index = self.form.study_file_tablewidget.rowCount()

                    self.form.study_file_tablewidget.insertRow(row_index)
                    self.form.study_file_tablewidget.setItem(row_index, 0, study)
                    self.form.study_file_tablewidget.setCellWidget(row_index, 1, analyse_widget)
                    self.form.study_file_tablewidget.setCellWidget(row_index, 2, check_box_widget)

                    header = self.form.study_file_tablewidget.horizontalHeader()
                    analyse_btn.setFixedWidth(30)
                    header.resizeSection(0, 180)
                    header.resizeSection(1, 40)
                    header.resizeSection(2, 30)

    def validate_directory(self, path):
        """
        Check the input directory exists.
        """
        try:
            valid = os.listdir(path)
            if not path or not valid:
                self.set_progress_text("Please select a valid directory.")
                return False
        except FileNotFoundError:
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
            self.set_loading_visible(True)
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
            GWASMiner.is_cancelled = True
            self.form.run_nlp_btn.setEnabled(False)
            self.set_progress_text("Cancelling...")
            self.form.run_nlp_btn.setText("Please \nWait...")
        else:
            if not self.validate_directory(self.form.study_directory_input.text()):
                return
            shortlist = self.get_selected_studies()
            if not shortlist:
                self.set_progress_text("No files selected.")
                return
            self.form.result_file_tablewidget.setRowCount(0)
            self.form.results_failed_listwidget.clear()
            self.form.result_tab_widget.setTabText(0, F"Succeeded")
            self.form.result_tab_widget.setTabText(1, F"No Results")
            self.is_running = True
            GWASMiner.is_cancelled = False
            self.form.run_nlp_btn.setText("Stop \nProcessing")
            self.run_worker(GWASMiner.process_studies, (self.form.study_directory_input.text(), None, shortlist),
                            self.update_results_files, False)

    def get_selected_studies(self):
        shortlist = []
        table = self.form.study_file_tablewidget
        for i in range(table.rowCount()):
            if table.cellWidget(i, 2).findChild(QCheckBox).isChecked():
                shortlist.append(table.item(i, 0).text())
        return shortlist

    @staticmethod
    def convert_svg_textpath(svg):
        import re
        skip = True
        coords = []
        result = svg
        # TODO Reduce the height of the SVG and scale all Y values to compensate
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
            result = re.sub(r"(<textPath[^>]+.)",
                            F"<tspan class=\"displacy-label\" fill=\"currentColor\" x=\"{x - 10}\" y=\"{y - 8}\">",
                            result, count=1)
        result = result.replace("</textPath>", "</tspan>")
        return result

    def convert_svg_to_png(self, svg):
        return_val = None
        rlg = svg2rlg(BytesIO(bytes(self.convert_svg_textpath(svg), encoding="utf-8")))
        if GWASMiner.theme() == "light":
            return_val = renderPM.drawToString(rlg, fmt="PNG")
        else:
            return_val = renderPM.drawToString(rlg, fmt="PNG", bg=0x19232D)
        return return_val

    def render_dependency_svg(self, index):
        qp = QPixmap()
        qp.loadFromData(self.convert_svg_to_png(self.dependency_svgs[index]))
        self.form.dependency_image_label.setPixmap(qp)
        self.update_dependency_index_text()
        self.form.dependency_scrollarea.horizontalScrollBar().setValue(
            self.form.dependency_scrollarea.horizontalScrollBar().maximum() / 2)
        self.form.dependency_scrollarea.verticalScrollBar().setValue(
            self.form.dependency_scrollarea.verticalScrollBar().maximum())

    def update_dependency_index_text(self):
        self.form.dependency_index_label.setText(F"{self.dependency_index + 1} of {len(self.dependency_svgs) + 1}")

    def reformat_html(self, html):
        result = html.replace("line-height: 2.5;", "line-height: 1.15;").replace("-apple-system", "calibri")
        result = result.replace("padding: 0.45em 0.6em;", "padding: 30px;")
        result = result.replace("border-radius: 0.35em;", "border-radius: 5px;")
        return result

    def visualisation_finished_callback(self, response):
        entity_html = self.reformat_html(response.data[0])
        self.form.entity_visualisation_browser.setHtml(entity_html)
        self.dependency_svgs.clear()
        self.dependency_index = 0
        self.form.dependency_image_label.setStyleSheet("background-color: transparent;")
        self.dependency_svgs = response.data[1]
        self.render_dependency_svg(0)
        self.form.visualise_stats_table.setRowCount(0)
        # Always disable sorting before changing contents to avoid missing data bug.
        self.form.visualise_stats_table.setSortingEnabled(False)
        for i in range(len(response.data[2])):
            phenotype = list(response.data[2].keys())[i]
            ontology = response.data[2][phenotype]["Ontology"]
            count = response.data[2][phenotype]["Count"]

            phenotype_cell = QTableWidgetItem()
            count_cell = QTableWidgetItem()
            ontology_cell = QTableWidgetItem()

            phenotype_cell.setTextAlignment(Qt.AlignCenter)
            count_cell.setTextAlignment(Qt.AlignCenter)
            ontology_cell.setTextAlignment(Qt.AlignCenter)

            phenotype_cell.setData(Qt.DisplayRole, phenotype)
            count_cell.setData(Qt.DisplayRole, count)
            ontology_cell.setData(Qt.DisplayRole, ontology)

            self.form.visualise_stats_table.insertRow(i)
            self.form.visualise_stats_table.setItem(i, 0, phenotype_cell)
            self.form.visualise_stats_table.setItem(i, 1, count_cell)
            self.form.visualise_stats_table.setItem(i, 2, ontology_cell)

        self.form.visualise_stats_table.resizeColumnsToContents()
        self.form.visualise_stats_table.setSortingEnabled(True)
        self.form.visualise_stats_table.sortByColumn(1, Qt.SortOrder(1))
        self.navigate_to_page(2)
        self.set_loading_visible(False)
        self.form.status_lbl.setText("Ready")

    def initial_loading_finished_callback(self, status):
        self.set_loading_visible(False)
        self.toggle_controls(True)
        self.form.run_nlp_btn.setText("Begin \nProcessing")
        self.set_progress_text("")
        self.navigate_to_page(1)

    def set_progress_text(self, progress):
        self.form.status_lbl.setText(progress)
        self.form.status_lbl.setHidden(False)

    def update_results_files(self, result):
        if result.data == 1:
            self.set_progress_text(result.text)
            self.is_running = False
            self.set_loading_visible(False)
            GWASMiner.is_cancelled = False
            self.form.run_nlp_btn.setText("Start \nProcessing")
            self.form.run_nlp_btn.setEnabled(True)
            return
        if result.status:
            self.form.result_tab_widget.setTabText(0, F"Succeeded ({self.form.result_file_tablewidget.rowCount() + 1})")

            view_btn = QPushButton()
            view_icon = QIcon("GWAS_Miner/res/icons/view_icon.png")
            view_btn.setIcon(view_icon)
            view_btn.clicked.connect(partial(self.view_result, result.text))
            view_widget = QWidget()
            view_layout = QHBoxLayout(view_widget)
            view_layout.setAlignment(Qt.AlignCenter)
            view_layout.setContentsMargins(0, 0, 0, 0)
            view_layout.addWidget(view_btn)

            study = QTableWidgetItem()
            study.setTextAlignment(Qt.AlignCenter)

            study.setData(Qt.DisplayRole, result.text)

            row_index = self.form.result_file_tablewidget.rowCount()

            self.form.result_file_tablewidget.insertRow(row_index)
            self.form.result_file_tablewidget.setItem(row_index, 0, study)
            self.form.result_file_tablewidget.setCellWidget(row_index, 1, view_widget)

            header = self.form.result_file_tablewidget.horizontalHeader()
            view_btn.setFixedWidth(30)
            header.resizeSection(0, 200)
            header.resizeSection(1, 30)
        else:
            self.form.results_failed_listwidget.addItem(result.text)
            self.form.result_tab_widget.setTabText(1, F"No Results ({self.form.results_failed_listwidget.count()})")

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
    def __init__(self, status, text, data=None):
        self.status = status  # Bool
        self.text = text  # String
        self.data = data


class Worker(QRunnable):
    """
    Worker thread

    Inherits from QRunnable to handler worker thread setup, signals and wrap-up.

    :param fn: The function callback to run on this worker thread. Supplied args and
                     kwargs will be passed through to the runner.
    :type fn: function
    :param args: Arguments to pass to the callback function

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
            import traceback
            traceback.print_exc()
            self.signals.error.emit(traceback.format_exc())
