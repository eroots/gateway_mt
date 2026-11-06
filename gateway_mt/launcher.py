"""Front window with a button for each main program (the ``mt_launcher`` command).

Each program is started as its own process through its normal command, so
nothing here imports them. A program that isn't installed in this environment
(mtinv and mtplot aren't part of gateway_mt) just has its button greyed out.

The programs each keep their own kind of project file. The window works on a
project folder: each button opens that program's project from the folder, and
starts the program there.
"""
import glob
import os
import shutil
import subprocess
import sys
import sysconfig

from PyQt5 import QtWidgets

#: The buttons: label, command, project files, the option that starts a blank
#: project (Gateway has none, it asks for a project itself) and a short description.
PROGRAMS = (('Gateway', 'gateway_mt', '*.pymt', None, 'Projects, plus Data Plot, Model Viewer and Mesh Designer'),
            ('mtinv', 'mtinv', '*.mtinv.json', '--new', 'Inversion setup'),
            ('mtplot', 'mtplot', '*.mtplot.json', '--new', 'Model plotting'))


def find_command(command):
    """Full path of ``command``, looking in this environment first.

    That finds the environment's programs even if it wasn't activated.

    Returns
    -------
    str or None
        None if it isn't found.
    """
    return shutil.which(command, path=sysconfig.get_path('scripts')) or shutil.which(command)


class Launcher(QtWidgets.QWidget):
    """A window with a project folder and one button per program.

    Closing it leaves any programs it started running.

    Parameters
    ----------
    folder : str, optional
        Project folder (default: the current folder).
    programs : sequence of tuple, optional
        ``(label, command, project files, new project option, description)``
        for each button.

    Attributes
    ----------
    folder : str
        The project folder.
    buttons : dict
        ``{label: QPushButton}``.
    processes : list of subprocess.Popen
        The programs started so far.
    """

    def __init__(self, folder=None, programs=PROGRAMS):
        super(Launcher, self).__init__()
        self.setWindowTitle('MT Launcher')
        self.folder = os.path.abspath(folder or os.getcwd())
        self.programs = programs
        self.buttons = {}
        self.processes = []
        layout = QtWidgets.QVBoxLayout(self)
        self.folderLabel = QtWidgets.QLabel()
        self.folderLabel.setWordWrap(True)
        layout.addWidget(self.folderLabel)
        choose = QtWidgets.QPushButton('Project Folder...')
        choose.clicked.connect(self.choose_folder)
        layout.addWidget(choose)
        for label, command, pattern, new_option, description in programs:
            button = QtWidgets.QPushButton()
            path = find_command(command)
            button.setEnabled(path is not None)
            button.setToolTip(description if path else "'{}' isn't installed in this environment".format(command))
            button.clicked.connect(lambda checked=False, label=label, path=path, pattern=pattern, new=new_option:
                                   self.start(label, path, pattern, new))
            layout.addWidget(button)
            self.buttons[label] = button
        self.show_folder()

    def project_files(self, pattern):
        """The files matching ``pattern`` in the project folder, sorted."""
        return sorted(glob.glob(os.path.join(glob.escape(self.folder), pattern)))

    def show_folder(self):
        """Show the folder, and on each button what it would open."""
        self.folderLabel.setText(self.folder)
        for label, command, pattern, new_option, description in self.programs:
            files = self.project_files(pattern)
            if len(files) == 1:
                project = os.path.basename(files[0])
            else:
                project = '{} projects'.format(len(files)) if files else 'new project'
            self.buttons[label].setText('{}: {}'.format(label, project))

    def choose_folder(self):
        """Ask for a new project folder."""
        folder = QtWidgets.QFileDialog.getExistingDirectory(self, 'Project Folder', self.folder)
        if folder:
            self.folder = os.path.abspath(folder)
            self.show_folder()

    def start(self, label, path, pattern, new_option=None):
        """Start a program on its project from the folder.

        If there are several projects, asks which to open. If there are none,
        ``new_option`` (if any) starts a blank project rather than whichever
        one the program would reopen. Shows a warning if the program can't be
        started.
        """
        files = self.project_files(pattern)
        if len(files) > 1:
            name, ok = QtWidgets.QInputDialog.getItem(self, label, 'Project:', [os.path.basename(f) for f in files],
                                                      0, False)
            if not ok:
                return
            files = [os.path.join(self.folder, name)]
        arguments = files or ([new_option] if new_option else [])
        try:
            self.processes.append(subprocess.Popen([path] + arguments, cwd=self.folder))
        except OSError as e:
            QtWidgets.QMessageBox.warning(self, 'Launcher', 'Could not start {}: {}'.format(label, e))


def main():
    """Console entry point (``mt_launcher [folder]``)."""
    app = QtWidgets.QApplication(sys.argv)
    window = Launcher(sys.argv[1] if len(sys.argv) > 1 else None)
    window.show()
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()
