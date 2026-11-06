"""Front window with a button for each main program (the ``mt_launcher`` command).

Each program is started as its own process through its normal command, so
nothing here imports them. A program that isn't installed in this environment
(mtinv and mtplot aren't part of gateway_mt) just has its button greyed out.
"""
import shutil
import subprocess
import sys
import sysconfig

from PyQt5 import QtWidgets

#: The buttons: label, command and a short description.
PROGRAMS = (('Gateway', 'gateway_mt', 'Projects, plus Data Plot, Model Viewer and Mesh Designer'),
            ('mtinv', 'mtinv', 'Inversion setup'),
            ('mtplot', 'mtplot', 'Model plotting'))


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
    """A window with one button per program.

    Closing it leaves any programs it started running.

    Parameters
    ----------
    programs : sequence of tuple, optional
        ``(label, command, description)`` for each button.

    Attributes
    ----------
    buttons : dict
        ``{label: QPushButton}``.
    processes : list of subprocess.Popen
        The programs started so far.
    """

    def __init__(self, programs=PROGRAMS):
        super(Launcher, self).__init__()
        self.setWindowTitle('MT Launcher')
        self.buttons = {}
        self.processes = []
        layout = QtWidgets.QVBoxLayout(self)
        for label, command, description in programs:
            button = QtWidgets.QPushButton(label)
            path = find_command(command)
            button.setEnabled(path is not None)
            button.setToolTip(description if path else "'{}' isn't installed in this environment".format(command))
            button.clicked.connect(lambda checked=False, label=label, path=path: self.start(label, path))
            layout.addWidget(button)
            self.buttons[label] = button

    def start(self, label, path):
        """Start a program, or show a warning if it can't be started."""
        try:
            self.processes.append(subprocess.Popen([path]))
        except OSError as e:
            QtWidgets.QMessageBox.warning(self, 'Launcher', 'Could not start {}: {}'.format(label, e))


def main():
    """Console entry point (``mt_launcher``)."""
    app = QtWidgets.QApplication(sys.argv)
    window = Launcher()
    window.show()
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()
