"""Helpers shared by the GUIs."""
from PyQt5 import QtCore


def check_key_presses(key_presses, verbose=False):
    """Decode the keyboard modifiers held during a click.

    In Data Plot's error multiplier tree, Shift / Ctrl / Alt extend an edit
    to all periods / components / sites; Mesh Designer uses Ctrl.

    Parameters
    ----------
    key_presses : Qt.KeyboardModifiers
    verbose : bool, optional
        Print what will be edited.

    Returns
    -------
    dict
        ``{'Control': bool, 'Shift': bool, 'Alt': bool}``
    """
    msg = ''
    retval = {'Control': False, 'Shift': False, 'Alt': False}
    if key_presses == QtCore.Qt.ShiftModifier:
        retval['Shift'] = True
    elif key_presses == QtCore.Qt.ControlModifier:
        retval['Control'] = True
    elif key_presses == QtCore.Qt.AltModifier:
        retval['Alt'] = True
    elif key_presses == (QtCore.Qt.ControlModifier |
                         QtCore.Qt.ShiftModifier):
        retval['Shift'] = True
        retval['Control'] = True
    elif key_presses == (QtCore.Qt.AltModifier |
                         QtCore.Qt.ShiftModifier):
        retval['Alt'] = True
        retval['Shift'] = True
    elif key_presses == (QtCore.Qt.ControlModifier |
                         QtCore.Qt.AltModifier):
        retval['Control'] = True
        retval['Alt'] = True
    elif key_presses == (QtCore.Qt.ControlModifier |
                         QtCore.Qt.ShiftModifier |
                         QtCore.Qt.AltModifier):
        retval['Shift'] = True
        retval['Alt'] = True
        retval['Control'] = True
    if any(retval.values()):
        msg = '+'.join([k for k in retval if retval[k]])
        msg = msg.replace('Shift', 'Periods')
        msg = msg.replace('Control', 'Components')
        msg = msg.replace('Alt', 'Sites')
    if verbose and msg:
        print(' '.join(['Editting all', msg]))
    return retval



