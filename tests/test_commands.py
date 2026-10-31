"""The installed command-line tools, run on the synthetic data in a temporary folder."""
import re

import numpy as np
import pytest

from gateway_mt.scripts import autogen_report, f2l, plot_dataset, plot_edi, to_vtk
import gateway_mt.IO as IO

from conftest import synthetic

DATA = synthetic('synthetic_modem.dat')
RESP = synthetic('synthetic_modem_resp.dat')
LIST = synthetic('edis', 'synthetic.lst')
MODEL = synthetic('synthetic.model')


# --------------------------------------------------------------------------- #
# f2l
# --------------------------------------------------------------------------- #
def test_f2l_lists_every_station(edi_folder):
    assert f2l.main(['-d', str(edi_folder)]) == 0
    names = IO.read_sites(str(edi_folder / 'all_sites.lst'))
    assert len(names) == 30 and names[0] == 'SYN001'
    assert all('.' not in n for n in names)          # written without extensions


def test_f2l_filter_and_no_overwrite(edi_folder):
    assert f2l.main(['syn00', '-d', str(edi_folder)]) == 0
    assert IO.read_sites(str(edi_folder / 'syn00_sites.lst')) == ['SYN00{}'.format(i) for i in range(1, 10)]
    assert f2l.main(['syn00', '-d', str(edi_folder)]) == 1                 # exists: refuses
    assert f2l.main(['syn00', '-d', str(edi_folder), '--force']) == 0
    assert f2l.main(['nothing_matches', '-d', str(edi_folder)]) == 1


# --------------------------------------------------------------------------- #
# plot_edi
# --------------------------------------------------------------------------- #
def test_plot_edi_saves_one_figure_per_file(edi_folder, tmp_path):
    out = tmp_path / 'figures'
    files = [str(edi_folder / 'SYN001.edi'), str(edi_folder / 'SYN030.edi')]
    assert plot_edi.main(files + ['--save', '-o', str(out), '--format', 'pdf']) == 0
    assert sorted(p.name for p in out.iterdir()) == ['SYN001.pdf', 'SYN030.pdf']


def test_plot_edi_reports_missing_files(tmp_path, capsys):
    assert plot_edi.main([str(tmp_path / 'nope.edi'), '--save', '-o', str(tmp_path)]) == 1
    assert 'not found' in capsys.readouterr().out


# --------------------------------------------------------------------------- #
# plot_dataset
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize('inputs', [['-l', LIST], ['-d', DATA], ['-r', RESP],
                                    ['-d', DATA, '-r', RESP], ['-l', LIST, '-d', DATA, '-r', RESP]],
                         ids=['raw', 'data', 'response', 'data+response', 'all'])
def test_plot_dataset_any_combination(tmp_path, inputs):
    out = tmp_path / 'plots.pdf'
    assert plot_dataset.main(inputs + ['--sites', 'SYN001', 'SYN030', '-o', str(out)]) == 0
    assert out.stat().st_size > 10000


def test_plot_dataset_options(tmp_path, capsys):
    out = tmp_path / 'z.pdf'
    assert plot_dataset.main(['-d', DATA, '--plot', 'impedance', '--sites', 'SYN005', 'SYN030',
                              '-o', str(out)]) == 0
    assert '1 page' in capsys.readouterr().out
    assert plot_dataset.main(['-d', DATA, '-o', str(out)]) == 1                    # exists
    assert plot_dataset.main(['-d', DATA, '--sites', 'NOPE', '-o', str(tmp_path / 'x.pdf')]) == 1


# --------------------------------------------------------------------------- #
# autogen_report
# --------------------------------------------------------------------------- #
def test_autogen_report(tmp_path):
    out = tmp_path / 'Report'
    assert autogen_report.main(['-d', DATA, '-r', RESP, '-l', LIST, '-m', MODEL, '-o', str(out)]) == 0
    summary = (out / 'summary.txt').read_text()
    rms = float(re.search(r'Final RMS: ([\d.]+)', summary).group(1))
    assert rms == pytest.approx(1.0, abs=0.1)
    assert 'Inversion type: 5 (full impedance, tipper)' in summary
    assert '32 x 32 x 35 cells' in summary
    assert (out / 'misfit_curves.pdf').stat().st_size > 10000
    # an existing report is only replaced with --force
    assert autogen_report.main(['-d', DATA, '-r', RESP, '-o', str(out)]) == 1


# --------------------------------------------------------------------------- #
# to_vtk
# --------------------------------------------------------------------------- #
def _vtk_point_count(path):
    match = re.search(r'POINTS (\d+)', path.read_text(errors='replace'))
    return int(match.group(1))


def test_to_vtk_model_and_sites(tmp_path):
    base = tmp_path / 'survey'
    assert to_vtk.main(['-o', str(base), '-m', MODEL, '-l', LIST, '--utm', '16U']) == 0
    assert _vtk_point_count(tmp_path / 'survey_sites.vtk') == 30
    assert (tmp_path / 'survey_model.vtk').stat().st_size > 100000


def test_to_vtk_sea_level_and_flat_sites(tmp_path):
    base = tmp_path / 'flat'
    assert to_vtk.main(['-o', str(base), '-d', DATA, '--origin', '0', '0',
                        '--no-elevation', '--sea-level', '50']) == 0
    lines = (tmp_path / 'flat_sites.vtk').read_text().splitlines()
    start = next(i for i, l in enumerate(lines) if l.startswith('POINTS')) + 1
    heights = {float(l.split()[2]) for l in lines[start:start + 30]}
    assert heights == {50.0}


def test_to_vtk_trims_the_model_and_its_resolution(tmp_path):
    """--trim cuts each side in x / y and the bottom in z, from the resolution as well.

    The model is used as its own resolution, which is written inverted, so
    the two fields multiply to 1 only if they line up.
    """
    base = tmp_path / 'trimmed'
    assert to_vtk.main(['-o', str(base), '-m', MODEL, '--origin', '0', '0',
                        '--trim', '3', '2', '4', '--resolution', MODEL]) == 0
    text = (tmp_path / 'trimmed_model.vtk').read_text(errors='replace')
    assert 'DIMENSIONS 29 27 32' in text          # (ny, nx, nz) + 1 for a 32 x 32 x 35 model trimmed to 26 x 28 x 31

    def field(name):
        return np.array(text.split('SCALARS ' + name)[1].split('SCALARS')[0].split()[3:], dtype=float)

    np.testing.assert_allclose(field('Resistivity') * field('Raw_Resolution'), 1, rtol=1e-4)


def test_to_vtk_rejects_gofem(tmp_path):
    assert to_vtk.main(['-o', str(tmp_path / 'x'), '-d', synthetic('synthetic_gofem.gdat'),
                        '--origin', '0', '0']) == 1
