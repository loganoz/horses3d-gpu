import importlib
import sys
import tempfile
import types
import unittest
from pathlib import Path
from unittest import mock


REPO_ROOT = Path(__file__).resolve().parents[1]
UTILS_DIR = REPO_ROOT / "utils" / "PythonUtilities"


def _install_stub_modules():
    matplotlib = types.ModuleType("matplotlib")
    matplotlib.use = lambda *_args, **_kwargs: None

    pyplot = types.ModuleType("matplotlib.pyplot")
    pyplot.figure = lambda *_args, **_kwargs: None
    pyplot.ion = lambda *_args, **_kwargs: None
    pyplot.pause = lambda *_args, **_kwargs: None
    pyplot.subplot = lambda *_args, **_kwargs: None
    pyplot.legend = lambda *_args, **_kwargs: None
    pyplot.rcParams = {}

    gridspec = types.ModuleType("matplotlib.gridspec")
    numpy = types.ModuleType("numpy")

    sys.modules["matplotlib"] = matplotlib
    sys.modules["matplotlib.pyplot"] = pyplot
    sys.modules["matplotlib.gridspec"] = gridspec
    sys.modules["numpy"] = numpy


def _import_module(name):
    if str(UTILS_DIR) not in sys.path:
        sys.path.insert(0, str(UTILS_DIR))
    _install_stub_modules()
    sys.modules.pop(name, None)
    return importlib.import_module(name)


class _FakeAxis:
    def __init__(self):
        self.plots = []

    def plot(self, _x, y, *_args, **_kwargs):
        self.plots.append(list(y))

    def get_legend_handles_labels(self):
        return [], []

    def legend(self, *_args, **_kwargs):
        return None


class TestMonitorsUtilities(unittest.TestCase):
    def _write_temp_file(self, content):
        tmp = tempfile.NamedTemporaryFile("w", delete=False)
        self.addCleanup(lambda: Path(tmp.name).unlink(missing_ok=True))
        tmp.write(content)
        tmp.close()
        return tmp.name

    def test_get_which_monitor_type_detects_residuals_and_scalar(self):
        monitors = _import_module("Monitors")

        self.assertEqual(
            monitors.getWhichMonitorType("case.residuals"), monitors.RESIDUALS_FILE
        )
        self.assertEqual(
            monitors.getWhichMonitorType("case.monitor"), monitors.MONITOR_FILE
        )

    def test_gather_new_residuals_values_reads_all_columns(self):
        residuals = _import_module("ResidualsMonitor")
        content = (
            "header one\n"
            "header two\n"
            "1 0.1 0 1.1 2.1 3.1 4.1 5.1\n"
            "2 0.2 0 1.2 2.2 3.2 4.2 5.2\n"
        )
        path = self._write_temp_file(content)

        counter, time, continuity, x_momentum, y_momentum, z_momentum, energy = (
            residuals.GatherNewResidualsValues(path, 0)
        )
        self.assertEqual(counter, 2)
        self.assertEqual(time, [0.1, 0.2])
        self.assertEqual(continuity, [1.1, 1.2])
        self.assertEqual(x_momentum, [2.1, 2.2])
        self.assertEqual(y_momentum, [3.1, 3.2])
        self.assertEqual(z_momentum, [4.1, 4.2])
        self.assertEqual(energy, [5.1, 5.2])

    def test_gather_new_residuals_values_respects_skip_data(self):
        residuals = _import_module("ResidualsMonitor")
        content = (
            "header one\n"
            "header two\n"
            "1 0.1 0 1.1 2.1 3.1 4.1 5.1\n"
            "2 0.2 0 1.2 2.2 3.2 4.2 5.2\n"
            "3 0.3 0 1.3 2.3 3.3 4.3 5.3\n"
        )
        path = self._write_temp_file(content)

        counter, time, *_ = residuals.GatherNewResidualsValues(path, 2)
        self.assertEqual(counter, 1)
        self.assertEqual(time, [0.3])

    def test_gather_new_scalar_values_parses_metadata_and_data(self):
        scalar = _import_module("ScalarMonitor")
        content = (
            "Monitor name LiftMonitor\n"
            "irrelevant heading\n"
            "Selected variable: kinetic energy rate\n"
            "Iteration Time Value\n"
            "1 0.1 4.0\n"
            "2 0.2 3.0\n"
        )
        path = self._write_temp_file(content)

        counter, time, value, monitor_name, variable = scalar.GatherNewScalarValues(
            path, 0
        )
        self.assertEqual(counter, 2)
        self.assertEqual(time, [0.1, 0.2])
        self.assertEqual(value, [4.0, 3.0])
        self.assertEqual(monitor_name, "LiftMonitor")
        self.assertIn("kinetic energy rate", variable)

    def test_gather_new_scalar_values_skip_data_works_in_python3(self):
        scalar = _import_module("ScalarMonitor")
        content = (
            "Monitor name LiftMonitor\n"
            "Selected variable: drag\n"
            "Iteration Time Value\n"
            "1 0.1 4.0\n"
            "2 0.2 3.0\n"
            "3 0.3 2.0\n"
        )
        path = self._write_temp_file(content)

        counter, time, value, *_ = scalar.GatherNewScalarValues(path, 2)
        self.assertEqual(counter, 1)
        self.assertEqual(time, [0.3])
        self.assertEqual(value, [2.0])

    def test_new_scalar_plot_inverts_kinetic_energy_rate(self):
        scalar = _import_module("ScalarMonitor")
        ax = _FakeAxis()

        scalar.NewScalarPlot([0.1, 0.2], [3.0, 4.0], ax, "m", "kinetic energy rate")

        self.assertEqual(ax.plots[0], [-3.0, -4.0])

    def test_update_scalar_plot_uses_new_on_first_pass(self):
        scalar = _import_module("ScalarMonitor")
        ax = object()
        with mock.patch.object(
            scalar,
            "GatherNewScalarValues",
            return_value=(1, [1.0], [2.0], "monitor", "var"),
        ), mock.patch.object(scalar, "NewScalarPlot") as new_plot, mock.patch.object(
            scalar, "AppendToScalarPlot"
        ) as append_plot:
            count = scalar.UpdateScalarPlot("unused", 0, ax)

        self.assertEqual(count, 1)
        new_plot.assert_called_once()
        append_plot.assert_not_called()

    def test_update_residuals_plot_uses_append_after_first_pass(self):
        residuals = _import_module("ResidualsMonitor")
        ax = object()
        with mock.patch.object(
            residuals,
            "GatherNewResidualsValues",
            return_value=(2, [1.0], [2.0], [3.0], [4.0], [5.0], [6.0]),
        ), mock.patch.object(residuals, "NewResidualsPlot") as new_plot, mock.patch.object(
            residuals, "AppendToResidualsPlot"
        ) as append_plot:
            count = residuals.UpdateResidualsPlot("unused", 1, ax)

        self.assertEqual(count, 2)
        new_plot.assert_not_called()
        append_plot.assert_called_once()


if __name__ == "__main__":
    unittest.main()
