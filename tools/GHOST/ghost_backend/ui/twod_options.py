"""Advanced 2-D monostatic controls; numerical policy stays in runs.controls."""
from PySide6.QtCore import Signal
from PySide6.QtWidgets import QCheckBox, QComboBox, QFormLayout, QWidget

from ghost_backend.runs.controls import FACTORIZATION_LABELS, run_control_options


class TwoDOptionsWidget(QWidget):
    changed = Signal()

    def __init__(self, parent=None):
        super().__init__(parent)
        form = QFormLayout(self)
        form.setContentsMargins(0, 0, 0, 0)
        self.factorization = QComboBox()
        for key, label in FACTORIZATION_LABELS.items():
            self.factorization.addItem(label, key)
        self.factorization.setToolTip(
            'Adaptive chooses dense or compressed for the problem and memory budget. '
            'Dense can use hierarchical factorization internally, with LU fallback. '
            'Hierarchical and Compressed force the corresponding backend.')
        form.addRow('Factorization', self.factorization)
        self.mesh = QComboBox()
        for key in ('adaptive', 'local', 'global'):
            self.mesh.addItem(key.capitalize(), key)
        self.mesh.setToolTip(
            'Adaptive chooses and checks polynomial meshes when certification is on. '
            'Local and Global retain manual polynomial-order control.')
        form.addRow('Mesh strategy', self.mesh)
        self.order = QComboBox()
        self.order.addItem('Automatic', 0)
        for degree in (1, 2, 3):
            self.order.addItem('P{}'.format(degree), degree)
        self.order.setToolTip('Manual order applies to Local or Global meshes. Adaptive controls it automatically.')
        form.addRow('Polynomial order', self.order)
        self.grading = QCheckBox('Grade corners, open ends, and material junctions')
        self.grading.setChecked(True)
        form.addRow(self.grading)
        self._manual_order = 1
        self.mesh.currentIndexChanged.connect(self._sync)
        self.factorization.currentIndexChanged.connect(lambda *_: self.changed.emit())
        self.order.currentIndexChanged.connect(lambda *_: self.changed.emit())
        self.grading.toggled.connect(lambda *_: self.changed.emit())
        self._sync()

    def _sync(self, *_):
        adaptive = self.mesh.currentData() == 'adaptive'
        current = self.order.currentData()
        if adaptive and current:
            self._manual_order = current
        blocked = self.order.blockSignals(True)
        self.order.setCurrentIndex(self.order.findData(0 if adaptive else (current or self._manual_order)))
        # Automatic is display-only: fixed meshes must select a real degree.
        self.order.model().item(0).setEnabled(adaptive)
        self.order.blockSignals(blocked)
        self.order.setEnabled(not adaptive)
        self.changed.emit()

    def value(self):
        return dict(factorization=self.factorization.currentData(), mesh_strategy=self.mesh.currentData(),
                    basis_order=self.order.currentData() or 1, mesh_grading=self.grading.isChecked())

    def set_value(self, controls):
        # Accept the legacy factorization alias without displaying a fifth mode.
        value = run_control_options(controls)
        self.factorization.setCurrentIndex(self.factorization.findData(value['factorization']))
        self.mesh.setCurrentIndex(self.mesh.findData(value['mesh_strategy']))
        self._manual_order = value['basis_order']
        self.order.setCurrentIndex(self.order.findData(value['basis_order']))
        self.grading.setChecked(value['mesh_grading'])
        self._sync()
