"""Interactive native-resolution image cuts, published as independent datasets."""
from __future__ import annotations

import numpy as np
from PySide6.QtCore import Qt, QItemSelection, QItemSelectionModel
from PySide6.QtWidgets import QWidget, QVBoxLayout, QLabel, QPushButton, QFileDialog, QInputDialog
from matplotlib.widgets import RectangleSelector, PolygonSelector
from matplotlib.backend_bases import MouseEvent

from GRIM_Backend.ui.responsive_controls import ResponsiveButtonPanel
from GRIM_Backend.ui.no_wheel import NoWheelComboBox
from GRIM_Backend.execution.plot_jobs import PlotJobRunner
from GRIM_Backend.isar.masking import (finish_cut, comparison_cuts,
                                     MaskImageRequest, prepare_mask_source)


class ClickPolygonSelector(PolygonSelector):
    """Place and close vertices even when Qt coalesces mouse-move events."""

    def _release(self, event):
        if not self._selection_completed and self._active_handle_idx < 0:
            self._onmove(event)
        super()._release(event)
        self._draw_polygon()


class MaskCutPanel(QWidget):
    def __init__(self, owner):
        super().__init__(owner.plot_canvas.parentWidget())
        self.owner = owner
        self.setObjectName("maskCutPanel")
        self.figure = owner.plot_figure
        self.canvas = owner.plot_canvas
        self.ax = None
        self._pick_connection = None
        self._first_point = None
        self.source_name = ""
        self.image = None
        self.vertices = None
        self.selector = None
        self._mask_definition = None
        self._shapes = {}
        self.runner = PlotJobRunner(self)
        self.runner.publication_failed.connect(self._publication_error)
        layout = QVBoxLayout(self)
        row = ResponsiveButtonPanel()
        self.shape = NoWheelComboBox()
        self.shape.addItems(["Polygon", "Rectangle"])
        self.operation = NoWheelComboBox()
        self.operation.addItems(["Cut Include", "Cut Exclude"])
        self.operation.setToolTip("Include reconstructs the selected region. Exclude subtracts that contribution from the source data and retains any fit residual.")
        row.addWidget(self.shape)
        row.addWidget(self.operation)
        layout.addWidget(row)
        self.instructions = QLabel("Draw on the plot: click polygon vertices, then click the first vertex to finish. "
                                   "Finish Cut adds the result to the dataset table.")
        self.instructions.setWordWrap(True)
        layout.addWidget(self.instructions)
        self.details = QLabel()
        self.details.setWordWrap(True)
        self.details.setTextFormat(Qt.PlainText)
        layout.addWidget(self.details)
        self.message = QLabel()
        self.message.hide()
        self.message.setWordWrap(True)
        self.message.setTextFormat(Qt.PlainText)
        layout.addWidget(self.message)
        self.finish = QPushButton("Finish Cut")
        self.finish.setEnabled(False)
        self.comparison = QPushButton("Create comparison set")
        self.comparison.setToolTip("Add both Include and Exclude datasets to the dataset table using the same boundary.")
        self.comparison.setEnabled(False)
        self.close_button = QPushButton("Close")
        self.load_button = QPushButton("Load Mask…")
        self.save_button = QPushButton("Save Mask…")
        self.load_button.setEnabled(False)
        self.save_button.setEnabled(False)
        self.load_button.setToolTip("Load an exact mask boundary, or import a single XY overlay polygon.")
        self.save_button.setToolTip("Save the boundary in physical coordinates, with units, frame, and Include/Exclude.")
        row.addWidget(self.load_button)
        row.addWidget(self.save_button)
        row.addWidget(self.finish)
        row.addWidget(self.comparison)
        row.addWidget(self.close_button)
        self.shape.currentIndexChanged.connect(self._set_selector)
        self.finish.clicked.connect(self._finish)
        self.comparison.clicked.connect(self._finish_comparison)
        self.close_button.clicked.connect(self.deactivate)
        self.load_button.clicked.connect(self._load_mask)
        self.save_button.clicked.connect(self._save_mask)
        self.hide()

    def _set_details(self, text):
        self.details.setText(text)
        self.details.setVisible(bool(text))

    def _set_message(self, text):
        self.message.setText(text)
        self.message.setVisible(bool(text))

    def activate(self):
        self._set_message("")
        self.show()
        self._load_plot_source()

    def _disconnect_selector(self):
        if self.selector is not None:
            self.selector.set_visible(False)
            self.selector.disconnect_events()
            for artist in self.selector.artists:
                if artist.axes is not None:
                    artist.remove()
            self.selector = None
        if self._pick_connection is not None:
            self.canvas.mpl_disconnect(self._pick_connection)
            self._pick_connection = None

    def deactivate(self):
        self.runner.cancel()
        self._disconnect_selector()
        self.canvas.unsetCursor()
        self.image = None
        self.vertices = None
        self._mask_definition = None
        self.ax = None
        self._shapes.clear()
        self._first_point = None
        self._busy(False)
        self.hide()
        self.canvas.draw_idle()

    def _pick_source(self, event):
        if event.button != 1 or event.inaxes is None:
            return
        if any(getattr(artist, '_grim_mask_source', None) is not None
               for artist in (*event.inaxes.images, *event.inaxes.collections)):
            import weakref
            self.figure._grim_mask_active_axes = weakref.ref(event.inaxes)
            self._first_point = (event.xdata, event.ydata)
            self._load_plot_source()

    def _load_plot_source(self):
        self.runner.cancel()
        self._disconnect_selector()
        self.ax = None
        self.image = None
        self.vertices = None
        self._mask_definition = None
        self._shapes.clear()
        figure = getattr(self.owner, "plot_figure", None)
        candidates = []
        if figure is not None:
            for axes in figure.axes:
                for artist in (*axes.images, *axes.collections):
                    if getattr(artist, "_grim_mask_source", None) is not None:
                        candidates.append(artist)
        if len(candidates) > 1:
            selected = getattr(figure, "_grim_mask_active_axes", None)
            selected = selected() if selected is not None else None
            candidates = [artist for artist in candidates if artist.axes is selected]
            if len(candidates) != 1:
                self._show_source()
                self._set_details("Click the image panel you want to cut, then draw your boundary on that plot.")
                self._pick_connection = self.canvas.mpl_connect('button_press_event', self._pick_source)
                return
        if candidates:
            artist = candidates[0]
            self.ax = artist.axes
            self.image = artist._grim_mask_source
            display = (self.image.display if isinstance(self.image, MaskImageRequest)
                       else self.image.metadata.setdefault('display', {}))
            display.update(cmap=artist.get_cmap().name, limits=list(artist.get_clim()))
        self._show_source()

    def _current(self):
        return self.image

    def _show_source(self, *_args):
        image = self._current()
        if image is None:
            self._set_details("Select a dataset in the table, plot ISAR, then click Mask Cut.")
            self.vertices = None
            self.finish.setEnabled(False)
            self.comparison.setEnabled(False)
            self.load_button.setEnabled(False)
            self.save_button.setEnabled(False)
            self.canvas.draw_idle()
            return
        if not isinstance(image, MaskImageRequest) and not self.runner.busy:
            self._source_ready(image, None)
            return
        self._busy(True)
        self._set_details("Preparing full-resolution cut source…")
        self.runner.submit(
            lambda cancel: prepare_mask_source(image, cancel) if isinstance(image, MaskImageRequest) else image,
            self._source_ready,
        )

    def _source_ready(self, result, error):
        self._busy(False)
        if error is not None:
            self.image = None
            self.vertices = None
            self.finish.setEnabled(False)
            self.comparison.setEnabled(False)
            self.load_button.setEnabled(False)
            self.save_button.setEnabled(False)
            self._set_details(f"Could not prepare image: {error}")
            return
        image = result
        self.image = image
        if image.inverse is None or 'native_samples' not in image.inverse or image.complex_image is None:
            self.image = None
            self.vertices = None
            self.finish.setEnabled(False)
            self.comparison.setEnabled(False)
            self.load_button.setEnabled(False)
            self.save_button.setEnabled(False)
            self._set_details("This image has no retained complex reconstruction source. "
                              "Plot ISAR again, then click Mask Cut. Composite magnitude images cannot produce radar cuts.")
            return
        from GRIM_Backend.ui.dataset_actions import DATASET_ID_ROLE
        source_id = image.metadata.get('source_dataset_id')
        self.source_name = image.metadata.get('source_dataset_name') or image.name
        if source_id:
            for row in range(self.owner.table.rowCount()):
                item = self.owner.table.item(row, 0)
                if item is not None and item.data(DATASET_ID_ROLE) == source_id:
                    self.source_name = item.text()
                    break
        self._set_details("")
        self.setToolTip(image.reconstruction_reason)
        self._set_selector()
        self.canvas.setCursor(Qt.CrossCursor)
        details = "Mask Cut is ready — click on the image to draw the boundary."
        details += " Cuts fit the native samples; Exclude retains any fit residual."
        self._set_details(details)
        self.load_button.setEnabled(True)
        if self._first_point is not None and self.shape.currentText() == 'Polygon':
            x, y = self.ax.transData.transform(self._first_point)
            self._first_point = None
            event = MouseEvent('button_press_event', self.canvas, x, y, button=1)
            self.selector.press(event)
            event = MouseEvent('button_release_event', self.canvas, x, y, button=1)
            self.selector.release(event)
        self.canvas.draw_idle()

    def _set_selector(self, *_args):
        self._disconnect_selector()
        self.vertices = None
        self._mask_definition = None
        self.save_button.setEnabled(False)
        self.finish.setEnabled(False)
        self.comparison.setEnabled(False)
        if self._current() is None or isinstance(self._current(), MaskImageRequest) or self.runner.busy:
            return
        if self.shape.currentText() == "Polygon":
            self.selector = ClickPolygonSelector(
                self.ax, self._polygon, useblit=False,
                props={'color': '#ffb000', 'linewidth': 2.0, 'zorder': 100},
                handle_props={'markerfacecolor': '#ffb000', 'markeredgecolor': '#111111',
                              'markersize': 7, 'zorder': 101},
            )
        else:
            self.selector = RectangleSelector(
                self.ax, self._rectangle, useblit=False, button=[1], interactive=True,
                props={'facecolor': '#ffb000', 'edgecolor': '#ffb000', 'alpha': .3,
                       'linewidth': 2.0, 'zorder': 100},
            )
        saved = self._shapes.get((id(self._current()), self.shape.currentText()))
        if saved is not None:
            self.vertices = list(saved)
            if self.shape.currentText() == 'Polygon':
                self.selector.verts = self.vertices
            else:
                points = np.asarray(self.vertices)
                self.selector.extents = (points[:,0].min(), points[:,0].max(), points[:,1].min(), points[:,1].max())
            self.finish.setEnabled(True)
            self.comparison.setEnabled(True)
            self.save_button.setEnabled(True)
        self.canvas.draw_idle()

    def _polygon(self, vertices):
        self.vertices = list(vertices)
        self._mask_definition = None
        self._shapes[(id(self._current()), self.shape.currentText())] = list(vertices)
        self.finish.setEnabled(len(vertices) >= 3 and not self.runner.busy)
        self.comparison.setEnabled(len(vertices) >= 3 and not self.runner.busy)
        self.save_button.setEnabled(len(vertices) >= 3 and not self.runner.busy)

    def _save_mask(self):
        from pathlib import Path
        from GRIM_Backend.isar.mask_geometry import save_mask
        if self.image is None or self.vertices is None or self.runner.busy:
            return
        filename, _ = QFileDialog.getSaveFileName(
            self, "Save reusable mask", self.source_name + '_mask.csv', 'Mask / overlay files (*.csv)')
        if not filename:
            return
        if not Path(filename).suffix:
            filename += '.csv'
        try:
            self._mask_definition = save_mask(filename, self.image, self.vertices,
                include=self.operation.currentIndex() == 0, definition=self._mask_definition)
            self._set_message('Mask saved. Reuse it on datasets with the same physical target origin.')
        except (OSError, ValueError, KeyError, TypeError) as exc:
            self._set_message(f'Could not save mask: {exc}')

    def _load_mask(self):
        from GRIM_Backend.isar.mask_geometry import load_mask, MissingMaskMetadata
        from GRIM_Backend.plotting.overlay_data import LENGTH_UNITS
        if self.image is None or self.runner.busy:
            return
        filename, _ = QFileDialog.getOpenFileName(self, 'Load reusable mask', '',
            'Mask / overlay files (*.csv *.xy *.txt *.dat *.pts);;All files (*)')
        if not filename:
            return
        try:
            try:
                vertices, definition = load_mask(filename, self.image)
            except MissingMaskMetadata:
                units = list(LENGTH_UNITS)
                current = self.image.metadata.get('unit', 'm')
                unit, accepted = QInputDialog.getItem(self, 'Overlay mask units',
                    'This overlay has no mask frame. XY will use the current plot frame. File units:',
                    units, units.index(current) if current in units else 0, False)
                if not accepted:
                    return
                vertices, definition = load_mask(filename, self.image, legacy_unit=unit,
                    include=self.operation.currentIndex() == 0)
            # A rectangle rotated into a different look frame is a polygon.
            self.shape.setCurrentText('Polygon')
            self._set_selector()
            self.selector.verts = vertices.tolist()
            self._polygon(vertices)
            self._mask_definition = definition
            self.operation.setCurrentIndex(0 if definition['operation'] == 'include' else 1)
            self._set_message('Mask loaded in physical coordinates. The target origin must match; different grids select different pixel centers.')
            self.canvas.draw_idle()
        except (OSError, UnicodeError, ValueError, KeyError, TypeError) as exc:
            self._set_message(f'Could not load mask: {exc}')

    def _rectangle(self, first, last):
        if None in (first.xdata, first.ydata, last.xdata, last.ydata):
            return
        self._polygon([(first.xdata, first.ydata), (last.xdata, first.ydata),
                       (last.xdata, last.ydata), (first.xdata, last.ydata)])

    def _busy(self, busy):
        for widget in (self.shape, self.operation):
            widget.setEnabled(not busy)
        ready = not busy and self.vertices is not None and len(self.vertices) >= 3
        self.finish.setEnabled(ready)
        self.comparison.setEnabled(ready)
        self.save_button.setEnabled(ready)
        self.load_button.setEnabled(not busy and self.image is not None and self.ax is not None)
        if self.selector is not None:
            self.selector.set_active(not busy)

    def _publication_error(self, error):
        self._busy(False)
        self._set_message(f"Could not publish Mask Cut result: {error}")

    def _finish(self):
        image = self._current()
        if image is None or self.vertices is None:
            return
        vertices = np.array(self.vertices, copy=True)
        include = self.operation.currentIndex() == 0
        name = self.source_name + ("_include" if include else "_exclude")
        self._busy(True)
        self._set_message("Creating the cut dataset…")
        definition = self._mask_definition
        self.runner.submit(lambda cancel: finish_cut(image, vertices, include=include, name=name, cancel=cancel,
                                                    mask_definition=definition),
                           self._publish_cut)

    def _publish_cut(self, result, error):
        self._busy(False)
        if error is not None:
            self._set_message(f"Mask Cut failed: {error}")
            return
        image, radar = result
        if radar is None:
            self._set_message("No dataset was created. This image does not support radar reconstruction.")
            return
        self._export_results([(image, radar)], replot=True)

    def _export_results(self, results, *, replot=False):
        from GRIM_Backend.ui.dataset_actions import DATASET_ID_ROLE
        identifiers = [self.owner._add_dataset_row(
            radar, image.name, "Mask Cut: reconstructed radar data", notify=False,
        ) for image, radar in results]
        self.deactivate()
        table = self.owner.table
        selection = QItemSelection()
        rows = [row for row in range(table.rowCount())
                if table.item(row, 0).data(DATASET_ID_ROLE) in identifiers]
        for row in rows:
            selection.select(table.model().index(row, 0),
                             table.model().index(row, table.columnCount() - 1))
        blocked = table.blockSignals(True)
        try:
            table.selectionModel().select(selection, QItemSelectionModel.ClearAndSelect | QItemSelectionModel.Rows)
            table.setCurrentItem(table.item(rows[0], 0), QItemSelectionModel.NoUpdate)
        finally:
            table.blockSignals(blocked)
        self.owner._on_dataset_selection_changed()
        self.owner._notify_dataset_catalog_changed()
        table.scrollToItem(table.item(rows[0], 0))
        self.owner.status.showMessage("Created and selected " + ", ".join(image.name for image, _ in results) + ".")
        if replot:
            self.owner._show_isar_workspace()
            for listing in (self.owner.list_pol, self.owner.list_freq,
                            self.owner.list_elev, self.owner.list_az):
                blocked = listing.blockSignals(True)
                try:
                    listing.selectAll()
                finally:
                    listing.blockSignals(blocked)
            timer = getattr(self.owner, '_autoplot_timer', None)
            if timer is not None:
                timer.stop()
            self.owner._plot_isar_image()

    def _finish_comparison(self):
        source = self._current()
        if source is None or self.vertices is None:
            return
        vertices = np.array(self.vertices, copy=True)
        name = self.source_name
        self._busy(True)
        self._set_message("Creating Include and Exclude datasets from the same source…")
        definition = self._mask_definition
        self.runner.submit(lambda cancel: comparison_cuts(source, vertices, name=name, cancel=cancel,
                                                         mask_definition=definition),
                           self._publish_comparison)

    def _publish_comparison(self, results, error):
        if error is not None:
            self._publish_cut(None, error)
            return
        if len(results) != 2 or any(radar is None for _image, radar in results):
            self._busy(False)
            self._set_message("No comparison datasets were created. This image does not support radar reconstruction.")
            return
        self._export_results(results)

    def closeEvent(self, event):
        self.deactivate()
        super().closeEvent(event)


def open_mask_cut(owner):
    bar = getattr(owner, 'plot_slider', None)
    if bar is not None:
        bar.stop_interaction()
    for name in ('btn_pan', 'btn_zoom_box', 'btn_markers', 'btn_overlays'):
        button = getattr(owner, name, None)
        if button is not None:
            button.setChecked(False)
    overlays = getattr(owner, 'spatial_overlays', None)
    if overlays is not None:
        overlays.stop_interaction()
        overlays.hide()
        button = getattr(overlays, 'toggle_button', None)
        if button is not None:
            button.setChecked(False)
    analysis = getattr(owner, 'inline_analysis', None)
    if analysis is not None:
        analysis.deactivate()
    owner._clear_zoom_box_drag()
    owner._clear_pan_drag()
    panel = getattr(owner, '_mask_cut_panel', None)
    if panel is not None and panel.canvas is not owner.plot_canvas:
        panel.deactivate()
        panel.deleteLater()
        panel = None
    if panel is None:
        panel = MaskCutPanel(owner)
        layout = owner.plot_canvas.parentWidget().layout()
        layout.insertWidget(layout.indexOf(owner.plot_canvas), panel)
        owner._mask_cut_panel = panel
    if not panel.runner.busy or panel.isHidden():
        panel.activate()
    return panel
