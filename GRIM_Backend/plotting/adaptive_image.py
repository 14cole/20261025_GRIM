"""Screen-sized ISAR previews backed by the complete computed intensity grid."""
from __future__ import annotations

import numpy as np
from matplotlib.image import AxesImage
from matplotlib.patches import Rectangle

from .modes.common import block_max_image, finite_data_limits
from GRIM_Backend.isar.geometry import image_extent


class AdaptiveImage(AxesImage):
    """Pool only the visible region; savefig always renders native pixels.

    Updating at make_image time also covers resize, shared axes and Matplotlib's
    multi-image compositing without callbacks that could outlive a plot tab.
    """

    def __init__(self, ax, band, display, *, cmap=None, vmin=None, vmax=None):
        extent = image_extent(band)
        super().__init__(ax, origin="lower", interpolation="nearest", extent=extent, cmap=cmap)
        self.native_values = np.asarray(display).T
        self.native_extent = tuple(extent)
        self.native_x = np.asarray(band["x_range"])
        self.native_y = np.asarray(band["y_range"])
        self._preview_key = None
        self._grim_native_probe = self.probe
        # Initialize at a bounded size, before layout has established the viewport.
        self._set_preview(full=True)
        limits = finite_data_limits([self.native_values])
        self.set_clim(vmin if vmin is not None else limits[0] if limits else 0.,
                      vmax if vmax is not None else limits[1] if limits else 1.)

    @staticmethod
    def _slice(size, edges, limits, pixels):
        step = (edges[1] - edges[0]) / size
        positions = sorted((np.asarray(limits) - edges[0]) / step)
        start = max(0, min(size - 1, int(np.floor(positions[0])) - 1))
        stop = max(start + 1, min(size, int(np.ceil(positions[1])) + 1))
        stride = max(1, int(np.ceil((stop - start) / max(1, pixels))))
        start = start // stride * stride
        stop = min(size, int(np.ceil(stop / stride)) * stride)
        return start, stop, stride

    def _set_preview(self, *, full=False):
        ax = self.axes
        x0, x1, y0, y1 = self.native_extent
        ny, nx = self.native_values.shape
        # One pooled cell per screen pixel, capped for very large windows.
        width = max(32, min(1536, int(ax.bbox.width)))
        height = max(32, min(1536, int(ax.bbox.height)))
        xs = self._slice(nx, (x0, x1), (x0, x1) if full else ax.get_xlim(), width)
        ys = self._slice(ny, (y0, y1), (y0, y1) if full else ax.get_ylim(), height)
        key = (*xs, *ys)
        if key == self._preview_key:
            return
        left, right, sx = xs
        bottom, top, sy = ys
        tile = self.native_values[bottom:top, left:right]
        if sx > 1 or sy > 1:
            tile = block_max_image(tile, np.arange(0, len(tile), sy),
                                   np.arange(0, tile.shape[1], sx))
        self.set_data(tile)
        # The last partial block retains the same pixel pitch. Clip it at the
        # native boundary during rendering, rather than stretching all cells.
        dx, dy = (x1 - x0) / nx, (y1 - y0) / ny
        self._extent = (x0 + left * dx, x0 + (left + tile.shape[1] * sx) * dx,
                        y0 + bottom * dy, y0 + (bottom + tile.shape[0] * sy) * dy)
        self._preview_key = key

    def make_image(self, renderer, magnification=1.0, unsampled=False):
        if self.figure.canvas.is_saving():
            preview, extent, cache = self._A, self._extent, self._imcache
            try:
                self.set_data(self.native_values)
                self._extent = self.native_extent
                return super().make_image(renderer, magnification, unsampled)
            finally:
                self._A, self._extent, self._imcache = preview, extent, cache
        self._set_preview()
        return super().make_image(renderer, magnification, unsampled)

    def probe(self, x, y):
        if x is None or y is None or not np.isfinite(x) or not np.isfinite(y):
            return None
        x0, x1, y0, y1 = self.native_extent
        if not min(x0, x1) <= x <= max(x0, x1) or not min(y0, y1) <= y <= max(y0, y1):
            return None
        ny, nx = self.native_values.shape
        ix = min(nx - 1, int((x - x0) / (x1 - x0) * nx))
        iy = min(ny - 1, int((y - y0) / (y1 - y0) * ny))
        return self.native_x[ix], self.native_y[iy], self.native_values[iy, ix]


def draw_adaptive_image(ax, band, display, **options):
    image = AdaptiveImage(ax, band, display, **options)
    ax.add_image(image)
    x0, x1, y0, y1 = image.native_extent
    boundary = Rectangle((min(x0, x1), min(y0, y1)), abs(x1 - x0), abs(y1 - y0),
                         transform=ax.transData)
    image.set_clip_path(boundary.get_path(), boundary.get_transform())
    image.set_clip_box(ax.bbox)
    # Keep the true outer bounds for autoscale/Home, independent of preview size.
    image.set_extent(image.native_extent)
    image._preview_key = None
    ax.set_aspect("auto")
    return image
