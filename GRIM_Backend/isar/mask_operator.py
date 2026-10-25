"""Adjoint-consistent, oversampled Fourier operator for native polar samples.

The interpolation is a bounded 16-tap Kaiser-windowed sinc on each axis.
Forward/adjoint use the same sparse stencil; no measured sample is fabricated.
"""
import numpy as np
from scipy.fft import fft2, ifft2, next_fast_len
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import LinearOperator
from .quality import C0
from GRIM_Backend.plotting.requests import PlotCancelled


class MaskPolarOperator(LinearOperator):
    def __init__(self, x, y, azimuth_degrees, frequency_hz, elevation_degrees=0.,
                 *, valid=None, cancel=None, maximum_bytes=384 * 1024**2):
        self.cancel = cancel
        self.x, self.y = np.asarray(x, float), np.asarray(y, float)
        for axis in (self.x, self.y):
            if (axis.ndim != 1 or len(axis) < 2 or not np.isfinite(axis).all()
                    or not np.allclose(np.diff(axis), np.diff(axis)[0], rtol=1e-7, atol=1e-12)
                    or np.diff(axis)[0] == 0):
                raise ValueError('Mask reconstruction needs a uniform spatial grid.')
        self.grid_shape = (len(x), len(y))
        self.fft_shape = tuple(next_fast_len(2 * n) for n in self.grid_shape)
        self.size_fft = int(np.prod(self.fft_shape))
        az_values, freq_values = np.asarray(azimuth_degrees, float), np.asarray(frequency_hz, float)
        if (az_values.ndim != 1 or freq_values.ndim != 1 or not len(az_values) or not len(freq_values)
                or not np.isfinite(az_values).all() or not np.isfinite(freq_values).all()
                or np.any(freq_values <= 0) or not np.isfinite(elevation_degrees)):
            raise ValueError('Mask reconstruction needs finite native angles and positive frequencies.')
        az = np.deg2rad(az_values)[:, None]
        freq = freq_values[None, :]
        native_shape = (len(az), freq.shape[1])
        self.valid = np.ones(native_shape, bool) if valid is None else np.asarray(valid, bool)
        if self.valid.shape != native_shape or not self.valid.any():
            raise ValueError('Mask reconstruction requires finite native samples.')
        count = int(self.valid.sum())
        # COO construction plus CSR stencil, FFT vectors, and solver vectors.
        estimate = count * 256 * 40 + self.size_fft * 64 + np.prod(self.grid_shape) * 128
        if estimate > maximum_bytes:
            raise ValueError('Mask reconstruction exceeds its memory budget. Select a smaller aperture, band, or scene.')
        factor = 2 / C0 * abs(np.cos(np.deg2rad(float(elevation_degrees))))
        kx = np.broadcast_to(factor * freq * np.sin(az), native_shape)[self.valid]
        ky = np.broadcast_to(factor * freq * np.cos(az), native_shape)[self.valid]
        self.phase = np.exp(-2j*np.pi*(kx*self.x[len(x)//2] + ky*self.y[len(y)//2]))
        self.ix = (np.arange(len(x)) - len(x)//2) % self.fft_shape[0]
        self.iy = (np.arange(len(y)) - len(y)//2) % self.fft_shape[1]
        coordinates = (kx*np.diff(self.x)[0]*self.fft_shape[0], ky*np.diff(self.y)[0]*self.fft_shape[1])
        indices, weights = [], []
        for q, size in zip(coordinates, self.fft_shape):
            base = np.floor(q).astype(np.int64)
            ids = base[:, None] + np.arange(-7, 9)
            delta = q[:, None]-ids
            window = np.i0(8.6*np.sqrt(np.maximum(0., 1-(delta/8)**2))) / np.i0(8.6)
            weight = np.sinc(delta)*window
            weight /= weight.sum(axis=1)[:, None]
            indices.append(ids % size)
            weights.append(weight)
        self._check()
        cols = (indices[0][:, :, None]*self.fft_shape[1]+indices[1][:, None, :]).reshape(-1)
        vals = (weights[0][:, :, None]*weights[1][:, None, :]).reshape(-1)
        rows = np.repeat(np.arange(count), 256)
        self.sampling = coo_matrix((vals, (rows, cols)), shape=(count, self.size_fft)).tocsr()
        self.sampling_transpose = self.sampling.T
        self.normalization = np.sqrt(count)
        super().__init__(dtype=np.dtype(np.complex128), shape=(count, int(np.prod(self.grid_shape))))

    def _check(self):
        if self.cancel is not None and self.cancel():
            raise PlotCancelled()

    def _matvec(self, value):
        self._check()
        padded = np.zeros(self.fft_shape, complex)
        padded[np.ix_(self.ix, self.iy)] = np.asarray(value).reshape(self.grid_shape)
        return self.phase * (self.sampling @ fft2(padded).ravel()) / self.normalization

    def _rmatvec(self, value):
        self._check()
        spectrum = self.sampling_transpose @ (self.phase.conj()*np.asarray(value).ravel())
        pixels = ifft2(spectrum.reshape(self.fft_shape))*self.size_fft
        return pixels[np.ix_(self.ix, self.iy)].ravel() / self.normalization
