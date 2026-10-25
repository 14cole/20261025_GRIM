"""Overlay-compatible XY masks with explicit physical coordinates and lineage."""
import hashlib
import json
import os
from pathlib import Path
import tempfile
import numpy as np
from GRIM_Backend.plotting.overlay_data import LENGTH_UNITS, read_overlay_points

PREFIX = '# GRIM-MASK '
MAX_VERTICES = 4096
BOUNDARY_TOLERANCE_METRES = 1e-10


class MissingMaskMetadata(ValueError):
    pass


def polygon_vertices(vertices):
    points = np.asarray(vertices, float)
    if points.ndim != 2 or points.shape[1] != 2 or len(points) < 3:
        raise ValueError(f'A mask needs 3 to {MAX_VERTICES} XY vertices in one closed boundary.')
    if not np.isfinite(points).all():
        raise ValueError('Mask vertices must be finite and form a single path.')
    if np.array_equal(points[0], points[-1]):
        points = points[:-1]
    if len(points) > MAX_VERTICES:
        raise ValueError(f'A mask can contain at most {MAX_VERTICES} vertices.')
    area = np.sum(points[:, 0]*np.roll(points[:, 1], -1)-points[:, 1]*np.roll(points[:, 0], -1))
    if len(points) < 3 or not np.isfinite(area) or area == 0:
        raise ValueError('The mask boundary has no area.')
    return points


def _frame(image):
    contract = image.metadata.get('contract', {})
    if contract.get('image_plane') != 'horizontal' or contract.get('phase_center') != 'dataset_fixed_origin':
        raise ValueError('This image does not declare a horizontal frame and fixed origin for reusable masks.')
    scale = float(contract['distance_scale_per_metre'])
    basis = np.asarray([contract['cross_range_basis'][:2], contract['range_basis'][:2]], float)
    signs = np.asarray(contract.get('coordinate_signs', [1, 1]), float)
    if (not np.isfinite(scale) or scale <= 0 or basis.shape != (2, 2)
            or not np.allclose(basis @ basis.T, np.eye(2), atol=1e-10)
            or signs.shape != (2,) or not np.isin(signs, [-1, 1]).all()):
        raise ValueError('The image coordinate frame is invalid.')
    return scale, basis, signs


def _hash(definition):
    physical = {key: definition[key] for key in ('format', 'version', 'units', 'frame', 'operation',
                                                'boundary_rule', 'boundary_tolerance_metres', 'vertices')}
    return hashlib.sha256(json.dumps(physical, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()


def with_operation(definition, include):
    result = dict(definition, operation='include' if include else 'exclude')
    result['sha256'] = _hash(result)
    return result


def validate_definition(image, vertices, definition):
    """Prevent cached provenance from disagreeing with the boundary being cut."""
    current = define_mask(image, vertices, definition.get('operation') == 'include')
    for key in ('format', 'version', 'units', 'frame', 'operation', 'boundary_rule', 'boundary_tolerance_metres'):
        if definition.get(key) != current[key]:
            raise ValueError('The saved mask settings do not match the current boundary.')
    expected = polygon_vertices(definition['vertices'])
    actual = np.asarray(current['vertices'])
    if (expected.shape != actual.shape or not np.allclose(expected, actual, rtol=1e-12, atol=1e-12)
            or definition.get('sha256') != _hash(definition)):
        raise ValueError('The saved mask coordinates do not match the current boundary.')
    return definition


def define_mask(image, vertices, include=True):
    scale, basis, signs = _frame(image)
    body = (polygon_vertices(vertices)*signs/scale) @ basis
    definition = dict(format='GRIM-mask', version=1, units='m', frame='body_xy_fixed_origin',
                      operation='include' if include else 'exclude', vertices=body.tolist(),
                      boundary_rule='pixel_centers_including_boundary',
                      boundary_tolerance_metres=BOUNDARY_TOLERANCE_METRES,
                      source_recipe=dict(reconstruction=image.metadata.get('reconstruction'),
                          window=image.metadata.get('window'), elevation_degrees=image.metadata.get('elevation_degrees'),
                          grid_shape=list(image.magnitude.shape),
                          x_extent=[float(image.x[0]), float(image.x[-1])],
                          y_extent=[float(image.y[0]), float(image.y[-1])],
                          display_units_per_metre=scale,
                          look_center_degrees=image.metadata.get('contract', {}).get('look_center_degrees')))
    definition['sha256'] = _hash(definition)
    return definition


def mask_header(path):
    path = Path(path)
    if path.stat().st_size > 2 * 1024**2:
        raise ValueError('Mask files must be smaller than 2 MiB.')
    with path.open(encoding='utf-8-sig') as stream:
        for raw in stream:
            if raw.startswith(PREFIX):
                header = json.loads(raw[len(PREFIX):])
                if not isinstance(header, dict):
                    raise ValueError('Invalid mask header.')
                return header
            if raw.strip() and not raw.startswith('#'):
                break
    return None


def save_mask(path, image, vertices, *, include=True, definition=None):
    definition = (validate_definition(image, vertices, definition) if definition is not None
                  else define_mask(image, vertices, include))
    definition = with_operation(definition, include)
    points = polygon_vertices(definition['vertices'])
    header = {key: value for key, value in definition.items() if key != 'vertices'}
    target, temporary = Path(path), None
    try:
        with tempfile.NamedTemporaryFile(mode='w', encoding='utf-8', newline='\n',
                dir=target.parent, prefix='.'+target.name+'.', suffix='.tmp', delete=False) as stream:
            temporary = Path(stream.name)
            stream.write(PREFIX + json.dumps(header, sort_keys=True, allow_nan=False) + '\n')
            stream.write('# Body X/Y in meters, relative to the common fixed target origin.\nx,y\n')
            for x, y in np.vstack((points, points[0])):
                stream.write(f'{x:.17g},{y:.17g}\n')
        os.replace(temporary, target)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)
    return definition


def load_mask(path, image, *, legacy_unit=None, include=True):
    header = mask_header(path)
    points = polygon_vertices(read_overlay_points(path))
    if header is None:
        if legacy_unit is None:
            raise MissingMaskMetadata('This overlay has no mask frame or units. Choose its units; XY will use the current plot frame.')
        if legacy_unit not in LENGTH_UNITS:
            raise ValueError('Unsupported overlay length unit.')
        scale, _, _ = _frame(image)
        displayed = points * LENGTH_UNITS[legacy_unit] * scale
        return displayed, define_mask(image, displayed, include)
    if (header.get('format') != 'GRIM-mask' or header.get('version') != 1
            or header.get('units') != 'm' or header.get('frame') != 'body_xy_fixed_origin'
            or header.get('operation') not in ('include', 'exclude')
            or header.get('boundary_tolerance_metres') != BOUNDARY_TOLERANCE_METRES
            or header.get('boundary_rule') != 'pixel_centers_including_boundary'):
        raise ValueError('Unsupported mask format, units, frame, or boundary rule.')
    definition = dict(header, vertices=points.tolist())
    if definition.get('sha256') != _hash(definition):
        raise ValueError('The mask coordinates or settings changed after saving (checksum mismatch).')
    scale, basis, signs = _frame(image)
    displayed = (points @ basis.T)*scale*signs
    return displayed, definition
