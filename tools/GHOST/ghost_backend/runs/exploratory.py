"""User-selected 2D meshes, without error-driven refinement or certification."""
from copy import deepcopy
import math


MESH_MODES = ('global', 'geometry',
              'certification_base', 'certification_fine')


def _integer(value, name, positive=False):
    if type(value) is not int or (positive and value < 1):
        raise ValueError(f'{name} must be a {"positive " if positive else ""}integer.')
    return value


def select_exploratory_mesh(snapshot, frequency_ghz, *, mesh_mode='certification_base',
                            geometry_units='inches', material_base_dir=None,
                            factorization='adaptive', ram_budget_gib=None,
                            panels_per_wavelength=None,
                            basis_order=None, mesh_grading=True,
                            execution_overrides=None):
    """Return an independent geometry, execution options, and selection record.

    Base/fine refer to the initial candidates, not an accepted or subsequently
    refined certification result. Eligible inputs use the same coarsened P2/P3
    pair as certification. Other inputs retain its linear base/fine pair.
    Global requires panels_per_wavelength.
    Geometry mode preserves the file's N values. Global applies one
    requested sizing value to all existing segments without changing primitives.
    Manual modes never coarsen or change the requested basis order. All modes
    disable the explicit-N density floor; input and resource checks remain.
    """
    from ghost_backend.execution.options import efficient_defaults, validate_options
    from ghost_backend.runs.quality import validate_mesh_convergence_policy
    from ghost_backend.twod.adaptive_geometry import candidate_meshes, eligible_snapshot
    from ghost_backend.twod.preparation import prepare_geometry

    if mesh_mode not in MESH_MODES:
        raise ValueError('MESH_MODE must be one of: ' + ', '.join(MESH_MODES))
    frequency = float(frequency_ghz)
    if not math.isfinite(frequency) or frequency <= 0:
        raise ValueError('Frequency must be finite and positive in GHz.')
    overrides = dict(execution_overrides or {})
    reserved = {'mesh_strategy', 'basis_order', 'mesh_grading', 'enforce_mesh_density',
                'factorization', 'ram_budget_gib'}
    if reserved & overrides.keys():
        raise ValueError('Use named mesh/factorization arguments instead of execution_overrides.')
    options = validate_options(dict(efficient_defaults(), **overrides))
    options = validate_options(dict(options, factorization=factorization,
                                    mesh_strategy='global', basis_order=1,
                                    mesh_grading=mesh_grading, enforce_mesh_density=False,
                                    ram_budget_gib=ram_budget_gib))
    candidate = deepcopy(snapshot)
    # A saved candidate must not smuggle previous refinement into a manual run.
    for key in ('_2d_hp_coarsening', '_2d_hp_refinements',
                '_2d_certification_refinement_factor', '_2d_certification_base_segment_n'):
        candidate.pop(key, None)
    if mesh_mode == 'global':
        n = -_integer(panels_per_wavelength, 'panels_per_wavelength', positive=True)
    else:
        n = None
    if panels_per_wavelength is not None and mesh_mode != 'global':
        raise ValueError('panels_per_wavelength applies only to global mode.')
    presets = mesh_mode.startswith('certification_')
    if presets and basis_order is not None:
        raise ValueError('Certification presets choose their own basis; use a manual mode to select basis_order.')
    segments = candidate.get('segments', [])
    for segment in segments if n is not None else ():
        props = list(segment.get('properties', []) or [])
        if not props:
            props.append(str(segment.get('seg_type', 2)))
        while len(props) < 2:
            props.append('0')
        props[1] = str(n)
        segment['properties'] = props
    if not presets:
        options = validate_options(dict(options, basis_order=1 if basis_order is None else basis_order))
    record = dict(requested=mesh_mode, candidate='geometry', basis_order=1,
                  mesh_strategy='global', coarsening_factor=1.,
                  certification_performed=False, convergence_owner='user',
                  mesh_grading=mesh_grading, enforce_mesh_density=False,
                  panels_per_wavelength=panels_per_wavelength)
    if presets:
        _, _, materials, scale = prepare_geometry(
            candidate, material_base_dir, geometry_units)
        factor = validate_mesh_convergence_policy()['fine_factor']
        candidates = candidate_meshes(candidate, materials, factor, True,
                                      [frequency], scale)
        phase = 'base' if mesh_mode == 'certification_base' else 'fine'
        _, candidate, degree = next(item for item in candidates if item[0] == phase)
        candidate = deepcopy(candidate)
        polynomial = '_2d_hp_coarsening' in candidate
        strategy = 'local' if polynomial else 'global'
        reason = ('Initial higher-order candidate used by mesh certification.' if polynomial else
                  eligible_snapshot(snapshot, materials, [frequency], scale)[1])
        options.update(mesh_strategy=strategy, basis_order=degree)
        record.update(candidate=phase, basis_order=degree, mesh_strategy=strategy,
                      coarsening_factor=float(candidate.get('_2d_hp_coarsening', 1.)),
                      fine_factor=factor, selection_reason=reason,
                      higher_order_candidate=polynomial)
    else:
        record.update(candidate=mesh_mode, basis_order=options['basis_order'],
                      mesh_strategy=options['mesh_strategy'])
    record['geometry_n'] = [(list(s.get('properties', []) or []) + [0, 0])[1]
                           for s in candidate.get('segments', [])]
    return candidate, options, record
