"""Shared SOFA scene helpers (controllers) across MMS 1D/2D/3D drivers."""

import Sofa
import Sofa.Core


class NodalForceAssembler(Sofa.Core.Controller):
    """Fill a placeholder ConstantForceField after Sofa init has run.

    SOFA topology components (RegularGridTopology,
    EdgeSetTopologyContainer, QuadSetTopologyContainer,
    HexahedronSetTopologyContainer, …) are only populated after init runs,
    not during Python scene-build time. This controller defers nodal-force
    assembly to onSimulationInitDoneEvent, where it reads rest positions
    off the MechanicalObject and hands them to `compute_forces` together
    with the topology container. The returned array is written into the
    force field.

    Parameters
    ----------
    dofs           : SOFA MechanicalObject with rest_position
    topology       : SOFA topology container (read by compute_forces)
    force_field    : SOFA ConstantForceField to populate
    compute_forces : callable (nodes, topology) -> ndarray sized like
                     force_field.forces
    """

    def __init__(self, dofs, topology, force_field, compute_forces,
                 *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.dofs           = dofs
        self.topology       = topology
        self.force_field    = force_field
        self.compute_forces = compute_forces

    def onSimulationInitDoneEvent(self, event):
        nodes = self.dofs.rest_position.array().copy()
        F = self.compute_forces(nodes, self.topology)
        with self.force_field.forces.writeableArray() as forces:
            forces[:] = F
