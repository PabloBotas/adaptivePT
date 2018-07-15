# Adaptive Proton Therapy

This code adapts the position, energy and weight of beamlets to keep planned dose levels to patients taking the specific geometry of a given fraction into account.

The code reads a patient directory structure as given by MCAuto-Astroid to get the necessary treatment parameters.

## List of tasks:
- [ ] None at the moment

### Some details about the implementation

#### Geometry
- To adapt to internal coordinates the CT's coordinates undergo: `x ↔ z` transformation
- Offsets are internally defined as the distance from the lower left corner (not the center of the voxel) to the isocenter. They are transformed as: `x ↔ z`, `x → -x`, `offset -= 0.5*voxn*voxd`
- Couch angle is reversed and gantry angle is taken as `ang = 270-ang`
- Positions and directions read from Tramp are transformed as in: `x → -y`, `y → -x`
- Same transformation is applied to treatment planes