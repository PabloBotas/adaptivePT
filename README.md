# Adaptive Proton Therapy

This code adapts the energy of the spots/layers in an attempt to keep planned dose levels to patients taking the specific geometry of a given fraction into account.

Patient geometry may change due to positioning or weight loss. Adaptive therapy is then desired to correct these effects in a daily basis. The code reads a patient directory structure as given by MCAuto-Astroid to get the necessary treatment parameters.

## List of tasks:
- [ ] Figure out what is messing up the results
- [ ] Study distribution of changes to spots and collapse them to meaningful clusters

### Some details about the implementation

#### Geometry
- To adapt to internal coordinates the CT's coordinates undergo: x <-> z transformation
- Offsets are internally defined as the distance from the lower left corner (not the center of the voxel) to the isocenter
- Positions and directions read from Tramp are transformed as in: x → -y, y → -x
- Same transformation is applied to treatment planes