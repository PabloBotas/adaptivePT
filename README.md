# Adaptive Proton Therapy

This code adapts the energy of the spots/layers in an attempt to keep planned dosages.

Patient geometry may change due to positioning or weight loss. Adaptive therapy is then desired to correct these effects in a daily basis. The code reads a patient directory structure as given by MCAuto-Astroid to get the necessary treatment parameters.

## List of tasks:
- [x] Read tramp files
- [x] Read patient files
    - [x] Read ctvolume.dat
    - [x] Create internal volume representation (float)
- [ ] Raytrace spots to final positions, output them
- [ ] Get vector field values at those positions:
    - [ ] Internally or use `plastimatch probe`?
- [ ] Apply vector field and get **intended positions**
- [ ] Save *x,y* as final positions and get WEPL of the correspondent depth. Get new desired energy.
- [ ] Save changes to log file for further study
- [ ] Write new tramps

## Additional tasks
- [ ] Study distribution of changes to spots and collapse them to meaningfull averages
