# Adaptive Proton Therapy

This code adapts the energy of the spots/layers in an attempt to keep planned dosages.

Patient geometry may change due to positioning or weight loss. Adaptive therapy is then desired to correct these effects in a daily basis. The code reads a patient directory structure as given by MCAuto-Astroid to get the necessary treatment parameters.

## List of tasks:
- [x] Read tramp files
- [x] Read patient structure
- [x] Read MHA files
- [ ] Read ctvolume.dat
- [ ] Create internal volume representation (float)
- [ ] Function to get distal coordinates for each spot in tramp file
  * [ ] Read physical data
  * [ ] Copy data to GPU
  * [ ] Create CUDA kernel
- [ ] Apply function to planning CT (a) and CBCT (b)
- [ ] Apply vector field to coordinates extracted from CBCT
- [ ] Get distance (geometrical and energy) between a and b and correct
- [ ] Write new tramp

## Additional tasks
- [ ] Study distribution of changes to spots and collapse them to meaningfull averages
