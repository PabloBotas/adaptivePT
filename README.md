# Adaptive Proton Therapy

This code adapts the energy of the spots/layers in an attempt to keep planned dosages.

Patient geometry may change due to positioning or weight loss. Adaptive therapy is then desired to correct these effects in a daily basis. The code reads a patient directory structure as given by MCAuto-Astroid to get the necessary treatment parameters.

## List of tasks:
- [x] Read tramp files
- [x] Read patient structure
- [x] Read MHA files
- [x] Read ctvolume.dat
- [x] Create internal volume representation (float)
- [ ] Translate spot map to mha. Warp it with vf and read it in associating it with given tramp files.
- [ ] Get distal coordinates for each spot in tramp file
  * [ ] Read physical data
  * [ ] Copy data to GPU
  * [ ] Create CUDA kernel
- [ ] Get physical distance between warped (plan) spot map (3D) and CBCT's spot map (3D)
- [ ] Save changes to log file for further study
- [ ] Write new tramps

## Additional tasks
- [ ] Study distribution of changes to spots and collapse them to meaningfull averages
