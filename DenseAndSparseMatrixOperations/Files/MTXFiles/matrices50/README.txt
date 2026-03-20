Mesh folder: 5x10
These four matrices all use the SAME mesh-derived connectivity pattern.
They are isotropic Laplace-type FVM matrices stored in Matrix Market general format,
with explicit mirrored off-diagonal entries.

Files:
  elliptic_iso_v1_general.mtx  -> uniform isotropic diffusivity scale 1.00
  elliptic_iso_v2_general.mtx  -> uniform isotropic diffusivity scale 1.20
  elliptic_iso_v3_general.mtx  -> uniform isotropic diffusivity scale 1.45
  elliptic_iso_v4_general.mtx  -> uniform isotropic diffusivity scale 1.75

All four are symmetric matrices written in GENERAL format.
Positive diagonal entries equal the sum of neighbouring coupling magnitudes.
Off-diagonals are negative face-neighbour couplings derived from face length / centre distance.
