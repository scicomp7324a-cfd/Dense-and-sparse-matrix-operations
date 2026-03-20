Mesh: structured 200x250 cells on [0,0.1]x[0,0.1]
Files: points.txt, faces.txt, cells.txt, boundary.txt follow the same format/orientation as the provided example.
Face ordering:
  - faces.txt first lists all horizontal edges row-by-row (j=0..ny, i=0..nx-1),
    then all vertical edges column-by-column (i=0..nx, j=0..ny-1).
Cell ordering:
  - cells are numbered row-major (i + nx*j), and each cell entry is 4(bottom right top left) face indices.
Boundary patches:
  - movingWall: top boundary horizontal faces (count=nx)
  - fixedWall: bottom horizontal + left vertical + right vertical faces (count=nx+2*ny)

Matrices (MatrixMarket, symmetric):
  - laplacian_1.mtx: isotropic diffusion (gamma_x=1, gamma_y=1), Dirichlet treatment on domain boundary (boundary coefficient added to diagonal).
  - laplacian_2.mtx: anisotropic diffusion (gamma_x=2, gamma_y=0.5), same stencil structure and boundary treatment.
