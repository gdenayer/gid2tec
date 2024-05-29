clc
clear all

% linear 2D membrane element/Triangle/MEMBRANE1
basename='test_air_inflated_structure_FEM_triangle_1926_steady'
gid2tec(basename)

% linear 2D membrane element/Triangle/MEMBRANE1
basename='test_air_inflated_structure_100000_FSI_medium_grid_FEM_1926'
gid2tec(basename)

% 2D shell element/Quadrilateral/SHELL8
basename='Tstruct_dynamic_quad_fine'
gid2tec(basename)

% linear 1D truss element/cable/TRUSS1
basename='hypar_mesh_fine_only_cable'
gid2tec(basename)

% linear 3D solid element/Hexahedra/SOLIDHEXA1
basename='carat_fsi3_EMPIRE_ruk3_hexa_32x2x1_z0pt02'
gid2tec(basename)
