clc
clear all

basename='test_air_inflated_structure_FEM_triangle_1926_steady'
gid2tec(basename)

basename='test_air_inflated_structure_100000_FSI_medium_grid_FEM_1926'
gid2tec(basename)

basename='Tstruct_dynamic_quad_fine'
gid2tec(basename)

basename='hypar_mesh_fine_only_cable'
gid2tec(basename)
