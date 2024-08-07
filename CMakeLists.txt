cmake_minimum_required(VERSION 3.16)

set(addon_name ng_mesh_cutting)

project(${addon_name})
include(ngsolve_addon.cmake)

add_ngsolve_addon(${addon_name} src/ng_mesh_cutting.cpp)

set(install_dir ${ADDON_INSTALL_DIR_PYTHON}/${addon_name})
install(TARGETS ${addon_name} DESTINATION ${install_dir})
install(FILES src/__init__.py DESTINATION ${install_dir})

ngsolve_generate_stub_files(${addon_name})
