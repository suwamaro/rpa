add_library(rpa_base OBJECT
parameters.cpp
array3.cpp
rpa_util.cpp
cuba_helper.cpp
find_root.cpp
optimization.cpp
Hartree_Fock.cpp
mat_elem.cpp
find_critical_U.cpp
find_critical_U_bilayer.cpp
find_critical_T.cpp
find_critical_T_bilayer.cpp
find_critical_point.cpp
find_critical_point_bilayer.cpp
#plot_chi0.cpp
#plot_chi0_AF.cpp
self_consistent_eq.cpp
self_consistent_eq_square.cpp
self_consistent_eq_bilayer.cpp
self_consistent_eq_cubic.cpp
self_consistent_eq_two_site.cpp
#plot_self_consistent_eq_square.cpp
calc_gap.cpp
calc_gap_square.cpp
calc_gap_cubic.cpp
calc_gap_bilayer.cpp
calc_gap_two_site.cpp
calc_dispersion_square.cpp
calc_dispersion_cubic.cpp
calc_velocity.cpp
BinarySearch.cpp
calc_spectrum_square.cpp
calc_spectrum_bilayer.cpp
calc_spectrum_cubic.cpp
calc_intensity_square.cpp
calc_intensity_bilayer.cpp
calc_intensity_cubic.cpp
calc_chemical_potential.cpp
calc_chemical_potential_bilayer.cpp
calc_single_particle_energy.cpp
hoppings.cpp
calc_wave_func_bilayer.cpp
calc_binding_energy_bilayer.cpp
calc_phase_boundary_bilayer.cpp
calc_current_bilayer.cpp
calc_Raman_bilayer.cpp
calc_mean_field.cpp
#calc_two_site.cpp
mean_field_bilayer.cpp
#particle_hole1.cpp
)

# rpa
add_executable(rpa rpa.cpp $<TARGET_OBJECTS:rpa_base>)
target_link_libraries(rpa ${LIBS})

# librpa_func
add_library(rpa_func SHARED $<TARGET_OBJECTS:rpa_base>)
target_link_libraries(rpa_func ${LIBS})
