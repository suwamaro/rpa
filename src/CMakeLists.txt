add_library(rpa_base OBJECT
parameters.cpp
rpa_util.cpp
plot_chi0.cpp
plot_chi0_AF.cpp
self_consistent_eq.cpp
self_consistent_eq_square.cpp
self_consistent_eq_bilayer.cpp
self_consistent_eq_cubic.cpp
plot_self_consistent_eq_square.cpp
calc_gap.cpp
calc_gap_square.cpp
calc_gap_cubic.cpp
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
calc_chemical_potential_bilayer.cpp
hoppings.cpp) 

add_executable(rpa rpa.cpp $<TARGET_OBJECTS:rpa_base>)
target_link_libraries(rpa ${LIBS})
