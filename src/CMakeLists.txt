add_library(rpa_base OBJECT rpa_util.cpp plot_chi0.cpp plot_chi0_AF.cpp self_consistent_eq_square.cpp self_consistent_eq_cubic.cpp plot_self_consistent_eq_square.cpp calc_gap.cpp calc_gap_square.cpp calc_gap_cubic.cpp calc_dispersion_square.cpp calc_dispersion_cubic.cpp calc_velocity.cpp BinarySearch.cpp) 

add_executable(rpa rpa.cpp $<TARGET_OBJECTS:rpa_base>)
target_link_libraries(rpa ${LIBS})
