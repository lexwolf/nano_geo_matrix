add_executable(ngm_smoke_qs_s smoke_qs_s.cpp)
target_link_libraries(ngm_smoke_qs_s PRIVATE nano_geo_matrix)
add_test(NAME ngm_smoke_qs_s COMMAND ngm_smoke_qs_s)

add_executable(ngm_smoke_qs_ns smoke_qs_ns.cpp)
target_link_libraries(ngm_smoke_qs_ns PRIVATE nano_geo_matrix)
add_test(NAME ngm_smoke_qs_ns COMMAND ngm_smoke_qs_ns)

add_executable(ngm_smoke_mie_s smoke_mie_s.cpp)
target_link_libraries(ngm_smoke_mie_s PRIVATE nano_geo_matrix)
add_test(NAME ngm_smoke_mie_s COMMAND ngm_smoke_mie_s)

add_executable(ngm_cup_qs_smoke cup_qs_smoke.cpp)
target_link_libraries(ngm_cup_qs_smoke PRIVATE nano_geo_matrix)
target_include_directories(ngm_cup_qs_smoke PRIVATE
    ${PROJECT_SOURCE_DIR}/modules/cup
    ${PROJECT_SOURCE_DIR}/include
)
add_test(NAME ngm_cup_qs_smoke COMMAND ngm_cup_qs_smoke)

add_executable(ngm_cup_qs_analytical_duration_test cup_qs_analytical_duration_test.cpp)
target_link_libraries(ngm_cup_qs_analytical_duration_test PRIVATE nano_geo_matrix)
target_include_directories(ngm_cup_qs_analytical_duration_test PRIVATE
    ${PROJECT_SOURCE_DIR}/modules/cup
)
add_test(NAME ngm_cup_qs_analytical_duration_test COMMAND ngm_cup_qs_analytical_duration_test)
set_tests_properties(ngm_cup_qs_analytical_duration_test PROPERTIES TIMEOUT 10)

add_executable(ngm_cup_stopping_arithmetic_test cup_stopping_arithmetic_test.cpp)
target_compile_features(ngm_cup_stopping_arithmetic_test PRIVATE cxx_std_17)
target_include_directories(ngm_cup_stopping_arithmetic_test PRIVATE
    ${PROJECT_SOURCE_DIR}/modules/cup
)
add_test(NAME ngm_cup_stopping_arithmetic_test COMMAND ngm_cup_stopping_arithmetic_test)

add_executable(ngm_cup_output_cadence_test cup_output_cadence_test.cpp)
target_compile_features(ngm_cup_output_cadence_test PRIVATE cxx_std_17)
target_include_directories(ngm_cup_output_cadence_test PRIVATE ${PROJECT_SOURCE_DIR}/modules/cup)
add_test(NAME ngm_cup_output_cadence_test COMMAND ngm_cup_output_cadence_test)

# Default sustained observations cover tens of ps; allow the Mie fixtures time to finish.
foreach(backend qs mie)
    add_executable(ngm_cup_${backend}_automatic_output_test cup_automatic_output_test.cpp)
    target_link_libraries(ngm_cup_${backend}_automatic_output_test PRIVATE nano_geo_matrix)
    target_include_directories(ngm_cup_${backend}_automatic_output_test PRIVATE ${PROJECT_SOURCE_DIR}/modules/cup)
    if(backend STREQUAL "mie")
        target_compile_definitions(ngm_cup_${backend}_automatic_output_test PRIVATE CUP_TEST_MIE)
    endif()
    add_test(NAME ngm_cup_${backend}_automatic_output_test COMMAND ngm_cup_${backend}_automatic_output_test)
    set_tests_properties(ngm_cup_${backend}_automatic_output_test PROPERTIES TIMEOUT 120)
endforeach()

add_executable(ngm_cup_mie_smoke cup_mie_smoke.cpp)
target_link_libraries(ngm_cup_mie_smoke PRIVATE nano_geo_matrix)
target_include_directories(ngm_cup_mie_smoke PRIVATE
    ${PROJECT_SOURCE_DIR}/modules/cup
    ${PROJECT_SOURCE_DIR}/include
)
add_test(NAME ngm_cup_mie_smoke COMMAND ngm_cup_mie_smoke)

add_executable(ngm_cup_qs_module_smoke cup_qs_module_smoke.cpp)
target_include_directories(ngm_cup_qs_module_smoke PRIVATE
    ${PROJECT_SOURCE_DIR}/modules/cup
    ${PROJECT_SOURCE_DIR}/include
)
target_link_libraries(ngm_cup_qs_module_smoke PRIVATE nano_geo_matrix)
add_test(NAME ngm_cup_qs_module_smoke COMMAND ngm_cup_qs_module_smoke)

add_executable(ngm_cup_mie_module_smoke cup_mie_module_smoke.cpp)
target_include_directories(ngm_cup_mie_module_smoke PRIVATE
    ${PROJECT_SOURCE_DIR}/modules/cup
    ${PROJECT_SOURCE_DIR}/include
)
target_link_libraries(ngm_cup_mie_module_smoke PRIVATE nano_geo_matrix)
add_test(NAME ngm_cup_mie_module_smoke COMMAND ngm_cup_mie_module_smoke)

add_executable(ngm_metal_database_test metal_database_test.cpp)
target_include_directories(ngm_metal_database_test PRIVATE
    ${PROJECT_SOURCE_DIR}/modules/cup
    ${PROJECT_SOURCE_DIR}/include
)
target_link_libraries(ngm_metal_database_test PRIVATE nano_geo_matrix)
add_test(NAME ngm_metal_database_test COMMAND ngm_metal_database_test)

add_executable(ngm_dielectric_database_test dielectric_database_test.cpp)
target_include_directories(ngm_dielectric_database_test PRIVATE
    ${PROJECT_SOURCE_DIR}/modules/cup
    ${PROJECT_SOURCE_DIR}/include
)
target_link_libraries(ngm_dielectric_database_test PRIVATE nano_geo_matrix)
add_test(NAME ngm_dielectric_database_test COMMAND ngm_dielectric_database_test)

add_executable(ngm_cup_convergence_test cup_convergence_test.cpp)
target_compile_features(ngm_cup_convergence_test PRIVATE cxx_std_17)
target_include_directories(ngm_cup_convergence_test PRIVATE ${PROJECT_SOURCE_DIR}/modules/cup)
add_test(NAME ngm_cup_convergence_test COMMAND ngm_cup_convergence_test)

add_executable(ngm_cup_nanoshell_convergence_test cup_nanoshell_convergence_test.cpp)
target_link_libraries(ngm_cup_nanoshell_convergence_test PRIVATE nano_geo_matrix)
target_include_directories(ngm_cup_nanoshell_convergence_test PRIVATE ${PROJECT_SOURCE_DIR}/modules/cup)
add_test(NAME ngm_cup_nanoshell_convergence_test COMMAND ngm_cup_nanoshell_convergence_test)
set_tests_properties(ngm_cup_nanoshell_convergence_test PROPERTIES TIMEOUT 60)

# Pure scale calculation; no geometry, solver, or simulation is needed.
add_executable(ngm_saturation_scale_test saturation_scale_test.cpp)
target_compile_features(ngm_saturation_scale_test PRIVATE cxx_std_17)
target_include_directories(ngm_saturation_scale_test PRIVATE ${PROJECT_SOURCE_DIR}/include)
add_test(NAME ngm_saturation_scale_test COMMAND ngm_saturation_scale_test)
set_tests_properties(ngm_saturation_scale_test PROPERTIES TIMEOUT 10)
