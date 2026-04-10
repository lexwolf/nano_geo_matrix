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
