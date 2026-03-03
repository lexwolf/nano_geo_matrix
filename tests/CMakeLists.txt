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
add_test(NAME ngm_cup_qs_smoke COMMAND ngm_cup_qs_smoke)
