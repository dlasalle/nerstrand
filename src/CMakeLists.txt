include_directories(.)
file(GLOB nerstrand_sources *.c)
file(GLOB domlib_sources ${CMAKE_SOURCE_DIR}/${DOMLIB_PATH}/*.c)
# libraries
add_library(nerstrand ${NERSTRAND_LIBRARY_TYPE} 
  ${nerstrand_sources}
  ${domlib_sources}
) 
install(TARGETS nerstrand 
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION lib
  ARCHIVE DESTINATION lib
)



# binary
add_executable(nerstrand_bin nerstrand_bin.c)
set_target_properties(nerstrand_bin PROPERTIES OUTPUT_NAME nerstrand)
target_link_libraries(nerstrand_bin nerstrand wildriver ${REQ_LIBS})
install(TARGETS nerstrand_bin
  RUNTIME DESTINATION bin
)
