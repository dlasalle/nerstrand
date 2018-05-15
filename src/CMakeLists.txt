include_directories(.)
file(GLOB nerstrand_sources *.c)
file(GLOB domlib_sources ${CMAKE_SOURCE_DIR}/${DOMLIB_PATH}/*.c)
# libraries
add_library(nerstrand ${NERSTRAND_LIBRARY_TYPE} 
  ${nerstrand_sources}
  ${domlib_sources}
) 

if (NOT WIN32)
  install(TARGETS nerstrand 
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  )
endif()



# binary
add_executable(nerstrand_bin nerstrand_bin.c)
set_target_properties(nerstrand_bin PROPERTIES OUTPUT_NAME nerstrand)
target_link_libraries(nerstrand_bin nerstrand wildriver ${REQ_LIBS})
if (NOT WIN32)
  install(TARGETS nerstrand_bin
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )
endif()
