set(LIB_DIR ${CMAKE_SOURCE_DIR}/code/lib)
set(CMAKE_SHARED_LIBRARY_SUFFIX ".so")
set(SRC_DIR ${CMAKE_CURRENT_LIST_DIR}/src/)
set(INC_DIR ${CMAKE_CURRENT_LIST_DIR}/include/)
add_library(snax SHARED ${SRC_DIR}/physics.cpp ${SRC_DIR}/alp_spectrum.cpp ${SRC_DIR}/alp_decays.cpp)
target_include_directories(snax PRIVATE ${INC_DIR})
set_target_properties(snax PROPERTIES VERSION ${PROJECT_VERSION} SOVERSION ${PROJECT_VERSION_MAJOR})
set_target_properties(snax PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${LIB_DIR})

pybind11_add_module(pysnax ${SRC_DIR}/python_wrapper.cpp)
target_include_directories(pysnax PRIVATE ${INC_DIR})
target_link_libraries(pysnax PUBLIC snax)
set_target_properties(pysnax PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${LIB_DIR})