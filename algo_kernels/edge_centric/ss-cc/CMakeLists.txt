set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS}  -Wall -Wextra -Wl,--no-as-needed -funroll-loops -march=native -flto -pedantic -std=c++14")
set(CMAKE_C_FLAGS  "${CMAKE_C_FLAGS}")

file(GLOB_RECURSE SRC_CC ${CMAKE_SOURCE_DIR}/FedGraph/algo_kernels/common_harness/*.cpp)
include_directories(${CMAKE_SOURCE_DIR}/FedGraph/algo_kernels/edge_centric/ss-cc)

# add_executable(cc-sgx ${SRC_CC})
# target_compile_definitions(cc-sgx PRIVATE -DSGXBACKEND)

# target_include_directories(cc-sgx PRIVATE "${CMAKE_SOURCE_DIR}/SGX-Worker/Common/")

# # add_untrusted_executable(cc SRCS ${SRC_CC} EDL ${CMAKE_SOURCE_DIR}/SGX-Worker/Enclave/Enclave.edl EDL_SEARCH_PATHS ${CMAKE_SOURCE_DIR}/SGX-Worker/Enclave/)

# #############################################
# #            Link our libraries             #
# #############################################
# target_link_libraries(cc-sgx cryptoTools TaskHandler)

# # add_dependencies(cc libenclave-sign)

add_executable(cc-ss ${SRC_CC})

set(Boost_USE_STATIC_LIBS        ON) # only find static libs
set(Boost_USE_MULTITHREADED      ON)
set(Boost_USE_STATIC_RUNTIME     ON)

find_package(Boost  COMPONENTS system thread serialization)

if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIR}) 
  #message( "Found Boost at ${Boost_LIBRARIES}")
else()
   message(FATAL_ERROR "Failed to find boost at " ${Boost_Lib_Dirs} " Need system thread")
endif()

target_compile_definitions(cc-ss PRIVATE -DSSHEBACKEND)

target_include_directories(cc-ss PRIVATE "${CMAKE_SOURCE_DIR}/SSHE-Worker/include/")

#############################################
#            Link our libraries             #
#############################################
target_link_libraries(cc-ss cryptoTools TaskHandler)
