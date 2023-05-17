cmake_minimum_required(VERSION 3.11.0)
project(comp_math_coursework)

find_package(Matplot++ REQUIRED)

set (EXEC_FILE "coursework")

set (SOURCES 
            srcs/main.cc
            srcs/shmemes.cc
    )

add_executable(${EXEC_FILE} ${SOURCES})
target_include_directories(${EXEC_FILE} PRIVATE libs)

target_compile_features(${EXEC_FILE} PRIVATE cxx_std_20)

target_link_libraries(${EXEC_FILE} PRIVATE
                        Matplot++::matplot
                    )

set (COMPILER_SETTINGS -O3)

if (DEBUG)
    target_compile_definitions(${EXEC_FILE} PRIVATE DEBUG)
endif()

target_compile_options(${EXEC_FILE} PRIVATE ${COMPILER_SETTINGS})
