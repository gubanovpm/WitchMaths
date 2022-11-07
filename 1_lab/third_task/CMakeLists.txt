cmake_minimum_required(VERSION 3.4.0)
project(third_task)

find_package(matplotplusplus)
find_package(Boost 1.40)

set (EXEC_FILE "third_task.exe")

set (SOURCES 
            srcs/main.cc
    )

add_executable(${EXEC_FILE} ${SOURCES})
target_include_directories(${EXEC_FILE} PRIVATE libs)

target_compile_features(${EXEC_FILE} PRIVATE cxx_std_20)

target_link_libraries(${EXEC_FILE} PRIVATE
                        matplot
                    )

set (COMPILER_SETTINGS
    -O2
)

target_compile_options(${EXEC_FILE} PRIVATE ${COMPILER_SETTINGS})