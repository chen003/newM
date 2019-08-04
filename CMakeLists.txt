cmake_minimum_required(VERSION 2.8)

project(newM)
set(CMAKE_CXX_FLAGS "-std=c++11 -Wall -Wextra")

# set(CMAKE_CXX_FLAGS "-std=c++11 -Wall -Wextra")

add_executable(${PROJECT_NAME} main.cpp eigen_calc.cpp bitconfig.c sp_shell.c setup_TBME.c readsp.c readint.c sigint.c)

add_executable(testRun testRun.c eigen_calc.cpp bitconfig.c sp_shell.c setup_TBME.c readsp.c readint.c sigint.c)

add_executable(t2pn t2pn.c sp_shell.c setup_TBME.c readsp.c readint.c sigint.c)

target_link_libraries(${PROJECT_NAME} m gsl cblas)
target_link_libraries(testRun m gsl cblas)
target_link_libraries(t2pn m gsl cblas)
