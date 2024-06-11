cmake_minimum_required(VERSION 3.5)#最低CMake版本

project (Test)# 工程名

find_package(ifopt)
# Formulate (ifopt:ifopt_core) and solve (ifopt::ifopt_ipopt) the problem
find_package(PkgConfig REQUIRED)
pkg_search_module(MYLIB REQUIRED cppad)

add_executable(Test ex_test_ipopt.cpp)

target_include_directories(Test
    PRIVATE 
        ${PROJECT_SOURCE_DIR}
)#设置这个可执行文件 Test 需要包含的库的路径

# 设置编译和链接选项
target_include_directories(Test PUBLIC ${MYLIB_INCLUDE_DIRS})
# Pull in include directories, libraries, ... 
target_link_libraries(Test PUBLIC ifopt::ifopt_ipopt ${MYLIB_LIBRARIES}) 
