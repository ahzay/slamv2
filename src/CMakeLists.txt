#cmake_minimum_required(VERSION 3.5)

project(mappernode LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(PYTHON_INCLUDE_DIRS ${PYTHON_INCLUDE_DIRS} /usr/lib/python3.10/site-packages/numpy/core/include)

add_executable(mappernode main.cpp vecext.hpp visualizer.hpp simulate.hpp model.hpp handler.hpp functions.hpp ekf.hpp fsm.hpp preproc.hpp model.cpp model.h data.cpp data.h aggregate.cpp aggregate.h ellipsemodel.cpp ellipsemodel.h entity.cpp entity.h tools.cpp tools.h scan.cpp scan.h cmdlineoptions.cpp cmdlineoptions.h)

# external libs
# find_package(FastAD CONFIG REQUIRED)
find_package(Boost COMPONENTS system filesystem REQUIRED)
find_package(Ceres REQUIRED)
find_package(autodiff REQUIRED)
# matplot
find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
find_package(Python3 COMPONENTS NumPy)
find_package(OpenMP REQUIRED)

if (Python3_NumPy_FOUND AND OpenMP_CXX_FOUND)
    target_link_libraries(mappernode INTERFACE
            Python3::NumPy OpenMP::OpenMP_CXX autodiff::autodiff
            )
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} ")
else ()
    target_compile_definitions(mappernode INTERFACE WITHOUT_NUMPY)
endif ()
install(
        TARGETS mappernode
        EXPORT install_targets
)

target_include_directories(mappernode
        PRIVATE
        ${CERES_INCLUDE_DIRS}
        ${PYTHON_INCLUDE_DIRS}
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}
        $<INSTALL_INTERFACE:include>
        )

target_link_libraries(mappernode
        PRIVATE
        ${Boost_FILESYSTEM_LIBRARY}
        ${Boost_SYSTEM_LIBRARY}
        ${CERES_LIBRARIES}
        Python3::Python
        Python3::Module
        )
