cmake_minimum_required(VERSION 3.5)

project(mappernode LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(PYTHON_INCLUDE_DIRS ${PYTHON_INCLUDE_DIRS} /usr/lib/python3.10/site-packages/numpy/core/include)

add_executable(mc-simulator mc-simulator/mc-simulator.cpp
        mc-simulator/cmdlineoptions.cpp
        mc-simulator/cmdlineoptions.h
        mc-simulator/simulator.cpp
        mc-simulator/simulator.h)
add_executable(mappernode mappernode/mappernode.cpp mappernode/mappernode.h
        mappernode/model.cpp mappernode/model.h
        mappernode/data.cpp mappernode/data.h
        mappernode/aggregate.cpp mappernode/aggregate.h
        mappernode/ellipsemodel.cpp mappernode/ellipsemodel.h
        mappernode/entity.cpp mappernode/entity.h
        mappernode/tools.cpp mappernode/tools.h
        mappernode/scan.cpp mappernode/scan.h
        mappernode/cmdlineoptions.cpp mappernode/cmdlineoptions.h
        mappernode/ekf.cpp mappernode/ekf.h mappernode/iekf.cpp mappernode/iekf.h
        mappernode/fsm.cpp mappernode/fsm.h
        mappernode/entitymap.cpp mappernode/entitymap.h
        mappernode/handler.cpp mappernode/handler.h mappernode/visualizer.cpp mappernode/visualizer.h
        mappernode/linemodel.cpp mappernode/linemodel.h)

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
target_include_directories(mc-simulator
        PRIVATE
        ${CERES_INCLUDE_DIRS}
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}
        $<INSTALL_INTERFACE:include>
)

target_link_libraries(mc-simulator
        PRIVATE
        ${Boost_FILESYSTEM_LIBRARY}
        ${Boost_SYSTEM_LIBRARY}
        ${CERES_LIBRARIES}
)
add_executable(evaluator evaluator/evaluator.cpp)
target_include_directories(evaluator
        PRIVATE
        ${CERES_INCLUDE_DIRS}
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}
        $<INSTALL_INTERFACE:include>
)

target_link_libraries(evaluator
        PRIVATE
        ${Boost_FILESYSTEM_LIBRARY}
        ${Boost_SYSTEM_LIBRARY}
        ${CERES_LIBRARIES}
)