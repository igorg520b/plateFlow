cmake_minimum_required(VERSION 3.20)

project(plateFlow LANGUAGES CXX C)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

find_package(Boost HINTS "C:\\libraries\\boost_1_79_0\\stage\\lib\\cmake\\Boost-1.79.0" REQUIRED)

if(UNIX)
    find_package(spdlog REQUIRED)
    find_package(HDF5 REQUIRED COMPONENTS CXX)
endif()

find_package(VTK 9.0 HINTS "C:\\libraries\\vtk-build\\lib\\cmake\\vtk-9.1" REQUIRED COMPONENTS
  CommonColor
  CommonCore
  FiltersSources
  FiltersGeometry
  InteractionStyle
  RenderingContextOpenGL2
  RenderingCore
  RenderingFreeType
  RenderingGL2PSOpenGL2
  RenderingOpenGL2
  ViewsQt
  IOGeometry
  RenderingAnnotation
  RenderingLabel
  IOImage
  QUIET
)

find_package(Qt5 COMPONENTS Widgets REQUIRED)
find_package(Qt5Charts)

add_compile_definitions(QT_NO_KEYWORDS)

add_executable(plateFlow
    gui/main.cpp

    gui/windows/mainwindow.cpp
    gui/windows/mainwindow.h
    gui/windows/mainwindow.ui
    icons.qrc

    gui/worker/backgroundworker.cpp
    gui/worker/backgroundworker.h
    gui/worker/modelcontrollerinterface.h

    geometry/ConcurrentPool.h


    geometry/node.cpp
    geometry/node.h

    geometry/elements/element.cpp
    geometry/elements/element.h

    geometry/floemesh.cpp
    geometry/floemesh.h

    model/parameters_sim.h
    model/parameters_sim.cpp

    model/paramswrapper.h
    model/paramswrapper.cpp


    model/modelstepinfo.h
    model/modelstepinfo.cpp
    model/model.cpp
    model/model.h

    gui/vtk/meshrepresentation.cpp
    gui/vtk/meshrepresentation.h

    solver/equationofmotionsolver.cpp
    solver/equationofmotionsolver.h
)

#target_include_directories(plateFlow PRIVATE ${Boost_INCLUDE_DIRS})
target_include_directories(plateFlow PRIVATE gui geometry geometry/elements solver model)
target_include_directories(plateFlow PRIVATE gui/windows gui/worker gui/vtk )


target_link_libraries(plateFlow PRIVATE Qt5::Widgets Qt5::Charts Qt5::Core)


# OpenMP
if(WIN32)
    target_compile_options(plateFlow PRIVATE /openmp)
else()
    target_compile_options(plateFlow PRIVATE -fopenmp -Wno-stringop-overflow -Wno-invalid-offsetof)
    target_link_libraries(plateFlow PRIVATE gomp)
    target_include_directories(plateFlow PRIVATE /usr/lib/gcc/x86_64-linux-gnu/11/include/)
endif()

# spdlog
if (WIN32)
    target_link_libraries(plateFlow PRIVATE spdlog)
    target_include_directories(plateFlow PRIVATE "C:\\libraries\\spdlog\\include")
    target_link_directories(plateFlow PRIVATE "C:\\libraries\\spdlog-build\\Release")
else()
    target_link_libraries(plateFlow PRIVATE spdlog::spdlog)
endif()

# VTK
target_link_libraries(plateFlow PRIVATE ${VTK_LIBRARIES})

# TBB and MKL
if (WIN32)

    target_include_directories(plateFlow PRIVATE "C:\\Program Files (x86)\\Intel\\oneAPI\\tbb\\2021.5.2\\include")
    target_link_directories(plateFlow PRIVATE "C:\\Program Files (x86)\\Intel\\oneAPI\\tbb\\2021.5.2\\lib\\intel64\\vc_mt")
    target_link_libraries(plateFlow PRIVATE tbb)

    target_include_directories(plateFlow PRIVATE "C:\\Program Files (x86)\\Intel\\oneAPI\\mkl\\2022.0.3\\include")
    target_link_directories(plateFlow PRIVATE "C:\\Program Files (x86)\\Intel\\oneAPI\\mkl\\2022.0.3\\lib\\intel64")
    target_link_libraries(plateFlow PRIVATE mkl_intel_lp64_dll mkl_tbb_thread_dll mkl_core_dll)

else()
    target_include_directories(plateFlow PRIVATE /usr/include/mkl)
    target_link_libraries(plateFlow PRIVATE tbb mkl_intel_lp64 mkl_gnu_thread mkl_core gomp m dl)
    target_link_libraries(plateFlow PRIVATE ${TBB_IMPORTED_TARGETS})
endif ()


# Gmsh
if(WIN32)
target_include_directories(plateFlow PRIVATE "C:\\libraries\\gmsh-4.9.5-source\\api")
target_link_directories(plateFlow PRIVATE "C:\\libraries\\gmsh-build\\Release")
endif()
target_link_libraries(plateFlow PRIVATE gmsh)

# Eigen
find_package(Eigen3 REQUIRED NO_MODULE HINTS "C:\\libraries\\eigen-build")
target_link_libraries(plateFlow PRIVATE Eigen3::Eigen)

# HDF5
if(WIN32)
    add_compile_definitions(H5_BUILT_AS_DYNAMIC_LIB)
    target_include_directories(plateFlow PRIVATE "C:\\libraries\\hdf5-build\\src")
    target_include_directories(plateFlow PRIVATE "C:\\libraries\\hdf5-develop\\src")
    target_include_directories(plateFlow PRIVATE "C:\\libraries\\hdf5-develop\\c++\\src")
    target_link_directories(plateFlow PRIVATE "C:\\libraries\\hdf5-build\\bin\\Release")
    target_link_libraries(plateFlow PRIVATE hdf5_cpp hdf5 )
else()
    target_link_libraries(plateFlow PRIVATE ${HDF5_LIBRARIES})
    target_include_directories(plateFlow PRIVATE ${HDF5_INCLUDE_DIRS})
endif()

#PBrowser
if(WIN32)
    target_link_directories(plateFlow PRIVATE "C:\\libraries\\pbrowser-build")
    target_include_directories(plateFlow PRIVATE "C:\\libraries\\PBrowser")
    target_include_directories(plateFlow PRIVATE "C:\\libraries\\PBrowser\\propertybrowser")
    target_include_directories(plateFlow PRIVATE "C:\\libraries\\PBrowser\\objectpropertybrowser")
else()
    target_link_directories(plateFlow PRIVATE "/home/s2/Projects/build-PBrowser-Desktop_Qt_5_15_2_GCC_64bit-Release")
    target_include_directories(plateFlow PRIVATE /home/s2/Projects/PBrowser/objectpropertybrowser)
    target_include_directories(plateFlow PRIVATE /home/s2/Projects/PBrowser/propertybrowser)
    target_include_directories(plateFlow PRIVATE /home/s2/Projects/PBrowser)
endif()
target_link_libraries(plateFlow PRIVATE PBrowser)

target_include_directories(plateFlow PRIVATE ${Boost_INCLUDE_DIRS})
