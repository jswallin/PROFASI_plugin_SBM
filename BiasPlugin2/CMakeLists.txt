
cmake_minimum_required(VERSION 3.30)
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
if(CMAKE_VERSION VERSION_GREATER_EQUAL "4.0.3")
	set(CMAKE_EXPERIMENTAL_CXX_IMPORT_STD "d0edc3af-4c50-42ea-a356-e2862fe7a444")
elseif(CMAKE_VERSION VERSION_GREATER_EQUAL "4.0.0")
    set(CMAKE_EXPERIMENTAL_CXX_IMPORT_STD "a9e1cf81-9932-4810-974b-6eccaf14e457")
elseif(CMAKE_VERSION VERSION_GREATER_EQUAL "3.29.20240416")
    set(CMAKE_EXPERIMENTAL_CXX_IMPORT_STD "0e5b6991-d74f-4b3d-a41c-cf096e0b2508")
endif()
set (CMAKE_POSITION_INDEPENDENT_CODE ON)
project(BiasPlugin2 CXX)
set(CMAKE_CXX_MODULE_STD 1)
set(plugin BiasPlugin2)
find_package(PROFASI REQUIRED COMPONENTS
    profasi_combined
)
add_library(BiasPlugin2 BiasPlugin2.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/BiasEnergy.cc
)
target_sources(BiasPlugin2 PUBLIC FILE_SET CXX_MODULES FILES BiasPlugin2.ixx     ${CMAKE_CURRENT_SOURCE_DIR}/BiasEnergy.ixx
)

add_library(${plugin}::${plugin} ALIAS ${plugin})
set_target_properties(${plugin} PROPERTIES
    version 1.0
    soversion 1.0
    clean_direct_output 1
)
target_include_directories(${plugin}
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
        $<INSTALL_INTERFACE:include>
    PRIVATE
        ${PROJECT_BINARY_DIR}
)
target_link_libraries(${plugin}
    PUBLIC
    profasi::profasi_combined
    profasi::io_nompi

    PRIVATE
)
install(TARGETS ${plugin} EXPORT ${plugin}
    FILE_SET CXX_MODULES
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/modules"
)
export(EXPORT ${plugin} NAMESPACE ${plugin}:: FILE ${plugin}Config.cmake )
install(EXPORT ${plugin} NAMESPACE ${plugin}:: DESTINATION share/cmake FILE ${plugin}Config.cmake )
#export(TARGETS ${plugin} NAMESPACE ${plugin}:: FILE ${plugin}Config.cmake)
        