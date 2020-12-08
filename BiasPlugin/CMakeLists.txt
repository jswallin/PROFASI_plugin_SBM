cmake_minimum_required(VERSION 3.1)

FILE (GLOB sources ./ *.cc)

FILE (GLOB Headers "${CMAKE_CURRENT_SOURCE_DIR}/*.hh")
install (FILES ${Headers} DESTINATION include/plugins)
add_library(bias_plugin ${sources})
add_library(BiasPlugin::BiasPlugin ALIAS bias_plugin)
target_link_libraries(bias_plugin PRIVATE prf_xml::prf_xml profasi::profasi)
target_include_directories(bias_plugin PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
set_target_properties(bias_plugin PROPERTIES
    version ${${PROJECT_NAME}_version}
    soversion ${${PROJECT_NAME}_soversion}
    clean_direct_output 1
    output_name ${PROJECT_NAME}
)

install( TARGETS bias_plugin DESTINATION lib PERMISSIONS
    OWNER_READ OWNER_WRITE OWNER_EXECUTE
    GROUP_READ GROUP_EXECUTE
    WORLD_READ WORLD_EXECUTE
    )

set (profasi_plugin_libs ${profasi_plugin_libs} bias_plugin PARENT_SCOPE)

