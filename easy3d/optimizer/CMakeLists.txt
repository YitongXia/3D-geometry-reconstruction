cmake_minimum_required(VERSION 3.1)

get_filename_component(MODULE_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
set(PROJECT_NAME "easy3d_${MODULE_NAME}")
project(${PROJECT_NAME})


set(${PROJECT_NAME}_HEADERS
        optimizer_lm.h
        )

set(${PROJECT_NAME}_SOURCES
        optimizer_lm.cpp
        )


add_library(${PROJECT_NAME} STATIC ${${PROJECT_NAME}_SOURCES} ${${PROJECT_NAME}_HEADERS})

set_target_properties(${PROJECT_NAME} PROPERTIES FOLDER "easy3d")

target_include_directories(${PROJECT_NAME} PRIVATE ${EASY3D_INCLUDE_DIR})

set(CMINPACK_INCLUDE_DIR ${EASY3D_THIRD_PARTY}/cminpack)
target_include_directories(${PROJECT_NAME} PRIVATE ${CMINPACK_INCLUDE_DIR})

target_link_libraries(${PROJECT_NAME} easy3d_util)

if (MSVC)
    target_compile_definitions(${PROJECT_NAME} PRIVATE _CRT_SECURE_NO_DEPRECATE)
endif ()
