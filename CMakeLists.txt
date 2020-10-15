cmake_minimum_required(VERSION 3.5.0)

project(iVol)
set(TARGET ${CMAKE_PROJECT_NAME})
set(CMAKE_BUILD_TYPE Debug)

set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTORCC ON)

# Dependencies
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})
find_package(OpenGL REQUIRED)

# Qt5 setup
find_package(Qt5Widgets REQUIRED)
find_package(Qt5Concurrent REQUIRED)
find_package(Qt5OpenGL REQUIRED)
find_package(Qt5PrintSupport REQUIRED)

#set(CMAKE_MACOSX_RPATH 1)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR})
set(SOURCES 
	"src/main.cpp"
	"src/ui/mainwindow.cpp"
	"src/algorithm/volumemanager.cpp"
    	"src/ui/transferfunctionwindow.cpp"
	)
set(HEADERS 
	"src/ui/mainwindow.h"
	"src/algorithm/volumemanager.h"
    	"src/ui/transferfunctionwindow.h"
if(APPLE)	
    	"src/osxhelper.h"
endif()	
	)
set(UI_SOURCES
	"src/ui/mainwindow.ui"
    	"src/ui/transferfunctionwindow.ui"
if(APPLE)	
    	"src/osxhelper.mm"
endif()	
	)
set(RESOURCES
	"res/icons.qrc"
	)

add_executable(${TARGET} ${SOURCES} ${HEADERS} ${UI_SOURCES} ${RESOURCES})
target_include_directories(${TARGET} PRIVATE
	${PROJECT_SOURCE_DIR} 
	${PROJECT_SOURCE_DIR}/src 
	${CMAKE_CURRENT_BINARY_DIR}
	)
qt5_use_modules(${TARGET} Widgets Concurrent OpenGL PrintSupport)
target_link_libraries(${TARGET} ${VTK_LIBRARIES} ${OPENGL_LIBRARIES})
