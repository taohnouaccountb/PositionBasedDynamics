set(SIMULATION_LINK_LIBRARIES PositionBasedDynamics Simulation Utils)
set(SIMULATION_DEPENDENCIES PositionBasedDynamics Simulation Utils)


set(SIMULATION_LINK_LIBRARIES
  ${SIMULATION_LINK_LIBRARIES}
)

############################################################
# GenericParameters
############################################################
include_directories(${GenericParameters_INCLUDE_DIR})
if(TARGET Ext_GenericParameters)
	set(SIMULATION_DEPENDENCIES ${SIMULATION_DEPENDENCIES} Ext_GenericParameters)
endif()

add_executable(BarDemo
	  main.cpp     
    # ${VIZ_FILES}
    ${PROJECT_PATH}/Common/Common.h
	  
	  CMakeLists.txt
)

add_definitions(-DTW_NO_LIB_PRAGMA -DTW_STATIC)

find_package( Eigen3 REQUIRED )
include_directories( ${EIGEN3_INCLUDE_DIR} )
# include_directories(${PROJECT_PATH}/extern/freeglut/include)
# include_directories(${PROJECT_PATH}/extern/glew/include)

set_target_properties(BarDemo PROPERTIES FOLDER "Demos")
set_target_properties(BarDemo PROPERTIES DEBUG_POSTFIX ${PBD_BINARY_DEBUG_POSTFIX})
set_target_properties(BarDemo PROPERTIES RELWITHDEBINFO_POSTFIX ${PBD_BINARY_RELWITHDEBINFO_POSTFIX})
set_target_properties(BarDemo PROPERTIES MINSIZEREL_POSTFIX ${PBD_BINARY_MINSIZEREL_POSTFIX})
add_dependencies(BarDemo ${SIMULATION_DEPENDENCIES})
target_link_libraries(BarDemo ${SIMULATION_LINK_LIBRARIES})
# VIS_SOURCE_GROUPS()
