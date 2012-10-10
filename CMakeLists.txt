CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

PROJECT(GetMatrix)  # name of the project

FIND_PACKAGE(ITK)
IF(ITK_FOUND)
	INCLUDE(${ITK_USE_FILE})
ELSE(ITK_FOUND)
	MESSAGE(FATAL_ERROR "ITK not found. Please set ITK_DIR")
ENDIF(ITK_FOUND)

FIND_PACKAGE(GenerateCLP REQUIRED)
IF(GenerateCLP_FOUND)
	INCLUDE(${GenerateCLP_USE_FILE})
ELSE(GenerateCLP_FOUND)
	MESSAGE(FATAL_ERROR "GenerateCLP not found. Please set GenerateCLP_DIR")
ENDIF(GenerateCLP_FOUND)

SET(sources GetMatrix.cxx) #define the variable "sources" that contains the names of files

GENERATECLP(sources GetMatrix.xml) # include the GCLP file to the project

ADD_EXECUTABLE(GetMatrix ${sources})  # add the files contained by "sources" to the project

TARGET_LINK_LIBRARIES(GetMatrix ${ITK_LIBRARIES})
