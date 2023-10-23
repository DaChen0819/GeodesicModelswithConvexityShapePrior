set(PathTransformation_Headers
	${CMAKE_CURRENT_SOURCE_DIR}/path_Simplified.h
	${CMAKE_CURRENT_SOURCE_DIR}/pathTransformation.h
	${CMAKE_CURRENT_SOURCE_DIR}/path_Simplified.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/pathTransformation.cpp
	PARENT_SCOPE
	)


set(PROJECT_DIR ${PROJECT_DIR} ${CMAKE_CURRENT_SOURCE_DIR} PARENT_SCOPE)