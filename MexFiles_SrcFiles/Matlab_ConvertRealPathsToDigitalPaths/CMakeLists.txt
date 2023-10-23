cmake_minimum_required(VERSION 3.6)
set (CMAKE_CXX_STANDARD 17)

project(ConvertRealPathsToDigitalPaths)

# Matlab setting
set(MATLAB_APP "/Applications/MATLAB_R2023b.app")
if(EXISTS "/Applications/MATLAB_R2023b.app/bin/maci64")
  set(MATLAB_LIB_SUBDIR "/bin/maci64")
else()
  set(MATLAB_LIB_SUBDIR "/bin/maca64")
endif()

set(MATLAB_INCLUDE_SUBDIR "/extern/include")

set(CMAKE_MACOSX_RPATH 0)

set(JMM_CPPLibs_dir "${CMAKE_CURRENT_SOURCE_DIR}/../JMM_CPPLibs" CACHE FILEPATH "JMM_CPPLibs directory")
set(DummyBinDir "${CMAKE_CURRENT_BINARY_DIR}/Dummy")
add_subdirectory("../JMM_CPPLibs" "${DummyBinDir}/JMM_CPPLibs")
add_subdirectory("AlgebraDataSturctures" "${DummyBinDir}/AlgebraDataSturctures")
add_subdirectory("PathTransformation" "${DummyBinDir}/PathTransformation")

include_directories(${PROJECT_DIR} ${CMAKE_CURRENT_SOURCE_DIR})
include_directories("${MATLAB_APP}${MATLAB_INCLUDE_SUBDIR}")
include_directories("../MatlabInterface")


set(MATLAB_LIB_DIR "${MATLAB_APP}${MATLAB_LIB_SUBDIR}")
set(MATLAB_LIBRARIES "${MATLAB_LIB_DIR}/libmat.dylib" "${MATLAB_LIB_DIR}/libmex.dylib" "${MATLAB_LIB_DIR}/libmx.dylib")


set(Project_Headers 
 	${DataStructures_Headers}
 	${LinearAlgebra_Headers}
 	${Output_Headers}
 	${Macros_Headers}
        ${AlgebraDataSturctures_Headers}
        ${PathTransformation_Headers}
)


 # ---- Headers IDE layout -----
source_group("DataStructures" FILES ${DataStructures_Headers})
source_group("Macros" FILES ${Macros_Headers})
source_group("LinearAlgebra" FILES ${LinearAlgebra_Headers})
source_group("Output" FILES ${Output_Headers})
source_group("AlgebraDataSturctures" FILES  ${AlgebraDataSturctures_Headers})
source_group("PathTransformation" FILES  ${PathTransformation_Headers})


add_library(ConvertRealPathsToDigitalPaths SHARED ConvertRealPathsToDigitalPaths.cxx  ${PathTransformation_Headers} ${AlgebraDataSturctures_Headers})

target_link_libraries(ConvertRealPathsToDigitalPaths ${MATLAB_LIBRARIES})









