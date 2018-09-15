# C++ interface source files
set(SRCCXXIF AlgebraicVector.cpp DA.cpp DACEException.cpp DAFormatter.cpp MathExtension.cpp Monomial.cpp compiledDA.cpp)

if(WITH_ALGEBRAICMATRIX)
  list(APPEND SRCCXXIF "AlgebraicMatrix.cpp")
endif(WITH_ALGEBRAICMATRIX)

# Build object files for static and dynamic libraries of C++ interface
add_library(dacecxx_s OBJECT ${SRCCXXIF})
add_library(dacecxx OBJECT ${SRCCXXIF})
set_target_properties(dacecxx PROPERTIES POSITION_INDEPENDENT_CODE True)
if(WIN32)
  target_compile_definitions(dacecxx PRIVATE "DACE_API=__declspec(dllexport)")
  target_compile_definitions(dacecxx_s PRIVATE "DACE_API=")
endif(WIN32)

# statically and dynamically linked C++ interface library with the corresponding c core object files compiled in
add_library(dace SHARED $<TARGET_OBJECTS:c_core> $<TARGET_OBJECTS:dacecxx>)
add_library(dace_s STATIC $<TARGET_OBJECTS:c_core_s> $<TARGET_OBJECTS:dacecxx_s>)

if(WITH_PTHREAD)
  target_link_libraries(dace PUBLIC pthread)
  target_link_libraries(dace_s PUBLIC pthread)
endif(WITH_PTHREAD)

# Build MacOS X framework
if(APPLE)
  add_library(daceFramework SHARED $<TARGET_OBJECTS:c_core> $<TARGET_OBJECTS:dacecxx>)
  set_target_properties(daceFramework PROPERTIES
                        FRAMEWORK TRUE
                        FRAMEWORK_VERSION A
                        MACOSX_FRAMEWORK_IDENTIFIER org.dacelib.dace2
#                        MACOSX_FRAMEWORK_INFO_PLIST Info.plist
                        VERSION ${DACE_MAJOR_VERSION}.${DACE_MINOR_VERSION}.${DACE_PATCH_VERSION}
                        SOVERSION 1.0.0
# open bug in CMake: PUBLIC_HEADER doesn't support include file hierarchies (https://gitlab.kitware.com/cmake/cmake/issues/16739)
#                        PUBLIC_HEADER ${PROJECT_SOURCE_DIR}/interfaces/cxx/include/dace/*.h
  )

# workaround for PUBLIC_HEADER limitation
  install(DIRECTORY ${PROJECT_SOURCE_DIR}/interfaces/cxx/include/dace DESTINATION Library/Frameworks/daceFramework.framework/Headers COMPONENT framework_headers FILES_MATCHING PATTERN "*.h")
  install(TARGETS daceFramework FRAMEWORK DESTINATION Library/Frameworks COMPONENT framework)
endif(APPLE)

install(TARGETS dace LIBRARY DESTINATION lib COMPONENT libraries
        RUNTIME DESTINATION lib COMPONENT libraries
        ARCHIVE DESTINATION lib COMPONENT libraries
        PUBLIC_HEADER DESTINATION include COMPONENT headers)
install(TARGETS dace_s LIBRARY DESTINATION lib COMPONENT libraries
        RUNTIME DESTINATION lib COMPONENT libraries
        ARCHIVE DESTINATION lib COMPONENT libraries
        PUBLIC_HEADER DESTINATION include COMPONENT headers)
install(DIRECTORY ${PROJECT_SOURCE_DIR}/interfaces/cxx/include/dace DESTINATION include COMPONENT headers FILES_MATCHING PATTERN "*.h")
