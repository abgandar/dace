cmake_minimum_required (VERSION 3.10)

# since MacOS 10.14 (XCode 10.0), default includes are no longer installed to /usr/include or /usr/local/include
# the SDK to use must be specified explicitly in CMAKE_OSX_SYSROOT as it seems at least cmake 3.12.4 does not find it automatically if left empty
# this must happen before project()
if(APPLE AND NOT DEFINED ENV{SDKROOT})
  execute_process(COMMAND xcodebuild -sdk macosx -version Path OUTPUT_VARIABLE RES OUTPUT_STRIP_TRAILING_WHITESPACE)
  set(CMAKE_OSX_SYSROOT "${RES}" CACHE PATH "Build against headers and libraries in this SDK.")
endif(APPLE AND NOT DEFINED ENV{SDKROOT})

# global and OS specific settings
cmake_policy(SET CMP0048 NEW)
if(WIN32)
  cmake_policy(SET CMP0054 NEW)
elseif(APPLE)
  cmake_policy(SET CMP0042 NEW)
  cmake_policy(SET CMP0068 NEW)
  # build libraries use full hard-coded binary path, install libraries replace it by @rpath
  set(CMAKE_MACOSX_RPATH OFF)
  set(CMAKE_INSTALL_NAME_DIR "@rpath")
endif(WIN32)


project(DACE
  VERSION 2.1.0
  DESCRIPTION "Differential Algebra Core Engine"
  HOMEPAGE_URL "https://www.dacelib.org"
  LANGUAGES C CXX
)

# CMake includes
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
include(CMakeDependentOption)
include(CheckSafeStrings)
include(CMakePackageConfigHelpers)
include(CheckIPOSupported)

# Set generally understood compiler flags (enforcing -O2 for release builds)
add_compile_options("$<$<CONFIG:Release>:-O2>" "$<$<CONFIG:RelWithDebInfo>:-O2>")

# Set compiler specific flags
if(CMAKE_C_COMPILER_ID STREQUAL "GNU")
  add_compile_options("-Wall" "$<$<BOOL:${WITH_FPSTRICT}>:-frounding-math>" "$<$<BOOL:${WITH_FPSTRICT}>:-fsignaling-nans>")
elseif(CMAKE_C_COMPILER_ID MATCHES ".*Clang")        # Also catch IDs like "AppleClang"
  add_compile_options("-Wall" "$<$<BOOL:${WITH_FPSTRICT}>:-ffp-model=strict>")
elseif(CMAKE_C_COMPILER_ID STREQUAL "MSVC")
  add_compile_options("/W4" "$<$<BOOL:${WITH_FPSTRICT}>:/fp:strict>")
endif(CMAKE_C_COMPILER_ID STREQUAL "GNU")

# Boolean options
option(WITH_PTHREAD "Enable threading via pthread library" OFF)
option(WITH_FPSTRICT "Enable strict IEEE-754 floating point handling" ON)
option(WITH_ALGEBRAICMATRIX "Include the AlgebraicMatrix type (experimental)" OFF)

# Tri-state option for memory model and the associated choices
set(DACE_MEMORY_MODEL DYNAMIC CACHE STRING "DACE memory model")
set_property(CACHE DACE_MEMORY_MODEL PROPERTY STRINGS DYNAMIC STATIC HYBRID)
string(TOUPPER "${DACE_MEMORY_MODEL}" DACE_MEMORY_MODEL_U)

# Default build type option if not set already
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Build type: Debug Release RelWithDebInfo MinSizeRel" FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release RelWithDebInfo MinSizeRel)
endif(NOT CMAKE_BUILD_TYPE)
string(TOUPPER "${CMAKE_BUILD_TYPE}" CMAKE_BUILD_TYPE_U)

# Debug flag
if(CMAKE_BUILD_TYPE_U STREQUAL "DEBUG")
  set(WITH_DEBUG 1)
endif(CMAKE_BUILD_TYPE_U STREQUAL "DEBUG")

# Static memory management default sizes
# If you do not know what these values are, you probably do not want to touch them!
if(DACE_MEMORY_MODEL_U STREQUAL "STATIC")
  # Pick any maximum NO, NV combination to support, then set these values like this:
  # NMMAX = (NO+NV)!/NO!/NV!         LIAMAX = (NO+1)^((NV+1)/2)
  # For example, for NO=10 and NV=10 set DACE_STATIC_NMMAX=184756 and DACE_STATIC_LIAMAX=161051
  # The following default values are sized for NO=10 and NV=10
  set(DACE_STATIC_NOMAX 10 CACHE STRING "(EXPERT) Maximum order.")
  set(DACE_STATIC_NVMAX 10 CACHE STRING "(EXPERT) Maximum number of variables.")
  set(DACE_STATIC_NMMAX 184756 CACHE STRING "(EXPERT) Maximum number of monomials [(NO+NV)!/NO!/NV!].")
  set(DACE_STATIC_LIAMAX 161051 CACHE STRING "(EXPERT) Maximum addressing array size [(NO+1)^((NV+1)/2)].")
  set(DACE_STATIC_VAR_SIZE 100 CACHE STRING "(EXPERT) Maximum number of available DA objects.")
  set(DACE_STATIC_MEM_SIZE 18475600 CACHE STRING "(EXPERT) Maximum memory size [NMMAX*VAR_SIZE].")
endif(DACE_MEMORY_MODEL_U STREQUAL "STATIC")

# Check platform specific features and pass options by configuring config.h file in core
check_ipo_supported(RESULT HAVE_IPO LANGUAGES C CXX)
check_safe_strings(HAVE_SAFE_STRINGS)
configure_file(core/include/dace/config.h.in core/include/dace/config.h)

# process the subdirectories with actual code
add_subdirectory(core/contrib)
add_subdirectory(core)
add_subdirectory(interfaces/cxx)

# statically and dynamically linked DACE library with C core and C++ interface
add_library(dace SHARED $<TARGET_OBJECTS:c_core> $<TARGET_OBJECTS:c_core_contrib> $<TARGET_OBJECTS:dacecxx>)
add_library(dace_s STATIC $<TARGET_OBJECTS:c_core_s> $<TARGET_OBJECTS:c_core_contrib_s> $<TARGET_OBJECTS:dacecxx_s>)

# versioning and target properties (c_std_99 seems broken on 3.31.7 on Mac with XCode 16.3)
set_target_properties(dace PROPERTIES
  POSITION_INDEPENDENT_CODE TRUE
  INTERPROCEDURAL_OPTIMIZATION ${HAVE_IPO}
#  INTERFACE_COMPILE_FEATURES "cxx_std_14;c_std_99"
  INTERFACE_COMPILE_FEATURES "cxx_std_14"
  VERSION "${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.${PROJECT_VERSION_PATCH}"
  SOVERSION ${PROJECT_VERSION_MAJOR}
  INTERFACE_PROJECT_VERSION_MAJOR ${PROJECT_VERSION_MAJOR})
set_target_properties(dace_s PROPERTIES
  INTERPROCEDURAL_OPTIMIZATION ${HAVE_IPO}
#  INTERFACE_COMPILE_FEATURES "cxx_std_14;c_std_99"
  INTERFACE_COMPILE_FEATURES "cxx_std_14"
  VERSION "${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.${PROJECT_VERSION_PATCH}"
  SOVERSION ${PROJECT_VERSION_MAJOR}
  INTERFACE_dace_s_MAJOR_VERSION ${PROJECT_VERSION_MAJOR})

# hard code all required include directories when building from source (C++ includes, core includes, config.h include)
target_include_directories(dace INTERFACE $<INSTALL_INTERFACE:include>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/interfaces/cxx/include>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/core/include>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/core/include>
)
target_include_directories(dace_s INTERFACE $<INSTALL_INTERFACE:include>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/interfaces/cxx/include>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/core/include>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/core/include>
)

if(WITH_PTHREAD)
  target_link_libraries(dace PUBLIC pthread)
  target_link_libraries(dace_s PUBLIC pthread)
endif(WITH_PTHREAD)

# install DACE libraries and required headers
install(TARGETS dace EXPORT dace
        LIBRARY DESTINATION lib COMPONENT libraries
        RUNTIME DESTINATION lib COMPONENT libraries
        ARCHIVE DESTINATION lib COMPONENT libraries
        PUBLIC_HEADER DESTINATION include COMPONENT headers)
install(TARGETS dace_s EXPORT dace
        LIBRARY DESTINATION lib COMPONENT libraries
        RUNTIME DESTINATION lib COMPONENT libraries
        ARCHIVE DESTINATION lib COMPONENT libraries
        PUBLIC_HEADER DESTINATION include COMPONENT headers)
install(DIRECTORY ${PROJECT_SOURCE_DIR}/interfaces/cxx/include/dace DESTINATION include
        COMPONENT headers FILES_MATCHING PATTERN "*.h")
install(DIRECTORY ${PROJECT_SOURCE_DIR}/core/include/ DESTINATION include/ COMPONENT headers FILES_MATCHING PATTERN "*.h")
install(FILES ${PROJECT_BINARY_DIR}/core/include/dace/config.h DESTINATION include/dace COMPONENT headers)

# Automatically create cmake package files to find the DACE
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/dace/daceConfigVersion.cmake"
  VERSION "${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.${PROJECT_VERSION_PATCH}"
  COMPATIBILITY AnyNewerVersion
)

export(EXPORT dace
  FILE "${CMAKE_CURRENT_BINARY_DIR}/dace/dace.cmake"
  NAMESPACE dace::
)

configure_file("cmake/daceConfig.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}/dace/daceConfig.cmake"
  COPYONLY
)

# install cmake package files
install(EXPORT dace NAMESPACE dace:: DESTINATION lib/cmake/dace COMPONENT cmake)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/dace/daceConfig.cmake" "${CMAKE_CURRENT_BINARY_DIR}/dace/daceConfigVersion.cmake" DESTINATION lib/cmake/dace COMPONENT cmake)

# Doxygen documentation
find_package(Doxygen OPTIONAL_COMPONENTS dot mscgen dia)

set(DOXYGEN_OUTPUT_DIRECTORY "docs")        # relative to CURRENT_BINARY_DIR
set(DOXYGEN_INCLUDE_PATH "${CMAKE_CURRENT_BINARY_DIR}/core/include")    # to find config.h with CMake settings
set(DOXYGEN_STRIP_FROM_INC_PATH "${CMAKE_CURRENT_SOURCE_DIR}/interfaces/cxx/include" "${CMAKE_CURRENT_SOURCE_DIR}/core/include" "${CMAKE_CURRENT_BINARY_DIR}/core/include")
set(DOXYGEN_STRIP_FROM_PATH "${CMAKE_CURRENT_SOURCE_DIR}" "${CMAKE_CURRENT_BINARY_DIR}")
set(DOXYGEN_EXCLUDE "core/contrib/netlib" )
if(NOT WITH_ALGEBRAICMATRIX)
  set(DOXYGEN_EXCLUDE_PATTERNS "*AlgebraicMatrix*" )
endif(NOT WITH_ALGEBRAICMATRIX)
#set(DOXYGEN_EXTRACT_ALL YES)    # for DACE developers: extract everything including private, internal APIs
set(DOXYGEN_JAVADOC_AUTOBRIEF YES)
set(DOXYGEN_BUILTIN_STL_SUPPORT YES)
set(DOXYGEN_MACRO_EXPANSION YES)
set(DOXYGEN_PREDEFINED "REF(x)=*x" "DACE_API=")   # help doxygen match declarations in header to definitions
set(DOXYGEN_HTML_COLORSTYLE "LIGHT")
set(DOXYGEN_HTML_DYNAMIC_SECTIONS YES)
set(DOXYGEN_HTML_PROJECT_COOKIE org.dacelib.dace)
set(DOXYGEN_GENERATE_TREEVIEW YES)
set(DOXYGEN_FULL_SIDEBAR YES)
set(DOXYGEN_TIMESTAMP YES)

if(DOXYGEN_FOUND)
  doxygen_add_docs(docs
    "docs"
    "interfaces/cxx"
    "core"
  )
endif(DOXYGEN_FOUND)