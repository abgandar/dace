# C++ interface source files
set(SRCCXXIF DA.cpp DACEException.cpp DAFormatter.cpp MathExtension.cpp Monomial.cpp compiledDA.cpp)
include_directories(include ${PROJECT_BINARY_DIR}/core/include ${PROJECT_SOURCE_DIR}/core/include)

if(WITH_ALGEBRAICMATRIX)
  list(APPEND SRCCXXIF "AlgebraicMatrix.cpp")
endif(WITH_ALGEBRAICMATRIX)

# Build object files for static and dynamic libraries of C++ interface
add_library(dacecxx OBJECT ${SRCCXXIF})
add_library(dacecxx_s OBJECT ${SRCCXXIF})

set_target_properties(dacecxx PROPERTIES
  POSITION_INDEPENDENT_CODE TRUE
  CXX_STANDARD 14
  CXX_STANDARD_REQUIRED ON
  INTERPROCEDURAL_OPTIMIZATION ${HAVE_IPO})
set_target_properties(dacecxx_s PROPERTIES
  CXX_STANDARD 14
  CXX_STANDARD_REQUIRED ON
  INTERPROCEDURAL_OPTIMIZATION ${HAVE_IPO})

if(WIN32)
  target_compile_definitions(dacecxx PRIVATE "DACE_API=__declspec(dllexport)")
  target_compile_definitions(dacecxx_s PRIVATE "DACE_API=")
endif(WIN32)
