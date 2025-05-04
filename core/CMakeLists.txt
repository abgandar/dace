set(SRCCCORE daceaux.c dacebasic.c daceerror.c daceeval.c daceinit.c daceio.c dacemath.c dacememory.c dacenorm.c)
include_directories(include ${PROJECT_BINARY_DIR}/core/include contrib/include)

add_library(c_core_s OBJECT ${SRCCCORE})
add_library(c_core OBJECT ${SRCCCORE})

set_target_properties(c_core PROPERTIES
  POSITION_INDEPENDENT_CODE TRUE
  C_STANDARD 99
  C_STANDARD_REQUIRED ON
  INTERPROCEDURAL_OPTIMIZATION ${HAVE_IPO})
set_target_properties(c_core_s PROPERTIES
  C_STANDARD 99
  C_STANDARD_REQUIRED ON
  INTERPROCEDURAL_OPTIMIZATION ${HAVE_IPO})

if(WIN32)
  target_compile_definitions(c_core PRIVATE "DACE_API=__declspec(dllexport)")
  target_compile_definitions(c_core_s PRIVATE "DACE_API=")
endif(WIN32)
