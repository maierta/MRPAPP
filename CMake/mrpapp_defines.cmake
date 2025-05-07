# // Copyright (C) 2024 UT-Battelle, LLC
# // All rights reserved.
# //
# // See LICENSE for terms of usage.
# //

################################################################################
# Adds a definition to the global 'have definitions' string.
function(mrpapp_add_define definition)
  if(ARGN)
    set_property(GLOBAL APPEND PROPERTY MRPAPP_CXX_DEFINITIONS "${definition} ${ARGN}")
  else()
    set_property(GLOBAL APPEND PROPERTY MRPAPP_CXX_DEFINITIONS "${definition}")
  endif()
endfunction()

################################################################################
# Generates in the build directory the haves_defines.hpp that contains all 'haves preprocessor
# definitions'.
function(mrpapp_write_definitions_file)
  get_property(MRPAPP_CXX_DEFINITIONS_VAR GLOBAL PROPERTY MRPAPP_CXX_DEFINITIONS)

  list(SORT MRPAPP_CXX_DEFINITIONS_VAR)
  list(REMOVE_DUPLICATES MRPAPP_CXX_DEFINITIONS_VAR)
  list(REMOVE_ITEM MRPAPP_CXX_DEFINITIONS_VAR "")

  set(mrpapp_cxx_defines_string "")
  message("cxx definitions: ${MRPAPP_CXX_DEFINITIONS_VAR}")
  foreach(def ${MRPAPP_CXX_DEFINITIONS_VAR})
    string(CONCAT mrpapp_cxx_defines_string
        "${mrpapp_cxx_defines_string}"
        "#ifndef ${def}\n"
        " #define ${def} ${${def}_define}\n"
        "#endif\n\n")
  endforeach()

  configure_file("${PROJECT_SOURCE_DIR}/defines.hpp.in"
    "${CMAKE_BINARY_DIR}/defines.hpp")
endfunction()
