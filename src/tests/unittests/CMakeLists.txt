# Add any new tests to this list or the one below!
set( SOURCES_UNITTESTS_TESTS
    unittest-gas-electrostatics-mbpol.cpp
    unittest-dispersion-gas-cutoff.cpp
    unittest-dispersion-pme.cpp
    unittest-gas-nbodyterms-mbpol.cpp
    unittest-systools.cpp
    unittest-system.cpp
    unittest-mathtools.cpp
    unittest-gas-nopolarization-findif.cpp
    unittest-gas-withpolarization-findif.cpp
    unittest-pme-nopolarization.cpp
    unittest-pme-nopolarization-findif.cpp
    unittest-pme-withpolarization.cpp
    unittest-pme-withpolarization-findif.cpp
    unittest-gamma.cpp
    unittest-pbc-1b-mbpol-findif.cpp
    unittest-pbc-2bpoly-mbpol-findif.cpp
    unittest-pbc-dispersion-mbpol-findif.cpp
    unittest-pbc-3bpoly-mbpol-findif.cpp
    unittest-pbc-electrostatics-mbpol-findif.cpp
    unittest-h2o-f-dimer.cpp
    unittest-h2o-cl-dimer.cpp
    unittest-h2o-br-dimer.cpp
    unittest-h2o-virial.cpp
    unittest-h2o-i-dimer.cpp
    unittest-h2o-li-dimer.cpp
    unittest-h2o-na-dimer.cpp
    unittest-h2o-k-dimer.cpp
    unittest-h2o-rb-dimer.cpp
    unittest-h2o-cs-dimer.cpp
    unittest-h2o-h2o-li-trimer.cpp
    unittest-h2o-h2o-na-trimer.cpp
    unittest-h2o-h2o-k-trimer.cpp
    unittest-h2o-h2o-rb-trimer.cpp
    unittest-h2o-h2o-cs-trimer.cpp
    unittest-co2-monomer.cpp
    unittest-ch4-monomer.cpp
    unittest-dummy-monomer.cpp
)

#
# Stuff below shouldn't need to be modified
#

include_directories(${PROJECT_SOURCE_DIR}
                    ${PROJECT_SOURCE_DIR}/../external/
                    ${PROJECT_SOURCE_DIR}/potential/electrostatics )

set( SOURCES_UNITTESTS_MAIN unittest-main.cpp )

# A function to check that the given hard-coded list contains all headers + sources in the given folder.
function(CheckFileList LIST_VAR FOLDER)
  set(MESSAGE " should be added to the variable ${LIST_VAR}")
  set(MESSAGE "${MESSAGE} in ${CMAKE_CURRENT_LIST_FILE}\n")
  file(GLOB GLOBBED_LIST "${FOLDER}/*.cpp"
                         "${FOLDER}/*.hpp")
  list(REMOVE_ITEM GLOBBED_LIST ${${LIST_VAR}})
  foreach(EXTRA_ITEM ${GLOBBED_LIST})
    string(REPLACE "${CATCH_DIR}/" "" RELATIVE_FILE_NAME "${EXTRA_ITEM}")
    message(AUTHOR_WARNING "The file \"${RELATIVE_FILE_NAME}\"${MESSAGE}")
  endforeach()
endfunction()

set( UNITTEST_DIR ${PROJECT_SOURCE_DIR}/tests/unittests )
set( HEADER_DIR ${PROJECT_SOURCE_DIR}/../external/Catch2/single_include )

# Check if all sources are listed; warn if not.
set( SOURCES_ALL
    ${SOURCES_UNITTESTS_MAIN}
    ${SOURCES_UNITTESTS_TESTS}
)
foreach( name ${SOURCES_ALL} )
    list( APPEND SOURCES_ALL_PATH ${UNITTEST_DIR}/${name} )
endforeach()
CheckFileList( SOURCES_ALL_PATH ${UNITTEST_DIR} )

# Create serial target names.
string( REPLACE ".cpp" "" BASENAMES_UNITTESTS_TESTS "${SOURCES_UNITTESTS_TESTS}" )
set( TARGETS_UNITTESTS_TESTS ${BASENAMES_UNITTESTS_TESTS} )

set( TARGETS_ALL ${TARGETS_UNITTESTS_TESTS} CatchMain )

# Define program targets:
add_library( CatchMain OBJECT ${UNITTEST_DIR}/${SOURCES_UNITTESTS_MAIN} ${HEADER_DIR}/catch.hpp )

foreach( name ${TARGETS_UNITTESTS_TESTS} )
    add_executable( ${name} ${name}.cpp $<TARGET_OBJECTS:CatchMain> ${HEADER_DIR}/catch.hpp )
    target_link_libraries( ${name} PRIVATE mbxlib fftw::fftw)
    add_dependencies( ${name} mbxlib )
    add_test(NAME ${name}
             COMMAND ${name})
endforeach()

foreach( name ${TARGETS_ALL} )
    target_include_directories( ${name} PRIVATE ${HEADER_DIR} )

    set_property(TARGET ${name} PROPERTY CXX_STANDARD 11)

    # GCC specific warning go here
    if ( CMAKE_CXX_COMPILER_ID MATCHES "GCC" )
        # Actually keep these
        target_compile_options( ${name}  PRIVATE -Wweak-vtables -Wexit-time-destructors -Wglobal-constructors -Wmissing-noreturn -Wno-int-in-bool-context)
    endif()
    # Clang specific warning go here
    if ( CMAKE_CXX_COMPILER_ID MATCHES "Clang" )
        # Actually keep these
        target_compile_options( ${name}  PRIVATE -Wweak-vtables -Wexit-time-destructors -Wglobal-constructors)
    endif()
    if ( CMAKE_CXX_COMPILER_ID MATCHES "MSVC" )
        target_compile_options( ${name}  PRIVATE /W4 /w44265 /WX )
    endif()
endforeach()
