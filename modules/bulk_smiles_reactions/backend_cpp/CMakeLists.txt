cmake_minimum_required(VERSION 3.18)

project(reaction_engine)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(RDBASE "./RDKit_git")
set(CMAKE_MODULE_PATH "${RDBASE}/Code/cmake/Modules")

set(CMAKE_PREFIX_PATH "$CONDA_PREFIX")
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED OFF)
set(Boost_USE_STATIC_RUNTIME OFF)
find_package( Boost COMPONENTS iostreams filesystem system)
# find_package( Cairo REQUIRED )

set(RDKit_LIBS RDKitChemReactions RDKitFileParsers RDKitSmilesParse RDKitDepictor
	       RDKitRDGeometryLib RDKitRDGeneral RDKitSubstructMatch RDKitSubgraphs
	       RDKitMolDraw2D RDKitGraphMol RDKitDistGeometry RDKitDistGeomHelpers
  	       RDKitMolAlign RDKitOptimizer RDKitForceField RDKitForceFieldHelpers
	       RDKitAlignment RDKitForceField  RDKitMolTransforms RDKitEigenSolvers )
set(Boost_DIR "$CONDA_PREFIX")
set(RDKit_DIR "$CONDA_PREFIX/lib/cmake/RDKit")
set(RDKit_INCLUDE_DIRS "$CONDA_PREFIX/include")
find_package(Boost REQUIRED COMPONENTS system serialization iostreams) 
find_package(rdkit REQUIRED COMPONENTS GraphMol)
include_directories ( ${RDBASE}/Code)
link_directories ( ${RDBASE}/lib )


set( LIBS ${RDKIT_LIBRARIES} Boost::iostreams  )

include_directories(${RDKIT_INCLUDE_DIR})


add_executable(reaction_engine reaction_engine.cpp)
target_link_libraries(reaction_engine ${LIBS} ${RDKit_LIBS})