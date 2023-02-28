cmake_minimum_required(VERSION 3.20)
project(gCore)

configure_file(
        "${PROJECT_SOURCE_DIR}/config.h.in"
        "${PROJECT_SOURCE_DIR}/config.h"
)

option(USE_MEMTRACK "Use memtrack library" OFF)

if (USE_MEMTRACK)
    include_directories("${PROJECT_SOURCE_DIR}/memtrack")
    add_subdirectory(memtrack)
    set(EXTRA_LIBS ${EXTRA_LIBS} memtrack)
endif(USE_MEMTRACK)

set(CMAKE_CXX_STANDARD 17)


set(CMAKE_CXX_FLAGS_DEBUG  "${CMAKE_C_FLAGS_DEBUG} -DMY_DEBUG")



add_executable(gCore main.cpp Graphs/Graph.cpp Graphs/Graph.h Header.h
        Graphs/MultilayerGraph.cpp Graphs/MultilayerGraph.h Utils/ArrayUtils.h Utils/Timer.h
        Graphs/BipartiteGraph.cpp Graphs/BipartiteGraph.h Graphs/UndirectedGraph.cpp Graphs/UndirectedGraph.h
        Structures/IntLinearHeap.h
        Structures/IntLinearHeap.cpp Core/CoreDec.cpp Core/CoreDec.h Core/KPCore.cpp Core/KPCore.h
        Core/PTreeBuilder.cpp Core/PTreeBuilder.h Core/PTree.h Utils/FracUtils.h Utils/Viz.h
        Utils/CollectionUtils.h Core/KPTree.h Core/RCD.cpp Core/RCD.h Utils/StringUtils.h
        Core/PTree.cpp Core/Frac2IntPri.cpp Core/Frac2IntPri.h Core/KPTree.cpp
        Core/KPTreeBuilder.cpp Core/KPTreeBuilder.h Utils/FileUtils.h Utils/MemoryUtils.h
        Core/ExtendedDCC.cpp Core/ExtendedDCC.h KPCEffTest.cpp KPCEffTest.h memtrack/memtrack.cpp memtrack/memtrack.h
        Structures/DataBuf.h Core/CoreIndex.h KPCEvalTest.h KPCEvalTest.cpp KPCCaseStudy.h CmdParser.hpp Params.h
        KPTreeUtils.cpp KPTreeUtils.h KPVecGenerator.h KPVecGenerator.cpp KPCCaseStudy.cpp
        )

target_link_libraries(gCore ${EXTRA_LIBS})
