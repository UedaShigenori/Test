# @author Martin Sjölund <martin.sjolund@liu.se>

QMAKE_CC  = clang
QMAKE_CXX = clang++
QMAKE_LINK = clang++

OPENMODELICAHOME = /home/uedashige/work/OM/OMBuild/20200101_2/OpenModelica/build
host_short = x86_64-linux-gnu

QMAKE_CXXFLAGS_RELEASE -= -O1
QMAKE_CXXFLAGS_RELEASE -= -O2

QMAKE_CFLAGS = -g -O2 
QMAKE_CXXFLAGS =  
QMAKE_LFLAGS +=  -Wl,-rpath-link,/home/uedashige/work/OM/OMBuild/20200101_2/OpenModelica/build/lib/x86_64-linux-gnu/omc

CONFIG += osg
