# @author Adeel Asghar <adeel.asghar@liu.se>

QMAKE_CC  = clang
QMAKE_CXX = clang++
QMAKE_LINK = clang++

OPENMODELICAHOME = /home/uedashige/work/OM/OMBuild/20200101_2/OpenModelica/build
host_short = x86_64-linux-gnu

LIBS += -L/home/uedashige/work/OM/OMBuild/20200101_2/OpenModelica/build/lib/$$host_short/omc -lOMPlot -lomqwt -lfmilib_shared -L$$OMEDIT_ROOT/OMEditLIB/Debugger/Parser -lGDBMIParser -lomantlr3 -Wl,-z,origin -Wl,-rpath,\'\$$ORIGIN/../lib/x86_64-linux-gnu/omc\' -Wl,-rpath,\'\$$ORIGIN/../lib\' -Wl,-rpath,\'\$$ORIGIN\' -lOpenModelicaCompiler -lOpenModelicaRuntimeC -lomcgc -L/home/uedashige/work/OM/OMBuild/20200101_2/OpenModelica/build/lib/$$host_short/omc -losg -losgViewer -losgDB -losgGA -lOpenThreads -lomopcua -L/home/uedashige/work/OM/OMBuild/20200101_2/OpenModelica/build/lib -lOMSimulator -lqjson 

QMAKE_CXXFLAGS_RELEASE -= -O1
QMAKE_CXXFLAGS_RELEASE -= -O2

QMAKE_CFLAGS = -g -O2 
QMAKE_CXXFLAGS =  
QMAKE_LFLAGS +=  -Wl,-rpath-link,/home/uedashige/work/OM/OMBuild/20200101_2/OpenModelica/build/lib/x86_64-linux-gnu/omc

# required for backtrace
# On unix we use backtrace of execinfo.h which requires -rdynamic
# The symbol names may be unavailable without the use of special linker
# options.  For systems using the GNU linker, it is necessary to use
# the -rdynamic linker option.  Note that names of "static" functions
# are not exposed, and won't be available in the backtrace.
CONFIG(release, debug|release) {
  QMAKE_CXXFLAGS += -g
  QMAKE_LFLAGS_RELEASE = -rdynamic
}
equals(QT_ARCH, i386)|equals(QT_ARCH, i486)|equals(QT_ARCH, i586)|equals(QT_ARCH, i686) { # 32-bit
  LIBS += -latomic -lboost_atomic
}
