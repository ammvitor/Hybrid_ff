TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt
QMAKE_CXXFLAGS += -std=c++11

SOURCES += main.cpp \
    aminoacid.cpp \
    amber_parm_parser.cpp \
    protein.cpp \
    dssp_parser.cpp \
    run_engine.cpp \
    writer.cpp

HEADERS += \
    aminoacid.h \
    amber_parm_parser.h \
    protein.h \
    dssp_parser.h \
    run_engine.h \
    writer.h

