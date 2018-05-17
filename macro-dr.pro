TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt
CONFIG += c++1z
QMAKE_CXXFLAGS += -std=c++17 -ftemplate-backtrace-limit=0
LIBS +=  -lblas  -llapack


SOURCES += \
    main.cpp

HEADERS += \
    myCommandManagement.h \
    myTuples.h \
    CommandManager.h \
    mySerializer.h \
    Matrix.h \
    myorderoperators.h \
    Markov.h \
    myDistributions.h \
    Experiment.h \
    simulation.h \
    qmodel.h \
    mySerializer.h \
    myreadwriter.h \
    mygrammar.h \
    mycompilation.h \
    mynewcommandmanager.h \
    myscriptmanager.h \
    myoptional.h \
    myfields.h \
    mytypetraits.h \
    mydataframe.h \
    newcommandmanager.h \
    measure_markov_process.h \
    mymath.h \
    mysmartpointerstools.h \
    commands.h
