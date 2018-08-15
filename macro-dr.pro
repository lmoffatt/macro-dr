TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt
CONFIG += c++1z
QMAKE_CXXFLAGS += -std=c++17 -Werror=return-type -ftemplate-backtrace-limit=0
LIBS +=  -lblas  -llapack

#QMAKE_CC = /usr/bin/clang-6.0
#QMAKE_CXX = /usr/bin/clang++-6.0
SOURCES += \
    main.cpp

HEADERS += \
    myTuples.h \
    mySerializer.h \
    Matrix.h \
    myDistributions.h \
    Experiment.h \
    qmodel.h \
    mySerializer.h \
    mygrammar.h \
    mycompilation.h \
    mynewcommandmanager.h \
    myscriptmanager.h \
    myoptional.h \
    myfields.h \
    mytypetraits.h \
    mydataframe.h \
    measure_markov_process.h \
    mymath.h \
    mysmartpointerstools.h \
    commands.h \
    myparameters.h \
    mycontainer.h \
    likelihood_markov_process.h \
    mydata.h \
    myevidence.h \
    qlikelihood.h

DISTFILES += \
    simulation.txt
