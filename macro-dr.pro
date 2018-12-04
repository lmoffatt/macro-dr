TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt
CONFIG += c++1z
#CONFIG += sanitizer sanitize_address
QMAKE_CXXFLAGS += -std=c++17 -Werror=return-type -ftemplate-backtrace-limit=0  -Wnon-virtual-dtor -Wnull-dereference
QMAKE_CXXFLAGS_RELEASE += -lpthread
QMAKE_CXXFLAGS_DEBUG += -lpthread


QMAKE_CXXFLAGS += -fopenmp
QMAKE_LFLAGS +=  -fopenmp
LIBS += -lpthread
LIBS +=  -lblas  -llapack

#QMAKE_CC = /usr/bin/clang-6.0
#QMAKE_CXX = /usr/bin/clang++-6.0

CONFIG(debug, debug|release) {

} else {
  DEFINES += NDEBUG
}
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
    qlikelihood.h \
    mytests.h \
    myoptimization.h \
    qsimulation.h \
    mylikelihood.h \
    myprobabilitytest.h \
    myoperators.h \
    myderivatives.h \
    matrixderivative.h \
    qmodel_derivative.h \
    likelihood_markov_process_derivative.h \
    qlikelihood_derivative.h \
    myparameters_derivative.h \
    mydistributions_derivative.h \
    mydynamicfunctions.h

DISTFILES += \
    simulation.txt \
    models.txt
