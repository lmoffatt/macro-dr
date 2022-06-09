TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt
CONFIG += c++17
#CONFIG += sanitizer sanitize_address
QMAKE_CXXFLAGS -= -std=gnu++11
QMAKE_CXXFLAGS +=  -std=c++17 -Werror=return-type -ftemplate-backtrace-limit=0  -Wnon-virtual-dtor -Wnull-dereference  -fdiagnostics-show-template-tree
QMAKE_CXXFLAGS_RELEASE += -lpthread
QMAKE_CXXFLAGS_DEBUG += -lpthread


QMAKE_CXXFLAGS += -fopenmp -lstdc++fs
QMAKE_LFLAGS +=  -fopenmp  -lstdc++fs
LIBS += -lpthread -lstdc++fs
LIBS +=  -lblas  -llapack

#QMAKE_CC = /usr/bin/clang-6.0
#QMAKE_CXX = /usr/bin/clang++-6.0

CONFIG(debug, debug|release) {

} else {
  DEFINES += NDEBUG
}
SOURCES += \
    commands_evidence.cpp \
    commands_evidence_derivative.cpp \
    commands_evidence_emcee.cpp \
commands_evidence_state.cpp \
mynewcommandmanager_insert_constructor.cpp \
mynewcommandmanager_insert_commands.cpp \
main.cpp \
    mynewcommandmanager.cpp \
commands.cpp

HEADERS += \
    myTuples.h \
    mySerializer.h \
    Matrix.h \
    myDistributions.h \
    Experiment.h \
    mydataindex.h \
    mymoments.h \
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
mydynamicfunctions.h \


DISTFILES += \
    Makefile \
    cecar/scripts/m1_EvidenceDProb02_MacroDMR.txt \
    cecar/scripts/m1_EvidenceDProb_MacroDMR.txt \
    cecar/scripts/m1_EvidenceD_MacroDMR.txt \
    cecar/scripts/m1_Evidence_MacroDMR.txt \
    cecar/scripts/m1_Evidence_emcee_MacroDMR.txt \
    cecar/scripts/m1_Evidence_emcee_MacroDVR.txt \
    cecar/scripts/m1_Evidence_prob_02_MacroDMR.txt \
    cecar/scripts/m1_Evidence_prob_0_MacroDMR.txt \
    cecar/scripts/s_m1_1 \
    cecar/scripts/s_m1_2 \
    cecar/scripts/s_m1_3 \
    cecar/scripts/s_m1_4 \
    cecar/scripts/s_m1_5 \
    cecar/scripts/s_m1_6 \
    cecar/scripts/s_m1_7 \
    cecar/scripts/srun \
    data/Moffatt_Hume_2007_ATP.txt \
    run/Evidence_cecar.txt \
    run/m_1/c_m_1 \
    run/s_cecar_1 \
    simulation.txt \
    models.txt \
    mode_1.txt \
    model_1_Evidence.txt
