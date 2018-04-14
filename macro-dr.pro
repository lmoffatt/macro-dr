TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt
CONFIG += c++1z

SOURCES += \
    main.cpp

HEADERS += \
    myCommandManagement.h \
    myTuples.h \
    CommandManager.h \
    myInputSerializer.h \
    myOutputSerializer.h \
    Matrix.h \
    myorderoperators.h \
    Markov.h \
    myDistributions.h \
    Experiment.h \
    simulation.h \
    qmodel.h
