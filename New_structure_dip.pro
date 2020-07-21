TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
        Grid.cpp \
        Matrix.cpp \
        Setting.cpp \
        Solver.cpp \
        Tensor.cpp \
        Vector.cpp \
        main.cpp

HEADERS += \
    Grid.h \
    Matrix.h \
    Setting.h \
    Solver.h \
    Tensor.h \
    Vector.h
