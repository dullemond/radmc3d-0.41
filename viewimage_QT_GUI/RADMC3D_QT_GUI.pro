SOURCES += \
    mainwindow.cpp \
    imagecontroller.cpp \
    main.cpp \
    logger.cpp \
    setup.cpp \
    physicalConstants.cpp \
    cartesian.cpp \
    polar.cpp \
    octtree.cpp \
    amrgrid.cpp \
    image.cpp \
    makeimage.cpp \
    imagelabelwidget.cpp \
    makespectrum.cpp \
    spectrum.cpp \
    spectrumlabelwidget.cpp \
    spectrumsetting.cpp \
    unformattedsetting.cpp
HEADERS += \
    mainwindow.h \
    imagecontroller.h \
    physicalConstants.h \
    logger.h \
    setup.h \
    cartesian.h \
    polar.h \
    octtree.h \
    amrgrid.h \
    image.h \
    makeimage.h \
    imagelabelwidget.h \
    makespectrum.h \
    spectrum.h \
    spectrumlabelwidget.h \
    spectrumsetting.h \
    unformattedsetting.h

RESOURCES += \
    images.qrc
TRANSLATIONS=languages/Translation_en_EN.ts \
             languages/Translation_de_DE.ts


CONFIG += release
CONFIG += c++11

CONFIG-=app_bundle

CONFIG(debug, debug|release) {

DESTDIR=$$PWD-build__Debug
}
CONFIG(release, debug|release) {
DESTDIR=$$PWD-build

}
equals(QT_MAJOR_VERSION, 5){
   QT += widgets
   QT += concurrent
  # LIBS += -lQt5Concurrent
    QT += printsupport
}
equals(QT_MAJOR_VERSION, 4) {
   QT += gui
}
TEMPLATE = app
OBJECTS_DIR = $$DESTDIR/.obj
MOC_DIR = $$DESTDIR/.moc

QMAKE_POST_LINK += perl $$PWD/install.perl



