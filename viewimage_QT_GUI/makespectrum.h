#ifndef MAKESPECTRUM_H
#define MAKESPECTRUM_H

#include <QObject>
#include "imagecontroller.h"
#include "setup.h"
#include "logger.h"
#include "physicalConstants.h"


class MakeSpectrum : public QObject
{
    Q_OBJECT
public:
    explicit MakeSpectrum(ImageControllerWidget *imageControllerWidget,QObject *parent = 0);
    QString getSpectrumCommand();
    void sendMakeSpectrumCommandToRadmc3d(QProcess *process);
    void sendWriteSpectrumCommandToRadmc3d(QProcess *process);
private:
    QFile *file;
    QFile *outputfile;
    PhysicalConstants *physicalConstants;
    Logger *logger;
    SetUp *setup;
    ImageControllerWidget *imageControllerWidget;
signals:
    
public slots:
    
};

#endif // MAKESPECTRUM_H
