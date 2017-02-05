#ifndef MAKEIMAGE_H
#define MAKEIMAGE_H

#include <QObject>
#include "imagecontroller.h"
#include "setup.h"
#include "logger.h"
#include "physicalConstants.h"

class MakeImage : public QObject
{
    Q_OBJECT
public:
    explicit MakeImage(ImageControllerWidget *imageControllerWidget,QObject *parent = 0);
    QString getImageCommand();
    void sendMakeImageCommandToRadmc3d(QProcess *process);
    void sendWriteImageCommandToRadmc3d(QProcess *process);


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

#endif // MAKEIMAGE_H
