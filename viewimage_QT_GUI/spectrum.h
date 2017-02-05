#ifndef SPECTRUM_H
#define SPECTRUM_H

#if QT_VERSION < 0x050000
      #include <QtGui>
#else
       #include <QtWidgets>
#endif
#include "logger.h"
#include "physicalConstants.h"
#include <algorithm>
#include <cstdio>
#include <csignal>

class Spectrum : public QObject
{
    Q_OBJECT
public:
    enum frequencyUnit{
        MICRON=0,HZ=1,EV=2,KEV=3
    };
    enum fluxUnit{
        NUFNU=0,NUFNUJY=1,FNU=2,FNUJY=3,NULNU=4,NULULSUN=5,LNU=6
    };
    explicit Spectrum(QObject *parent = 0);
    bool readSpectrumFormatted(QString spectrumFilename,QString spectrumFileEnding,Qt::CaseSensitivity cs=Qt::CaseSensitive);
    bool findSpectrumFile(QString *currentSpectrumFileName,QString spectrumFilename,QString spectrumFileEnding,Qt::CaseSensitivity cs, QStringList *list);
    QVector<double> getFrequencies(){return this->frequencies;}
    QVector<double> getFlux(){return this->flux;}
    int getNumberOfFrequencies(){return this->numberOfFrequencies;}
    void readSpectrumFormattedFromRadmc3d(QProcess *radmc3dProcess,bool *exit);
    void setXcoord(QVector<double> coord){this->xcoord=coord;}
    void setYcoord(QVector<double> coord){this->ycoord=coord;}
    QVector<double> getXcoord(){return this->xcoord;}
    QVector<double> getYcoord(){return this->ycoord;}
    void setXStart(double nr){this->xStart=nr;}
    void setYStart(double nr){this->yStart=nr;}
    void setXEnd(double nr){this->xEnd=nr;}
    void setYEnd(double nr){this->yEnd=nr;}
    double getXStart(){return this->xStart;}
    double getYStart(){return this->yStart;}
    double getYEnd(){return this->yEnd;}
    double getXEnd(){return this->xEnd;}
    bool isError(){return this->error;}
    QPixmap getPixmap(){return this->pixmap.copy();}
    void setPixmap( QPixmap pixmap){this->pixmap=pixmap;}
    int getPixmapWidth(){return this->pixmap.width();}
    int getPixmapHeight(){return this->pixmap.height();}
private:
    QPixmap pixmap;
    bool error;
    QVector<double> xcoord,ycoord;
    int format;
    double xStart;
    double yStart;
    double xEnd;
    double yEnd;
    int numberOfFrequencies;
    PhysicalConstants *physicalConstants;
    Logger *logger;
    QVector<double> frequencies;
    QVector<double> flux;
signals:
    void finished();
    void setCancelable(bool);
    
public slots:
    
};

#endif // SPECTRUM_H
