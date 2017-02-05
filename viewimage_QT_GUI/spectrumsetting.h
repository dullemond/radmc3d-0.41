#ifndef SPECTRUMSETTING_H
#define SPECTRUMSETTING_H

#include <QObject>

class spectrumSetting : public QObject
{
public:
    explicit spectrumSetting(QObject *parent = 0);
    void setXStart(double number){this->xStart=number;}
    void setYStart(double number){this->yStart=number;}
    void setXEnd(double number){this->xEnd=number;}
    void setYEnd(double number){this->yEnd=number;}
    double getXStart(){return this->xStart;}
    double getYStart(){return this->yStart;}
    double getXEnd(){return this->xEnd;}
    double getYEnd(){return this->yEnd;}
private:
    double xStart,xEnd;
    double yStart,yEnd;
signals:

public slots:

};

#endif // SPECTRUMSETTING_H
