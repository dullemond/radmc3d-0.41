#ifndef SPECTRUMLABELWIDGET_H
#define SPECTRUMLABELWIDGET_H

#include <QLabel>
#include "spectrum.h"
#include <QPrinter>
#include <QPrintDialog>
class SpectrumLabelWidget : public QLabel
{
    Q_OBJECT
public:
    explicit SpectrumLabelWidget(int coordDistanceToImage,QWidget *parent = 0);
    void setSpectrum(Spectrum *spectrum){this->spectrum=spectrum;}
    Spectrum::frequencyUnit getfrequencyUnit(){return this->unit;}
    void setFrequencyUnit(Spectrum::frequencyUnit unit){this->unit=unit;}
    Spectrum::fluxUnit getFluxUnit(){return this->fluxunit;}
    void setFluxUnit(Spectrum::fluxUnit unit){this->fluxunit=unit;}
    bool isXLogOn(){return this->xLog;}
    bool isYLogOn(){return this->yLog;}
    void setXLogOn(bool log){this->xLog=log;}
    void setYLogOn(bool log){this->yLog=log;}
    void setZoomActive(bool active){this->zoomActive=active;this->setZoomOutVisible(active);}
    bool isZoomActive(){return this->zoomActive;}
    QToolButton *zoomOutButton;

private:
    bool mouseMoved;
    bool zoomActive;
    bool xLog;
    bool yLog;
    bool mousePressed;
    Spectrum::frequencyUnit unit;
    Spectrum::fluxUnit fluxunit;
    Spectrum *spectrum;
    QPoint first,last;
    int coordDistanceToImage;
protected:
    void mouseMoveEvent(QMouseEvent *event);
    void mousePressEvent(QMouseEvent *event);
    void mouseReleaseEvent(QMouseEvent *event);
signals:
    void drawCircle(int,int);
    void pixelValueChange(QString);
    void drawRect(QPoint,QPoint);
    void calculateZoom(QPoint,QPoint);
    void removeRect();
    void removeZoom();
private slots:
    void savePixmap();
    void setZoomOutVisible(bool);
    void printPixmap();
};

#endif // SPECTRUMLABELWIDGET_H
