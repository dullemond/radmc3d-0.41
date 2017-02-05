#ifndef IMAGELABELWIDGET_H
#define IMAGELABELWIDGET_H
#if QT_VERSION < 0x050000
      #include <QtGui>
#else
       #include <QtWidgets>
#endif
#include <QPrinter>
#include <QPrintDialog>
#include "image.h"
class ImageLabelWidget : public QLabel
{
    Q_OBJECT
public:
    explicit ImageLabelWidget(int coordDistanceToImag,QWidget *parent = 0);
    void setImage(Image *image){this->image=image;}
    void setCoordDistanceToImage(int coordDistanceToImage){this->coordDistanceToImage=coordDistanceToImage;}
    void setPhiAndInclination(double phi,double inclination);
    bool isPlot(){return this->plot;}
private:
    void setPlot(bool active){this->plot=active;}
    QPoint first,last;
    bool mousePressed;
    Image *image;
    int coordDistanceToImage;
    bool plot;
    bool zoomActive;
    bool mouseMoved;
    double inclination;
    QAction *zoomAction;
    QAction *rotateAction;
    double phi;
protected:
        void mouseMoveEvent(QMouseEvent *event);
        void mousePressEvent(QMouseEvent *event);
        void mouseReleaseEvent(QMouseEvent *event);
signals:
        void pixelValueChange(QString);
        void drawRect(QPoint,QPoint);
        void removeRect();
        void calculateZoom(QPoint,QPoint);
        void phiAndinclinactionValueChanged(double phi,double inclination);
        void reloadImage();
private slots:
        void zoomClicked(bool);
        void plotClicked(bool);
        void savePixmap();
        void printPixmap();
        void rotateClicked(bool);

public slots:
    
};

#endif // IMAGELABELWIDGET_H
