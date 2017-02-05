#ifndef IMAGE_H
#define IMAGE_H

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
#include "unformattedsetting.h"
using namespace  std;

class Image : public QObject
{
    Q_OBJECT
public:
    explicit Image(QObject *parent = 0);
    bool readImageFormatted(QString imageFilename,QString imageFileEnding,Qt::CaseSensitivity cs=Qt::CaseSensitive);
    bool readImageUnformatted(QString imageFilename,QString imageFileEnding,Qt::CaseSensitivity cs=Qt::CaseSensitive);
    bool isLocalObserver(){return this->localObserver;}
    void setLocalObserver(bool local){this->localObserver=local;}
    bool findImageFile(QString* currentImageFileName,QString imageFilename,QString imageFileEnding,Qt::CaseSensitivity cs, QStringList *list);
    void logImage();
    void setByteOrder(QDataStream::ByteOrder order){this->byteOrder=order;}
    QDataStream::ByteOrder getByteOrder(){return this->byteOrder;}
    void setNumberOfPixelX(int number){this->numberOfPixelX=number;}
    void setNumberOfPixelY(int number){this->numberOfPixelY=number;}
    int getNumberOfPixelX(){return this->numberOfPixelX;}
    int getNumberOfPixelY(){return this->numberOfPixelY;}
    void setNumberOfImages(int number){this->numberOfImages=number;}
    int getNumberOfImages(){return this->numberOfImages;}
    void setSizeOfPixelX(double size){this->sizeOfPixelX=size;}
    void setSizeOfPixelY(double size){this->sizeOfPixelY=size;}
    double getSizeOfPixelX(){return this->sizeOfPixelX;}
    double getSizeOfPixelY(){return this->sizeOfPixelY;}
    QVector<double> getLambda(){return this->lambda;}
    void setLambda(QVector<double>vec){this->lambda=vec;}
    double getLambdaAt(int index){return this->lambda.at(index);}
    void setRecordLengthSize(qint32 size){this->recordLengthSize=size;}
    QVector< QVector< QVector<double>  > > getImages(){return this->images;}
    QVector< QVector<double>  > getImageAt(int index){return this->images.at(index);}
    QVector <double > flux;
    qint32 getRecordLengthSize(){return this->recordLengthSize;}
    double getMax(int);
    double getMin(int);
    QVector<double> getX(){return this->x;}
    QVector<double>  createX(int size);
    QVector<double>  createY(int size);
    static  Image *image;
    QVector<double> getY(){return this->y;}
    bool isError(){return this->error;}
    QPixmap getPixmap(){return this->pixmap.copy();}
    void setPixmap( QPixmap pixmap){this->pixmap=pixmap;}
    int getPixmapWidth(){return pixmap.width();}
    int getPixmapHeight(){return pixmap.height();}
private:
    QPixmap pixmap;
    int numberOfPixelX;
    PhysicalConstants *physicalConstants;
    Logger *logger;
    QDataStream::ByteOrder byteOrder;
    qint32 recordLengthSize;
    QVector<double> lambda;
    QVector< QVector< QVector<double>  > > images;
    int numberOfImages;
    int numberOfPixelY;
    double sizeOfPixelX;
    double sizeOfPixelY;
    QVector<double> x,y;
    bool error;
    bool localObserver;
    QVector<double> dindgen(int size);
signals:
    void finished();
    void setCancelable(bool);
public slots:
    void readImageFormattedFromRadmc3d(QProcess *process,bool *exit);
    
};

#endif // IMAGE_H
