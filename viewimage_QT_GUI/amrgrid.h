#ifndef AMRGRID_H
#define AMRGRID_H
#if QT_VERSION < 0x050000
    #include <QtGui>
#else
    #include <QtWidgets>
#endif

#include "logger.h"
#include "octtree.h"
#include "polar.h"
#include "physicalConstants.h"
#include "cartesian.h"
#include <algorithm>
#include "unformattedsetting.h"
using namespace  std;
class AmrGrid : public QObject
{
    Q_OBJECT
public:

    explicit AmrGrid(QObject *parent = 0);
    void setGridBasicInfo(bool mirror=false);
    void setAmrGridFilenames(QStringList filenames);
    static AmrGrid *getInstanece();
    QStringList getAmrGridFilenames ();
    void setIparentAt(int index,int number);
    void setByteOrder(QDataStream::ByteOrder order){this->byteOrder=order;}
    QDataStream::ByteOrder getByteOrder(){return this->byteOrder;}
    void setRecordLengthSize(qint32 size){this->recordLengthSize=size;}
    qint32 getRecordLengthSize(){return this->recordLengthSize;}
    void setAmrGrid(bool basic=true,bool mirror=false);
    void setIFormat(qint32 number){this->iformat=number;}
    qint32 getIFormat(){return this->iformat;}
    void setAmrStyle(qint32 number){this->amrstyle=number;}
    qint32 getAmrStyle(){return this->amrstyle;}
    void setCoordSystem(qint32 number){this->coordsystem=number;}
    void setGridInfo(qint32 number){this->gridinfo=number;}
    qint32 getGridInfo(){return this->gridinfo;}
    qint32 iformat,amrstyle,coordsystem,gridinfo;
    qint32 getCoordSystem(){return this->coordsystem;}
    bool getDimensionActivityX(){return this->dimensionActivityX;}
    bool getDimensionActivityY(){return this->dimensionActivityY;}
    bool getDimensionActivityZ(){return this->dimensionActivityZ;}
    void setDimensionSizeX(qint32 size){this->dimensionSizeX=size;}
    qint32 getDimensionSizeX(){return this->dimensionSizeX;}
    void setDimensionSizeY(qint32 size){this->dimensionSizeY=size;}
    qint32 getDimensionSizeY(){return this->dimensionSizeY;}
    void setDimensionSizeZ(qint32 size){this->dimensionSizeZ=size;}
    qint32 getDimensionSizeZ(){return this->dimensionSizeZ;}
    void setDimensionSizeXMax(qint32 size){this->dimensionSizeXMax=size;}
    qint32 getDimensionSizeXMax(){return this->dimensionSizeXMax;}
    void setDimensionSizeYMax(qint32 size){this->dimensionSizeYMax=size;}
    qint32 getDimensionSizeYMax(){return this->dimensionSizeYMax;}
    void setDimensionSizeZMax(qint32 size){this->dimensionSizeZMax=size;}
    qint32 getDimensionSizeZMax(){return this->dimensionSizeZMax;}
    void setLevelMax(qint32 number){this->levelMax=number;}
    void setNumberOfLeafsMax(qint32 number){this->numberOfLeafsMax=number;}
    void setNumberOfBranchMax(qint32 number){this->numberOfBranchMax=number;}
    qint32 getLevelMax(){return this->levelMax;}
    qint32 getNumberOfLeafsMax(){return this->numberOfLeafsMax;}
    qint32 getNumberOfBranchMax(){return this->numberOfBranchMax;}

    QVector< QVector<double> > getLayerXi(){return layerXi;}
    QVector< QVector<double> > getLayerYi(){return layerYi;}
    QVector< QVector<double> > getLayerZi(){return layerZi;}
    QVector< QVector<double> > getLayerX(){return layerX;}
    QVector< QVector<double> > getLayerY(){return layerY;}
    QVector< QVector<double> > getLayerZ(){return layerZ;}
    qint32 getNumberOfLayer(){return this->numberOfLayer;}
    void setNumberOfLayer(qint32 number){this->numberOfLayer=number;}
    void logAmrGrid(QString amrGridFilename);
    void setDimensionActivities(bool activityX,bool activityY,bool activityZ){this->dimensionActivityX=activityX;
                                                                             this->dimensionActivityY=activityY;
                                                                            this->dimensionActivityZ=activityZ;}
    void setMirror();
    void setNCell(qint32 number);
    qint32 getNCell();
    void setNCellInp(qint32 number);
    qint32 getNCellInp();
    void setRThetaPhi();
    bool isAmrGridRead(){return this->amrGridRead;}
    void setAmrGridRead(bool read){this->amrGridRead=read;}
private:
    PhysicalConstants *physicalConstants;
    QVector< QVector<double> > layerXi,layerYi,layerZi,layerX,layerY,layerZ;
    Cartesian *cartesian;
    static AmrGrid *instance;
    qint32 dimensionSizeZMax;
    qint32 dimensionSizeZ;
    qint32 numberOfLayer;
    QDataStream::ByteOrder byteOrder;
    Logger *logger;
    bool readFormattedAmrGrid(QFile *file,bool basic=false,bool mirror=false);
    bool readUnformattedAmrGrid(QFile *file,bool basic=false,bool mirror=false);
    bool readBinaryAmrGrid(QFile *file,bool basic=false,bool mirror=false);
    bool mirror;
    bool dimensionActivityX;
    bool dimensionActivityY;
    bool dimensionActivityZ;
    bool amrGridRead;
     QStringList amrGridFilenames;
    qint32 recordLengthSize;
     Polar *polar;
    QVector<int> iparent;
    QVector< QVector<int> > ixyz,nxyz,nnxyz;
    qint32 nCell,nCellInp;
    void setNnxyzAt(int index,QVector<int>vec);
    void setIxyzAt(int index,QVector<int>vec);
    void setNxyzAt(int index,QVector<int>vec);
    Octtree *octtree;
    qint32 dimensionSizeX;
    bool readAmrGrid(bool basic=false,bool mirror=false);
    void setRiPhiiThetai();
    qint32 dimensionSizeY;
    qint32 dimensionSizeXMax;
    qint32 dimensionSizeYMax;
    qint32 numberOfBranchMax;
    qint32 levelMax;
    qint32 numberOfLeafsMax;

public slots:
    
};

#endif // AMRGRID_H
