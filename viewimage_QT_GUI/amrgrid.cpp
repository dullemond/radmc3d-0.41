#include "amrgrid.h"
AmrGrid* AmrGrid::instance=0;
AmrGrid* AmrGrid::getInstanece(){
    if(!instance){
        instance=new AmrGrid();
    }
    return instance;
}

AmrGrid::AmrGrid(QObject *parent) :
    QObject(parent)
{
    this->mirror=false;
    this->dimensionActivityX=false;
    this->dimensionActivityY=false;
    this->dimensionActivityZ=false;
    this->amrGridRead=false;
    this->polar=Polar::getInstance();
    this->cartesian=Cartesian::getInstance();
    this->physicalConstants=PhysicalConstants::getInstance();
    this->logger=Logger::getInstance();
    this->octtree=Octtree::getInstance();
    if(QSysInfo::ByteOrder==QSysInfo::BigEndian)
        this->setByteOrder(QDataStream::BigEndian);
    else
        this->setByteOrder(QDataStream::LittleEndian);
    this->setRecordLengthSize(4);
}
/**
 * @brief This method sets AMR-GRID-filenames appropriate to the given filenames.
 * @param filenames as QStringlist.
 */
void AmrGrid::setAmrGridFilenames(QStringList filenames){
    this->amrGridFilenames=filenames;
    this->logger->writeToLogFile(trUtf8("AMR-grid's filenames are %1").arg(this->amrGridFilenames.join(" ")));

}


/**
 * @brief This method returns AMR-GRID's filenames
 * @return amrGridFilenames as QStringList.
 */
QStringList AmrGrid::getAmrGridFilenames(){
    return this->amrGridFilenames;
}
/**
 * @brief This method sets AMR-Grid.
 *
 */
void AmrGrid::setAmrGrid(bool basic,bool mirror){
    if(readAmrGrid(basic, mirror)){
        amrGridRead=true;
        if(!basic){
            int nx=getDimensionSizeX();
            int ny=getDimensionSizeY();
            int nz=getDimensionSizeZ();
            int nlayer=getNumberOfLayer();
            if(getAmrStyle()==10){
                this->layerXi.resize(nlayer+1);
                this->layerYi.resize(nlayer+1);
                this->layerZi.resize(nlayer+1);
                this->layerX.resize(nlayer+1);
                this->layerY.resize(nlayer+1);
                this->layerZ.resize(nlayer+1);
                for (int index = 0; index<=nlayer; ++index) {
                    this->layerXi[index].resize(getDimensionSizeXMax()+1);
                    this->layerXi[index].fill(0);
                    this->layerYi[index].resize(getDimensionSizeYMax()+1);
                    this->layerYi[index].fill(0);
                    this->layerZi[index].resize(getDimensionSizeZMax()+1);
                    this->layerZi[index].fill(0);
                    this->layerX[index].resize(getDimensionSizeXMax());
                    this->layerX[index].fill(0);
                    this->layerY[index].resize(getDimensionSizeYMax());
                    this->layerY[index].fill(0);
                    this->layerZ[index].resize(getDimensionSizeZMax());
                    this->layerZ[index].fill(0);
                }
                QVector<double> xi,yi,zi;
                xi=this->cartesian->getXi();
                yi=this->cartesian->getYi();
                zi=this->cartesian->getZi();
                this->layerXi.replace(0,xi);
                this->layerYi.replace(0,yi);
                this->layerZi.replace(0,zi);
                QVector<double> layer;
                double value;
                for (int index = 0; index < nx; ++index) {
                    value=0.5*(xi.at(index)+xi.at(index+1));
                    layer.append(value);
                }
                this->layerX.replace(0,layer);
                layer.clear();
                for (int index = 0; index < ny; ++index) {
                    value=0.5*(yi.at(index)+yi.at(index+1));
                    layer.append(value);
                }
                this->layerY.replace(0,layer);
                layer.clear();
                for (int index = 0; index < nz; ++index) {
                    value=0.5*(zi.at(index)+zi.at(index+1));
                    layer.append(value);
                }
                this->layerZ.replace(0,layer);
                layer.clear();
                int nlayer=getNumberOfLayer();
                for (int index = 1; index <= nlayer; ++index) {
                    int total=this->nnxyz.at(index).at(0);
                    for (int k = 0; k <=total; k=k+2) {
                        this->layerXi[index].replace(k,this->layerXi.at(iparent.at(index)).at(ixyz.at(index).at(0)-1+k/2));
                    }
                    total=this->nnxyz.at(index).at(1);
                    for (int k = 0; k <=total; k=k+2) {
                        this->layerYi[index].replace(k,this->layerYi.at(iparent.at(index)).at(ixyz.at(index).at(1)-1+k/2));
                    }
                    total=this->nnxyz.at(index).at(2);
                    for (int k = 0; k <=total; k=k+2) {
                        this->layerZi[index].replace(k,this->layerZi.at(iparent.at(index)).at(ixyz.at(index).at(2)-1+k/2));
                    }
                    total=this->nnxyz.at(index).at(0)-1;
                    double value;
                    if(getCoordSystem()<100){
                        for (int k = 1; k <=total; k=k+2) {
                            value=0.5*(this->layerXi.at(iparent.at(index)).at(ixyz.at(index).at(0)-1+(k-1)/2)+
                                       this->layerXi.at(iparent.at(index)).at(ixyz.at(index).at(0)-1+(k+1)/2));
                            this->layerXi[index].replace(k,value);
                        }
                    }else{
                        for (int k = 1; k <=total; k=k+2) {
                            value=sqrt(this->layerXi.at(iparent.at(index)).at(ixyz.at(index).at(0)-1+(k-1)/2)+
                                       this->layerXi.at(iparent.at(index)).at(ixyz.at(index).at(0)-1+(k+1)/2));
                            this->layerXi[index].replace(k,value);
                        }
                    }
                    total=this->nnxyz.at(index).at(1)-1;
                    for (int k = 1; k <=total; k=k+2) {
                        value=0.5*(this->layerYi.at(iparent.at(index)).at(ixyz.at(index).at(1)-1+(k-1)/2)+
                                   this->layerYi.at(iparent.at(index)).at(ixyz.at(index).at(1)-1+(k+1)/2));

                        this->layerYi[index].replace(k,value);
                    }
                    total=this->nnxyz.at(index).at(2)-1;
                    for (int k = 1; k <=total; k=k+2) {
                        value=0.5*(this->layerZi.at(iparent.at(index)).at(ixyz.at(index).at(2)-1+(k-1)/2)+
                                   this->layerZi.at(iparent.at(index)).at(ixyz.at(index).at(2)-1+(k+1)/2));
                        this->layerZi[index].replace(k,value);
                    }
                    if(getCoordSystem()<100){
                        total=this->nnxyz.at(index).at(0)-1;
                        for (int k = 0; k <=total; ++k) {
                            this->layerX[index].replace(k,0.5*(this->layerXi.at(index).at(k)+this->layerXi.at(index).at(k+1)));
                        }
                    }else{
                        total=this->nnxyz.at(index).at(0)-1;
                        for (int k = 0; k <=total; ++k) {
                            this->layerX[index].replace(k,sqrt(this->layerXi.at(index).at(k)+this->layerXi.at(index).at(k+1)));
                        }
                    }
                    total=this->nnxyz.at(index).at(1)-1;
                    for (int k = 0; k <=total; ++k) {
                        this->layerY[index].replace(k,0.5*(this->layerYi.at(index).at(k)+this->layerYi.at(index).at(k+1)));
                    }
                    total=this->nnxyz.at(index).at(2)-1;
                    for (int k = 0; k <=total; ++k) {
                        this->layerZ[index].replace(k,0.5*(this->layerZi.at(index).at(k)+this->layerZi.at(index).at(k+1)));
                    }
                }
            }
            setGridBasicInfo(mirror);
        }
    }else
        amrGridRead=false;

}

/**
  * @brief This method reads one of the AMR-Grid file, that exists in current directory.
  *
  * @details Just one of the 3 possible file should be existed. The possible file are formatted
  * /unformatted/binary file.
  *
  */
bool AmrGrid::readAmrGrid(bool basic,bool mirror){
    QDir dir(QDir::currentPath());
    QStringList AmrGridFileList=dir.entryList(this->amrGridFilenames);
    qint32 numberOfAmrGridFiles=0;
    numberOfAmrGridFiles=AmrGridFileList.size();
    if(numberOfAmrGridFiles==0){
        this->logger->writeToLogFile(trUtf8("Error: None of these files exists: %1").arg(this->amrGridFilenames.join(" ")));
        return false;
    }else if(numberOfAmrGridFiles>1){
        this->logger->writeToLogFile(trUtf8("Error:Just one of these files can be existed: %1").arg(AmrGridFileList.join(" ")));
        return false;
    }
    QFile file(AmrGridFileList.at(0));
    if(file.fileName().contains(".b")){
        readBinaryAmrGrid(&file,basic,mirror);
    }else if(file.fileName().contains(".u")){
        readUnformattedAmrGrid(&file,basic,mirror);
    }else{
        readFormattedAmrGrid(&file,basic,mirror);
    }
    logAmrGrid(file.fileName());
    return true;
}
/**
 * @brief This method logs the AMR-grid's info.
 *
 */
void AmrGrid::logAmrGrid(QString AmrGridFilename){
    this->logger->writeToLogFile(trUtf8("iformat = %1").arg(QString::number(this->getIFormat())));
    QString amrstyleString="unkown";

    switch(this->getAmrStyle()){
    case 0:
        amrstyleString="regular grid (No mesh refinement)";
        break;
    case 1:
        amrstyleString="octree-style AMR(Adaptive Mesh Refinement)";
        break;
    case 10:
        amrstyleString="layer-style AMR(Adaptive Mesh Refinement)";
    }

    QString coordsystemString="unkown";

    if(this->getCoordSystem()<100){
        coordsystemString="cartesian";
    }else if(this->getCoordSystem()<200){
        coordsystemString="spherical (polar)";
    }else if(this->getCoordSystem()<300){
        coordsystemString="cylindrical";
    }
    this->logger->writeToLogFile(trUtf8("AMR-Style = %1").arg(amrstyleString));
    this->logger->writeToLogFile(trUtf8("coordinate system = %1").arg(coordsystemString));
    if(this->getGridInfo()==0){
        this->logger->writeToLogFile(trUtf8("there will not be abundant grid information written into %1").arg(AmrGridFilename));
    }else{
        this->logger->writeToLogFile(trUtf8("there will be abundant grid information written into %1").arg(AmrGridFilename));
    }
    bool x=this->dimensionActivityX!=0,y=this->dimensionActivityY!=0,z=this->dimensionActivityZ!=0;
    QString dimensionActiveStrig="";
    if(x)dimensionActiveStrig+="X ";
    if(y)dimensionActiveStrig+="Y ";
    if(z)dimensionActiveStrig+="Z ";
    this->logger->writeToLogFile(trUtf8("Active dimension: %1").arg(dimensionActiveStrig));
    QString dimensionSizeString=trUtf8("Dimension's size: ");
    if(x)dimensionSizeString+=trUtf8(" X= %1").arg(QString::number(this->getDimensionSizeX()));
    if(y)dimensionSizeString+=trUtf8(" Y= %1").arg(QString::number(this->getDimensionSizeY()));
    if(z)dimensionSizeString+=trUtf8(" Z= %1").arg(QString::number(this->getDimensionSizeZ()));
    this->logger->writeToLogFile(dimensionSizeString);
    if(this->getAmrStyle()==1){
        this->logger->writeToLogFile(trUtf8("Max. level=%1").arg(QString::number(this->getLevelMax())));
        this->logger->writeToLogFile(trUtf8("Max. number of leafs=%1").arg(QString::number(this->getNumberOfLeafsMax())));
        this->logger->writeToLogFile(trUtf8("Max. number of branch=%1").arg(QString::number(this->getNumberOfBranchMax())));
    }else if(this->amrstyle==10){
        this->logger->writeToLogFile(trUtf8("Top refinement level of the grid = %1").arg(QString::number(this->getLevelMax())));
        this->logger->writeToLogFile(trUtf8("Number of layer=%1").arg(QString::number(this->getNumberOfLayer())));
    }
}

/**
  * @brief This method reads the binary AMR-Grid file, that exists in current directory.
  *
  * @details This method can just read the amr_grid.bnip, if the following conditions are given:
  * The saved number are either  integer(8)Bytes or double(8)Bytes.
  * You can additonally change the byteOrder with setByteOrder(default value is your system endian).
  * The byteOrder is very important for right reading the amr file.
  *
  * @attention qint64 is 64bit
  * @param basic read just basic,if it is set;
  * @param mirror
  *
  * @see setByteOrder(QDataStream::ByteOrder)
  * @Kees warum ist qint32 8 bytes?
  */
bool AmrGrid::readBinaryAmrGrid(QFile *file, bool basic, bool mirror){
    if(!file->exists())return false;
    UnformattedSetting *unformattedSetting;
    unformattedSetting=UnformattedSetting::getInstance();
    QDataStream in(file);
    if(!file->open(QIODevice::ReadOnly)){

        this->logger->writeToLogFile(trUtf8("Error: Can not open %1").arg(file->fileName()));
        this->logger->stop(EXIT_FAILURE);
    }
    in.setByteOrder(unformattedSetting->getByteOrder());



    qint64 format;
    in>>format;
    this->setIFormat((qint32)format);

    qint64 style;
    in>>style;
    this->setAmrStyle((qint32)style);

    qint64 coord;
    in>>coord;
    this->setCoordSystem((qint32)coord);

    qint64 info;
    in>>info;
    this->setGridInfo((qint32)info);

    qint64 dimensionActivity[3];


    for (int index = 0; index < 3; ++index) {
        in>>dimensionActivity[index];
    }
    this->setDimensionActivities(dimensionActivity[0]==1,dimensionActivity[1]==1,dimensionActivity[2]==1);


    qint64 sizes[3];

    for (int index = 0; index < 3; ++index) {
        in>>sizes[index];
    }
    this->setDimensionSizeX((qint32)sizes[0]);
    this->setDimensionSizeY((qint32)sizes[1]);
    this->setDimensionSizeZ((qint32)sizes[2]);

    if((style>=1) && (style<10)){
        qint64 lMax;
        in>>lMax;
        setLevelMax((qint32)lMax);

        qint64 nLM;
        in>>nLM;
        setNumberOfLeafsMax((qint32)nLM);

        qint64 nBM;
        in>>nBM;
        setNumberOfBranchMax((qint32)nBM);

    }else if((style>=10) && (style<20)){
        qint64 lMax;
        in>>lMax;
        setLevelMax((qint32)lMax);

        qint64 nL;
        in>>nL;
        setNumberOfLayer((qint32)nL);

    }else if(style!=0){
        this->logger->writeToLogFile(trUtf8("Unknown AMR-Style: %1").arg(QString::number(this->getAmrStyle())));

    }

    this->setDimensionSizeXMax(sizes[0]);
    this->setDimensionSizeYMax(sizes[1]);
    this->setDimensionSizeZMax(sizes[2]);
    QVector<double> xi,yi,zi;
    xi.clear();
    yi.clear();
    zi.clear();
    qint32 nxWall=sizes[0]+1;
    qint32 nyWall=sizes[1]+1;
    qint32 nzWall=sizes[2]+1;
    xi.resize(nxWall);
    yi.resize(nyWall);
    zi.resize(nzWall);

    double value;

    for (int index = 0; index < nxWall; ++index) {
        in>>value;
        xi.append(value);

    }

    for (int index = 0; index < nyWall; ++index) {
        in>>value;
        yi.append(value);

    }

    for (int index = 0; index < nzWall; ++index) {
        in>>value;
        zi.append(value);

    }
    this->cartesian->setXi(xi);
    this->cartesian->setYi(yi);
    this->cartesian->setZi(zi);
    if(basic){
        setGridBasicInfo(mirror);
        file->close();
        return true;
    }

    if(style==0){
        setNCell(getDimensionSizeX()*getDimensionSizeY()*getDimensionSizeZ());
        setNCellInp(getNCell());
    }else if( (style>=1) && (style<10)){
        qint8 value;
        QVector<qint32>values;
        qint32 ncell=0;
        for (qint32 index = 0; index < this->numberOfBranchMax; ++index) {
            in>>value;
            if(value==0)ncell++;
            values.append(value);
        }
        this->logger->writeToLogFile(trUtf8("Octtree Number of leaf equal zero= %1").arg(QString::number(ncell)));
        this->octtree->setOcttree(values);
        this->setNCell(ncell);
        this->setNCellInp(ncell);

    } else if( (style>=10) && (style<20)){

        qint32 nx,ny,nz;
        nx=getDimensionSizeX();
        ny=getDimensionSizeY();
        nz=getDimensionSizeZ();
        QVector<qint32> sizes;
        sizes.append(nx);
        sizes.append(ny);
        sizes.append(nz);
        qint32 ncell=nx*ny*nz;
        qint32 ncellinp=nx*ny*nz;
        qint32 nlayer=getNumberOfLayer();

        QVector<int> nxyzVector,nnxyzVector,vector;
        QVector<int> idat;
        int dat;
        this->iparent.resize(nlayer+1);
        this->nnxyz.resize(nlayer+1);
        this->nxyz.resize(nlayer+1);
        this->ixyz.resize(nlayer+1);
        for (int index = 0; index <= nlayer; ++index) {
            this->nnxyz[index].resize(3);
            this->nnxyz[index].fill(0);
            this->nxyz[index].resize(3);
            this->nxyz[index].fill(0);
            this->ixyz[index].resize(3);
            this->ixyz[index].fill(0);
        }
        this->setNnxyzAt(0,sizes);
        this->iparent.resize(nlayer+1);
        for (int index = 1; index <= nlayer; ++index) {
            idat.clear();
            for(int i=0;i<7;i++){
                in>>dat;
                idat.append((int)dat);
            }
            this->setIparentAt(index,idat[0]);
            vector.clear();
            vector<<idat[1]<<idat[2]<<idat[3];
            this->setIxyzAt(index,vector);
            nxyzVector.clear();
            nnxyzVector.clear();
            nxyzVector<<idat[4]<<idat[5]<<idat[6];
            nnxyzVector=nxyzVector;
            this->setNxyzAt(index,nxyzVector);
            if(getDimensionActivityX())nnxyzVector[0]=nnxyzVector[0]*2;
            if(getDimensionActivityY())nnxyzVector[1]=nnxyzVector[1]*2;
            if(getDimensionActivityZ())nnxyzVector[2]=nnxyzVector[2]*2;
            this->setNnxyzAt(index,nnxyzVector);
            ncell=ncell+nnxyzVector[0]*nnxyzVector[1]*nnxyzVector[2]-nxyzVector[0]*nxyzVector[1]*nxyzVector[2];
            ncellinp=ncellinp+nnxyzVector[0]*nnxyzVector[1]*nnxyzVector[2];
            setDimensionSizeXMax(max(getDimensionSizeX(),nnxyzVector[0]));

        }
        this->setNCell(ncell);
        this->setNCellInp(ncellinp);

    }


    file->close();

    return true;
}


/**
  * @brief This method reads the unformatted AMR-Grid file, that exists in current directory.
  *
  * @details The right record length's size and byteOrder is very important for reading the amr_grid.unip. Because there
  * is no standart implementation for unfomratted sequentiel file.
  * This method can just read the amr_grid.unip, if the following conditions are given:
  * Fortran write's function writes every time, that it called, the record length  at the start and end of each record.
  * The size of record length is given in bytes, usually 4 bytes. You can change the size of record length by calling
  * setRecordLengthSize-method. You can additonally change the byteOrder with setByteOrder(default value is your system endian).
  * The saved number are either  qint32(4)Bytes or double(8)Bytes.
  *
  * @attention qint32 is 4 byte(32 bit)
  *
  * @param basic read just basic,if it is set;
  * @param mirror
  *
  * @see setByteOrder(QDataStream::ByteOrder)
  * @see setRecordLengthSize(qint32)
  *
  * @return a false, if file does not exist, otherwise true.
  *
  */
bool AmrGrid::readUnformattedAmrGrid(QFile *file,bool basic, bool mirror){
    if(!file->exists())return false;
    UnformattedSetting *unformattedSetting;
    unformattedSetting=UnformattedSetting::getInstance();
    QDataStream in(file);
    if(!file->open(QIODevice::ReadOnly)){
        this->logger->writeToLogFile(trUtf8("Error: Can not open %1").arg(file->fileName()));
        this->logger->stop(EXIT_FAILURE);
    }

    qint32 recordlengthsize=unformattedSetting->getRecordLengthSize();
    in.setByteOrder(unformattedSetting->getByteOrder());



    qint32 format;

    in.skipRawData(recordlengthsize);

    in>>format;
    this->setIFormat(format);

    in.skipRawData(recordlengthsize*2);

    qint32 style;
    in>>style;
    this->setAmrStyle(style);

    in.skipRawData(recordlengthsize*2);

    qint32 coord;
    in>>coord;
    this->setCoordSystem(coord);

    in.skipRawData(recordlengthsize*2);

    qint32 info;
    in>>info;
    this->setGridInfo(info);

    in.skipRawData(recordlengthsize*2);

    qint32 dimensionActivity[3];

    for (qint32 index = 0; index < 3; ++index) {
        in>>dimensionActivity[index];
    }
    this->setDimensionActivities(dimensionActivity[0]==1,dimensionActivity[1]==1,dimensionActivity[2]==1);

    in.skipRawData(recordlengthsize*2);

    qint32 sizes[3];

    for (qint32 index = 0; index < 3; ++index) {
        in>>sizes[index];
    }
    this->setDimensionSizeX(sizes[0]);
    this->setDimensionSizeY(sizes[1]);
    this->setDimensionSizeZ(sizes[2]);

    in.skipRawData(recordlengthsize*2);


    if((style>=1) && (style<10)){
        qint32 lMax;
        in>>lMax;
        setLevelMax(lMax);

        qint32 nLM;
        in>>nLM;
        setNumberOfLeafsMax(nLM);

        qint32 nBM;
        in>>nBM;
        setNumberOfBranchMax(nBM);

        in.skipRawData(recordlengthsize*2);
    }else if((style>=10) && (style<20)){
        qint32 lMax;
        in>>lMax;
        setLevelMax(lMax);

        qint32 nL;
        in>>nL;
        setNumberOfLayer(nL);

        in.skipRawData(recordlengthsize*2);

    }else if(style!=0){
        this->logger->writeToLogFile(trUtf8("Unknown AMR-Style: %1").arg(QString::number(this->getAmrStyle())));

    }
    this->setDimensionSizeXMax(sizes[0]);
    this->setDimensionSizeYMax(sizes[1]);
    this->setDimensionSizeZMax(sizes[2]);
    QVector<double>xi,yi,zi;
    xi.clear();
    yi.clear();
    zi.clear();
    qint32 nxWall=sizes[0]+1;
    qint32 nyWall=sizes[1]+1;
    qint32 nzWall=sizes[2]+1;
    xi.resize(nxWall);
    yi.resize(nyWall);
    zi.resize(nzWall);

    double value;



    for (int index = 0; index < nxWall; ++index) {
        in>>value;
        xi.append(value);
    }

    in.skipRawData(recordlengthsize*2);

    for (int index = 0; index < nyWall; ++index) {
        in>>value;
        yi.append(value);
    }
    in.skipRawData(recordlengthsize*2);

    for (int index = 0; index < nzWall; ++index) {
        in>>value;
        zi.append(value);
    }
    this->cartesian->setXi(xi);
    this->cartesian->setYi(yi);
    this->cartesian->setZi(zi);
    if(basic){
        setGridBasicInfo(mirror);
        file->close();
        return true;
    }



    if(style==0){
        setNCell(getDimensionSizeX()*getDimensionSizeY()*getDimensionSizeZ());
        setNCellInp(getNCell());
    }else if((style>=1) && (style<10)){
        in.skipRawData(recordlengthsize*2);

        qint32 value;
        QVector<qint32>values;
        qint32 ncell=0;
        for (qint32 index = 0; index < this->numberOfBranchMax; ++index) {
            in>>value;
            in.skipRawData(recordlengthsize*2);
            if(value==0)ncell++;
            values.append(value);
        }
        this->logger->writeToLogFile(trUtf8("Octtree Number of leaf equal zero= %1").arg(QString::number(ncell)));
        this->octtree->setOcttree(values);
        this->setNCell(ncell);
        this->setNCellInp(ncell);

    } else if( (style>=10) && (style<20)){
        in.skipRawData(recordlengthsize*2);
        qint32 nx,ny,nz;
        nx=getDimensionSizeX();
        ny=getDimensionSizeY();
        nz=getDimensionSizeZ();
        QVector<qint32> sizes;
        sizes.append(nx);
        sizes.append(ny);
        sizes.append(nz);
        qint32 ncell=nx*ny*nz;
        qint32 ncellinp=nx*ny*nz;
        qint32 nlayer=getNumberOfLayer();

        QVector<int> nxyzVector,nnxyzVector,vector;
        QVector<int> idat;
        int dat;
        this->iparent.resize(nlayer+1);
        this->nnxyz.resize(nlayer+1);
        this->nxyz.resize(nlayer+1);
        this->ixyz.resize(nlayer+1);
        for (int index = 0; index <= nlayer; ++index) {
            this->nnxyz[index].resize(3);
            this->nnxyz[index].fill(0);
            this->nxyz[index].resize(3);
            this->nxyz[index].fill(0);
            this->ixyz[index].resize(3);
            this->ixyz[index].fill(0);
        }
        this->setNnxyzAt(0,sizes);

        this->iparent.resize(nlayer+1);
        for (int index = 1; index <= nlayer; ++index) {
            idat.clear();
            for(int i=0;i<7;i++){
                in>>dat;
                idat.append(dat);
            }
            this->setIparentAt(index,idat[0]);
            vector.clear();
            vector<<idat[1]<<idat[2]<<idat[3];
            this->setIxyzAt(index,vector);
            nxyzVector.clear();
            nnxyzVector.clear();
            nxyzVector<<idat[4]<<idat[5]<<idat[6];
            nnxyzVector=nxyzVector;
            this->setNxyzAt(index,nxyzVector);
            if(getDimensionActivityX())nnxyzVector[0]=nnxyzVector[0]*2;
            if(getDimensionActivityY())nnxyzVector[1]=nnxyzVector[1]*2;
            if(getDimensionActivityZ())nnxyzVector[2]=nnxyzVector[2]*2;
            this->setNnxyzAt(index,nnxyzVector);
            ncell=ncell+nnxyzVector[0]*nnxyzVector[1]*nnxyzVector[2]-nxyzVector[0]*nxyzVector[1]*nxyzVector[2];
            ncellinp=ncellinp+nnxyzVector[0]*nnxyzVector[1]*nnxyzVector[2];
            setDimensionSizeXMax(max(getDimensionSizeX(),nnxyzVector[0]));

        }
        this->setNCell(ncell);
        this->setNCellInp(ncellinp);

    }
    file->close();
    return true;
}
/**
  * @brief This method reads the formatted AMR-Grid file, that exists in current directory
  * @param basic read just basic,if it is set;
  * @param mirror
  * @return a false, if file does not exist, otherwise true.
  */

bool AmrGrid::readFormattedAmrGrid(QFile *file,bool basic,bool mirror){
    if(!file->exists())return false;
    QTextStream in(file);
    if(!file->open(QIODevice::Text|QIODevice::ReadOnly)){
        this->logger->writeToLogFile(trUtf8("Error:Can not open %1").arg(file->fileName()));
        this->logger->stop(EXIT_FAILURE);
    }
    qint32 format;
    in>>format;
    this->setIFormat(format);
    qint32 style;
    in>>style;
    this->setAmrStyle(style);
    qint32 coord;
    in>>coord;
    this->setCoordSystem(coord);
    qint32 info;
    in>>info;
    this->setGridInfo(info);
    qint32 dimensionActivity[3];
    in>>dimensionActivity[0];
    in>>dimensionActivity[1];
    in>>dimensionActivity[2];
    this->setDimensionActivities(dimensionActivity[0]==1,dimensionActivity[1]==1,dimensionActivity[2]==1);
    qint32 sizes[3];
    in>>sizes[0];
    in>>sizes[1];
    in>>sizes[2];
    this->setDimensionSizeX(sizes[0]);
    this->setDimensionSizeY(sizes[1]);
    this->setDimensionSizeZ(sizes[2]);

    if((style>=1) && (style<10)){
        qint32 lMax;
        in>>lMax;
        setLevelMax(lMax);
        qint32 nLM;
        in>>nLM;
        setNumberOfLeafsMax(nLM);
        qint32 nBM;
        in>>nBM;
        setNumberOfBranchMax(nBM);
    }else if((style>=10) && (style<20)){
        qint32 lMax;
        in>>lMax;
        setLevelMax(lMax);
        qint32 nL;
        in>>nL;
        setNumberOfLayer(nL);
    }else if(style!=0){
        this->logger->writeToLogFile(trUtf8("Unknown AMR-Style: %1").arg(QString::number(this->getAmrStyle())));
    }
    this->setDimensionSizeXMax(sizes[0]);
    this->setDimensionSizeYMax(sizes[1]);
    this->setDimensionSizeZMax(sizes[2]);
    QVector<double>xi,yi,zi;
    xi.clear();
    yi.clear();
    zi.clear();
    qint32 nxWall=getDimensionSizeX()+1;
    qint32 nyWall=getDimensionSizeY()+1;
    qint32 nzWall=getDimensionSizeZ()+1;
    double value;
    for (int index = 0; index < nxWall; ++index) {
        in>>value;
        xi.append(value);
    }
    for (int index = 0; index < nyWall; ++index) {
        in>>value;
        yi.append(value);
    }
    for (int index = 0; index < nzWall; ++index) {
        in>>value;
        zi.append(value);
    }
    this->cartesian->setXi(xi);
    this->cartesian->setYi(yi);
    this->cartesian->setZi(zi);
    if(basic){
        setGridBasicInfo(mirror);
        file->close();
        return true;
    }
    if(style==0){
        setNCell(getDimensionSizeX()*getDimensionSizeY()*getDimensionSizeZ());
        setNCellInp(getNCell());

    }  else if((style>=1) && (style<10)){
        qint32 value;
        QVector<qint32>values;
        qint32 ncell=0;
        for (qint32 index = 0; index < this->numberOfBranchMax; ++index) {
            in>>value;
            if(value==0)ncell++;
            values.append(value);
        }
        this->logger->writeToLogFile(trUtf8("Octtree Number of leaf equal zero= %1").arg(QString::number(ncell)));
        this->octtree->setOcttree(values);
        this->setNCell(ncell);
        this->setNCellInp(ncell);
    }else if((style>=10) && (style<20)){
        qint32 nx,ny,nz;
        nx=getDimensionSizeX();
        ny=getDimensionSizeY();
        nz=getDimensionSizeZ();
        QVector<qint32> sizes;
        sizes.append(nx);
        sizes.append(ny);
        sizes.append(nz);

        qint32 ncell=nx*ny*nz;
        qint32 ncellinp=nx*ny*nz;
        qint32 nlayer=getNumberOfLayer();

        QVector<int> nxyzVector,nnxyzVector,vector;
        QVector<int> idat;
        int dat;
        this->iparent.resize(nlayer+1);
        this->nnxyz.resize(nlayer+1);
        this->nxyz.resize(nlayer+1);
        this->ixyz.resize(nlayer+1);
        for (int index = 0; index <= nlayer; ++index) {
            this->nnxyz[index].resize(3);
            this->nnxyz[index].fill(0);
            this->nxyz[index].resize(3);
            this->nxyz[index].fill(0);
            this->ixyz[index].resize(3);
            this->ixyz[index].fill(0);
        }
        this->setNnxyzAt(0,sizes);
        for (int index = 1; index <= nlayer; ++index) {
            idat.clear();
            for(int i=0;i<7;i++){
                in>>dat;
                idat.append(dat);
            }
            this->setIparentAt(index,idat[0]);
            vector.clear();
            vector<<idat[1]<<idat[2]<<idat[3];
            this->setIxyzAt(index,vector);
            nxyzVector.clear();
            nnxyzVector.clear();
            nxyzVector<<idat[4]<<idat[5]<<idat[6];
            nnxyzVector=nxyzVector;
            this->setNxyzAt(index,nxyzVector);
            if(getDimensionActivityX())nnxyzVector[0]=nnxyzVector[0]*2;
            if(getDimensionActivityY())nnxyzVector[1]=nnxyzVector[1]*2;
            if(getDimensionActivityZ())nnxyzVector[2]=nnxyzVector[2]*2;
            this->setNnxyzAt(index,nnxyzVector);
            ncell=ncell+nnxyzVector[0]*nnxyzVector[1]*nnxyzVector[2]-nxyzVector[0]*nxyzVector[1]*nxyzVector[2];
            ncellinp=ncellinp+nnxyzVector[0]*nnxyzVector[1]*nnxyzVector[2];
            setDimensionSizeXMax(max(getDimensionSizeX(),nnxyzVector[0]));

        }
        this->setNCell(ncell);
        this->setNCellInp(ncellinp);
    }

    file->close();
    return true;
}
void AmrGrid::setIparentAt(int index,int number){
    this->iparent.replace(index,number);
}

void AmrGrid::setNnxyzAt(int index,QVector<qint32> vec){
    this->nnxyz.replace(index,vec);
}
void AmrGrid::setNxyzAt(int index,QVector<qint32> vec){
    this->nxyz.replace(index,vec);
}
void AmrGrid::setIxyzAt(int index,QVector<qint32> vec){
    this->ixyz.replace(index,vec);
}


void AmrGrid::setNCell(qint32 number){
    this->nCell=number;
}
void AmrGrid::setNCellInp(qint32 number){
    this->nCellInp=number;
}
qint32 AmrGrid::getNCell(){
    return this->nCell;
}
qint32 AmrGrid::getNCellInp(){
    return this->nCellInp;
}

void AmrGrid::setGridBasicInfo(bool mirror){
    QVector<double> xi,yi,zi;
    if(getCoordSystem()<100){
        if(getDimensionActivityX()){
            xi=this->cartesian->getXi();
            QVector<double> x;
            for (int var = 0; var < this->dimensionSizeX; ++var) {
                x.append(0.5*(xi.at(var)+xi.at(var+1)));
            }
            this->cartesian->setX(x);

        }
        if(getDimensionActivityY()){
            yi=this->cartesian->getYi();
            QVector<double> y;
            for (int var = 0; var < this->dimensionSizeY; ++var) {
                y.append(0.5*(yi.at(var)+yi.at(var+1)));
            }
            this->cartesian->setY(y);

        }
        if(getDimensionActivityZ()){
            QVector<double> z;
            zi=this->cartesian->getZi();
            for (int var = 0; var < this->dimensionSizeZ; ++var) {
                z.append(0.5*(zi.at(var)+zi.at(var+1)));
            }
            this->cartesian->setZ(z);

        }
    }else{
        setRiPhiiThetai();
        if(mirror){
            setMirror();
        }
        setRThetaPhi();

    }
}
void AmrGrid::setRThetaPhi(){
    if(getDimensionActivityX()){
        QVector<double> r,ri;
        ri= this->cartesian->getXi();;
        for (int index = 0; index < this->polar->getNumberOfR(); ++index) {
            r.append(sqrt(ri.at(index)*ri.at(index+1)));
        }
        this->polar->setR(r);
    }

    if(getDimensionActivityY()){
        QVector<double> theta,thetai;
        thetai=this->cartesian->getYi();
        for (int index = 0; index < this->polar->getNumberOfTheta(); ++index) {
            theta.append(0.5*(thetai.at(index)+thetai.at(index+1)));
        }
        this->polar->setTheta(theta);
    }
    if(getDimensionActivityZ()){
        QVector<double> phi,phii;
        phii=this->cartesian->getZi();
        for (int index = 0; index <this->polar->getNumberOfPhi(); ++index) {
            phi.append(0.5*(phii.at(index)+phii.at(index+1)));
        }
        this->polar->setPhi(phi);
    }
}

void AmrGrid::setMirror(){

    QVector<double> thetai=this->polar->getThetai();
    QVector<double> vector;
    double pi=physicalConstants->getPi();
    for (int index = this->dimensionSizeY-1; index > -1; --index) {
        thetai.append(pi-thetai.at(index));
    }
    this->polar->setThetai(thetai);
    this->polar->setNumberOfTheta(this->polar->getNumberOfTheta()*2);
}

void AmrGrid::setRiPhiiThetai(){
    if(getDimensionActivityX()){
        this->polar->setRi(this->cartesian->getXi());
        this->polar->setNumberOfR(this->getDimensionSizeX());
    }
    if(getDimensionActivityY()){
        this->polar->setThetai(this->cartesian->getYi());
        this->polar->setNumberOfTheta(this->getDimensionSizeY());
    }

    if(getDimensionActivityZ()){
        this->polar->setPhii(this->cartesian->getZi());
        this->polar->setNumberOfPhi(this->getDimensionSizeZ());
    }
}
