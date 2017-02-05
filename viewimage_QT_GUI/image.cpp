#include "image.h"
Image * Image::image=NULL;
Image::Image(QObject *parent) :
    QObject(parent)
{
    this->localObserver=false;
    this->error=false;
    this->logger=Logger::getInstance();
    this->physicalConstants=PhysicalConstants::getInstance();
    if(QSysInfo::ByteOrder==QSysInfo::BigEndian)
        this->setByteOrder(QDataStream::BigEndian);
    else
        this->setByteOrder(QDataStream::LittleEndian);
    this->setRecordLengthSize(4);
}
double Image::getMin(int index){
    double minValue=0;
    QVector< QVector<double>  > intensity=this->images[index];
    for (int var = 0; var < numberOfPixelY; ++var){
        double* min=min_element(intensity[var].begin(),intensity[var].end());
        if(var==0||(*min)<minValue)minValue=*min;
    }
    return minValue;

}
double Image::getMax(int index){
    double maxValue=0;
    QVector< QVector<double> > intensity=this->images[index];
    for (int var = 0; var < numberOfPixelY; ++var){
        double* max=max_element(intensity[var].begin(),intensity[var].end());
        if(var==0||(*max)>maxValue)maxValue=*max;
    }
    return maxValue;
}
void Image::readImageFormattedFromRadmc3d(QProcess *radmc3dProcess,bool *exit){
    QTextStream in(radmc3dProcess);
    QString str;
    bool calculate=false;
    int nx=0,ny=0;
    int nlam=0;
    this->error=true;
    int total=-999;
    QStringList list;
    bool done=false;
    while(!done){
        QCoreApplication::processEvents();
        radmc3dProcess->waitForReadyRead(200);
        str+=QString::fromLatin1(radmc3dProcess->readAll());
        str.remove("\n");
        list.clear();
        list=str.split(QRegExp("\\s+"),QString::SkipEmptyParts);
        if((list.size()>5) & (!calculate)){
            calculate=true;
            nx=list.at(1).toInt();
            ny=list.at(2).toInt();
            nlam=list.at(3).toInt();
            total=(nx*ny*nlam)+(nlam+6);
        }
        if(calculate){
            if(list.size()>=total){
                done=true;
                error=false;
                emit setCancelable(false);
            }
        }
        if(radmc3dProcess->state()==QProcess::NotRunning||(*exit)){
            if(*exit)this->logger->writeToLogFile(trUtf8("Killing RADMC3D is caused from RADMC3D-GUI"));
            if(calculate)
                this->logger->writeToLogFile(trUtf8("%1 of %2 values are received, bevor radmc3d is crashed. %3\% was completed.")
                                             .arg(QString::number(list.size()),QString::number(total),QString::number((list.size()*100./total))));
            else
                this->logger->writeToLogFile(trUtf8("%1 values are received, bevor radmc3d is crashed").arg(QString::number(list.size())));
            error=true;
            break;
        }
    }

    if(!error){
        list=str.split(QRegExp("\\s+"),QString::SkipEmptyParts);
        int iformat=list.at(0).toInt();
        this->setLocalObserver(iformat==2);



        this->setNumberOfPixelX(nx);
        this->setNumberOfPixelY(ny);



        this->setNumberOfImages(nlam);

        double sizex, sizey;
        sizex=list.at(4).toDouble();
        sizey=list.at(5).toDouble();
        this->setSizeOfPixelX(sizex);
        this->setSizeOfPixelY(sizey);
        QVector<double> lambda;
        int pos=nlam+6;
        for (int var = 6; var < pos; ++var) {
            lambda.append(list.at(var).toDouble());
        }


        this->setLambda(lambda);
        QVector< QVector<double> > imageVec;
        QVector<double> xvec;
        int positon;
        for (int index = 0; index < nlam; ++index) {
            imageVec.clear();
            imageVec.resize(ny);
            for (int iy = 0; iy < ny; ++iy){
                xvec.clear();
                for (int ix =0 ; ix < nx; ++ix){
                    positon=pos+(index*ny*nx)+(iy*nx)+ix;
                    double value=list.at(positon).toDouble();
                    xvec.append(value);
                }
                imageVec[ny-iy-1]=xvec;
            }

            this->images.append(imageVec);
        }
        flux.resize(nlam);
        flux.fill(0);
        for(int n=0;n<nlam;++n){
            for (int iy = 0; iy < ny; ++iy){
                for (int ix = 0; ix < nx; ++ix){
                    flux[n]=flux[n]+images[n][iy][ix];
                }
            }
            flux[n]=flux[n]*sizex*sizey;
            if(!localObserver) flux[n]=flux[n]/(pow(physicalConstants->getParsec(),2));
        }

        this->x =this->createX(nx);
        this->y=this->createY(ny);


        this->logImage();
        error=false;
    }
    emit finished();
}

/** @brief This method read imagefile, if the file is formatted. This method calls findImageFile to look after image file.
 *
 *  @param cs CaseSensitivit, default value is CaseSensitive
 *  @param imageFilename (for example image)
 *  @param imageFileEnding (for example .out)
 *
 *  @returns true, if there is just one image file,otherwise returns false.
 *
 *
 *  @see  bool findImageFile(QString* currentImageFileName,QString imageFilename,QString imageFileEnding,Qt::CaseSensitivity cs);
**/
bool Image::readImageFormatted(QString imageFilename,QString imageFileEnding,Qt::CaseSensitivity cs){
    QString filename;
    QTextStream in;
    QFile file;
    QStringList list;
    if(imageFileEnding.trimmed()!=""){
        if(!this->findImageFile(&filename,imageFilename,imageFileEnding,cs,&list)){
            if(list.size()>1){
                this->logger->writeToLogFile(trUtf8("There is more then 1 image file in this directory(%1)").arg(list.join(" ")));
            }
            return false;
        }
    }else{
        filename=imageFilename;
    }
    file.setFileName(filename);
    in.setDevice(&file);
    if(!file.open(QIODevice::ReadOnly|QIODevice::Text)){
        this->logger->writeToLogFile(trUtf8("Erorr: can not open %1").arg(file.fileName()));
        return false;
    }
    qint32 iformat;
    in>>iformat;
    this->setLocalObserver(iformat==2);
    qint32 nx,ny;
    in>>nx;
    in>>ny;

    this->setNumberOfPixelX(nx);
    this->setNumberOfPixelY(ny);
    qint32 nlam;

    in>>nlam;

    this->setNumberOfImages(nlam);
    double sizex,sizey;
    in>>sizex;
    in>>sizey;
    this->setSizeOfPixelX(sizex);
    this->setSizeOfPixelY(sizey);
    QVector<double>lambdavec;
    lambdavec.clear();
    double value;
    for (int index = 0; index < nlam; ++index) {

        in>>value;
        lambdavec.append(value);
    }
    this->setLambda(lambdavec);
    this->images.clear();
    QVector< QVector<double> > imageVec;
    QVector<double> xvec;
    for (int index = 0; index < nlam; ++index) {
        imageVec.clear();
        imageVec.resize(ny);
        for (int iy = ny-1; 0 <= iy; --iy){
            xvec.clear();
            for (int ix =0 ; ix < nx; ++ix){
                in>>value;
                xvec.append(value);
            }
            imageVec[iy]=xvec;
        }
        this->images.append(imageVec);
    }

    if(file.isOpen())file.close();
    flux.resize(this->numberOfImages);
    flux.fill(0);
    for(int n=0;n<numberOfImages;++n){
        for (int iy = 0; iy < ny; ++iy){
            for (int ix = 0; ix < nx; ++ix){
                flux[n]=flux[n]+images[n][iy][ix];
            }
        }
        flux[n]=flux[n]*sizex*sizey;
        if(!localObserver) flux[n]=flux[n]/(pow(physicalConstants->getParsec(),2));
    }

    this->x =this->createX(nx);
    this->y=this->createY(ny);


    this->logImage();
    return true;
}

/** @brief This method read imagefile, if the file is formatted. This method calls findImageFile to look after image file.
 *
 *  @param cs CaseSensitivit, default value is CaseSensitive
 *  @param imageFilename (for example image)
 *  @param imageFileEnding (for example .out)
 *
 *  @returns true, if there is just one image file,otherwise returns false.
 *
 *
 *  @see  bool findImageFile(QString* currentImageFileName,QString imageFilename,QString imageFileEnding,Qt::CaseSensitivity cs);
**/
bool Image::readImageUnformatted( QString imageFilename,QString imageFileEnding,Qt::CaseSensitivity cs){
    UnformattedSetting *unformattedSetting;
    unformattedSetting=UnformattedSetting::getInstance();
    QString filename;
    QStringList list;
    if(imageFileEnding.trimmed()!=""){
        if(!this->findImageFile(&filename,imageFilename,imageFileEnding,cs,&list)){
            return false;
        }
    }else{
        filename=imageFilename;
    }
    QFile file(filename);
    QDataStream in(&file);
    if(!file.open(QIODevice::ReadOnly)){
        this->logger->writeToLogFile(trUtf8("Erorr: can not open %1").arg(file.fileName()));
        return false;
    }

    qint32 recordlengthsize=unformattedSetting->getRecordLengthSize();
    in.setByteOrder(unformattedSetting->getByteOrder());

    in.skipRawData(recordlengthsize);
    qint32 iformat;
    in>>iformat;
    this->setLocalObserver(iformat==2);
    in.skipRawData(2*recordlengthsize);

    qint32 nx,ny;
    in>>nx;
    in>>ny;
    this->setNumberOfPixelX(nx);
    this->setNumberOfPixelY(ny);

    in.skipRawData(2*recordlengthsize);
    qint32 nlam;
    in>>nlam;

    this->setNumberOfImages(nlam);
    double sizex,sizey;
    in.skipRawData(2*recordlengthsize);
    in>>sizex;
    in>>sizey;
    this->setSizeOfPixelX(sizex);
    this->setSizeOfPixelY(sizey);


    QVector<double>lambdavec;
    lambdavec.clear();
    double value;
    for (int index = 0; index < nlam; ++index) {
        in.skipRawData(2*recordlengthsize);
        in>>value;
        lambdavec.append(value);
    }
    this->setLambda(lambdavec);
    this->images.clear();
    QVector< QVector<double> > imageVec;
    QVector<double> xvec;
    for (int index = 0; index < nlam; ++index) {
        imageVec.clear();
        in.skipRawData(2*recordlengthsize);
        imageVec.resize(ny);
        for (int iy = ny-1; 0 <= iy; --iy){
            xvec.clear();
            for (int ix =0 ; ix < nx; ++ix){
                in>>value;
                xvec.append(value);
            }
            imageVec[iy]=xvec;
        }
        this->images.append(imageVec);
    }

    file.close();
    flux.resize(this->numberOfImages);
    flux.fill(0);
    for(int n=0;n<numberOfImages;++n){
        for (int iy = 0; iy < ny; ++iy){
            for (int ix = 0; ix < nx; ++ix){
                flux[n]=flux[n]+images[n][iy][ix];
            }
        }
        flux[n]=flux[n]*sizex*sizey;
        if(!localObserver) flux[n]=flux[n]/(pow(physicalConstants->getParsec(),2));
    }

    this->x =this->createX(nx);
    this->y=this->createY(ny);


    this->logImage();
    return true;
}

QVector<double> Image::createX(int size){
    QVector<double> vec=this->dindgen(size);
    QVector<double> result;
    for (int var = 0; var < size; ++var) {
        result.append(((vec.at(var)+0.5)/(size*1.)-0.5)*this->sizeOfPixelX*size);
    }
    return result;
}
QVector<double> Image::createY(int size){
    QVector<double> result;
    QVector<double> vec=this->dindgen(size);
    for (int var = 0; var < size; ++var) {
        result.append(((vec.at(var)+0.5)/(size*1.)-0.5)*this->sizeOfPixelY*size);
    }
    return result;
}


void  Image::logImage(){
    QString unit;
    if(this->isLocalObserver()){
        this->logger->writeToLogFile(trUtf8("The Observer is local"));
        unit="radian";
    }else{
        this->logger->writeToLogFile(trUtf8("The Observer is at infinity"));
        unit="cm";
    }
    this->logger->writeToLogFile(trUtf8("Number of Pixel in x direction %1").arg(QString::number(this->getNumberOfPixelX())));

    this->logger->writeToLogFile(trUtf8("Number of Pixel in y direction %1").arg(QString::number(this->getNumberOfPixelY())));

    this->logger->writeToLogFile(trUtf8("Number of images = %1.").arg(QString::number(this->getNumberOfImages())));

    this->logger->writeToLogFile(trUtf8("Size of Pixel in direction x = %1 %2").arg(QString::number(this->getSizeOfPixelX()),unit));
    this->logger->writeToLogFile(trUtf8("Size of Pixel in direction y = %1 %2").arg(QString::number(this->getSizeOfPixelY()),unit));
    QString lambdaString;
    QString fluxString;
    for (int index = 0; index < this->numberOfImages; ++index) {
        lambdaString.append(QString::number(getLambdaAt(index),'e')+" ");
        if(!this->localObserver){
            fluxString.append(QString::number((this->flux.at(index)*pow(10.,23)),'e')+" ");
        }else{
            fluxString.append(QString::number(this->flux.at(index),'e')+" ");
        }
    }
    if(!this->localObserver) fluxString.append("Jy @ 1pc");
    this->logger->writeToLogFile(trUtf8("lambda(s) = %1").arg(lambdaString));
    this->logger->writeToLogFile(trUtf8("Flux(es) = %1").arg(fluxString));

}

QVector<double> Image::dindgen(int size){
    QVector<double> vector;
    for (int var = 0; var < size; ++var) {
        vector.append(var);
    }
    return vector;
}

/**
 *  @brief this method finds the image file in current directory and give that back.
 *  @param cs CaseSensitivit, default value is CaseSensitive
 *  @param imageFilename (for example image)
 *  @param imageFileEnding (for example .out)
 *
 *  @returns true, if there is a image file such imageFilenameimageFileEnding(for example image.out) or such
 *   imageFilename_***imageFileEnding (for example image_***.out)  in the current directory. If
 *   there is no image file or there are more than 1 file, returns false.
 *  @returns currentImageFileName, if there is just one image file, Otherwise returns a empty name.
 *
 */
bool Image::findImageFile(QString *currentImageFileName,QString imageFilename,QString imageFileEnding,Qt::CaseSensitivity cs, QStringList *list){
    QDirIterator iterator(QDir::currentPath());
    bool found=false;
    list->clear();
    QString fileName;
    *currentImageFileName="";
    while(!found && iterator.hasNext()){
        iterator.next();
        fileName=iterator.fileName();
        //  if(((fileName.startsWith(imageFilename+"_",cs)) && (fileName.endsWith(imageFileEnding,cs)))
        if(((fileName.startsWith(imageFilename,cs)) && (fileName.endsWith(imageFileEnding,cs))) && fileName.trimmed()!=""){
            list->append(iterator.fileName());
        }
    }
    if(list->size()==0){
        return false;
    }else if(list->size()>1){
        return false;
    }
    *currentImageFileName=list->at(0);
    return true;
}
