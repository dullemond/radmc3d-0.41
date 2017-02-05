#include "spectrum.h"

Spectrum::Spectrum(QObject *parent) :
    QObject(parent)
{
    numberOfFrequencies=0;
    this->logger=Logger::getInstance();
    this->physicalConstants=PhysicalConstants::getInstance();
}
/** @brief This method read spectrumfile, if the file is formatted. This method calls findSpectrumFile to look after spectrum file.
 *
 *  @param cs CaseSensitivit, default value is CaseSensitive
 *  @param spectrumFilename (for example spectrum)
 *  @param spectrumFileEnding (for example .out)
 *
 *  @returns true, if there is just one spectrum file,otherwise returns false.
 *
 *
 *  @see  bool findSpectrumFile(QString* currentSpectrumFileName,QString spectrumFilename,QString spectrumFileEnding,Qt::CaseSensitivity cs);
**/
bool Spectrum::readSpectrumFormatted(QString spectrumFilename,QString spectrumFileEnding,Qt::CaseSensitivity cs){
    QString filename;
    QTextStream in;
    QFile file;
    QStringList list;
    if(spectrumFileEnding!=""){
    if(!this->findSpectrumFile(&filename,spectrumFilename,spectrumFileEnding,cs,&list)){
        if(list.size()>1)
            this->logger->writeToLogFile(trUtf8("There is more then 1 spectrum file in this directory(%1)").arg(list.join(" ")));
        return false;
    }
    }else{
        filename=spectrumFilename;
    }
    file.setFileName(filename);
    in.setDevice(&file);
    if(!file.open(QIODevice::ReadOnly|QIODevice::Text)){
        this->logger->writeToLogFile(trUtf8("Erorr: can not open %1").arg(file.fileName()));
        return false;
    }

    qint32 iformat;
    in>>iformat;
    if(iformat!=1){
        this->logger->writeToLogFile(trUtf8("ERROR: File format of %1 not recognized.(File format 1 is just valid)").arg(file.fileName()));
        file.close();
        return false;
    }
    this->format=iformat;

    qint32 nrfr;
    in>>nrfr;


    this->numberOfFrequencies=nrfr;
    double lambda;
    double flx;
    this->flux.clear();
    this->frequencies.clear();
    this->frequencies.resize(nrfr);
    this->flux.resize(nrfr);
    double velocityDivideByMicron=this->physicalConstants->getLightSpeed()/this->physicalConstants->getMicron();
    for (int var = 0; var < nrfr; ++var) {
        in>>lambda;
        in>>flx;
        frequencies[var]=velocityDivideByMicron/lambda;
        flux[var]=flx;
    }
    file.close();
    return true;
}
/**
 *  @brief this method finds the spectrum file in current directory and give that back.
 *  @param cs CaseSensitivit, default value is CaseSensitive
 *  @param spectrumFilename (for example spectrum)
 *  @param spectrumFileEnding (for example .out)
 *
 *  @returns true, if there is a spectrum file such spectrumFilenamespectrumFileEnding(for example spectrum.out) or such
 *   spectrumFilename_***spectrumFileEnding (for example spectrum_***.out)  in the current directory. If
 *   there is no spectrum file or there are more than 1 file, returns false.
 *  @returns currentSpectrumFileName, if there is just one spectrum file, Otherwise returns a empty name.
 *
 */
bool Spectrum::findSpectrumFile(QString *currentSpectrumFileName,QString spectrumFilename,QString spectrumFileEnding,Qt::CaseSensitivity cs, QStringList *list){
    QDirIterator iterator(QDir::currentPath());
    bool found=false;
    list->clear();
    QString fileName;
    *currentSpectrumFileName="";
    while(!found && iterator.hasNext()){
        iterator.next();
        fileName=iterator.fileName();
        //  if(((fileName.startsWith(spectrumFilename+"_",cs)) && (fileName.endsWith(spectrumFileEnding,cs)))
            if(((fileName.startsWith(spectrumFilename,cs)) && (fileName.endsWith(spectrumFileEnding,cs)))){
                list->append(iterator.fileName());
            }
    }
    if(list->size()==0){
        return false;
    }else if(list->size()>1){

        return false;
    }
    *currentSpectrumFileName=list->at(0);
    return true;
}

void Spectrum::readSpectrumFormattedFromRadmc3d(QProcess *radmc3dProcess,bool *exit){
    QTextStream in(radmc3dProcess);
    QString str;
    bool calculate=false;
    int nrfr=0;
    int iformat=0;
    error=true;
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
        if((list.size()>2) & (!calculate)){
            iformat=list.at(0).toInt();
            nrfr=list.at(1).toInt();
            if(iformat!=1){
                this->logger->writeToLogFile(trUtf8("Spectrum format 1 is just valid(current:%1)").arg(list.at(0)));
                break;
                error=true;
            }
            calculate=true;
            total=nrfr*2+2;
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
        this->numberOfFrequencies=nrfr;
        list=str.split(QRegExp("\\s+"),QString::SkipEmptyParts);
        this->flux.clear();
        this->frequencies.clear();
        this->frequencies.resize(nrfr);
        this->flux.resize(nrfr);
        int index=0;
        double lambda;
        double velocityDivideByMicron=this->physicalConstants->getLightSpeed()/this->physicalConstants->getMicron();
        for (int var = 2; var < total-1; var=var+2) {
            lambda=list.at(var).toDouble();
            frequencies[index]=velocityDivideByMicron/lambda;
            flux[index]=list.at(var+1).toDouble();
            index++;
        }

    }
    emit finished();
}
