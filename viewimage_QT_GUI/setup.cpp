#include "setup.h"
#include <algorithm>
#include <cmath>
using namespace  std;
SetUp *SetUp::instance=0;
/**
 *  @brief This method gives a uniq instance of SetUp class.
 *
 *
 */
SetUp  *SetUp::getInstance(QObject *parent){
    if(!instance){
        instance=new SetUp(parent);
    }
    return instance;
}
double SetUp::calculateImageSize(){
    double imageSize;
    if(this->amrGrid->getCoordSystem()<100){
        QVector<double> XiYiZi,xi,yi,zi;
        xi=this->cartesian->getXi();
        yi=this->cartesian->getYi();
        zi=this->cartesian->getZi();
        int xiSize=xi.size();
        int yiSize=yi.size();
        int ziSize=zi.size();

        for (int index = 0; index < xiSize; ++index) {
            XiYiZi.append(abs(xi.at(index)));
        }
        for (int index = 0; index < yiSize; ++index) {
            XiYiZi.append(abs(yi.at(index)));
        }
        for (int index = 0; index < ziSize; ++index) {
            XiYiZi.append(abs(zi.at(index)));
        }
        double *max=max_element(XiYiZi.begin(),XiYiZi.end());
        imageSize=3*(*max);
    }else{
        QVector<double> ri;
        ri=this->polar->getRi();
        double *max=max_element(ri.begin(),ri.end());
        imageSize=2.2*(*max);
    }

    if(this->imageunit==PC){
        imageSize=imageSize/ this->physicalConstants->getParsec();
    }else if(this->imageunit==AU){
        imageSize=imageSize/ this->physicalConstants->getAstronomicalUnit();
    }
    return imageSize;
}

/**
 *  @brief This constructor call sets the basic properties of RADMC3DGUI such filenames.
 *
 * @attention This constructor should not be called directly. instead of that you should call the getInstance, to have a uniq instance
 * of SetUp class
 * @see getInstance
 *
*/

SetUp::SetUp(QObject *parent):
    QObject(parent)
{
    this->lambdaRead=false;
    this->pipe=false;
    this->radmc3dLocal=false;
    this->logger=Logger::getInstance();
    this->physicalConstants=PhysicalConstants::getInstance();
    this->amrGrid=AmrGrid::getInstanece();
    this->cartesian=Cartesian::getInstance();
    this->polar=Polar::getInstance();
    this->logger->setLogFilename("gui_inout.log");
    imageunit=CM;
    this->logger->createLogFile();
    this->logger->redirectToLogFile(true);
    this->logger->writeToLogFile(trUtf8("Current directory is %1").arg(QDir::currentPath()));

    this->setStandardFilenames();




}

void SetUp::setLambda(){
    if(this->readLambda()){
        this->setBlueLambda(1.0);
        this->setGreenLambda(10.0);
        this->setRedLambda(1000.0);
        this->setImageUnit(CM);
        double *max=max_element((this->lambdaArray.begin()),(this->lambdaArray.end()));
        double *min=min_element((this->lambdaArray.begin()),(this->lambdaArray.end()));
        this->blueLambdaIndex=getIndexOfLambda(blueLambda,*min,*max);
        this->greenLambdaIndex=getIndexOfLambda(greenLambda,*min,*max);
        this->redLambdaIndex=getIndexOfLambda(redLambda,*min,*max);
        this->logger->writeToLogFile("First index is 0");
        this->logger->writeToLogFile(trUtf8("Index of Blue Lambda (%1 µM) is %2").
                                     arg(QString::number(this->blueLambda),QString::number(this->blueLambdaIndex)));
        this->logger->writeToLogFile(trUtf8("Index of Green Lambda (%1 µM) is %2").
                                     arg(QString::number(this->greenLambda),QString::number(this->greenLambdaIndex)));
        this->logger->writeToLogFile(trUtf8("Index of Red Lambda (%1 µM) is %2").arg(
                                         QString::number(this->redLambda),QString::number(this->redLambdaIndex)));
        this->lambdaRead=true;
    }else{
        this->lambdaRead=false;
    }
}

void SetUp::setAmrGrid(bool basic,bool mirror){

    this->amrGrid->setAmrGrid(basic,mirror);
}

/**
 * @brief This method returns index of the given lambda.
 * @param lambda.
 *
 */
int SetUp::getIndexOfLambda(double lambda, double min, double max){
    int indexOfLambda=0;
    if(lambda<=max && lambda>=min){
        for (int index = 0; index < this->lambdaSize-2; ++index) {
            if((lambda-this->lambdaArray.at(index))*(lambda-this->lambdaArray.at(index+1))<=0){
                indexOfLambda=index;
            }
        }
    }else{
        if(this->lambdaArray.at(1)>this->lambdaArray.at(0)){
            if(lambda>max)indexOfLambda=this->lambdaSize-1;
        }else{
            if(lambda<=max)indexOfLambda=this->lambdaSize-1;
        }

    }
    return indexOfLambda;
}

/**
 * @brief This method sets blue lambda.
 * @param lambda.
 *
 */
void SetUp::setBlueLambda(double lambda){
    this->blueLambda=lambda;
}
/**
 * @brief This method sets green lambda.
 * @param lambda.
 *
 */
void SetUp::setGreenLambda(double lambda){
    this->greenLambda=lambda;
}
/**
 * @brief This method sets red lambda.
 * @param lambda.
 *
 */
void SetUp::setRedLambda(double lambda){
    this->redLambda=lambda;
}

/**
*
* @brief This method reads  wavelengths to the lambda array,
* If there is a wavelength file on current directory,
* If there is not any wavelengths's file, the method will looking after
* a frequency file.
*
* @attention The number of wavelengths/freuencies should at least one,
*/
bool SetUp::readLambda(){
    QFile *file=new QFile(this->waveLengthFilename);
    int nrOfLambdas=0;
    if(file->exists()){
        this->logger->writeToLogFile(trUtf8("The wavelengths will be read from %1").arg(this->waveLengthFilename));
        nrOfLambdas=this->readWaveLengths(file);
    }else{
        file=new QFile(this->frequencyFilename);
        if(file->exists()){
            this->logger->writeToLogFile(trUtf8("The frequencies will be read from %1").arg(this->waveLengthFilename));
            nrOfLambdas=this->readFrequencies(file);
        }else{
            this->logger->writeToLogFile(trUtf8("Erorr:There is either %1 or %2.").arg(this->waveLengthFilename,this->frequencyFilename));
            return false;
        }
    }
    if(nrOfLambdas<1){
        this->logger->writeToLogFile(trUtf8("Erorr:The number of wavelengths/freuencies should at least one."));
        return false;
    }
    this->lambdaSize=nrOfLambdas;
    return true;
}
/**
 * @return lambda's array.
 */
QVector<double> SetUp::getLambdaArray(){
    return this->lambdaArray;
}

/**
 * @return size of lambda's array.
 */
int SetUp::getLambdaSize(){
    return this->lambdaSize;
}

/**
*  @brief This method reads the frequencies from frequency's file.
*
* @details This method reads the first line of the frequency's file.
* The first line contains the number of frequencies in this file. This
* method reads frequencies appropriate to this number.
* If there is a problem with the opening/reading the files, the method
* protocols that. This method calculates the wavelength and saves the results in lambda's array.
*
* @param The frequency's file, that should be read.
*
* @return The number of frequencies.
*/
int SetUp::readFrequencies(QFile *file){
    int nrOfLambdas=0;
    bool ok=false;
    if(file->exists()){
        QTextStream *in=new QTextStream(file);
        file->open(QIODevice::ReadOnly);
        nrOfLambdas=in->readLine().toInt(&ok);
        if(!ok){
            this->logger->writeToLogFile(trUtf8("Erorr: Can not reads the number of frequencies correctly in %1.")
                                         .arg(this->frequencyFilename));
            this->logger->stop(EXIT_FAILURE);
        }
        if(nrOfLambdas<1){
            this->logger->writeToLogFile(trUtf8("Erorr: The number of frequencies should be greater than zero in %1. Current value is %2.")
                                         .arg(this->frequencyFilename,QString::number(nrOfLambdas)));
            this->logger->stop(EXIT_FAILURE);
        }
        this->lambdaArray.clear();
        double erg=PhysicalConstants::getInstance()->getLightSpeed()*1e4;
        int frequency;
        for (int lambdaIndex = 0; lambdaIndex < nrOfLambdas; ++lambdaIndex) {
            frequency=in->readLine().toDouble(&ok);
            if(!ok){
                this->logger->writeToLogFile(trUtf8("Erorr: Can not reads frequency at index %1 correctly in %2.")
                                             .arg(QString::number(lambdaIndex),this->frequencyFilename));
                this->logger->stop(EXIT_FAILURE);
            }
            this->lambdaArray.append(erg/frequency);
        }

    } else{
        this->logger->writeToLogFile(trUtf8("Erorr: The %1 does not exists in curren Directory %2.")
                                     .arg(this->frequencyFilename,QDir::currentPath()));
        this->logger->stop(EXIT_FAILURE);
    }
    file->close();
    this->logger->writeToLogFile(trUtf8("%1 wavelengths are saved in lambda's Array from %2").arg(QString::number(nrOfLambdas),frequencyFilename));
    return nrOfLambdas;
}


/**
*  @brief This method reads the wavelengths from wavelength's file.
*
* @details This method reads the first line of the wavelength's file.
* The first line contains the number of wavelengths in this file. This
* method reads wavelengths appropriate to this number and saves them in lambda's array
* If there is a problem with the opening/reading the files, the method
* protocols that.
*
* @param The Wavelength's file, that should be read.
*
* @return The number of wavelengths.
*/
int SetUp::readWaveLengths(QFile *file){
    int nrOfLambdas=0;
    bool ok=false;
    if(file->exists()){
        QTextStream *in=new QTextStream(file);
        file->open(QIODevice::ReadOnly);
        nrOfLambdas=in->readLine().toInt(&ok);
        if(!ok){
            this->logger->writeToLogFile(trUtf8("Erorr: Can not reads the number of lambda correctly in %1.")
                                         .arg(this->waveLengthFilename));
            this->logger->stop(EXIT_FAILURE);
        }
        if(nrOfLambdas<1){
            this->logger->writeToLogFile(trUtf8("Erorr: The number of lambda should be greater than zero in %1. Current value is %2.")
                                         .arg(this->waveLengthFilename,QString::number(nrOfLambdas)));
            this->logger->stop(EXIT_FAILURE);
        }
        this->lambdaArray.clear();
        for (int lambdaIndex = 0; lambdaIndex < nrOfLambdas; ++lambdaIndex) {
            this->lambdaArray.append(in->readLine().toDouble(&ok));
            if(!ok){
                this->logger->writeToLogFile(trUtf8("Erorr: Can not reads lambda at index %1 correctly in %2.")
                                             .arg(QString::number(lambdaIndex),this->waveLengthFilename));
                this->logger->stop(EXIT_FAILURE);
            }
        }

    } else{
        this->logger->writeToLogFile(trUtf8("Erorr:The %1 does not exists in curren Directory %2")
                                     .arg(this->waveLengthFilename,QDir::currentPath()));
        this->logger->stop(EXIT_FAILURE);
    }
    file->close();
    this->logger->writeToLogFile(trUtf8("%1 wavelengths are saved in lambda's Array from %2")
                                 .arg(QString::number(nrOfLambdas),this->waveLengthFilename));
    return nrOfLambdas;
}

/**
 *  @brief This method returns local setting's file.
 *  @return QSettings
 */
QSettings *SetUp::getLocalSettings(){
    this->localSettings=new QSettings(localSettingsFilename,QSettings::IniFormat);
    return this->localSettings;
}
/**
 *  @brief This method returns global setting's file.
 *  @return QSettings
 */
QSettings *SetUp::getGlobalSettings(){
    this->globalSettings=new QSettings(globalSettingsFilename,QSettings::IniFormat);
    return this->globalSettings;
}

/**
 * @brief This method checks the existence of RADMC3D's executable file in current working directory
 *  and sets radmc3dLocal's value appropriate to it's result.
 */
void SetUp::checkRadmc3dlocal(){
    if(!radmc3dFilename.isEmpty()){
        QFile file(radmc3dFilename);
    }else{
        logger->writeToLogFile(trUtf8("Erorr: while finding the location of radmc3d's executable file: "
                                      "the RADMC3D's filename is empty"));
        logger->stop(EXIT_FAILURE);
    }
}

/**
 * @brief This method returns boolean value about existence of RADMC3D's executable file in current working directory.
 * @return  a boolean value. If RADMC3D's executable file exists, it will return true, otherwise false.
 */
bool SetUp::isRadmc3dLocal(){
    return radmc3dLocal;
}

/**
 * @brief This method returns temporary's filename
 * @return temporaryFolder as QString.
 */
QString SetUp::getTemporaryFolder(){
    return this->temporaryFolder;
}

/**
 * @brief This method returns local setting's filename
 * @return settingsFilename as QString.
 */

QString SetUp::getLocalSettingsFilename(){
    return this->localSettingsFilename;
}
/**
 * @brief This method returns global setting's filename
 * @return settingsFilename as QString.
 */

QString SetUp::getGlobalSettingsFilename(){
    return this->globalSettingsFilename;
}
/**
 * @brief This method returns camera-wavelength's filename
 * @return camerWavelengthFilename as QString.
 */

QString SetUp::getCamerWavelengthFilename(){
    return this->camerWavelengthFilename;
}
/**
 * @brief This method returns color's filename
 * @return colorInusFilename as QString.
 */
QString SetUp::getColorInusFilename(){
    return this->colorInusFilename;
}

/**
 * @brief This method returns frequency's filename
 * @return frequencyFilename as QString.
 */
QString SetUp::getFrequencyFilename(){
    return this->frequencyFilename;
}

/**
 * @brief This method returns image's fileending
 * @return imageEnding as QString.
 */

QString SetUp::getImageFileEnding(){
    return this->imageEnding;
}
/**
 * @brief This method returns image's filename
 * @return imageFilename as QString.
 */

QString SetUp::getImageFilename(){
    return this->imageFilename;
}
/**
 * @brief This method returns movieinput's filename
 * @return movieInputFilename as QString.
 */

QString SetUp::getMovieInputFilename(){
    return this->movieInputFilename;
}
/**
 * @brief This method returns RADMC3D's logfilename
 * @return radmc3dLogFilename as QString.
 */

QString SetUp::getRadmc3dLogFilename(){
    return this->radmc3dLogFilename;
}
/**
 * @brief This method returns RADMC3D's executable filename
 * @return radmc3dFilename as QString.
 */

QString SetUp::getradmc3dFilename(){
    return this->radmc3dFilename;
}
/**
 * @brief This method returns spectrum's fileending
 * @return spectrumEnding as QString.
 */

QString SetUp::getSpectrumEnding(){
    return this->spectrumEnding;
}
/**
 * @brief This method returns spectrum's filename
 * @return spectrumFilename as QString.
 */

QString SetUp::getSpectrumFilename(){
    return this->spectrumFilename;
}

/**
 * @brief This method returns transfer's filename
 * @return transferFilename as QString.
 */

QString SetUp::getTransferDefaultFilename(){
    return this->transferDefaultFilename;
}
/**
 * @brief This method returns wavelength's filename
 * @return waveLengthFilename as QString.
 */

QString SetUp::getWaveLengthFilename(){
    return this->waveLengthFilename;
}

/**
 * @brief Thismethod sets the standard filenames.
 */
void SetUp::setStandardFilenames(){
    QStringList amrGridList;
    amrGridList<<"amr_grid.inp"<<"amr_grid.uinp"<<"amr_grid.binp";
    this->setAmrGridFilenames(amrGridList);
    this->setCamerWavelengthFilename("camera_wavelength_micron.inp");
    this->setColorInusFilename("color_inus.inp");
    this->setFrequencyFilename("frequency.inp");
    this->setLinesFilename("lines.inp");
    this->setImageFilename("image",".out");
    this->setMovieInputFilename("movie.inp");
    this->setRadmc3dLogFilename("radmc3d.out");
    this->setradmc3dFilename("radmc3d");
    this->setSpectrumFilename("spectrum",".out");
    this->setTransferDefaultFilename("transfer.inp.default");
    this->setTransferFilename("transfer.inp");
    this->setWaveLengthFilename("wavelength_micron.inp");
    this->setLocalSettingsFilename("setting.ini");
    QString filename=QCoreApplication::applicationDirPath();
    if(!filename.endsWith(QDir::separator ()))filename=filename+QDir::separator ();
    filename=filename+"setting.ini";
    this->setGlobalSettingsFilename(filename);
    this->setHistoryFileName("history");

}



/**
 * @brief This method returns AMR-GRID's filenames
 * @return amrGridFilenames as QStringList.
 */
QStringList SetUp::getAmrGridFilenames(){
    return this->amrGrid->getAmrGridFilenames();
}
/**
 * @brief This method sets temporary's filename
 * @param path as QString.
 */

void SetUp::setTemporaryFolder(QString path){
    this->temporaryFolder=path;
    this->logger->writeToLogFile(trUtf8("The Temporary folder is %1").arg(this->temporaryFolder));
}
/**
 * @brief This method sets AMR-GRID-filenames appropriate to the given filenames.
 * @param filenames as QStringlist.
 */
void SetUp::setAmrGridFilenames(QStringList filenames){

    this->amrGrid->setAmrGridFilenames(filenames);
}
/**
 * @brief This method sets local settingfilename appropriate to the given filename.
 * @param filename as QString.
 */
void SetUp::setLocalSettingsFilename(QString filename){
    this->localSettingsFilename=filename;
    this->logger->writeToLogFile(trUtf8("local Setting's filename is %1").arg(this->localSettingsFilename));

}
/**
 * @brief This method sets gloabl settingfilename appropriate to the given filename.
 * @param filename as QString.
 */
void SetUp::setGlobalSettingsFilename(QString filename){
    this->globalSettingsFilename=filename;
    this->logger->writeToLogFile(trUtf8("global Setting's filename is %1").arg(this->localSettingsFilename));

}
/**
 * @brief This method sets camera-wavelength-filename appropriate to the given filename.
 * @param filename as QString.
 */
void SetUp::setCamerWavelengthFilename(QString filename){
    this->camerWavelengthFilename=filename;
    this->logger->writeToLogFile(trUtf8("Camre-Wavelength's filename is %1").arg(this->camerWavelengthFilename));

}

/**
 * @brief This method sets color-filename appropriate to the given filename.
 * @param filename as QString.
 */

void SetUp::setColorInusFilename(QString filename){
    this->colorInusFilename=filename;
    this->logger->writeToLogFile(trUtf8("Color's filename is %1").arg(this->colorInusFilename ));

}

/**
 * @brief This method sets image-filename appropriate to the given filename and fileending.
 * @param filename as QString.
 * @param fileending as QString.
 */

void SetUp::setImageFilename(QString filename, QString fileEnding){
    this->imageFilename=filename;
    this->imageEnding=fileEnding;
    this->logger->writeToLogFile(trUtf8("image's filename is %1%2").arg(this->imageFilename,imageEnding));
}

/**
 * @brief This method sets  Movie input's filename appropriate to the given filename.
 * @param filename as QString.
 */
void SetUp::setMovieInputFilename(QString filename){
    this->movieInputFilename=filename;
    this->logger->writeToLogFile(trUtf8("Movie's input's filename is %1 ").arg(this->movieInputFilename));

}
/**
 * @brief This method sets  RADMC3D's logfilename appropriate to the given filename.
 * @param filename as QString.
 */
void SetUp::setRadmc3dLogFilename(QString filename){
    this->radmc3dLogFilename=filename;
    this->logger->writeToLogFile(trUtf8("RADMC3D's logfilename is %1").arg(this->radmc3dLogFilename));
}
/**
 * @brief This method sets  RADMC3D's executable filename appropriate to the given filename.
 * @param filename as QString.
 */
void SetUp::setradmc3dFilename(QString filename){
    this->radmc3dFilename=filename;
    this->logger->writeToLogFile(trUtf8("RADMC3D's filename is %1").arg(this->radmc3dFilename));

}

void SetUp::setTransferFilename (QString filename){
    this->transferFilename=filename;
    this->logger->writeToLogFile(trUtf8("Transfer's  filename is %1").arg( this->transferFilename));
}

/**
 * @brief This method sets  spectrum's filename appropriate to the given filename and fileending.
 * @param filename as QString.
 * @param fileending as QString.
 */
void SetUp::setSpectrumFilename(QString filename, QString fileEnding){
    this->spectrumFilename=filename;
    this->spectrumEnding=fileEnding;
    this->logger->writeToLogFile(trUtf8("spectrum's filename is %1%2").arg(this->spectrumFilename,spectrumEnding));

}
/**
 * @brief This method sets  transfer's filename appropriate to the given filename.
 * @param filename as QString.
 */
void SetUp::setTransferDefaultFilename(QString filename){
    this->transferDefaultFilename=filename;
    this->logger->writeToLogFile(trUtf8("Transfer's default filename is %1").arg( this->transferDefaultFilename));

}
/**
 * @brief This method sets  frequency's filename appropriate to the given filename.
 * @param filename as QString.
 */
void SetUp::setFrequencyFilename(QString filename){
    this->frequencyFilename=filename;
    this->logger->writeToLogFile(trUtf8("frequency's filename is %1").arg( this->frequencyFilename));

}
/**
 * @brief This method sets  wavelength's filename appropriate to the given filename.
 * @param filename as QString.
 */
void SetUp::setWaveLengthFilename(QString filename){
    this->waveLengthFilename=filename;
    this->logger->writeToLogFile(trUtf8("waveLength's filename is %1").arg(this->waveLengthFilename));

}
