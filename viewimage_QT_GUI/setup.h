#ifndef SETUP_H
#define SETUP_H
#include <QtCore>
#include "logger.h"

#include "amrgrid.h"
class SetUp:public QObject
{
    Q_OBJECT
public:
    explicit SetUp(QObject *parent = 0);
    static SetUp* getInstance(QObject *parent = 0);
    void setStandardFilenames();
    void setWaveLengthFilename (QString filename);
    QString getWaveLengthFilename ();
    void setFrequencyFilename (QString filename);
    QString getLinesFilname(){return this->linesFilename;}
    void setLinesFilename(QString filename ){this->linesFilename=filename;}
    QString getFrequencyFilename ();
    void setradmc3dFilename (QString filename);
    QString getradmc3dFilename ();
    void setAmrGridFilenames (QStringList filenames);
    QStringList getAmrGridFilenames ();
    void setTransferDefaultFilename (QString filename);
    QString getTransferDefaultFilename ();
    void setTransferFilename (QString filename);
    QString getTransferFilename (){return this->transferFilename;}
    void setImageFilename (QString filename,QString fileEnding);
    QString getImageFilename ();
    QString getImageFileEnding ();
    void setSpectrumFilename (QString filename,QString fileEnding);
    QString getSpectrumFilename ();
    QString getSpectrumEnding ();
    void setCamerWavelengthFilename (QString filename);
    QString getCamerWavelengthFilename ();
    void setColorInusFilename (QString filename);
    QString getColorInusFilename ();
    void setMovieInputFilename (QString filename);
    QString getMovieInputFilename ();
    void setRadmc3dLogFilename (QString filename);
    QString getRadmc3dLogFilename ();
    void setTemporaryFolder(QString path);
    QString getTemporaryFolder();
    void setLambda();
    bool isRadmc3dLocal();
    QString getLocalSettingsFilename();
    QString getGlobalSettingsFilename();
    void setLocalSettingsFilename(QString filename);
    void setGlobalSettingsFilename(QString filename);
    QSettings* getLocalSettings();
    QSettings* getGlobalSettings();
    bool readLambda();
    int readFrequencies(QFile* file);
    bool IsLambaRead(){return this->lambdaRead;}
    void setLambdaRead(bool read){this->lambdaRead=read;}
    int readWaveLengths(QFile* file);
    QVector<double> getLambdaArray();
    int getLambdaSize();
    void setGreenLambda(double lambda);
    void setRedLambda(double lambda);
    void setBlueLambda(double lambda);
    int getGreenLambda(){return this->greenLambda;}
    int getBlueLambda(){return this->blueLambda;}
    int getRedLambda(){return this->redLambda;}
    int getBlueLambdaIndex(){return this->blueLambdaIndex;}
    int getRedLambdaIndex(){return this->redLambdaIndex;}
    int getGreenLambdaIndex(){return this->greenLambdaIndex;}
    int getIndexOfLambda(double lambda,double min,double max);
    enum unit{
        CM=0,PC=1,AU=2
    };
    double calculateImageSize();
    unit getImageUnit(){return this->imageunit;}
    void setImageUnit(unit iUnit){this->imageunit=iUnit;}
    void setAmrGrid(bool basic=true,bool mirror=false);
    void setAmrGridRead(bool read){this->amrGrid->setAmrGridRead(read);}
    bool isAmrGridRead(){return this->amrGrid->isAmrGridRead();}
    QString getHistoryFileName(){return this->historyFileName;}
    void setHistoryFileName(QString filename){this->historyFileName=filename;}
    void checkRadmc3dlocal();
private:
    QString historyFileName;
    QString linesFilename;
    Cartesian *cartesian;
    Polar *polar;
    unit imageunit;
    bool lambdaRead;
    bool pipe;
    bool radmc3dLocal;
    bool isPipeMode(){return this->pipe;}
    void setPipeMode(bool on){this->pipe=on;}
    AmrGrid *amrGrid;
    double blueLambda;
    double redLambda;
    PhysicalConstants *physicalConstants;
    double greenLambda;
    int blueLambdaIndex;
    int redLambdaIndex;
    int greenLambdaIndex;
    QSettings *localSettings;
    QSettings *globalSettings;
    static SetUp* instance;
    QString temporaryFolder;
    Logger *logger;
    QVector<double> lambdaArray;
    int lambdaSize;
    QString standardOutputDuplicator;
    QString localSettingsFilename;
    QString globalSettingsFilename;
    QString waveLengthFilename;
    QString guiInput;
    QString guiOutput;
    QString frequencyFilename;
    QString radmc3dFilename;
    QString transferDefaultFilename;
    QString transferFilename;
    QString imageFilename;
    QString spectrumFilename;
    QString spectrumEnding;
    QString imageEnding;
    QString camerWavelengthFilename;
    QString colorInusFilename;
    QString movieInputFilename;
    QString radmc3dLogFilename;
};

#endif // SETUP_H
