#ifndef IMAGECONTROLLER_H
#define IMAGECONTROLLER_H
#include "logger.h"
#include "setup.h"

#if QT_VERSION < 0x050000
      #include <QtGui>
#else
       #include <QtWidgets>
#endif

class ImageControllerWidget : public QWidget
{
    Q_OBJECT
public:
    explicit ImageControllerWidget(QWidget *parent = 0);
    bool isLocalModus(){return this->localModus;}
    void reset();
    void setSizeLineEdit();
    void setZoomBox(QVector<double> vector){this->zoomBox=vector;zoomActive=true;}
    QVector<double> getZoomBox(){return this->zoomBox;}
    bool isZoomActive(){return this->zoomActive;}
    void setZoomActivity(bool active){this->zoomActive=active;}
    QMap<QString,QString> *labelColorMap;
    static const QString labelFontColor;
    static const QString CheckBoxInactiveFontColor;
    void save();
    void saveAsDefaults();
    void load();
    void loadAsDefaults();
    const static int whatsThisIconHeight=64;
    const static int whatsThisIconWidth=64;
    QVector <QRgb> getCurrentColorLookUpTable();
    QVector <QRgb> getGrayColorLookUpTable(){return this->colorLookUpTables.at(0);}
    QVector <QRgb> getRedColorLookUpTable(){return this->colorLookUpTables.at(4);}
    void loadColorLookUpTables(QString path,bool saveGlobal);
    QString getColorTableSliderObjectName();
    QString getAntiAliasingCheckBoxObjectName(){return this->antiAliasingCheckBox->objectName();}
    QString getInvertColorsCheckBoxObjectName(){return this->invertColorsCheckBox->objectName();}
    QString getContourCheckBoxObjectName(){return this->contourCheckBox->objectName();}
    QString getNumberOfContoursObjectName(){return this->numberOfContoursSpinBox->objectName();}
    QString getRenderImagePushButtonObjectName(){return this->renderImagePushButton->objectName();}
    QString getRenderSpectrumPushButtonObjectName(){return this->renderSpectrumPushButton->objectName();}
    QString getPositionComboxBoxObjectName(){return this->positionComboxBox->objectName();}
    QString getAbsoluteScaleCheckBoxObjectName(){return this->absoluteScaleCheckBox->objectName();}
    QString getSizeLineEditObjectName(){return this->sizeLineEdit->objectName();}
    int getPosition(){return this->positionComboxBox->currentIndex();}
    bool isColorInvert(){return this->invertColorsCheckBox->isChecked();}
    bool isSecondOrder(){return this->secondOrderCheckBox->isChecked();}
    bool isPreview(){return this->previewCheckBox->isChecked();}
    bool isStar(){return this->starCheckBox->isChecked();}
    bool isAntiAliasing(){return this->antiAliasingCheckBox->isChecked();}
    bool isTau(){return this->tauCheckBox->isChecked();}
    double getSaturate(){return this->saturationSpinBox->value()/100;}
    bool isLinear(){return this->linearCheckBox->isChecked();}
    bool isContour(){return this->contourCheckBox->isChecked();}
    bool isColorBar(){return this->colorbarCheckBox->isChecked();}
    QString getSizeUnit();
    QString getContrastObjectName(){return this->contrastSpinBox->objectName();}
    QString getLinearObjectName(){return this->linearCheckBox->objectName();}
    QString getColorBarObjectName(){return this->colorbarCheckBox->objectName();}
    QString getSaturateObjectName(){return this->saturationSpinBox->objectName();}
    int getNumberOfContours(){return this->numberOfContoursSpinBox->value();}
    QVector <bool> appendCommandVector;
    void setNumberOfPixelsSpinBox(int number){
        this->numberOfPixelsSpinBox->blockSignals(true);
        this->numberOfPixelsSpinBox->setValue(number);
        this->numberOfPixels=number;
        this->numberOfPixelsSpinBox->blockSignals(false);
    }
    double getMaxLog(){return this->contrastSpinBox->value();}
    void setSizeUnit(SetUp::unit);
    int getNumberOfPixels(){return this->numberOfPixelsSpinBox->value();}
    double getPhi(){return this->phiSpinBox->value();}
    double getInclination(){return this->inclinationSpinBox->value();}
    double getLambdaRed(){return this->redLambdaValueLineEdit->text().remove(lambdaUnit).toDouble();}
    double getLambdaBlue(){return this->blueLambdaValueLineEdit->text().remove(lambdaUnit).toDouble();}
    double getLambdaGreen(){return this->greenLambdaValueLineEdit->text().remove(lambdaUnit).toDouble();}
    double getSize(){
        return this->sizeLineEdit->text().remove(sizeUnit).toDouble();}
    double getPositionAngle(){return this->positionAngleSpinBox->value();}
    QString getShowRadmc3dPushButtonObjectName(){return this->showRadmc3dPushButton->objectName();}
    bool isPipe(){return this->pipeCheckBox->isChecked();}
    QString getPipeCheckBoxObjectName(){return this->pipeCheckBox->objectName();}
    QString  getUnzoomImagePushButtonObjectName(){return this->unzoomImagePushButton->objectName();}
    QString getInclinationObjectName(){return this->inclinationSpinBox->objectName();}
    QString getPhiObjectName(){return this->phiSpinBox->objectName();}
    QString getRedColorTuneValueLineEditObjectName(){return this->redColorTuneValueLineEdit->objectName();}
    QString getBlueColorTuneValueLineEditObjectName(){return this->greenColorTuneValueLineEdit->objectName();}
    QString getGreenColorTuneValueLineEditObjectName(){return this->blueColorTuneValueLineEdit->objectName();}
    bool isRGBModus(){return this->rgbModus;}
    void writeLambdasToCamerWavelengthFilename();
    void setLastColorMaximumArray(QVector<double> vec);
    void setLastColorMinimumArray(QVector<double> vec);
    double* getLastColorMaximumArray(){return this->lastColorMaximumArray;}
    double* getLastColorMinimumArray(){return this->lastColorMinimumArray;}
    bool isAbsoluteScale(){return this->absoluteScaleCheckBox->isChecked();}
    void updateImageMenuActionVisibilties();
    QVector<double> lastColorTune;
    QString getRedColorTuneValueLineEditText(){return this->redColorTuneValueLineEdit->text();}
    QString getBlueColorTuneValueLineEditText(){return this->blueColorTuneValueLineEdit->text();}
    QString getGreenColorTuneValueLineEditText(){return this->greenColorTuneValueLineEdit->text();}
    void setUpLambdaIndexSlider();
    bool isDopplerCatchingAvailable(){return this->dopplerCatchingAvailable;}
    void setDoLine(bool doLine){this->doLine=doLine;}
    bool isDoLine(){return this->doLine;}
    bool isDopplerCatchingActive(){return this->dopplerCatchCheckBox->isChecked();}
    void blockMoleculeSpinBoxSignal(bool block){this->moleculeComboBox->blockSignals(block);this->moleculeComboBox->setEnabled(!block);}
    void blockLineSpinBoxSignal(bool block){this->lineSpinBox->blockSignals(block);this->lineSpinBox->setEnabled(!block);}
    void blockVelocitySpinBoxSignal(bool block){this->velocitySpinBox->blockSignals(block);this->velocitySpinBox->setEnabled(!block);}
    int getIMolecule(){return this->moleculeComboBox->currentIndex()+1;}
    int getILine(){return this->lineSpinBox->value();}
    double getVelocity(){return this->velocitySpinBox->value();}
    bool isLineModus(){return this->lineModus;}
    void setLineModus(bool active);
    void setRedLambdaLineEdit(double);
    void setUpColorTableSlider();
    void setUpColorTableComboBoxValues();
    void setUpLocalTab();
    bool isRelativeObserver(){return this->relobCheckbox->isChecked();}
    void setLocalModus(bool active){this->localModus=active;}
    double getPointXLineEditValue(){return this->pointXLineEdit->text().remove(QRegExp(sizeRegExp)).toDouble();}
    double getPointYLineEditValue(){return this->pointYLineEdit->text().remove(QRegExp(sizeRegExp)).toDouble();}
    double getPointZLineEditValue(){return this->pointZLineEdit->text().remove(QRegExp(sizeRegExp)).toDouble();}
    double getObserverposXLineEditValue(){return this->observerposXLineEdit->text().remove(QRegExp(sizeRegExp)).toDouble();}
    double getObserverposYLineEditValue(){return this->observerposYLineEdit->text().remove(QRegExp(sizeRegExp)).toDouble();}
    double getObserverposZLineEditValue(){return this->observerposZLineEdit->text().remove(QRegExp(sizeRegExp)).toDouble();}
    double getViewangleSpinBoxValue(){return this->viewangleSpinBox->value();}
    double getScaleSpinBoxValue(){return this->scaleSpinBox->value();}
    void setRenderSpectrumActive(bool active){this->renderSpectrumActive=active;}
    bool isRenderSpectrumActive(){return this->renderSpectrumActive;}
    bool isObjectNameBelongsToturnToRenderImageList(QString objectname){return turnToRenderImageList.contains(objectname);}
    QString getUpPushButtonObjectName(){return this->upPushButton->objectName();}
    QString getDownPushButtonObjectName(){return this->downPushButton->objectName();}
    QString getLeftPushButtonObjectName(){return this->leftPushButton->objectName();}
    QString getRightPushButtonObjectName(){return this->rightPushButton->objectName();}
    void setSendCommandLineList(bool send){this->sendCommandLineList=send;}
    bool getsendCommandLineList(){ return this->sendCommandLineList;}
    bool isUsertransferOn(){return this->usertransferOn;}
    QMap<QString,QString> getAppendCommand();
private:
    bool dopplerCatchingAvailable;
    bool vertical;
    bool zoomActive;
    bool renderSpectrumActive;
    bool sendCommandLineList;
    bool lineModus;
    bool doLine;
    bool rgbModus;
    bool localModus;
    QSet<QString> turnToRenderImageList;
    QLabel *pointingLabel;
    QLabel *viewangleLabel;
    QCheckBox *relobCheckbox;
    QDoubleSpinBox *scaleSpinBox;
    QDoubleSpinBox *viewangleSpinBox;
    QLineEdit *pointXLineEdit;
    QLineEdit *pointYLineEdit;
    QLineEdit *pointZLineEdit;
    bool usertransferOn;
    double lastPointX;
    double lastPointY;
    double lastPointZ;
    double lastobsX;
    double lastobsY;
    double lastobsZ;
    QLabel *observerposLabel;
    QLineEdit *observerposXLineEdit;
    QLineEdit *observerposYLineEdit;
    QLineEdit *observerposZLineEdit;
    QLineEdit *redColorTuneValueLineEdit;
    QLineEdit *blueColorTuneValueLineEdit;
    QLineEdit *greenColorTuneValueLineEdit;
    double lastVelocityValue;
    int lastLineValue;
    int lastMoleculeValue;
    QTableWidget *table;
    QMap<QString,QString> map;
    int numberOfContours;
    int numberOfPixels;
    double phi,inclination,saturation,posAngle,maxlog;
    QVector <QRgb>  createLookUpTableLinear(int rStart,int rEnd,int gStart,int gEnd,int bStart,int bEnd, int numberOfColor );
    QListView *debugger;
    QGroupBox *debuggerGroupBox;
    QPushButton *clearDebuggerPushButton;
    QPushButton *showRadmc3dPushButton;
    QPushButton *scollDebuggerToBottomPushButton;
    QPushButton *scollDebuggerToTopPushButton;
    QVBoxLayout *debuggerVBox;
    QVBoxLayout *usertransferLayout;
    QGroupBox *usertransferGroupBox;
    QPushButton *resetUserTransferPushButton;
    QTabWidget *middleTabWidget;
    QCheckBox *starCheckBox;
    QCheckBox *localCheckBox;
    QCheckBox *secondOrderCheckBox;
    QCheckBox *tauCheckBox;
    QCheckBox *previewCheckBox;
    QCheckBox *dopplerCatchCheckBox;
    QCheckBox *contourCheckBox;
    QCheckBox *antiAliasingCheckBox;
    QCheckBox *linearCheckBox;
    QCheckBox *colorbarCheckBox;
    QCheckBox *pipeCheckBox;
    QCheckBox *absoluteScaleCheckBox;
    QCheckBox *invertColorsCheckBox;
    QGroupBox *checkBoxGroupBox;
    QSlider *colorTableSlider;
    QLabel *colorTableLabel;
    QLabel *positionLabel;
    QComboBox *colorTableComboxBox;
    QComboBox *positionComboxBox;
    QLabel * colorTableMinLabel;
    QLabel * colorTableMaxLabel;
    QGridLayout *mainLayout;
    QLabel *contrastLabel;
    QDoubleSpinBox *contrastSpinBox;
    QLabel *saturationLabel;
    QDoubleSpinBox *saturationSpinBox;
    QLabel *numberOfContoursLabel;
    QSpinBox *numberOfContoursSpinBox;
    QLabel *positionAngleLabel;
    QDoubleSpinBox *positionAngleSpinBox;
    QLabel *numberOfPixelsLabel;
    QSpinBox *numberOfPixelsSpinBox;
    QLabel *phiLabel;
    QDoubleSpinBox *phiSpinBox;
    QLabel *inclinationLabel;
    QDoubleSpinBox *inclinationSpinBox;
    QLabel *sizeLabel;
    QLabel *moleculeLabel;
    QComboBox *moleculeComboBox;
    QLabel *lineLabel;
    QSpinBox *lineSpinBox;
    QLabel *velocityLabel;
    QDoubleSpinBox *velocitySpinBox;
    QSlider *redLambdaIndexSlider;
    QLineEdit *redLambdaValueLineEdit;
    QLabel *redLambdaLabel;
    QSlider *blueLambdaIndexSlider;
    QLineEdit *blueLambdaValueLineEdit;
    QLabel *blueLambdaLabel;
    QSlider *greenLambdaIndexSlider;
    QLineEdit *greenLambdaValueLineEdit;
    QLabel *greenLambdaLabel;
    QLineEdit *sizeLineEdit;
    QGridLayout * pushButtonLayout;
    QHBoxLayout * allpushButtonLayout;
    QGridLayout *checkBoxLayout;
    QGridLayout *spinBoxLayout;
    QGroupBox * pushButtonGroupBox;
    QGroupBox * specpushButtonGroupBox;
    QGroupBox *spinBoxGroupBox;
    QShortcut *phiShortCut;
    QShortcut *inclinationShortCut;
    QShortcut *redLambdaShortCut;
    QShortcut *blueLambdaShortCut;
    QShortcut *greenLambdaShortCut;
    QPushButton *renderImagePushButton;
    QPushButton *renderSpectrumPushButton;
    QVector<double> zoomBox;
    QString sizeRegExp;
    QString lambdaRegExp;
    QPushButton *leftPushButton;
    QPushButton *rightPushButton;
    QPushButton *upPushButton;
    QPushButton *downPushButton;
    QPushButton *unzoomImagePushButton;
    QGroupBox *pipeAntiGroupBox;
    SetUp *setup;
    double lastSizeValue;
    double lastRedColorTuneValue;
    double lastBlueColorTuneValue;
    double lastGreenColorTuneValue;
    const static int CheckBoxIconHeight;
    const static int CheckBoxIconWidth;
    const static int pushButtonIconWidth;
    const static int pushButtonIconheight;
    const static int labelIconHeight;
    const static int labelIconWidth;
    const static int comboBoxIconHeight;
    const static int comboBoxIconWidth;
    const static double contrastSpinBoxMaximum;
    const static double contrastSpinBoxMinimum;
    const static double contrastSpinBoxStep;
    const static double saturationSpinBoxMaximum;
    const static double saturationSpinBoxMinimum;
    const static double saturationSpinBoxStep;
    const static double positionAngleSpinBoxMaximum;
    const static double positionAngleSpinBoxMinimum;
    const static double positionAngleSpinBoxStep;
    const static double inclinationSpinBoxMaximum;
    const static double inclinationSpinBoxMinimum;
    const static double inclinationSpinBoxStep;
    const static double inclinationSpinBoxDefaultValue;

    const static double scaleSpinBoxMaximum;
    const static double scaleSpinBoxMinimum;
    const static double scaleSpinBoxStep;
    const static double scaleSpinBoxDefaultValue;

    const static double viewangleSpinBoxMaximum;
    const static double viewangleSpinBoxMinimum;
    const static double viewangleSpinBoxStep;
    const static double viewangleSpinBoxDefaultValue;

    const static double phiSpinBoxMaximum;
    const static double phiSpinBoxMinimum;
    const static double phiSpinBoxStep;
    const static double phiSpinBoxDefaultValue;

    const static double velocitySpinBoxMaximum;
    const static double velocitySpinBoxMinimum;
    const static double velocitySpinBoxStep;
    const static double velocitySpinBoxDefaultValue;


    const static int numberOfContourSpinBoxMaximum;
    const static int numberOfContourSpinBoxMinimum;
    const static int numberOfContourSpinBoxStep;
    const static int numberOfPixelsSpinBoxMaximum;
    const static int numberOfPixelsSpinBoxMinimum;
    const static int numberOfPixelsSpinBoxStep;
    const static int numberOfPixelsSpinBoxDefaultValue;
    static const QString tauCheckBoxActiveFontColor;
    static const QString secondOrderCheckBoxActiveFontColor;
    static const QString contourCheckBoxActiveFontColor;
    static const QString starCheckBoxActiveFontColor;
    static const QString previewCheckBoxActiveFontColor;
    static const QString linearCheckBoxActiveFontColor;
    static const QString colorbarCheckBoxActiveFontColor;
    static const QString invertColorsCheckBoxActiveFontColor;
    static const QString pipeCheckBoxActiveFontColor;
    static const QString northWestGroupBoxBorderColor;
    Logger *logger;
    void setColorLookUpTables();
private:
    double lastColorMinimumArray[3];
    double lastColorMaximumArray[3];
    void addToColorLookUpTables(QVector <QRgb> table,QString name);
    QVector<QRgb> readColorLookUpTable(QString filename,QString path,bool saveGlobal);
    QVector< QVector <QRgb> > colorLookUpTables;
    QVector< QString > colorLookUpTableNames;
    int lambdaSize;
    QVector<double> lambdaArray;
    QString sizeUnit;
    QString lambdaUnit;
    void addWidgetsToMainLayout(bool);
    void translate();
    void setLambdaUnit(QString);
    void setLambdaSize(int);
    int getLambdaSize();
    QVector<double> getLambdaArray();
    void setLambdaArray(QVector<double>);
    QString getLambdaUnit();
    void setUp();
    void setLambdaLabelText();
    void setUpWhatsThis();
    void setUpContrastSpinBox();
    void setUpNumberOfContoursSpinBox();
    void setUpNumberOfPixelsSpinBox();
    void setUpSaturationSpinBbox();
    void setUpPositionAngleSpinBoxSpinBbox();
    void setUpInclinationSpinBoxSpinBbox();
    void setUpPhiSpinBoxSpinBbox();
    void setUpLineSpinBox();
    void setUpVelocityBoxSpinBbox();
    void setUpMoleculeSpinBox();
    void setConstrastSpinBoxEnabled(bool active);
    void setUpPositionComboBoxValues();
    void setUpScaleSpinBox();
    void setUpViewangleSpinBox();
    void setSlots();
    void setStyleSheets();
    void setShortcuts();
    void setIcons();
    void setToolTips();
    void setStarToolTip(bool active);
    void setAntiAliasingToolTip(bool active);
    void setDopplerCatchingToolTip(bool active);
    void setContrastToolTip(double);
    void setSizeToolTip(QString);
    void setRedLambdaToolTip(QString);
    void setGreenLambdaToolTip(QString);
    void setBlueLambdaToolTip(QString);
    void setRelobToolTip(bool);
    void setNumberOfContoursToolTip(int);
    void setNumberOfPixelsToolTip(int);
    void setPositionAngleToolTip(double);
    void setPhiToolTip(double);
    void setInclinationToolTip(double);
    void setSaturationToolTip(double);
    void setTauToolTip(bool active);
    void setColorbarToolTip(bool active);
    void setSecondOrderToolTip(bool active);
    void setContourToolTip(bool active);
    void setLinearToolTip(bool active);
    void setPreviewToolTip(bool active);
    void setLocalToolTip(bool active);
    void setPipeToolTip(bool active);
    void setColorTuneEnabled(bool active);
    void setInvertColorsToolTip(bool active);
public slots:
    void positionComboxBoxValueChanged(int);
private slots:
    void showRadmc3dPushButtonClicked();
    void itemChanged(QTableWidgetItem*);
    void activateStar(bool active);
    void activateRelob(bool active);
    void activateAntiAliasing(bool active);
    void activateDopplerCatching(bool active);
    void activatePipe(bool active);
    void activateColorbar(bool active);
    void activateSecondOrder(bool active);
    void activateInvertColors(bool active);
    void activateAbsoluteScale(bool active);
    void activateTau(bool active);
    void activatePreview(bool active);
    void activateLocal(bool active);
    void activateContour(bool active);
    void activateLinear(bool active);
    void colorTableSliderValueChanged(int);
    void redLambdaIndexSliderValueChanged(int);
    void blueLambdaIndexSliderValueChanged(int);
    void greenLambdaIndexSliderValueChanged(int);
    void colorTableComboxBoxValueChanged(int);
    void contrastSpinBoxValueChanged();
    void saturationSpinBoxValueChanged();
    void positionAngleSpinBoxValueChanged();
    void phiSpinBoxValueChanged();
    void inclinationSpinBoxValueChanged();
    void sizeLineEditValueChanged();
    void redLambdaValueLineEditValueChanged();
    void greenLambdaValueLineEditValueChanged();
    void blueLambdaValueLineEditValueChanged();
    void redColorTuneValueLineEditValueChanged();
    void moleculeComboBoxValueChanged(int);
    void lineSpinBoxValueChanged();
    void viewangleSpinBoxValueChanged();
    void scaleSpinBoxValueChanged();
    void velocitySpinBoxValueChanged();
    void greenColorTuneValueLineEditValueChanged();
    void blueColorTuneValueLineEditValueChanged();
    void numberOfContoursSpinBoxValueChanged();
    void numberOfPixelsSpinBoxValueChanged();
    void pointXLineEditValueChanged();
    void pointYLineEditValueChanged();
    void pointZLineEditValueChanged();
    void observerposXLineEditValueChanged();
    void observerposYLineEditValueChanged();
    void observerposZLineEditValueChanged();
    void unzoomImagePushButtonClicked();
    void renderImagePushButtonClicked();
    void renderSpectrumPushButtonClicked();
    void upPushButtonClicked();
    void downPushButtonClicked();
    void rightPushButtonClicked();
    void leftPushButtonClicked();
    void middleTabWidgetCurrentIndex(int);
    void scrollDebuggerToTop();
    void scrollDebuggerToBottom();
    void clearDebugger();
    void resetTransferFile();
signals:
    void sendCommand(const QString objectName);
    void checktStateChanged(QString);
    void readMyAction();
    void removeRect();
    void setActionVisible(const QString objectName,bool visibility);
public slots:
    void addColorModusWidgets();
    void removeColorModusWidgets();
    void activateLocalObserver(bool);
    void sizeUnitChanged(QAction*);
    void removeTransferTab();
    void readTransferFile(bool load=true);
    void selectDebugger();
    void changePhiAndInclination(double phi,double inclination);
    
};

#endif // IMAGECONTROLLER_H
