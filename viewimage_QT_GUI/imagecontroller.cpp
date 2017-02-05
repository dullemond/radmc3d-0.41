#include "imagecontroller.h"
/**
*   @file mainwindow.cpp
*   @brief The class imageController is derived from Qwidget. This class will
*   contain elements to release image's events.
*   @author Farzin Sereshti
*   @version 1.0
*/
const QString ImageControllerWidget::CheckBoxInactiveFontColor="C0C0C0";
const QString ImageControllerWidget::contourCheckBoxActiveFontColor="05B513";
const QString ImageControllerWidget::tauCheckBoxActiveFontColor="F8083C";
const QString ImageControllerWidget::starCheckBoxActiveFontColor="C8870F";
const QString ImageControllerWidget::secondOrderCheckBoxActiveFontColor="0F60FF";
const QString ImageControllerWidget::previewCheckBoxActiveFontColor="45C4F6";
const QString ImageControllerWidget::linearCheckBoxActiveFontColor="1C7013";
const QString ImageControllerWidget::colorbarCheckBoxActiveFontColor="20DD30";
const QString ImageControllerWidget::invertColorsCheckBoxActiveFontColor="C52FC2";
const QString ImageControllerWidget::pipeCheckBoxActiveFontColor=contourCheckBoxActiveFontColor;
const QString ImageControllerWidget::northWestGroupBoxBorderColor="A09C9C";
const QString ImageControllerWidget::labelFontColor="1C7013";
const  int ImageControllerWidget::CheckBoxIconHeight=22;
const  int ImageControllerWidget::CheckBoxIconWidth=22;
const  int ImageControllerWidget::pushButtonIconWidth=48;
const  int ImageControllerWidget::pushButtonIconheight=48;
const  int ImageControllerWidget::labelIconHeight=14;
const  int ImageControllerWidget::labelIconWidth=14;
const  int ImageControllerWidget::comboBoxIconHeight=18;
const  int ImageControllerWidget::comboBoxIconWidth=18;
const  double ImageControllerWidget::contrastSpinBoxMaximum=255.00000000000000;
const  double ImageControllerWidget::contrastSpinBoxMinimum=0.00000000000000;
const  double ImageControllerWidget::contrastSpinBoxStep=1.00000000000000;
const  double ImageControllerWidget::saturationSpinBoxMaximum=100.0000000000000;
const  double ImageControllerWidget::saturationSpinBoxMinimum=0.00000000000000;
const  double ImageControllerWidget::saturationSpinBoxStep=1.000000000000000;
const  double ImageControllerWidget::positionAngleSpinBoxMaximum=36000.00000000000000;
const  double ImageControllerWidget::positionAngleSpinBoxMinimum=-36000.00000000000000;
const  double ImageControllerWidget:: positionAngleSpinBoxStep=1.000000000000000;
const  double ImageControllerWidget::inclinationSpinBoxMaximum=36000.00000000000000;
const  double ImageControllerWidget::inclinationSpinBoxMinimum=-36000.00000000000000;
const  double ImageControllerWidget::inclinationSpinBoxStep=30.000000000000000;
const  double ImageControllerWidget::inclinationSpinBoxDefaultValue=60.000000000000000;

const  double ImageControllerWidget::scaleSpinBoxMaximum=36000.00000000000000;
const  double ImageControllerWidget::scaleSpinBoxMinimum=-36000.00000000000000;
const  double ImageControllerWidget::scaleSpinBoxStep=1.000000000000000;
const  double ImageControllerWidget::scaleSpinBoxDefaultValue=0.00000000000000;

const  double ImageControllerWidget::viewangleSpinBoxMaximum=36000.00000000000000;
const  double ImageControllerWidget::viewangleSpinBoxMinimum=-36000.00000000000000;
const  double ImageControllerWidget::viewangleSpinBoxStep=1.000000000000000;
const  double ImageControllerWidget::viewangleSpinBoxDefaultValue=0.80000000000000;

const  double ImageControllerWidget::phiSpinBoxMaximum=36000.00000000000000;
const  double ImageControllerWidget::phiSpinBoxMinimum=-36000.00000000000000;
const  double ImageControllerWidget::phiSpinBoxStep=30.000000000000000;
const  double ImageControllerWidget::phiSpinBoxDefaultValue=30.000000000000000;

const  double ImageControllerWidget::velocitySpinBoxMaximum=300000.00000000000000;
const  double ImageControllerWidget::velocitySpinBoxMinimum=-300000.00000000000000;
const  double ImageControllerWidget::velocitySpinBoxStep=1.000000000000000;
const  double ImageControllerWidget::velocitySpinBoxDefaultValue=0.000000000000000;


const  int ImageControllerWidget::numberOfContourSpinBoxMaximum=1000;
const  int ImageControllerWidget::numberOfContourSpinBoxMinimum=1.0;
const  int ImageControllerWidget:: numberOfContourSpinBoxStep=1;
const  int ImageControllerWidget:: numberOfPixelsSpinBoxMaximum=800;
const  int ImageControllerWidget::numberOfPixelsSpinBoxMinimum=1.0;
const  int ImageControllerWidget::numberOfPixelsSpinBoxStep=100;
const  int ImageControllerWidget::numberOfPixelsSpinBoxDefaultValue=100;

/**
 *  @brief This is the contsructur.
 *
 *  @param a widget parent can be given. This a is optional parameter.
 *  @details This method initializes the image-controller's widget.
 *  @see setUp()
 *
 */
ImageControllerWidget::ImageControllerWidget(QWidget *parent) :
    QWidget(parent)
{
    renderSpectrumActive=false;
    zoomActive=false;
    doLine=false;

    this->dopplerCatchingAvailable=false;
    this->vertical=false;
    this->sendCommandLineList=false;
    this->lineModus=false;
    this->rgbModus=false;
    this->localModus=false;
    this->setup=SetUp::getInstance();
    this->logger=Logger::getInstance();
    this->setObjectName("image");
    this->starCheckBox=new QCheckBox();
    this->starCheckBox->setObjectName("star");



    this->previewCheckBox=new QCheckBox();
    this->previewCheckBox->setObjectName("preview");


    this->secondOrderCheckBox=new QCheckBox();
    this->secondOrderCheckBox->setObjectName("secondorder");


    this->contourCheckBox=new QCheckBox();
    this->contourCheckBox->setObjectName("contour");

    this->tauCheckBox=new QCheckBox();
    this->tauCheckBox->setObjectName("tau");

    this->linearCheckBox=new QCheckBox();
    this->linearCheckBox->setObjectName("linear");


    this->colorbarCheckBox=new QCheckBox();
    this->colorbarCheckBox->setObjectName("colorbar");

    this->colorTableSlider=new QSlider(Qt::Horizontal);
    this->colorTableSlider->setObjectName("colortable");


    this->colorTableLabel=new QLabel();
    this->colorTableComboxBox=new QComboBox();

    this->invertColorsCheckBox=new QCheckBox();
    this->invertColorsCheckBox->setObjectName("invertcolors");

    this->pipeCheckBox=new QCheckBox();
    this->pipeCheckBox->setObjectName("pipe");

    this->antiAliasingCheckBox=new QCheckBox();
    this->antiAliasingCheckBox->setObjectName("antiAliasing");

    this->absoluteScaleCheckBox=new QCheckBox();
    this->absoluteScaleCheckBox->setObjectName("scale");
    this->absoluteScaleCheckBox->setVisible(false);

    this->dopplerCatchCheckBox=new QCheckBox();
    this->dopplerCatchCheckBox->setObjectName("dc");
    this->dopplerCatchCheckBox->setVisible(true);

    this->localCheckBox=new QCheckBox();
    this->localCheckBox->setObjectName(trUtf8("local"));



    this->mainLayout=new QGridLayout();


    QGridLayout *box=new QGridLayout();
    box->addWidget(this->colorTableLabel,0,0,Qt::AlignLeft);
    box->addWidget(colorTableComboxBox,0,1);
    this->colorTableLabel->setBuddy(this->colorTableComboxBox);
    //box->addWidget(this->colorTableSlider,3,0,1,2);
    this->setColorLookUpTables();
    //    this->colorTableMinLabel=new QLabel(QString::number(1));
    //    this->colorTableMaxLabel=new QLabel(QString::number(this->colorLookUpTables.size()));
    //    this->colorTableMinLabel->setAlignment(Qt::AlignLeft);
    //    box->addWidget(this->colorTableMinLabel,4,0,1,1,Qt::AlignLeft);
    //    this->colorTableMaxLabel->setAlignment(Qt::AlignRight);
    //    box->addWidget(this->colorTableMaxLabel,4,1,1,1,Qt::AlignRight);


    QVBoxLayout *vDirectionBox=new QVBoxLayout();
    this->positionComboxBox=new QComboBox();
    this->positionComboxBox->setObjectName("position");
    this->positionLabel=new QLabel();
    this->positionLabel->setBuddy(this->positionComboxBox);
    this->positionLabel->setMargin(1);
    vDirectionBox->addWidget(this->positionLabel);
    vDirectionBox->addWidget(this->positionComboxBox);



    this->checkBoxLayout=new QGridLayout();
    this->checkBoxLayout->setMargin(0);
    this->checkBoxLayout->setContentsMargins(0,0,0,0);

    this->checkBoxLayout->addWidget(this->starCheckBox,0,0);
    this->checkBoxLayout->addWidget(this->secondOrderCheckBox,0,1);
    this->checkBoxLayout->addWidget(this->tauCheckBox,0,2);
    this->checkBoxLayout->addWidget(this->pipeCheckBox,1,0);
    this->checkBoxLayout->addWidget(this->previewCheckBox,1,1);
    this->checkBoxLayout->addWidget(this->linearCheckBox,1,2);
    QVBoxLayout *pipeAnti=new QVBoxLayout();
    pipeAnti->setAlignment(Qt::AlignLeft);
    pipeAnti->setContentsMargins(0,0,0,0);
    pipeAntiGroupBox=new QGroupBox();
    pipeAntiGroupBox->setLayout(pipeAnti);
    pipeAnti->addWidget(this->colorbarCheckBox,Qt::AlignLeft);
    pipeAnti->addWidget(this->antiAliasingCheckBox,Qt::AlignLeft);
    this->checkBoxLayout->addWidget(this->invertColorsCheckBox,2,0);
    this->checkBoxLayout->addWidget(this->contourCheckBox,2,1);
    this->checkBoxLayout->addLayout(vDirectionBox,2,2);
    this->checkBoxLayout->addWidget(pipeAntiGroupBox,3,0,1,0);
    this->checkBoxLayout->addWidget(this->absoluteScaleCheckBox,3,1,1,3);
    this->checkBoxLayout->addLayout(box,3,1,1,2);
    this->checkBoxLayout->addWidget(this->dopplerCatchCheckBox,4,0,1,2);
    this->checkBoxLayout->addWidget(this->localCheckBox,4,1,1,1);


    this->checkBoxGroupBox=new  QGroupBox();
    this->checkBoxGroupBox->setLayout(this->checkBoxLayout);

    this->spinBoxLayout=new QGridLayout();
    this->spinBoxLayout->setSpacing(0);
    this->spinBoxGroupBox=new QGroupBox();
    this->spinBoxGroupBox->setLayout(this->spinBoxLayout);

    this->middleTabWidget=new QTabWidget();
    middleTabWidget->setContentsMargins(0,0,0,0);
    this->middleTabWidget->addTab(this->checkBoxGroupBox,QIcon(tr(":/images/%1.png").arg(this->objectName())),tr("CheckBoxes"));
    this->middleTabWidget->addTab(this->spinBoxGroupBox,QIcon(tr(":/images/spin.png")),trUtf8("SpinBoxes"));
    this->debuggerGroupBox=new QGroupBox();
    this->debuggerVBox=new QVBoxLayout();
    this->debuggerGroupBox->setLayout(this->debuggerVBox);
    QHBoxLayout *debuggerButtonLayout=new QHBoxLayout();
    this->clearDebuggerPushButton=new QPushButton();
    this->clearDebuggerPushButton->setObjectName("clearDebugger");
    this->scollDebuggerToBottomPushButton=new QPushButton();
    this->scollDebuggerToBottomPushButton->setObjectName("scrollDebuggerToBottom");
    this->scollDebuggerToTopPushButton=new QPushButton();
    this->scollDebuggerToTopPushButton->setObjectName("scrollDebuggerToTop");
    this->showRadmc3dPushButton=new QPushButton();
    this->showRadmc3dPushButton->setObjectName("showradmc3d");
    debuggerButtonLayout->setMargin(0);
    this->debuggerVBox->setMargin(0);
    debuggerButtonLayout->setContentsMargins(0,0,0,0);
    this->debuggerVBox->setContentsMargins(0,0,0,0);
    debuggerButtonLayout->addWidget(this->scollDebuggerToTopPushButton);
    debuggerButtonLayout->addWidget(this->scollDebuggerToBottomPushButton);
    debuggerButtonLayout->addWidget(this->clearDebuggerPushButton);
    this->debugger=new QListView();
    this->debuggerVBox->addLayout(debuggerButtonLayout);
    this->debuggerVBox->addWidget(debugger);
    this->debuggerVBox->addWidget(this->showRadmc3dPushButton);
    this->debugger->setModel( this->logger->getMessageModel());
    this->middleTabWidget->addTab(this->debuggerGroupBox,QIcon(tr(":/images/bug.png")),trUtf8("Debugger"));
    this->middleTabWidget->setDocumentMode(true);


    this->contrastLabel=new QLabel();
    this->contrastSpinBox=new QDoubleSpinBox();
    this->contrastLabel->setBuddy(this->contrastSpinBox);
    this->contrastSpinBox->setObjectName("contrast");


    this->saturationLabel=new QLabel();
    this->saturationSpinBox=new QDoubleSpinBox();
    this->saturationLabel->setBuddy(this->saturationSpinBox);
    this->saturationSpinBox->setObjectName("saturation");


    this->numberOfContoursLabel=new QLabel();
    this->numberOfContoursSpinBox=new QSpinBox();
    this->numberOfContoursLabel->setBuddy(this->numberOfContoursSpinBox);
    this->numberOfContoursSpinBox->setObjectName("numberOfContours");

    this->positionAngleLabel=new QLabel();
    this->positionAngleSpinBox=new QDoubleSpinBox();
    this->positionAngleLabel->setBuddy(this->positionAngleSpinBox);
    this->positionAngleSpinBox->setObjectName("positionAngle");

    this->numberOfPixelsLabel=new QLabel();
    this->numberOfPixelsSpinBox=new QSpinBox();
    this->numberOfPixelsLabel->setBuddy(this->numberOfPixelsSpinBox);
    this->numberOfPixelsSpinBox->setObjectName("numberOfPixels");

    this->sizeLabel=new QLabel();
    this->sizeLineEdit=new QLineEdit();
    this->sizeLabel->setBuddy(sizeLineEdit);
    this->sizeLineEdit->setObjectName("size");

    this->inclinationLabel=new QLabel();
    this->inclinationSpinBox=new QDoubleSpinBox();
    this->inclinationSpinBox->setObjectName("inclination");

    this->phiLabel=new QLabel();
    this->phiSpinBox=new QDoubleSpinBox();
    this->phiSpinBox->setObjectName("phi");

    this->redLambdaLabel=new QLabel();
    this->redLambdaIndexSlider=new QSlider(Qt::Horizontal);
    this->redLambdaValueLineEdit=new QLineEdit();
    this->redLambdaIndexSlider->setObjectName("redLambdaIndex");
    this->redLambdaValueLineEdit->setObjectName("redLambda");

    QVBoxLayout *redLambdaLayout=new QVBoxLayout();
    redLambdaLayout->addWidget(this->redLambdaLabel);
    this->redLambdaLabel->setBuddy(this->redLambdaValueLineEdit);
    redLambdaLayout->addWidget(this->redLambdaIndexSlider);

    this->blueLambdaLabel=new QLabel();
    this->blueLambdaIndexSlider=new QSlider(Qt::Horizontal);
    this->blueLambdaValueLineEdit=new QLineEdit();
    this->blueLambdaIndexSlider->setObjectName("blueLambdaIndex");
    this->blueLambdaValueLineEdit->setObjectName("blueLambda");

    QVBoxLayout *blueLambdaLayout=new QVBoxLayout();
    blueLambdaLayout->addWidget(this->blueLambdaLabel);
    this->blueLambdaLabel->setBuddy(this->blueLambdaValueLineEdit);
    blueLambdaLayout->addWidget(this->blueLambdaIndexSlider);



    this->greenLambdaLabel=new QLabel();
    this->greenLambdaIndexSlider=new QSlider(Qt::Horizontal);
    this->greenLambdaValueLineEdit=new QLineEdit();
    this->greenLambdaIndexSlider->setObjectName("greenLambdaIndex");
    this->greenLambdaValueLineEdit->setObjectName("greenLambda");

    QVBoxLayout *greenLambdaLayout=new QVBoxLayout();
    greenLambdaLayout->addWidget(this->greenLambdaLabel);
    this->greenLambdaLabel->setBuddy(this->greenLambdaValueLineEdit);
    greenLambdaLayout->addWidget(this->greenLambdaIndexSlider);


    this->redColorTuneValueLineEdit=new QLineEdit();
    this->redColorTuneValueLineEdit->setObjectName("redColorTune");
    this->blueColorTuneValueLineEdit=new QLineEdit();
    this->blueColorTuneValueLineEdit->setObjectName("blueColorTune");
    this->greenColorTuneValueLineEdit=new QLineEdit();
    this->greenColorTuneValueLineEdit->setObjectName("greenColorTune");

    this->moleculeLabel=new QLabel();
    this->moleculeComboBox=new QComboBox();
    this->moleculeComboBox->setObjectName("imol");


    this->lineLabel=new QLabel();
    this->lineSpinBox=new QSpinBox();
    this->lineSpinBox->setObjectName("line");
    this->lineSpinBox->setKeyboardTracking(false);

    this->velocityLabel=new QLabel();
    this->velocitySpinBox=new QDoubleSpinBox();
    this->velocitySpinBox->setObjectName("velocity");
    this->velocitySpinBox->setKeyboardTracking(false);

    this->spinBoxLayout->addWidget(this->contrastLabel,0,0,1,1);
    this->spinBoxLayout->addWidget(this->contrastSpinBox,0,1,1,1);
    this->spinBoxLayout->addWidget(this->saturationLabel,0,2,1,1);
    this->spinBoxLayout->addWidget(this->saturationSpinBox,0,3,1,2);
    this->spinBoxLayout->addWidget(this->positionAngleLabel,1,0,1,1);
    this->spinBoxLayout->addWidget(this->positionAngleSpinBox,1,1,1,1);
    this->spinBoxLayout->addWidget(this->numberOfContoursLabel,1,2,1,1);
    this->spinBoxLayout->addWidget(this->numberOfContoursSpinBox,1,3,1,2);
    this->spinBoxLayout->addWidget(this->numberOfPixelsLabel,2,0,1,1);
    this->spinBoxLayout->addWidget(this->numberOfPixelsSpinBox,2,1,1,1);
    this->spinBoxLayout->addWidget(this->sizeLabel,2,2,1,1);
    this->spinBoxLayout->addWidget(this->sizeLineEdit,2,3,1,2);
    this->spinBoxLayout->addWidget(this->inclinationLabel,3,0,1,1);
    this->spinBoxLayout->addWidget(this->inclinationSpinBox,3,1,1,1);
    this->spinBoxLayout->addWidget(this->phiLabel,3,2,1,1);
    this->spinBoxLayout->addWidget(this->phiSpinBox,3,3,1,2);
    this->spinBoxLayout->addLayout(redLambdaLayout,4,0,1,3);
    this->spinBoxLayout->addWidget(this->redLambdaValueLineEdit,4,3,1,1);
    this->spinBoxLayout->addWidget(this->redColorTuneValueLineEdit,4,4,1,1);
    this->spinBoxLayout->addLayout(greenLambdaLayout,5,0,1,3);
    this->spinBoxLayout->addWidget(this->greenLambdaValueLineEdit,5,3,1,1);
    this->spinBoxLayout->addWidget(this->greenColorTuneValueLineEdit,5,4,1,1);
    this->spinBoxLayout->addLayout(blueLambdaLayout,6,0,1,3);
    this->spinBoxLayout->addWidget(this->blueLambdaValueLineEdit,6,3,1,1);
    this->spinBoxLayout->addWidget(this->blueColorTuneValueLineEdit,6,4,1,1);
    this->spinBoxLayout->addWidget(this->moleculeLabel,7,0,1,1);
    this->spinBoxLayout->addWidget(this->moleculeComboBox,7,1,1,1);
    this->spinBoxLayout->addWidget(this->lineLabel,7,2,1,1);
    this->spinBoxLayout->addWidget(this->lineSpinBox,7,3,1,1);
    this->spinBoxLayout->addWidget(this->velocityLabel,8,0,1,1);
    this->spinBoxLayout->addWidget(this->velocitySpinBox,8,1,1,1);



    this->pushButtonLayout=new QGridLayout();
    this->pushButtonGroupBox=new QGroupBox();

    this->pushButtonGroupBox->setLayout(this->pushButtonLayout);
    this->renderImagePushButton=new QPushButton();
    this->renderImagePushButton->setObjectName("renderImage");

    this->unzoomImagePushButton=new QPushButton();
    this->unzoomImagePushButton->setObjectName("zoomOut");

    this->upPushButton=new QPushButton();
    this->upPushButton->setObjectName("up");

    this->downPushButton=new QPushButton();
    this->downPushButton->setObjectName("down");

    this->leftPushButton=new QPushButton();
    this->leftPushButton->setObjectName("left");

    this->rightPushButton=new QPushButton();
    this->rightPushButton->setObjectName("right");


    this->renderSpectrumPushButton=new QPushButton();
    this->renderSpectrumPushButton->setObjectName("renderSpectrum");

    this->pushButtonLayout->addWidget(this->renderImagePushButton,0,0);
    this->pushButtonLayout->addWidget(this->unzoomImagePushButton,1,0);
    this->pushButtonLayout->setSpacing(0);

    this->pushButtonLayout->setSizeConstraint(QLayout::SetFixedSize);
    this->pushButtonLayout->addWidget(this->upPushButton,0,2);
    this->pushButtonLayout->addWidget(this->leftPushButton,1,1);
    this->pushButtonLayout->addWidget(this->downPushButton,1,2);
    this->pushButtonLayout->addWidget(this->rightPushButton,1,3);

    this->specpushButtonGroupBox=new QGroupBox();
    QVBoxLayout *vbx=new QVBoxLayout();
    vbx->addWidget(this->renderSpectrumPushButton);
    this->specpushButtonGroupBox->setLayout(vbx);



    this->addWidgetsToMainLayout(true);

    this->labelColorMap=new QMap<QString,QString>();
    this->labelColorMap->insert(this->starCheckBox->objectName(),starCheckBoxActiveFontColor);
    this->labelColorMap->insert(this->secondOrderCheckBox->objectName(),this->secondOrderCheckBoxActiveFontColor);
    this->labelColorMap->insert(this->tauCheckBox->objectName(),this->tauCheckBoxActiveFontColor);
    this->labelColorMap->insert(this->contourCheckBox->objectName(),this->contourCheckBoxActiveFontColor);
    this->labelColorMap->insert(this->previewCheckBox->objectName(),this->previewCheckBoxActiveFontColor);
    this->labelColorMap->insert(linearCheckBox->objectName(),this->linearCheckBoxActiveFontColor);
    this->labelColorMap->insert(this->colorbarCheckBox->objectName(),this->colorbarCheckBoxActiveFontColor);
    this->labelColorMap->insert(this->invertColorsCheckBox->objectName(),this->invertColorsCheckBoxActiveFontColor);
    this->labelColorMap->insert(this->pipeCheckBox->objectName(),this->pipeCheckBoxActiveFontColor);
    this->labelColorMap->insert(this->antiAliasingCheckBox->objectName(),this->pipeCheckBoxActiveFontColor);
    this->labelColorMap->insert(this->absoluteScaleCheckBox->objectName(),this->colorbarCheckBoxActiveFontColor);
    this->labelColorMap->insert(this->dopplerCatchCheckBox->objectName(),this->invertColorsCheckBoxActiveFontColor);
    this->labelColorMap->insert(this->localCheckBox->objectName(),this->previewCheckBoxActiveFontColor);



    this->setLayout(this->mainLayout);
    this->setUp();
    this->numberOfContours=this->numberOfContoursSpinBox->value();
    this->numberOfPixels=this->numberOfPixelsSpinBox->value();
    this->saturation=saturationSpinBox->value();
    this->phi=this->phiSpinBox->value();
    this->inclination=this->inclinationSpinBox->value();
    this->maxlog=this->contrastSpinBox->value();
    this->posAngle=this->positionAngleSpinBox->value();
    turnToRenderImageList.insert(this->starCheckBox->objectName());
    turnToRenderImageList.insert(this->inclinationSpinBox->objectName());
    turnToRenderImageList.insert(this->phiSpinBox->objectName());
    turnToRenderImageList.insert(this->renderSpectrumPushButton->objectName());
    turnToRenderImageList.insert(this->sizeLineEdit->objectName());
    turnToRenderImageList.insert(this->contrastSpinBox->objectName());
    turnToRenderImageList.insert(this->pipeCheckBox->objectName());
    turnToRenderImageList.insert(this->positionComboxBox->objectName());
    turnToRenderImageList.insert(this->upPushButton->objectName());
    turnToRenderImageList.insert(this->downPushButton->objectName());
    turnToRenderImageList.insert(this->leftPushButton->objectName());
    turnToRenderImageList.insert(this->rightPushButton->objectName());

}
void ImageControllerWidget::setLastColorMaximumArray(QVector<double> vec){
    this->lastColorMaximumArray[0]=vec.at(0);
    this->lastColorMaximumArray[1]=vec.at(1);
    this->lastColorMaximumArray[2]=vec.at(2);
}

void ImageControllerWidget::setLastColorMinimumArray(QVector<double> vec){
    this->lastColorMinimumArray[0]=vec.at(0);
    this->lastColorMinimumArray[1]=vec.at(1);
    this->lastColorMinimumArray[2]=vec.at(2);
}
void ImageControllerWidget::updateImageMenuActionVisibilties(){
    emit setActionVisible(this->relobCheckbox->objectName(),isLocalModus());
    emit setActionVisible(this->previewCheckBox->objectName(),!isLocalModus());
    emit setActionVisible(this->secondOrderCheckBox->objectName(),!isLocalModus());
    emit setActionVisible(this->localCheckBox->objectName(),isLocalModus());
    if(QFile(this->setup->getLinesFilname()).exists()){
        this->dopplerCatchCheckBox->setVisible(true);
        this->dopplerCatchingAvailable=true;
        emit setActionVisible(this->dopplerCatchCheckBox->objectName(),true);
    }else{
        this->dopplerCatchCheckBox->setVisible(false);
        this->dopplerCatchingAvailable=false;
        emit setActionVisible(this->dopplerCatchCheckBox->objectName(),false);
    }

    emit setActionVisible(this->absoluteScaleCheckBox->objectName(),rgbModus);
    emit setActionVisible(this->contourCheckBox->objectName(),!rgbModus);
    emit setActionVisible(this->colorbarCheckBox->objectName(),!rgbModus);
}

void ImageControllerWidget::writeLambdasToCamerWavelengthFilename(){
    QFile file(this->setup->getCamerWavelengthFilename());
    file.open(QFile::WriteOnly|QFile::Text);
    QTextStream in(&file);
    in<<"3"<<endl;
    in<<this->getLambdaRed()<<endl;
    in<<this->getLambdaGreen()<<endl;
    in<<this->getLambdaBlue()<<endl;
    file.close();
}
void ImageControllerWidget::addColorModusWidgets(){
    rgbModus=true;
    this->absoluteScaleCheckBox->setVisible(rgbModus);
    this->blueLambdaLabel->setVisible(rgbModus);
    this->blueLambdaValueLineEdit->setVisible(rgbModus);
    this->blueLambdaIndexSlider->setVisible(rgbModus);
    this->greenLambdaLabel->setVisible(rgbModus);
    this->greenLambdaValueLineEdit->setVisible(rgbModus);
    this->greenLambdaIndexSlider->setVisible(rgbModus);
    this->redColorTuneValueLineEdit->setVisible(rgbModus);
    this->blueColorTuneValueLineEdit->setVisible(rgbModus);
    this->greenColorTuneValueLineEdit->setVisible(rgbModus);
    //this->colorbarCheckBox->setVisible(false);
    this->colorTableComboxBox->setVisible(false);
    //    this->colorTableMaxLabel->setVisible(false);
    //    this->colorTableMinLabel->setVisible(false);
    this->colorTableLabel->setVisible(false);
    //    this->colorTableSlider->setVisible(false);
    this->updateImageMenuActionVisibilties();

}
void ImageControllerWidget::removeColorModusWidgets(){
    rgbModus=false;
    this->absoluteScaleCheckBox->setVisible(rgbModus);
    this->blueLambdaLabel->setVisible(rgbModus);
    this->blueLambdaValueLineEdit->setVisible(rgbModus);
    this->blueLambdaIndexSlider->setVisible(rgbModus);
    this->greenLambdaLabel->setVisible(rgbModus);
    this->greenLambdaValueLineEdit->setVisible(rgbModus);
    this->greenLambdaIndexSlider->setVisible(rgbModus);
    this->redColorTuneValueLineEdit->setVisible(rgbModus);
    this->blueColorTuneValueLineEdit->setVisible(rgbModus);
    this->greenColorTuneValueLineEdit->setVisible(rgbModus);
    // this->colorbarCheckBox->setVisible(true);
    this->colorTableComboxBox->setVisible(true);
    //    this->colorTableMaxLabel->setVisible(true);
    //    this->colorTableMinLabel->setVisible(true);
    //    this->colorTableSlider->setVisible(true);
    this->colorTableLabel->setVisible(true);
    this->updateImageMenuActionVisibilties();
}

void ImageControllerWidget::addWidgetsToMainLayout(bool vertical){
    this->vertical=vertical;

    this->mainLayout->removeWidget(this->middleTabWidget);
    allpushButtonLayout=new QHBoxLayout();
    allpushButtonLayout->addWidget(this->pushButtonGroupBox);
    allpushButtonLayout->addWidget(this->specpushButtonGroupBox);
    if(vertical){
        this->mainLayout->addWidget(this->middleTabWidget,0,0,1,1);
        this->mainLayout->addLayout(this->allpushButtonLayout,1,0,1,1);
    }else{
        this->mainLayout->addWidget(this->middleTabWidget,0,0,1,1);
        this->mainLayout->addLayout(this->allpushButtonLayout,1,0,1,1);
    }

}

/**
  * @brief This method sets colortable's slider
  *
  */
void ImageControllerWidget::setUpColorTableSlider(){
    this->colorTableSlider->setRange(1,this->colorLookUpTables.size());
    //this->colorTableMaxLabel->setText(QString::number(this->colorLookUpTables.size()));
    this->colorTableSlider->setSingleStep(1);
    this->colorTableSlider->setTickPosition(QSlider::TicksBelow);
    this->colorTableSlider->setTickInterval(1);

}
/**
  * @brief This method sets lambdaindex's slider
  *
  */
void ImageControllerWidget::setUpLambdaIndexSlider(){
    this->setLambdaArray(this->setup->getLambdaArray());
    this->setLambdaSize(this->setup->getLambdaSize());
    this->redLambdaIndexSlider->setRange(1,this->lambdaSize);
    this->redLambdaIndexSlider->setSingleStep(1);
    //because the first slider position is 1
    this->redLambdaIndexSlider->setValue(this->setup->getRedLambdaIndex()+1);
    this->redLambdaIndexSlider->setTickInterval(1);

    this->greenLambdaIndexSlider->setRange(1,this->lambdaSize);
    this->greenLambdaIndexSlider->setSingleStep(1);
    //because the first slider position is 1
    this->greenLambdaIndexSlider->setValue(this->setup->getGreenLambdaIndex()+1);
    this->greenLambdaIndexSlider->setTickInterval(1);

    this->blueLambdaIndexSlider->setRange(1,this->lambdaSize);

    this->blueLambdaIndexSlider->setSingleStep(1);
    //because the first slider position is 1
    this->blueLambdaIndexSlider->setValue(this->setup->getBlueLambdaIndex()+1);
    this->blueLambdaIndexSlider->setTickInterval(1);


}
/**  @brief this method return a CLUT.
  *  @return a vector, that contains a QRgb
  *
  */
QVector<QRgb> ImageControllerWidget::getCurrentColorLookUpTable(){
    int index=this->colorTableSlider->value()-1;
    if(index==-1)index=0;
    return colorLookUpTables.at(index);
}

void ImageControllerWidget::setColorLookUpTables(){
    QVector<QRgb>colorTable;
    for (int var = 0; var < 256; ++var) {
        colorTable.push_back(qRgb(var,var,var));
    }
    this->addToColorLookUpTables(colorTable,"Gray");


    this->addToColorLookUpTables(createLookUpTableLinear(191,255,95,255,0,188,256),"Blue/White linear");

    this->addToColorLookUpTables(createLookUpTableLinear(97,255,0,255,181,255,256),"Green/White linear");

    colorTable.clear();

    colorTable.append(qRgb(0, 0, 0));
    colorTable.append(qRgb(0, 0, 2));
    colorTable.append(qRgb(0, 0, 4));
    colorTable.append(qRgb(0, 0, 6));
    colorTable.append(qRgb(0, 0, 8));
    colorTable.append(qRgb(0, 0, 10));
    colorTable.append(qRgb(0, 0, 12));
    colorTable.append(qRgb(0, 0, 14));
    colorTable.append(qRgb(0, 0, 16));
    colorTable.append(qRgb(0, 0, 18));
    colorTable.append(qRgb(0, 0, 20));
    colorTable.append(qRgb(0, 0, 22));
    colorTable.append(qRgb(0, 0, 25));
    colorTable.append(qRgb(0, 0, 27));
    colorTable.append(qRgb(0, 0, 29));
    colorTable.append(qRgb(0, 0, 31));
    colorTable.append(qRgb(0, 0, 33));
    colorTable.append(qRgb(0, 0, 35));
    colorTable.append(qRgb(0, 0, 37));
    colorTable.append(qRgb(0, 0, 39));
    colorTable.append(qRgb(0, 0, 41));
    colorTable.append(qRgb(0, 0, 43));
    colorTable.append(qRgb(0, 0, 45));
    colorTable.append(qRgb(0, 0, 47));
    colorTable.append(qRgb(0, 0, 50));
    colorTable.append(qRgb(0, 0, 52));
    colorTable.append(qRgb(0, 0, 54));
    colorTable.append(qRgb(0, 0, 56));
    colorTable.append(qRgb(0, 0, 58));
    colorTable.append(qRgb(0, 0, 60));
    colorTable.append(qRgb(0, 0, 62));
    colorTable.append(qRgb(0, 0, 64));
    colorTable.append(qRgb(0, 0, 66));
    colorTable.append(qRgb(0, 3, 68));
    colorTable.append(qRgb(0, 6, 70));
    colorTable.append(qRgb(0, 9, 72));
    colorTable.append(qRgb(0, 12, 75));
    colorTable.append(qRgb(0, 15, 77));
    colorTable.append(qRgb(0, 18, 79));
    colorTable.append(qRgb(0, 21, 81));
    colorTable.append(qRgb(0, 25, 83));
    colorTable.append(qRgb(0, 28, 85));
    colorTable.append(qRgb(0, 31, 87));
    colorTable.append(qRgb(0, 34, 89));
    colorTable.append(qRgb(0, 37, 91));
    colorTable.append(qRgb(0, 40, 93));
    colorTable.append(qRgb(0, 43, 95));
    colorTable.append(qRgb(0, 46, 97));
    colorTable.append(qRgb(0, 50, 100));
    colorTable.append(qRgb(0, 53, 100));
    colorTable.append(qRgb(0, 56, 100));
    colorTable.append(qRgb(0, 59, 100));
    colorTable.append(qRgb(0, 62, 100));
    colorTable.append(qRgb(0, 65, 100));
    colorTable.append(qRgb(0, 68, 100));
    colorTable.append(qRgb(0, 71, 100));
    colorTable.append(qRgb(0, 75, 100));
    colorTable.append(qRgb(0, 78, 100));
    colorTable.append(qRgb(0, 81, 100));
    colorTable.append(qRgb(0, 84, 100));
    colorTable.append(qRgb(0, 87, 100));
    colorTable.append(qRgb(0, 90, 100));
    colorTable.append(qRgb(0, 93, 100));
    colorTable.append(qRgb(0, 96, 100));
    colorTable.append(qRgb(0, 100, 100));
    colorTable.append(qRgb(0, 103, 100));
    colorTable.append(qRgb(0, 106, 100));
    colorTable.append(qRgb(0, 109, 100));
    colorTable.append(qRgb(0, 112, 100));
    colorTable.append(qRgb(0, 115, 100));
    colorTable.append(qRgb(0, 118, 100));
    colorTable.append(qRgb(0, 121, 100));
    colorTable.append(qRgb(0, 125, 100));
    colorTable.append(qRgb(0, 128, 100));
    colorTable.append(qRgb(0, 131, 100));
    colorTable.append(qRgb(0, 134, 100));
    colorTable.append(qRgb(0, 137, 100));
    colorTable.append(qRgb(0, 140, 100));
    colorTable.append(qRgb(0, 143, 100));
    colorTable.append(qRgb(0, 146, 100));
    colorTable.append(qRgb(0, 150, 100));
    colorTable.append(qRgb(0, 150, 96));
    colorTable.append(qRgb(0, 150, 93));
    colorTable.append(qRgb(0, 150, 90));
    colorTable.append(qRgb(0, 150, 87));
    colorTable.append(qRgb(0, 150, 84));
    colorTable.append(qRgb(0, 150, 81));
    colorTable.append(qRgb(0, 150, 78));
    colorTable.append(qRgb(0, 150, 75));
    colorTable.append(qRgb(0, 150, 71));
    colorTable.append(qRgb(0, 150, 68));
    colorTable.append(qRgb(0, 150, 65));
    colorTable.append(qRgb(0, 150, 62));
    colorTable.append(qRgb(0, 150, 59));
    colorTable.append(qRgb(0, 150, 56));
    colorTable.append(qRgb(0, 150, 53));
    colorTable.append(qRgb(0, 150, 50));
    colorTable.append(qRgb(0, 149, 46));
    colorTable.append(qRgb(0, 148, 43));
    colorTable.append(qRgb(0, 148, 40));
    colorTable.append(qRgb(0, 147, 37));
    colorTable.append(qRgb(0, 146, 34));
    colorTable.append(qRgb(0, 146, 31));
    colorTable.append(qRgb(0, 145, 28));
    colorTable.append(qRgb(0, 145, 25));
    colorTable.append(qRgb(0, 144, 21));
    colorTable.append(qRgb(0, 143, 18));
    colorTable.append(qRgb(0, 143, 15));
    colorTable.append(qRgb(0, 142, 12));
    colorTable.append(qRgb(0, 141, 9));
    colorTable.append(qRgb(0, 141, 6));
    colorTable.append(qRgb(0, 140, 3));
    colorTable.append(qRgb(0, 140, 0));
    colorTable.append(qRgb(7, 137, 0));
    colorTable.append(qRgb(15, 135, 0));
    colorTable.append(qRgb(22, 132, 0));
    colorTable.append(qRgb(30, 130, 0));
    colorTable.append(qRgb(37, 127, 0));
    colorTable.append(qRgb(45, 125, 0));
    colorTable.append(qRgb(52, 122, 0));
    colorTable.append(qRgb(60, 120, 0));
    colorTable.append(qRgb(67, 117, 0));
    colorTable.append(qRgb(75, 115, 0));
    colorTable.append(qRgb(82, 112, 0));
    colorTable.append(qRgb(90, 110, 0));
    colorTable.append(qRgb(97, 107, 0));
    colorTable.append(qRgb(105, 105, 0));
    colorTable.append(qRgb(112, 102, 0));
    colorTable.append(qRgb(120, 100, 0));
    colorTable.append(qRgb(125, 93, 0));
    colorTable.append(qRgb(130, 87, 0));
    colorTable.append(qRgb(135, 81, 0));
    colorTable.append(qRgb(140, 75, 0));
    colorTable.append(qRgb(145, 68, 0));
    colorTable.append(qRgb(150, 62, 0));
    colorTable.append(qRgb(155, 56, 0));
    colorTable.append(qRgb(160, 50, 0));
    colorTable.append(qRgb(165, 43, 0));
    colorTable.append(qRgb(170, 37, 0));
    colorTable.append(qRgb(175, 31, 0));
    colorTable.append(qRgb(180, 25, 0));
    colorTable.append(qRgb(185, 18, 0));
    colorTable.append(qRgb(190, 12, 0));
    colorTable.append(qRgb(195, 6, 0));
    colorTable.append(qRgb(200, 0, 0));
    colorTable.append(qRgb(200, 2, 0));
    colorTable.append(qRgb(201, 4, 0));
    colorTable.append(qRgb(201, 6, 0));
    colorTable.append(qRgb(202, 9, 0));
    colorTable.append(qRgb(202, 11, 0));
    colorTable.append(qRgb(203, 13, 0));
    colorTable.append(qRgb(203, 16, 0));
    colorTable.append(qRgb(204, 18, 0));
    colorTable.append(qRgb(204, 20, 0));
    colorTable.append(qRgb(205, 23, 0));
    colorTable.append(qRgb(205, 25, 0));
    colorTable.append(qRgb(206, 27, 0));
    colorTable.append(qRgb(206, 29, 0));
    colorTable.append(qRgb(207, 32, 0));
    colorTable.append(qRgb(207, 34, 0));
    colorTable.append(qRgb(208, 36, 0));
    colorTable.append(qRgb(208, 39, 0));
    colorTable.append(qRgb(209, 41, 0));
    colorTable.append(qRgb(209, 43, 0));
    colorTable.append(qRgb(210, 46, 0));
    colorTable.append(qRgb(210, 48, 0));
    colorTable.append(qRgb(211, 50, 0));
    colorTable.append(qRgb(211, 53, 0));
    colorTable.append(qRgb(212, 55, 0));
    colorTable.append(qRgb(212, 57, 0));
    colorTable.append(qRgb(213, 59, 0));
    colorTable.append(qRgb(213, 62, 0));
    colorTable.append(qRgb(214, 64, 0));
    colorTable.append(qRgb(214, 66, 0));
    colorTable.append(qRgb(215, 69, 0));
    colorTable.append(qRgb(215, 71, 0));
    colorTable.append(qRgb(216, 73, 0));
    colorTable.append(qRgb(216, 76, 0));
    colorTable.append(qRgb(217, 78, 0));
    colorTable.append(qRgb(217, 80, 0));
    colorTable.append(qRgb(218, 83, 0));
    colorTable.append(qRgb(218, 85, 0));
    colorTable.append(qRgb(219, 87, 0));
    colorTable.append(qRgb(219, 89, 0));
    colorTable.append(qRgb(220, 92, 0));
    colorTable.append(qRgb(220, 94, 0));
    colorTable.append(qRgb(221, 96, 0));
    colorTable.append(qRgb(221, 99, 0));
    colorTable.append(qRgb(222, 101, 0));
    colorTable.append(qRgb(222, 103, 0));
    colorTable.append(qRgb(223, 106, 0));
    colorTable.append(qRgb(223, 108, 0));
    colorTable.append(qRgb(224, 110, 0));
    colorTable.append(qRgb(224, 113, 0));
    colorTable.append(qRgb(225, 115, 0));
    colorTable.append(qRgb(225, 117, 0));
    colorTable.append(qRgb(226, 119, 0));
    colorTable.append(qRgb(226, 122, 0));
    colorTable.append(qRgb(227, 124, 0));
    colorTable.append(qRgb(227, 126, 0));
    colorTable.append(qRgb(228, 129, 0));
    colorTable.append(qRgb(228, 131, 0));
    colorTable.append(qRgb(229, 133, 0));
    colorTable.append(qRgb(229, 136, 0));
    colorTable.append(qRgb(230, 138, 0));
    colorTable.append(qRgb(230, 140, 0));
    colorTable.append(qRgb(231, 142, 0));
    colorTable.append(qRgb(231, 145, 0));
    colorTable.append(qRgb(232, 147, 0));
    colorTable.append(qRgb(232, 149, 0));
    colorTable.append(qRgb(233, 152, 0));
    colorTable.append(qRgb(233, 154, 0));
    colorTable.append(qRgb(234, 156, 0));
    colorTable.append(qRgb(234, 159, 0));
    colorTable.append(qRgb(235, 161, 0));
    colorTable.append(qRgb(235, 163, 0));
    colorTable.append(qRgb(236, 166, 0));
    colorTable.append(qRgb(236, 168, 0));
    colorTable.append(qRgb(237, 170, 0));
    colorTable.append(qRgb(237, 172, 0));
    colorTable.append(qRgb(238, 175, 0));
    colorTable.append(qRgb(238, 177, 0));
    colorTable.append(qRgb(239, 179, 0));
    colorTable.append(qRgb(239, 182, 0));
    colorTable.append(qRgb(240, 184, 0));
    colorTable.append(qRgb(240, 186, 0));
    colorTable.append(qRgb(241, 189, 0));
    colorTable.append(qRgb(241, 191, 0));
    colorTable.append(qRgb(242, 193, 0));
    colorTable.append(qRgb(242, 196, 0));
    colorTable.append(qRgb(243, 198, 0));
    colorTable.append(qRgb(243, 200, 0));
    colorTable.append(qRgb(244, 202, 0));
    colorTable.append(qRgb(244, 205, 0));
    colorTable.append(qRgb(245, 207, 0));
    colorTable.append(qRgb(245, 209, 0));
    colorTable.append(qRgb(246, 212, 0));
    colorTable.append(qRgb(246, 214, 0));
    colorTable.append(qRgb(247, 216, 0));
    colorTable.append(qRgb(247, 219, 0));
    colorTable.append(qRgb(248, 221, 0));
    colorTable.append(qRgb(248, 223, 0));
    colorTable.append(qRgb(249, 226, 0));
    colorTable.append(qRgb(249, 228, 0));
    colorTable.append(qRgb(250, 230, 0));
    colorTable.append(qRgb(250, 232, 0));
    colorTable.append(qRgb(251, 235, 0));
    colorTable.append(qRgb(251, 237, 0));
    colorTable.append(qRgb(252, 239, 0));
    colorTable.append(qRgb(252, 242, 0));
    colorTable.append(qRgb(253, 244, 0));
    colorTable.append(qRgb(253, 246, 0));
    colorTable.append(qRgb(254, 249, 0));
    colorTable.append(qRgb(254, 251, 0));
    colorTable.append(qRgb(255, 253, 0));
    colorTable.append(qRgb(255, 255, 0));


    this->addToColorLookUpTables(colorTable,"Green Red Blue White");
    colorTable.clear();
    colorTable.append(qRgb(0, 0, 0));
    colorTable.append(qRgb(1, 0, 0));
    colorTable.append(qRgb(2, 0, 0));
    colorTable.append(qRgb(4, 0, 0));
    colorTable.append(qRgb(5, 0, 0));
    colorTable.append(qRgb(7, 0, 0));
    colorTable.append(qRgb(8, 0, 0));
    colorTable.append(qRgb(10, 0, 0));
    colorTable.append(qRgb(11, 0, 0));
    colorTable.append(qRgb(13, 0, 0));
    colorTable.append(qRgb(14, 0, 0));
    colorTable.append(qRgb(15, 0, 0));
    colorTable.append(qRgb(17, 0, 0));
    colorTable.append(qRgb(18, 0, 0));
    colorTable.append(qRgb(20, 0, 0));
    colorTable.append(qRgb(21, 0, 0));
    colorTable.append(qRgb(23, 0, 0));
    colorTable.append(qRgb(24, 0, 0));
    colorTable.append(qRgb(26, 0, 0));
    colorTable.append(qRgb(27, 0, 0));
    colorTable.append(qRgb(28, 0, 0));
    colorTable.append(qRgb(30, 0, 0));
    colorTable.append(qRgb(31, 0, 0));
    colorTable.append(qRgb(33, 0, 0));
    colorTable.append(qRgb(34, 0, 0));
    colorTable.append(qRgb(36, 0, 0));
    colorTable.append(qRgb(37, 0, 0));
    colorTable.append(qRgb(39, 0, 0));
    colorTable.append(qRgb(40, 0, 0));
    colorTable.append(qRgb(42, 0, 0));
    colorTable.append(qRgb(43, 0, 0));
    colorTable.append(qRgb(44, 0, 0));
    colorTable.append(qRgb(46, 0, 0));
    colorTable.append(qRgb(47, 0, 0));
    colorTable.append(qRgb(49, 0, 0));
    colorTable.append(qRgb(50, 0, 0));
    colorTable.append(qRgb(52, 0, 0));
    colorTable.append(qRgb(53, 0, 0));
    colorTable.append(qRgb(55, 0, 0));
    colorTable.append(qRgb(56, 0, 0));
    colorTable.append(qRgb(57, 0, 0));
    colorTable.append(qRgb(59, 0, 0));
    colorTable.append(qRgb(60, 0, 0));
    colorTable.append(qRgb(62, 0, 0));
    colorTable.append(qRgb(63, 0, 0));
    colorTable.append(qRgb(65, 0, 0));
    colorTable.append(qRgb(66, 0, 0));
    colorTable.append(qRgb(68, 0, 0));
    colorTable.append(qRgb(69, 0, 0));
    colorTable.append(qRgb(70, 0, 0));
    colorTable.append(qRgb(72, 0, 0));
    colorTable.append(qRgb(73, 0, 0));
    colorTable.append(qRgb(75, 0, 0));
    colorTable.append(qRgb(76, 0, 0));
    colorTable.append(qRgb(78, 0, 0));
    colorTable.append(qRgb(79, 0, 0));
    colorTable.append(qRgb(81, 0, 0));
    colorTable.append(qRgb(82, 0, 0));
    colorTable.append(qRgb(84, 0, 0));
    colorTable.append(qRgb(85, 0, 0));
    colorTable.append(qRgb(86, 0, 0));
    colorTable.append(qRgb(88, 0, 0));
    colorTable.append(qRgb(89, 0, 0));
    colorTable.append(qRgb(91, 0, 0));
    colorTable.append(qRgb(92, 0, 0));
    colorTable.append(qRgb(94, 0, 0));
    colorTable.append(qRgb(95, 0, 0));
    colorTable.append(qRgb(97, 0, 0));
    colorTable.append(qRgb(98, 0, 0));
    colorTable.append(qRgb(99, 0, 0));
    colorTable.append(qRgb(101, 0, 0));
    colorTable.append(qRgb(102, 0, 0));
    colorTable.append(qRgb(104, 0, 0));
    colorTable.append(qRgb(105, 0, 0));
    colorTable.append(qRgb(107, 0, 0));
    colorTable.append(qRgb(108, 0, 0));
    colorTable.append(qRgb(110, 0, 0));
    colorTable.append(qRgb(111, 0, 0));
    colorTable.append(qRgb(113, 0, 0));
    colorTable.append(qRgb(114, 0, 0));
    colorTable.append(qRgb(115, 0, 0));
    colorTable.append(qRgb(117, 0, 0));
    colorTable.append(qRgb(118, 0, 0));
    colorTable.append(qRgb(120, 0, 0));
    colorTable.append(qRgb(121, 0, 0));
    colorTable.append(qRgb(123, 0, 0));
    colorTable.append(qRgb(124, 0, 0));
    colorTable.append(qRgb(126, 0, 0));
    colorTable.append(qRgb(127, 0, 0));
    colorTable.append(qRgb(128, 0, 0));
    colorTable.append(qRgb(130, 0, 0));
    colorTable.append(qRgb(131, 0, 0));
    colorTable.append(qRgb(133, 0, 0));
    colorTable.append(qRgb(134, 0, 0));
    colorTable.append(qRgb(136, 0, 0));
    colorTable.append(qRgb(137, 0, 0));
    colorTable.append(qRgb(139, 0, 0));
    colorTable.append(qRgb(140, 0, 0));
    colorTable.append(qRgb(141, 0, 0));
    colorTable.append(qRgb(143, 0, 0));
    colorTable.append(qRgb(144, 0, 0));
    colorTable.append(qRgb(146, 0, 0));
    colorTable.append(qRgb(147, 0, 0));
    colorTable.append(qRgb(149, 0, 0));
    colorTable.append(qRgb(150, 0, 0));
    colorTable.append(qRgb(152, 0, 0));
    colorTable.append(qRgb(153, 0, 0));
    colorTable.append(qRgb(155, 0, 0));
    colorTable.append(qRgb(156, 0, 0));
    colorTable.append(qRgb(157, 0, 0));
    colorTable.append(qRgb(159, 0, 0));
    colorTable.append(qRgb(160, 0, 0));
    colorTable.append(qRgb(162, 0, 0));
    colorTable.append(qRgb(163, 0, 0));
    colorTable.append(qRgb(165, 0, 0));
    colorTable.append(qRgb(166, 0, 0));
    colorTable.append(qRgb(168, 0, 0));
    colorTable.append(qRgb(169, 0, 0));
    colorTable.append(qRgb(170, 0, 0));
    colorTable.append(qRgb(172, 0, 0));
    colorTable.append(qRgb(173, 0, 0));
    colorTable.append(qRgb(175, 1, 0));
    colorTable.append(qRgb(176, 3, 0));
    colorTable.append(qRgb(178, 5, 0));
    colorTable.append(qRgb(179, 7, 0));
    colorTable.append(qRgb(181, 9, 0));
    colorTable.append(qRgb(182, 11, 0));
    colorTable.append(qRgb(184, 13, 0));
    colorTable.append(qRgb(185, 15, 0));
    colorTable.append(qRgb(186, 17, 0));
    colorTable.append(qRgb(188, 18, 0));
    colorTable.append(qRgb(189, 20, 0));
    colorTable.append(qRgb(191, 22, 0));
    colorTable.append(qRgb(192, 24, 0));
    colorTable.append(qRgb(194, 26, 0));
    colorTable.append(qRgb(195, 28, 0));
    colorTable.append(qRgb(197, 30, 0));
    colorTable.append(qRgb(198, 32, 0));
    colorTable.append(qRgb(199, 34, 0));
    colorTable.append(qRgb(201, 35, 0));
    colorTable.append(qRgb(202, 37, 0));
    colorTable.append(qRgb(204, 39, 0));
    colorTable.append(qRgb(205, 41, 0));
    colorTable.append(qRgb(207, 43, 0));
    colorTable.append(qRgb(208, 45, 0));
    colorTable.append(qRgb(210, 47, 0));
    colorTable.append(qRgb(211, 49, 0));
    colorTable.append(qRgb(212, 51, 0));
    colorTable.append(qRgb(214, 52, 0));
    colorTable.append(qRgb(215, 54, 0));
    colorTable.append(qRgb(217, 56, 0));
    colorTable.append(qRgb(218, 58, 0));
    colorTable.append(qRgb(220, 60, 0));
    colorTable.append(qRgb(221, 62, 0));
    colorTable.append(qRgb(223, 64, 0));
    colorTable.append(qRgb(224, 66, 0));
    colorTable.append(qRgb(226, 68, 0));
    colorTable.append(qRgb(227, 69, 0));
    colorTable.append(qRgb(228, 71, 0));
    colorTable.append(qRgb(230, 73, 0));
    colorTable.append(qRgb(231, 75, 0));
    colorTable.append(qRgb(233, 77, 0));
    colorTable.append(qRgb(234, 79, 0));
    colorTable.append(qRgb(236, 81, 0));
    colorTable.append(qRgb(237, 83, 0));
    colorTable.append(qRgb(239, 85, 0));
    colorTable.append(qRgb(240, 86, 0));
    colorTable.append(qRgb(241, 88, 0));
    colorTable.append(qRgb(243, 90, 0));
    colorTable.append(qRgb(244, 92, 0));
    colorTable.append(qRgb(246, 94, 0));
    colorTable.append(qRgb(247, 96, 0));
    colorTable.append(qRgb(249, 98, 0));
    colorTable.append(qRgb(250, 100, 0));
    colorTable.append(qRgb(252, 102, 0));
    colorTable.append(qRgb(253, 103, 0));
    colorTable.append(qRgb(255, 105, 0));
    colorTable.append(qRgb(255, 107, 0));
    colorTable.append(qRgb(255, 109, 0));
    colorTable.append(qRgb(255, 111, 0));
    colorTable.append(qRgb(255, 113, 0));
    colorTable.append(qRgb(255, 115, 0));
    colorTable.append(qRgb(255, 117, 0));
    colorTable.append(qRgb(255, 119, 0));
    colorTable.append(qRgb(255, 120, 0));
    colorTable.append(qRgb(255, 122, 0));
    colorTable.append(qRgb(255, 124, 0));
    colorTable.append(qRgb(255, 126, 0));
    colorTable.append(qRgb(255, 128, 0));
    colorTable.append(qRgb(255, 130, 0));
    colorTable.append(qRgb(255, 132, 0));
    colorTable.append(qRgb(255, 134, 3));
    colorTable.append(qRgb(255, 136, 7));
    colorTable.append(qRgb(255, 137, 11));
    colorTable.append(qRgb(255, 139, 15));
    colorTable.append(qRgb(255, 141, 19));
    colorTable.append(qRgb(255, 143, 23));
    colorTable.append(qRgb(255, 145, 27));
    colorTable.append(qRgb(255, 147, 31));
    colorTable.append(qRgb(255, 149, 35));
    colorTable.append(qRgb(255, 151, 39));
    colorTable.append(qRgb(255, 153, 43));
    colorTable.append(qRgb(255, 154, 47));
    colorTable.append(qRgb(255, 156, 51));
    colorTable.append(qRgb(255, 158, 54));
    colorTable.append(qRgb(255, 160, 58));
    colorTable.append(qRgb(255, 162, 62));
    colorTable.append(qRgb(255, 164, 66));
    colorTable.append(qRgb(255, 166, 70));
    colorTable.append(qRgb(255, 168, 74));
    colorTable.append(qRgb(255, 170, 78));
    colorTable.append(qRgb(255, 171, 82));
    colorTable.append(qRgb(255, 173, 86));
    colorTable.append(qRgb(255, 175, 90));
    colorTable.append(qRgb(255, 177, 94));
    colorTable.append(qRgb(255, 179, 98));
    colorTable.append(qRgb(255, 181, 102));
    colorTable.append(qRgb(255, 183, 105));
    colorTable.append(qRgb(255, 185, 109));
    colorTable.append(qRgb(255, 187, 113));
    colorTable.append(qRgb(255, 188, 117));
    colorTable.append(qRgb(255, 190, 121));
    colorTable.append(qRgb(255, 192, 125));
    colorTable.append(qRgb(255, 194, 129));
    colorTable.append(qRgb(255, 196, 133));
    colorTable.append(qRgb(255, 198, 137));
    colorTable.append(qRgb(255, 200, 141));
    colorTable.append(qRgb(255, 202, 145));
    colorTable.append(qRgb(255, 204, 149));
    colorTable.append(qRgb(255, 205, 153));
    colorTable.append(qRgb(255, 207, 156));
    colorTable.append(qRgb(255, 209, 160));
    colorTable.append(qRgb(255, 211, 164));
    colorTable.append(qRgb(255, 213, 168));
    colorTable.append(qRgb(255, 215, 172));
    colorTable.append(qRgb(255, 217, 176));
    colorTable.append(qRgb(255, 219, 180));
    colorTable.append(qRgb(255, 221, 184));
    colorTable.append(qRgb(255, 222, 188));
    colorTable.append(qRgb(255, 224, 192));
    colorTable.append(qRgb(255, 226, 196));
    colorTable.append(qRgb(255, 228, 200));
    colorTable.append(qRgb(255, 230, 204));
    colorTable.append(qRgb(255, 232, 207));
    colorTable.append(qRgb(255, 234, 211));
    colorTable.append(qRgb(255, 236, 215));
    colorTable.append(qRgb(255, 238, 219));
    colorTable.append(qRgb(255, 239, 223));
    colorTable.append(qRgb(255, 241, 227));
    colorTable.append(qRgb(255, 243, 231));
    colorTable.append(qRgb(255, 245, 235));
    colorTable.append(qRgb(255, 247, 239));
    colorTable.append(qRgb(255, 249, 243));
    colorTable.append(qRgb(255, 251, 247));
    colorTable.append(qRgb(255, 253, 251));
    colorTable.append(qRgb(255, 255, 255));


    this->addToColorLookUpTables(colorTable,"Red Temperature");
    colorTable.clear();

    colorTable.append(qRgb(0, 0, 0));
    colorTable.append(qRgb(18, 1, 1));
    colorTable.append(qRgb(36, 2, 3));
    colorTable.append(qRgb(54, 3, 5));
    colorTable.append(qRgb(72, 4, 7));
    colorTable.append(qRgb(90, 5, 9));
    colorTable.append(qRgb(108, 6, 11));
    colorTable.append(qRgb(127, 7, 13));
    colorTable.append(qRgb(145, 8, 15));
    colorTable.append(qRgb(163, 9, 17));
    colorTable.append(qRgb(181, 10, 19));
    colorTable.append(qRgb(199, 11, 21));
    colorTable.append(qRgb(217, 12, 23));
    colorTable.append(qRgb(235, 13, 25));
    colorTable.append(qRgb(254, 14, 27));
    colorTable.append(qRgb(249, 15, 29));
    colorTable.append(qRgb(244, 16, 31));
    colorTable.append(qRgb(239, 17, 33));
    colorTable.append(qRgb(234, 18, 35));
    colorTable.append(qRgb(229, 19, 37));
    colorTable.append(qRgb(223, 20, 39));
    colorTable.append(qRgb(218, 21, 41));
    colorTable.append(qRgb(213, 22, 43));
    colorTable.append(qRgb(208, 23, 45));
    colorTable.append(qRgb(203, 24, 47));
    colorTable.append(qRgb(197, 25, 49));
    colorTable.append(qRgb(192, 26, 51));
    colorTable.append(qRgb(187, 27, 53));
    colorTable.append(qRgb(182, 28, 55));
    colorTable.append(qRgb(177, 29, 57));
    colorTable.append(qRgb(172, 30, 59));
    colorTable.append(qRgb(166, 31, 61));
    colorTable.append(qRgb(161, 32, 63));
    colorTable.append(qRgb(156, 33, 65));
    colorTable.append(qRgb(151, 34, 67));
    colorTable.append(qRgb(146, 35, 69));
    colorTable.append(qRgb(140, 36, 71));
    colorTable.append(qRgb(135, 37, 73));
    colorTable.append(qRgb(130, 38, 75));
    colorTable.append(qRgb(125, 39, 77));
    colorTable.append(qRgb(120, 40, 79));
    colorTable.append(qRgb(115, 41, 81));
    colorTable.append(qRgb(109, 42, 83));
    colorTable.append(qRgb(104, 43, 85));
    colorTable.append(qRgb(99, 44, 87));
    colorTable.append(qRgb(94, 45, 89));
    colorTable.append(qRgb(89, 46, 91));
    colorTable.append(qRgb(83, 47, 93));
    colorTable.append(qRgb(78, 48, 95));
    colorTable.append(qRgb(73, 49, 97));
    colorTable.append(qRgb(68, 50, 99));
    colorTable.append(qRgb(63, 51, 101));
    colorTable.append(qRgb(58, 52, 103));
    colorTable.append(qRgb(52, 53, 105));
    colorTable.append(qRgb(47, 54, 107));
    colorTable.append(qRgb(42, 55, 109));
    colorTable.append(qRgb(37, 56, 111));
    colorTable.append(qRgb(32, 57, 113));
    colorTable.append(qRgb(26, 58, 115));
    colorTable.append(qRgb(21, 59, 117));
    colorTable.append(qRgb(16, 60, 119));
    colorTable.append(qRgb(11, 61, 121));
    colorTable.append(qRgb(6, 62, 123));
    colorTable.append(qRgb(0, 63, 125));
    colorTable.append(qRgb(64, 64, 127));
    colorTable.append(qRgb(65, 65, 129));
    colorTable.append(qRgb(66, 66, 131));
    colorTable.append(qRgb(67, 67, 133));
    colorTable.append(qRgb(68, 68, 135));
    colorTable.append(qRgb(69, 69, 137));
    colorTable.append(qRgb(70, 70, 139));
    colorTable.append(qRgb(71, 71, 141));
    colorTable.append(qRgb(72, 72, 143));
    colorTable.append(qRgb(73, 73, 145));
    colorTable.append(qRgb(74, 74, 147));
    colorTable.append(qRgb(75, 75, 149));
    colorTable.append(qRgb(76, 76, 151));
    colorTable.append(qRgb(77, 77, 153));
    colorTable.append(qRgb(78, 78, 155));
    colorTable.append(qRgb(79, 79, 157));
    colorTable.append(qRgb(80, 80, 159));
    colorTable.append(qRgb(81, 81, 161));
    colorTable.append(qRgb(82, 82, 163));
    colorTable.append(qRgb(83, 83, 165));
    colorTable.append(qRgb(84, 84, 167));
    colorTable.append(qRgb(85, 85, 169));
    colorTable.append(qRgb(86, 86, 171));
    colorTable.append(qRgb(87, 87, 173));
    colorTable.append(qRgb(88, 88, 175));
    colorTable.append(qRgb(89, 89, 177));
    colorTable.append(qRgb(90, 90, 179));
    colorTable.append(qRgb(91, 91, 181));
    colorTable.append(qRgb(92, 92, 183));
    colorTable.append(qRgb(93, 93, 185));
    colorTable.append(qRgb(94, 94, 187));
    colorTable.append(qRgb(95, 95, 189));
    colorTable.append(qRgb(96, 96, 191));
    colorTable.append(qRgb(97, 97, 193));
    colorTable.append(qRgb(98, 98, 195));
    colorTable.append(qRgb(99, 99, 197));
    colorTable.append(qRgb(100, 100, 199));
    colorTable.append(qRgb(101, 101, 201));
    colorTable.append(qRgb(102, 102, 203));
    colorTable.append(qRgb(103, 103, 205));
    colorTable.append(qRgb(104, 104, 207));
    colorTable.append(qRgb(105, 105, 209));
    colorTable.append(qRgb(106, 106, 211));
    colorTable.append(qRgb(107, 107, 213));
    colorTable.append(qRgb(108, 108, 215));
    colorTable.append(qRgb(109, 109, 217));
    colorTable.append(qRgb(110, 110, 219));
    colorTable.append(qRgb(111, 111, 221));
    colorTable.append(qRgb(112, 112, 223));
    colorTable.append(qRgb(113, 113, 225));
    colorTable.append(qRgb(114, 114, 227));
    colorTable.append(qRgb(115, 115, 229));
    colorTable.append(qRgb(116, 116, 231));
    colorTable.append(qRgb(117, 117, 233));
    colorTable.append(qRgb(118, 118, 235));
    colorTable.append(qRgb(119, 119, 237));
    colorTable.append(qRgb(120, 120, 239));
    colorTable.append(qRgb(121, 121, 241));
    colorTable.append(qRgb(122, 122, 243));
    colorTable.append(qRgb(123, 123, 245));
    colorTable.append(qRgb(124, 124, 247));
    colorTable.append(qRgb(125, 125, 249));
    colorTable.append(qRgb(126, 126, 251));
    colorTable.append(qRgb(127, 127, 253));
    colorTable.append(qRgb(128, 128, 255));
    colorTable.append(qRgb(129, 129, 251));
    colorTable.append(qRgb(130, 130, 247));
    colorTable.append(qRgb(131, 131, 243));
    colorTable.append(qRgb(132, 132, 238));
    colorTable.append(qRgb(133, 133, 234));
    colorTable.append(qRgb(134, 134, 230));
    colorTable.append(qRgb(135, 135, 226));
    colorTable.append(qRgb(136, 136, 221));
    colorTable.append(qRgb(137, 137, 217));
    colorTable.append(qRgb(138, 138, 213));
    colorTable.append(qRgb(139, 139, 209));
    colorTable.append(qRgb(140, 140, 204));
    colorTable.append(qRgb(141, 141, 200));
    colorTable.append(qRgb(142, 142, 196));
    colorTable.append(qRgb(143, 143, 192));
    colorTable.append(qRgb(144, 144, 187));
    colorTable.append(qRgb(145, 145, 183));
    colorTable.append(qRgb(146, 146, 179));
    colorTable.append(qRgb(147, 147, 175));
    colorTable.append(qRgb(148, 148, 170));
    colorTable.append(qRgb(149, 149, 166));
    colorTable.append(qRgb(150, 150, 162));
    colorTable.append(qRgb(151, 151, 158));
    colorTable.append(qRgb(152, 152, 153));
    colorTable.append(qRgb(153, 153, 149));
    colorTable.append(qRgb(154, 154, 145));
    colorTable.append(qRgb(155, 155, 141));
    colorTable.append(qRgb(156, 156, 136));
    colorTable.append(qRgb(157, 157, 132));
    colorTable.append(qRgb(158, 158, 128));
    colorTable.append(qRgb(159, 159, 124));
    colorTable.append(qRgb(160, 160, 119));
    colorTable.append(qRgb(161, 161, 115));
    colorTable.append(qRgb(162, 162, 111));
    colorTable.append(qRgb(163, 163, 107));
    colorTable.append(qRgb(164, 164, 102));
    colorTable.append(qRgb(165, 165, 98));
    colorTable.append(qRgb(166, 166, 94));
    colorTable.append(qRgb(167, 167, 90));
    colorTable.append(qRgb(168, 168, 85));
    colorTable.append(qRgb(169, 169, 81));
    colorTable.append(qRgb(170, 170, 77));
    colorTable.append(qRgb(171, 171, 73));
    colorTable.append(qRgb(172, 172, 68));
    colorTable.append(qRgb(173, 173, 64));
    colorTable.append(qRgb(174, 174, 60));
    colorTable.append(qRgb(175, 175, 56));
    colorTable.append(qRgb(176, 176, 51));
    colorTable.append(qRgb(177, 177, 47));
    colorTable.append(qRgb(178, 178, 43));
    colorTable.append(qRgb(179, 179, 39));
    colorTable.append(qRgb(180, 180, 34));
    colorTable.append(qRgb(181, 181, 30));
    colorTable.append(qRgb(182, 182, 26));
    colorTable.append(qRgb(183, 183, 22));
    colorTable.append(qRgb(184, 184, 17));
    colorTable.append(qRgb(185, 185, 13));
    colorTable.append(qRgb(186, 186, 9));
    colorTable.append(qRgb(187, 187, 5));
    colorTable.append(qRgb(188, 188, 0));
    colorTable.append(qRgb(189, 189, 3));
    colorTable.append(qRgb(190, 190, 7));
    colorTable.append(qRgb(191, 191, 11));
    colorTable.append(qRgb(192, 192, 15));
    colorTable.append(qRgb(193, 193, 19));
    colorTable.append(qRgb(194, 194, 22));
    colorTable.append(qRgb(195, 195, 26));
    colorTable.append(qRgb(196, 196, 30));
    colorTable.append(qRgb(197, 197, 34));
    colorTable.append(qRgb(198, 198, 38));
    colorTable.append(qRgb(199, 199, 41));
    colorTable.append(qRgb(200, 200, 45));
    colorTable.append(qRgb(201, 201, 49));
    colorTable.append(qRgb(202, 202, 53));
    colorTable.append(qRgb(203, 203, 57));
    colorTable.append(qRgb(204, 204, 60));
    colorTable.append(qRgb(205, 205, 64));
    colorTable.append(qRgb(206, 206, 68));
    colorTable.append(qRgb(207, 207, 72));
    colorTable.append(qRgb(208, 208, 76));
    colorTable.append(qRgb(209, 209, 79));
    colorTable.append(qRgb(210, 210, 83));
    colorTable.append(qRgb(211, 211, 87));
    colorTable.append(qRgb(212, 212, 91));
    colorTable.append(qRgb(213, 213, 95));
    colorTable.append(qRgb(214, 214, 98));
    colorTable.append(qRgb(215, 215, 102));
    colorTable.append(qRgb(216, 216, 106));
    colorTable.append(qRgb(217, 217, 110));
    colorTable.append(qRgb(218, 218, 114));
    colorTable.append(qRgb(219, 219, 117));
    colorTable.append(qRgb(220, 220, 121));
    colorTable.append(qRgb(221, 221, 125));
    colorTable.append(qRgb(222, 222, 129));
    colorTable.append(qRgb(223, 223, 133));
    colorTable.append(qRgb(224, 224, 137));
    colorTable.append(qRgb(225, 225, 140));
    colorTable.append(qRgb(226, 226, 144));
    colorTable.append(qRgb(227, 227, 148));
    colorTable.append(qRgb(228, 228, 152));
    colorTable.append(qRgb(229, 229, 156));
    colorTable.append(qRgb(230, 230, 159));
    colorTable.append(qRgb(231, 231, 163));
    colorTable.append(qRgb(232, 232, 167));
    colorTable.append(qRgb(233, 233, 171));
    colorTable.append(qRgb(234, 234, 175));
    colorTable.append(qRgb(235, 235, 178));
    colorTable.append(qRgb(236, 236, 182));
    colorTable.append(qRgb(237, 237, 186));
    colorTable.append(qRgb(238, 238, 190));
    colorTable.append(qRgb(239, 239, 194));
    colorTable.append(qRgb(240, 240, 197));
    colorTable.append(qRgb(241, 241, 201));
    colorTable.append(qRgb(242, 242, 205));
    colorTable.append(qRgb(243, 243, 209));
    colorTable.append(qRgb(244, 244, 213));
    colorTable.append(qRgb(245, 245, 216));
    colorTable.append(qRgb(246, 246, 220));
    colorTable.append(qRgb(247, 247, 224));
    colorTable.append(qRgb(248, 248, 228));
    colorTable.append(qRgb(249, 249, 232));
    colorTable.append(qRgb(250, 250, 235));
    colorTable.append(qRgb(251, 251, 239));
    colorTable.append(qRgb(252, 252, 243));
    colorTable.append(qRgb(253, 253, 247));
    colorTable.append(qRgb(254, 254, 251));
    colorTable.append(qRgb(255, 255, 255));


    this->addToColorLookUpTables(colorTable,"Stern special");
    colorTable.clear();


    colorTable.append(qRgb(0, 0, 0));
    colorTable.append(qRgb(4, 0, 3));
    colorTable.append(qRgb(9, 0, 7));
    colorTable.append(qRgb(13, 0, 10));
    colorTable.append(qRgb(18, 0, 14));
    colorTable.append(qRgb(22, 0, 19));
    colorTable.append(qRgb(27, 0, 23));
    colorTable.append(qRgb(31, 0, 28));
    colorTable.append(qRgb(36, 0, 32));
    colorTable.append(qRgb(40, 0, 38));
    colorTable.append(qRgb(45, 0, 43));
    colorTable.append(qRgb(50, 0, 48));
    colorTable.append(qRgb(54, 0, 53));
    colorTable.append(qRgb(58, 0, 59));
    colorTable.append(qRgb(61, 0, 63));
    colorTable.append(qRgb(64, 0, 68));
    colorTable.append(qRgb(68, 0, 72));
    colorTable.append(qRgb(69, 0, 77));
    colorTable.append(qRgb(72, 0, 81));
    colorTable.append(qRgb(74, 0, 86));
    colorTable.append(qRgb(77, 0, 91));
    colorTable.append(qRgb(79, 0, 95));
    colorTable.append(qRgb(80, 0, 100));
    colorTable.append(qRgb(82, 0, 104));
    colorTable.append(qRgb(83, 0, 109));
    colorTable.append(qRgb(85, 0, 113));
    colorTable.append(qRgb(84, 0, 118));
    colorTable.append(qRgb(86, 0, 122));
    colorTable.append(qRgb(87, 0, 127));
    colorTable.append(qRgb(88, 0, 132));
    colorTable.append(qRgb(86, 0, 136));
    colorTable.append(qRgb(87, 0, 141));
    colorTable.append(qRgb(87, 0, 145));
    colorTable.append(qRgb(87, 0, 150));
    colorTable.append(qRgb(85, 0, 154));
    colorTable.append(qRgb(84, 0, 159));
    colorTable.append(qRgb(84, 0, 163));
    colorTable.append(qRgb(84, 0, 168));
    colorTable.append(qRgb(83, 0, 173));
    colorTable.append(qRgb(79, 0, 177));
    colorTable.append(qRgb(78, 0, 182));
    colorTable.append(qRgb(77, 0, 186));
    colorTable.append(qRgb(76, 0, 191));
    colorTable.append(qRgb(71, 0, 195));
    colorTable.append(qRgb(70, 0, 200));
    colorTable.append(qRgb(68, 0, 204));
    colorTable.append(qRgb(66, 0, 209));
    colorTable.append(qRgb(60, 0, 214));
    colorTable.append(qRgb(58, 0, 218));
    colorTable.append(qRgb(55, 0, 223));
    colorTable.append(qRgb(53, 0, 227));
    colorTable.append(qRgb(46, 0, 232));
    colorTable.append(qRgb(43, 0, 236));
    colorTable.append(qRgb(40, 0, 241));
    colorTable.append(qRgb(36, 0, 245));
    colorTable.append(qRgb(33, 0, 250));
    colorTable.append(qRgb(25, 0, 255));
    colorTable.append(qRgb(21, 0, 255));
    colorTable.append(qRgb(16, 0, 255));
    colorTable.append(qRgb(12, 0, 255));
    colorTable.append(qRgb(4, 0, 255));
    colorTable.append(qRgb(0, 0, 255));
    colorTable.append(qRgb(0, 4, 255));
    colorTable.append(qRgb(0, 8, 255));
    colorTable.append(qRgb(0, 16, 255));
    colorTable.append(qRgb(0, 21, 255));
    colorTable.append(qRgb(0, 25, 255));
    colorTable.append(qRgb(0, 29, 255));
    colorTable.append(qRgb(0, 38, 255));
    colorTable.append(qRgb(0, 42, 255));
    colorTable.append(qRgb(0, 46, 255));
    colorTable.append(qRgb(0, 51, 255));
    colorTable.append(qRgb(0, 55, 255));
    colorTable.append(qRgb(0, 63, 255));
    colorTable.append(qRgb(0, 67, 255));
    colorTable.append(qRgb(0, 72, 255));
    colorTable.append(qRgb(0, 76, 255));
    colorTable.append(qRgb(0, 84, 255));
    colorTable.append(qRgb(0, 89, 255));
    colorTable.append(qRgb(0, 93, 255));
    colorTable.append(qRgb(0, 97, 255));
    colorTable.append(qRgb(0, 106, 255));
    colorTable.append(qRgb(0, 110, 255));
    colorTable.append(qRgb(0, 114, 255));
    colorTable.append(qRgb(0, 119, 255));
    colorTable.append(qRgb(0, 127, 255));
    colorTable.append(qRgb(0, 131, 255));
    colorTable.append(qRgb(0, 135, 255));
    colorTable.append(qRgb(0, 140, 255));
    colorTable.append(qRgb(0, 144, 255));
    colorTable.append(qRgb(0, 152, 255));
    colorTable.append(qRgb(0, 157, 255));
    colorTable.append(qRgb(0, 161, 255));
    colorTable.append(qRgb(0, 165, 255));
    colorTable.append(qRgb(0, 174, 255));
    colorTable.append(qRgb(0, 178, 255));
    colorTable.append(qRgb(0, 182, 255));
    colorTable.append(qRgb(0, 187, 255));
    colorTable.append(qRgb(0, 195, 255));
    colorTable.append(qRgb(0, 199, 255));
    colorTable.append(qRgb(0, 203, 255));
    colorTable.append(qRgb(0, 208, 255));
    colorTable.append(qRgb(0, 216, 255));
    colorTable.append(qRgb(0, 220, 255));
    colorTable.append(qRgb(0, 225, 255));
    colorTable.append(qRgb(0, 229, 255));
    colorTable.append(qRgb(0, 233, 255));
    colorTable.append(qRgb(0, 242, 255));
    colorTable.append(qRgb(0, 246, 255));
    colorTable.append(qRgb(0, 250, 255));
    colorTable.append(qRgb(0, 255, 255));
    colorTable.append(qRgb(0, 255, 246));
    colorTable.append(qRgb(0, 255, 242));
    colorTable.append(qRgb(0, 255, 238));
    colorTable.append(qRgb(0, 255, 233));
    colorTable.append(qRgb(0, 255, 225));
    colorTable.append(qRgb(0, 255, 220));
    colorTable.append(qRgb(0, 255, 216));
    colorTable.append(qRgb(0, 255, 212));
    colorTable.append(qRgb(0, 255, 203));
    colorTable.append(qRgb(0, 255, 199));
    colorTable.append(qRgb(0, 255, 195));
    colorTable.append(qRgb(0, 255, 191));
    colorTable.append(qRgb(0, 255, 187));
    colorTable.append(qRgb(0, 255, 178));
    colorTable.append(qRgb(0, 255, 174));
    colorTable.append(qRgb(0, 255, 170));
    colorTable.append(qRgb(0, 255, 165));
    colorTable.append(qRgb(0, 255, 157));
    colorTable.append(qRgb(0, 255, 152));
    colorTable.append(qRgb(0, 255, 148));
    colorTable.append(qRgb(0, 255, 144));
    colorTable.append(qRgb(0, 255, 135));
    colorTable.append(qRgb(0, 255, 131));
    colorTable.append(qRgb(0, 255, 127));
    colorTable.append(qRgb(0, 255, 123));
    colorTable.append(qRgb(0, 255, 114));
    colorTable.append(qRgb(0, 255, 110));
    colorTable.append(qRgb(0, 255, 106));
    colorTable.append(qRgb(0, 255, 102));
    colorTable.append(qRgb(0, 255, 97));
    colorTable.append(qRgb(0, 255, 89));
    colorTable.append(qRgb(0, 255, 84));
    colorTable.append(qRgb(0, 255, 80));
    colorTable.append(qRgb(0, 255, 76));
    colorTable.append(qRgb(0, 255, 67));
    colorTable.append(qRgb(0, 255, 63));
    colorTable.append(qRgb(0, 255, 59));
    colorTable.append(qRgb(0, 255, 55));
    colorTable.append(qRgb(0, 255, 46));
    colorTable.append(qRgb(0, 255, 42));
    colorTable.append(qRgb(0, 255, 38));
    colorTable.append(qRgb(0, 255, 34));
    colorTable.append(qRgb(0, 255, 25));
    colorTable.append(qRgb(0, 255, 21));
    colorTable.append(qRgb(0, 255, 16));
    colorTable.append(qRgb(0, 255, 12));
    colorTable.append(qRgb(0, 255, 8));
    colorTable.append(qRgb(0, 255, 0));
    colorTable.append(qRgb(4, 255, 0));
    colorTable.append(qRgb(8, 255, 0));
    colorTable.append(qRgb(12, 255, 0));
    colorTable.append(qRgb(21, 255, 0));
    colorTable.append(qRgb(25, 255, 0));
    colorTable.append(qRgb(29, 255, 0));
    colorTable.append(qRgb(33, 255, 0));
    colorTable.append(qRgb(42, 255, 0));
    colorTable.append(qRgb(46, 255, 0));
    colorTable.append(qRgb(51, 255, 0));
    colorTable.append(qRgb(55, 255, 0));
    colorTable.append(qRgb(63, 255, 0));
    colorTable.append(qRgb(67, 255, 0));
    colorTable.append(qRgb(72, 255, 0));
    colorTable.append(qRgb(76, 255, 0));
    colorTable.append(qRgb(80, 255, 0));
    colorTable.append(qRgb(89, 255, 0));
    colorTable.append(qRgb(93, 255, 0));
    colorTable.append(qRgb(97, 255, 0));
    colorTable.append(qRgb(101, 255, 0));
    colorTable.append(qRgb(110, 255, 0));
    colorTable.append(qRgb(114, 255, 0));
    colorTable.append(qRgb(119, 255, 0));
    colorTable.append(qRgb(123, 255, 0));
    colorTable.append(qRgb(131, 255, 0));
    colorTable.append(qRgb(135, 255, 0));
    colorTable.append(qRgb(140, 255, 0));
    colorTable.append(qRgb(144, 255, 0));
    colorTable.append(qRgb(153, 255, 0));
    colorTable.append(qRgb(157, 255, 0));
    colorTable.append(qRgb(161, 255, 0));
    colorTable.append(qRgb(165, 255, 0));
    colorTable.append(qRgb(169, 255, 0));
    colorTable.append(qRgb(178, 255, 0));
    colorTable.append(qRgb(182, 255, 0));
    colorTable.append(qRgb(187, 255, 0));
    colorTable.append(qRgb(191, 255, 0));
    colorTable.append(qRgb(199, 255, 0));
    colorTable.append(qRgb(203, 255, 0));
    colorTable.append(qRgb(208, 255, 0));
    colorTable.append(qRgb(212, 255, 0));
    colorTable.append(qRgb(221, 255, 0));
    colorTable.append(qRgb(225, 255, 0));
    colorTable.append(qRgb(229, 255, 0));
    colorTable.append(qRgb(233, 255, 0));
    colorTable.append(qRgb(242, 255, 0));
    colorTable.append(qRgb(246, 255, 0));
    colorTable.append(qRgb(250, 255, 0));
    colorTable.append(qRgb(255, 255, 0));
    colorTable.append(qRgb(255, 250, 0));
    colorTable.append(qRgb(255, 242, 0));
    colorTable.append(qRgb(255, 238, 0));
    colorTable.append(qRgb(255, 233, 0));
    colorTable.append(qRgb(255, 229, 0));
    colorTable.append(qRgb(255, 221, 0));
    colorTable.append(qRgb(255, 216, 0));
    colorTable.append(qRgb(255, 212, 0));
    colorTable.append(qRgb(255, 208, 0));
    colorTable.append(qRgb(255, 199, 0));
    colorTable.append(qRgb(255, 195, 0));
    colorTable.append(qRgb(255, 191, 0));
    colorTable.append(qRgb(255, 187, 0));
    colorTable.append(qRgb(255, 178, 0));
    colorTable.append(qRgb(255, 174, 0));
    colorTable.append(qRgb(255, 170, 0));
    colorTable.append(qRgb(255, 165, 0));
    colorTable.append(qRgb(255, 161, 0));
    colorTable.append(qRgb(255, 153, 0));
    colorTable.append(qRgb(255, 148, 0));
    colorTable.append(qRgb(255, 144, 0));
    colorTable.append(qRgb(255, 140, 0));
    colorTable.append(qRgb(255, 131, 0));
    colorTable.append(qRgb(255, 127, 0));
    colorTable.append(qRgb(255, 123, 0));
    colorTable.append(qRgb(255, 119, 0));
    colorTable.append(qRgb(255, 110, 0));
    colorTable.append(qRgb(255, 106, 0));
    colorTable.append(qRgb(255, 102, 0));
    colorTable.append(qRgb(255, 97, 0));
    colorTable.append(qRgb(255, 89, 0));
    colorTable.append(qRgb(255, 85, 0));
    colorTable.append(qRgb(255, 80, 0));
    colorTable.append(qRgb(255, 76, 0));
    colorTable.append(qRgb(255, 72, 0));
    colorTable.append(qRgb(255, 63, 0));
    colorTable.append(qRgb(255, 59, 0));
    colorTable.append(qRgb(255, 55, 0));
    colorTable.append(qRgb(255, 51, 0));
    colorTable.append(qRgb(255, 42, 0));
    colorTable.append(qRgb(255, 38, 0));
    colorTable.append(qRgb(255, 34, 0));
    colorTable.append(qRgb(255, 29, 0));
    colorTable.append(qRgb(255, 21, 0));
    colorTable.append(qRgb(255, 17, 0));
    colorTable.append(qRgb(255, 12, 0));
    colorTable.append(qRgb(255, 8, 0));
    colorTable.append(qRgb(255, 0, 0));

    this->addToColorLookUpTables(colorTable,"Rainbow ");
    colorTable.clear();

    colorTable.append(qRgb(0, 0, 0));
    colorTable.append(qRgb(0, 0, 5));
    colorTable.append(qRgb(0, 0, 10));
    colorTable.append(qRgb(0, 0, 15));
    colorTable.append(qRgb(0, 0, 20));
    colorTable.append(qRgb(0, 0, 26));
    colorTable.append(qRgb(0, 0, 31));
    colorTable.append(qRgb(0, 0, 36));
    colorTable.append(qRgb(0, 0, 41));
    colorTable.append(qRgb(0, 0, 46));
    colorTable.append(qRgb(0, 0, 52));
    colorTable.append(qRgb(0, 0, 57));
    colorTable.append(qRgb(0, 0, 62));
    colorTable.append(qRgb(0, 0, 67));
    colorTable.append(qRgb(0, 0, 72));
    colorTable.append(qRgb(0, 0, 78));
    colorTable.append(qRgb(0, 0, 83));
    colorTable.append(qRgb(0, 0, 88));
    colorTable.append(qRgb(0, 0, 93));
    colorTable.append(qRgb(0, 0, 98));
    colorTable.append(qRgb(0, 0, 104));
    colorTable.append(qRgb(0, 0, 109));
    colorTable.append(qRgb(0, 0, 114));
    colorTable.append(qRgb(0, 0, 119));
    colorTable.append(qRgb(0, 0, 124));
    colorTable.append(qRgb(0, 0, 130));
    colorTable.append(qRgb(0, 0, 135));
    colorTable.append(qRgb(0, 0, 140));
    colorTable.append(qRgb(0, 0, 145));
    colorTable.append(qRgb(0, 0, 150));
    colorTable.append(qRgb(0, 0, 156));
    colorTable.append(qRgb(0, 0, 161));
    colorTable.append(qRgb(0, 0, 166));
    colorTable.append(qRgb(0, 0, 171));
    colorTable.append(qRgb(0, 0, 176));
    colorTable.append(qRgb(0, 0, 182));
    colorTable.append(qRgb(0, 0, 187));
    colorTable.append(qRgb(0, 0, 192));
    colorTable.append(qRgb(0, 0, 197));
    colorTable.append(qRgb(0, 0, 202));
    colorTable.append(qRgb(0, 0, 208));
    colorTable.append(qRgb(0, 0, 213));
    colorTable.append(qRgb(0, 0, 218));
    colorTable.append(qRgb(0, 0, 223));
    colorTable.append(qRgb(0, 0, 228));
    colorTable.append(qRgb(0, 0, 234));
    colorTable.append(qRgb(0, 0, 239));
    colorTable.append(qRgb(0, 0, 244));
    colorTable.append(qRgb(4, 0, 249));
    colorTable.append(qRgb(9, 0, 255));
    colorTable.append(qRgb(14, 0, 250));
    colorTable.append(qRgb(19, 0, 245));
    colorTable.append(qRgb(23, 0, 239));
    colorTable.append(qRgb(28, 0, 234));
    colorTable.append(qRgb(33, 0, 228));
    colorTable.append(qRgb(38, 0, 223));
    colorTable.append(qRgb(42, 0, 218));
    colorTable.append(qRgb(47, 0, 212));
    colorTable.append(qRgb(52, 0, 207));
    colorTable.append(qRgb(57, 0, 201));
    colorTable.append(qRgb(61, 0, 196));
    colorTable.append(qRgb(66, 0, 190));
    colorTable.append(qRgb(71, 0, 185));
    colorTable.append(qRgb(76, 0, 180));
    colorTable.append(qRgb(81, 0, 174));
    colorTable.append(qRgb(81, 0, 169));
    colorTable.append(qRgb(81, 0, 163));
    colorTable.append(qRgb(81, 0, 158));
    colorTable.append(qRgb(81, 0, 152));
    colorTable.append(qRgb(81, 0, 147));
    colorTable.append(qRgb(81, 0, 142));
    colorTable.append(qRgb(81, 0, 136));
    colorTable.append(qRgb(80, 0, 131));
    colorTable.append(qRgb(80, 0, 125));
    colorTable.append(qRgb(80, 0, 120));
    colorTable.append(qRgb(80, 0, 114));
    colorTable.append(qRgb(80, 0, 109));
    colorTable.append(qRgb(80, 0, 104));
    colorTable.append(qRgb(80, 0, 98));
    colorTable.append(qRgb(79, 0, 93));
    colorTable.append(qRgb(84, 0, 87));
    colorTable.append(qRgb(89, 0, 82));
    colorTable.append(qRgb(94, 0, 76));
    colorTable.append(qRgb(99, 0, 71));
    colorTable.append(qRgb(104, 0, 66));
    colorTable.append(qRgb(109, 0, 60));
    colorTable.append(qRgb(114, 0, 55));
    colorTable.append(qRgb(119, 0, 49));
    colorTable.append(qRgb(124, 0, 44));
    colorTable.append(qRgb(129, 0, 38));
    colorTable.append(qRgb(134, 0, 33));
    colorTable.append(qRgb(139, 0, 28));
    colorTable.append(qRgb(144, 0, 22));
    colorTable.append(qRgb(149, 0, 17));
    colorTable.append(qRgb(154, 0, 11));
    colorTable.append(qRgb(159, 0, 6));
    colorTable.append(qRgb(164, 0, 0));
    colorTable.append(qRgb(169, 0, 0));
    colorTable.append(qRgb(174, 0, 0));
    colorTable.append(qRgb(180, 0, 0));
    colorTable.append(qRgb(185, 0, 0));
    colorTable.append(qRgb(190, 0, 0));
    colorTable.append(qRgb(196, 0, 0));
    colorTable.append(qRgb(201, 0, 0));
    colorTable.append(qRgb(206, 0, 0));
    colorTable.append(qRgb(212, 0, 0));
    colorTable.append(qRgb(217, 0, 0));
    colorTable.append(qRgb(222, 0, 0));
    colorTable.append(qRgb(228, 0, 0));
    colorTable.append(qRgb(233, 0, 0));
    colorTable.append(qRgb(255, 0, 0));
    colorTable.append(qRgb(255, 0, 0));
    colorTable.append(qRgb(255, 0, 0));
    colorTable.append(qRgb(255, 0, 0));
    colorTable.append(qRgb(255, 5, 0));
    colorTable.append(qRgb(255, 10, 0));
    colorTable.append(qRgb(255, 16, 0));
    colorTable.append(qRgb(255, 21, 0));
    colorTable.append(qRgb(255, 27, 0));
    colorTable.append(qRgb(255, 32, 0));
    colorTable.append(qRgb(255, 37, 0));
    colorTable.append(qRgb(255, 43, 0));
    colorTable.append(qRgb(255, 48, 0));
    colorTable.append(qRgb(255, 54, 0));
    colorTable.append(qRgb(255, 59, 0));
    colorTable.append(qRgb(255, 64, 0));
    colorTable.append(qRgb(255, 70, 0));
    colorTable.append(qRgb(255, 75, 0));
    colorTable.append(qRgb(255, 81, 0));
    colorTable.append(qRgb(255, 85, 4));
    colorTable.append(qRgb(255, 90, 9));
    colorTable.append(qRgb(255, 95, 14));
    colorTable.append(qRgb(255, 100, 19));
    colorTable.append(qRgb(255, 105, 24));
    colorTable.append(qRgb(255, 109, 28));
    colorTable.append(qRgb(255, 114, 33));
    colorTable.append(qRgb(255, 119, 38));
    colorTable.append(qRgb(255, 124, 43));
    colorTable.append(qRgb(255, 129, 48));
    colorTable.append(qRgb(255, 134, 53));
    colorTable.append(qRgb(255, 138, 57));
    colorTable.append(qRgb(255, 143, 62));
    colorTable.append(qRgb(255, 148, 67));
    colorTable.append(qRgb(255, 153, 72));
    colorTable.append(qRgb(255, 158, 77));
    colorTable.append(qRgb(255, 163, 82));
    colorTable.append(qRgb(255, 163, 77));
    colorTable.append(qRgb(255, 163, 71));
    colorTable.append(qRgb(255, 163, 65));
    colorTable.append(qRgb(255, 163, 59));
    colorTable.append(qRgb(255, 163, 53));
    colorTable.append(qRgb(255, 163, 47));
    colorTable.append(qRgb(255, 163, 41));
    colorTable.append(qRgb(255, 163, 36));
    colorTable.append(qRgb(255, 163, 30));
    colorTable.append(qRgb(255, 163, 24));
    colorTable.append(qRgb(255, 163, 18));
    colorTable.append(qRgb(255, 163, 12));
    colorTable.append(qRgb(255, 163, 6));
    colorTable.append(qRgb(255, 163, 0));
    colorTable.append(qRgb(255, 163, 0));
    colorTable.append(qRgb(255, 163, 0));
    colorTable.append(qRgb(255, 163, 0));
    colorTable.append(qRgb(248, 163, 0));
    colorTable.append(qRgb(240, 163, 0));
    colorTable.append(qRgb(232, 163, 0));
    colorTable.append(qRgb(225, 163, 0));
    colorTable.append(qRgb(217, 163, 0));
    colorTable.append(qRgb(209, 163, 0));
    colorTable.append(qRgb(202, 163, 0));
    colorTable.append(qRgb(194, 163, 0));
    colorTable.append(qRgb(186, 163, 0));
    colorTable.append(qRgb(179, 163, 0));
    colorTable.append(qRgb(171, 163, 0));
    colorTable.append(qRgb(163, 163, 0));
    colorTable.append(qRgb(168, 163, 0));
    colorTable.append(qRgb(173, 163, 0));
    colorTable.append(qRgb(178, 169, 3));
    colorTable.append(qRgb(183, 175, 6));
    colorTable.append(qRgb(188, 181, 9));
    colorTable.append(qRgb(193, 187, 12));
    colorTable.append(qRgb(198, 193, 16));
    colorTable.append(qRgb(203, 199, 19));
    colorTable.append(qRgb(209, 205, 22));
    colorTable.append(qRgb(214, 212, 25));
    colorTable.append(qRgb(219, 218, 29));
    colorTable.append(qRgb(224, 224, 32));
    colorTable.append(qRgb(229, 230, 35));
    colorTable.append(qRgb(234, 236, 38));
    colorTable.append(qRgb(239, 242, 41));
    colorTable.append(qRgb(244, 248, 45));
    colorTable.append(qRgb(249, 255, 48));
    colorTable.append(qRgb(255, 255, 51));
    colorTable.append(qRgb(255, 255, 54));
    colorTable.append(qRgb(255, 255, 58));
    colorTable.append(qRgb(255, 255, 61));
    colorTable.append(qRgb(255, 255, 64));
    colorTable.append(qRgb(255, 255, 67));
    colorTable.append(qRgb(255, 255, 71));
    colorTable.append(qRgb(255, 255, 74));
    colorTable.append(qRgb(255, 255, 77));
    colorTable.append(qRgb(255, 255, 80));
    colorTable.append(qRgb(255, 255, 83));
    colorTable.append(qRgb(255, 255, 87));
    colorTable.append(qRgb(255, 255, 90));
    colorTable.append(qRgb(255, 255, 93));
    colorTable.append(qRgb(255, 255, 96));
    colorTable.append(qRgb(255, 255, 100));
    colorTable.append(qRgb(255, 255, 103));
    colorTable.append(qRgb(255, 255, 106));
    colorTable.append(qRgb(255, 255, 109));
    colorTable.append(qRgb(255, 255, 112));
    colorTable.append(qRgb(255, 255, 116));
    colorTable.append(qRgb(255, 255, 119));
    colorTable.append(qRgb(255, 255, 122));
    colorTable.append(qRgb(255, 255, 125));
    colorTable.append(qRgb(255, 255, 129));
    colorTable.append(qRgb(255, 255, 132));
    colorTable.append(qRgb(255, 255, 135));
    colorTable.append(qRgb(255, 255, 138));
    colorTable.append(qRgb(255, 255, 142));
    colorTable.append(qRgb(255, 255, 145));
    colorTable.append(qRgb(255, 255, 148));
    colorTable.append(qRgb(255, 255, 151));
    colorTable.append(qRgb(255, 255, 154));
    colorTable.append(qRgb(255, 255, 158));
    colorTable.append(qRgb(255, 255, 161));
    colorTable.append(qRgb(255, 255, 164));
    colorTable.append(qRgb(255, 255, 167));
    colorTable.append(qRgb(255, 255, 171));
    colorTable.append(qRgb(255, 255, 174));
    colorTable.append(qRgb(255, 255, 177));
    colorTable.append(qRgb(255, 255, 180));
    colorTable.append(qRgb(255, 255, 183));
    colorTable.append(qRgb(255, 255, 187));
    colorTable.append(qRgb(255, 255, 190));
    colorTable.append(qRgb(255, 255, 193));
    colorTable.append(qRgb(255, 255, 196));
    colorTable.append(qRgb(255, 255, 200));
    colorTable.append(qRgb(255, 255, 203));
    colorTable.append(qRgb(255, 255, 206));
    colorTable.append(qRgb(255, 255, 209));
    colorTable.append(qRgb(255, 255, 213));
    colorTable.append(qRgb(255, 255, 216));
    colorTable.append(qRgb(255, 255, 219));
    colorTable.append(qRgb(255, 255, 222));
    colorTable.append(qRgb(255, 255, 225));
    colorTable.append(qRgb(255, 255, 229));
    colorTable.append(qRgb(255, 255, 232));
    colorTable.append(qRgb(255, 255, 235));
    colorTable.append(qRgb(255, 255, 238));
    colorTable.append(qRgb(255, 255, 242));
    colorTable.append(qRgb(255, 255, 245));
    colorTable.append(qRgb(255, 255, 248));
    colorTable.append(qRgb(255, 255, 251));
    colorTable.append(qRgb(255, 255, 255));
    this->addToColorLookUpTables(colorTable,"Standard Gamma II");

    colorTable.clear();
    colorTable.append(qRgb(2, 19, 7));
    colorTable.append(qRgb(2, 19, 7));
    colorTable.append(qRgb(3, 20, 11));
    colorTable.append(qRgb(4, 22, 15));
    colorTable.append(qRgb(5, 25, 19));
    colorTable.append(qRgb(6, 29, 22));
    colorTable.append(qRgb(7, 33, 26));
    colorTable.append(qRgb(8, 38, 30));
    colorTable.append(qRgb(9, 44, 34));
    colorTable.append(qRgb(13, 51, 39));
    colorTable.append(qRgb(11, 58, 43));
    colorTable.append(qRgb(13, 65, 47));
    colorTable.append(qRgb(14, 73, 51));
    colorTable.append(qRgb(15, 81, 55));
    colorTable.append(qRgb(16, 89, 59));
    colorTable.append(qRgb(17, 98, 63));
    colorTable.append(qRgb(18, 106, 67));
    colorTable.append(qRgb(19, 115, 71));
    colorTable.append(qRgb(21, 124, 75));
    colorTable.append(qRgb(22, 132, 78));
    colorTable.append(qRgb(23, 140, 82));
    colorTable.append(qRgb(24, 148, 86));
    colorTable.append(qRgb(25, 155, 90));
    colorTable.append(qRgb(26, 162, 93));
    colorTable.append(qRgb(27, 169, 97));
    colorTable.append(qRgb(28, 175, 100));
    colorTable.append(qRgb(30, 180, 103));
    colorTable.append(qRgb(31, 184, 107));
    colorTable.append(qRgb(32, 188, 110));
    colorTable.append(qRgb(33, 191, 113));
    colorTable.append(qRgb(34, 193, 116));
    colorTable.append(qRgb(35, 194, 119));
    colorTable.append(qRgb(36, 195, 122));
    colorTable.append(qRgb(37, 194, 125));
    colorTable.append(qRgb(38, 193, 128));
    colorTable.append(qRgb(39, 191, 130));
    colorTable.append(qRgb(39, 188, 133));
    colorTable.append(qRgb(40, 184, 136));
    colorTable.append(qRgb(41, 180, 139));
    colorTable.append(qRgb(42, 175, 141));
    colorTable.append(qRgb(43, 169, 144));
    colorTable.append(qRgb(44, 162, 147));
    colorTable.append(qRgb(45, 155, 150));
    colorTable.append(qRgb(46, 148, 153));
    colorTable.append(qRgb(46, 140, 156));
    colorTable.append(qRgb(47, 132, 159));
    colorTable.append(qRgb(48, 124, 162));
    colorTable.append(qRgb(49, 115, 165));
    colorTable.append(qRgb(50, 107, 169));
    colorTable.append(qRgb(51, 98, 172));
    colorTable.append(qRgb(52, 89, 176));
    colorTable.append(qRgb(53, 81, 180));
    colorTable.append(qRgb(54, 73, 184));
    colorTable.append(qRgb(55, 65, 188));
    colorTable.append(qRgb(56, 58, 192));
    colorTable.append(qRgb(57, 51, 196));
    colorTable.append(qRgb(58, 44, 201));
    colorTable.append(qRgb(59, 38, 205));
    colorTable.append(qRgb(60, 33, 210));
    colorTable.append(qRgb(61, 29, 215));
    colorTable.append(qRgb(62, 25, 220));
    colorTable.append(qRgb(63, 22, 225));
    colorTable.append(qRgb(65, 20, 230));
    colorTable.append(qRgb(66, 19, 235));
    colorTable.append(qRgb(67, 19, 241));
    colorTable.append(qRgb(68, 19, 246));
    colorTable.append(qRgb(70, 20, 251));
    colorTable.append(qRgb(71, 22, 253));
    colorTable.append(qRgb(72, 25, 248));
    colorTable.append(qRgb(74, 29, 243));
    colorTable.append(qRgb(75, 33, 237));
    colorTable.append(qRgb(77, 38, 232));
    colorTable.append(qRgb(78, 44, 227));
    colorTable.append(qRgb(79, 51, 222));
    colorTable.append(qRgb(81, 58, 217));
    colorTable.append(qRgb(82, 65, 212));
    colorTable.append(qRgb(84, 73, 207));
    colorTable.append(qRgb(85, 81, 202));
    colorTable.append(qRgb(86, 89, 198));
    colorTable.append(qRgb(88, 98, 194));
    colorTable.append(qRgb(89, 107, 190));
    colorTable.append(qRgb(90, 115, 186));
    colorTable.append(qRgb(92, 124, 182));
    colorTable.append(qRgb(93, 132, 179));
    colorTable.append(qRgb(94, 140, 176));
    colorTable.append(qRgb(95, 148, 172));
    colorTable.append(qRgb(97, 155, 170));
    colorTable.append(qRgb(98, 162, 167));
    colorTable.append(qRgb(99, 169, 165));
    colorTable.append(qRgb(100, 175, 162));
    colorTable.append(qRgb(101, 180, 160));
    colorTable.append(qRgb(102, 184, 158));
    colorTable.append(qRgb(103, 188, 157));
    colorTable.append(qRgb(103, 191, 155));
    colorTable.append(qRgb(104, 193, 154));
    colorTable.append(qRgb(105, 194, 152));
    colorTable.append(qRgb(106, 195, 151));
    colorTable.append(qRgb(107, 194, 150));
    colorTable.append(qRgb(107, 193, 148));
    colorTable.append(qRgb(108, 191, 147));
    colorTable.append(qRgb(108, 188, 146));
    colorTable.append(qRgb(109, 184, 144));
    colorTable.append(qRgb(110, 180, 143));
    colorTable.append(qRgb(110, 175, 141));
    colorTable.append(qRgb(111, 169, 140));
    colorTable.append(qRgb(111, 162, 138));
    colorTable.append(qRgb(112, 155, 136));
    colorTable.append(qRgb(113, 148, 134));
    colorTable.append(qRgb(113, 140, 132));
    colorTable.append(qRgb(114, 132, 129));
    colorTable.append(qRgb(114, 124, 126));
    colorTable.append(qRgb(115, 115, 123));
    colorTable.append(qRgb(116, 106, 120));
    colorTable.append(qRgb(117, 98, 116));
    colorTable.append(qRgb(117, 89, 113));
    colorTable.append(qRgb(118, 81, 108));
    colorTable.append(qRgb(119, 73, 104));
    colorTable.append(qRgb(120, 65, 99));
    colorTable.append(qRgb(121, 58, 94));
    colorTable.append(qRgb(122, 51, 89));
    colorTable.append(qRgb(123, 44, 83));
    colorTable.append(qRgb(124, 38, 78));
    colorTable.append(qRgb(125, 33, 72));
    colorTable.append(qRgb(127, 29, 65));
    colorTable.append(qRgb(128, 25, 59));
    colorTable.append(qRgb(129, 22, 53));
    colorTable.append(qRgb(131, 20, 46));
    colorTable.append(qRgb(132, 19, 39));
    colorTable.append(qRgb(134, 19, 32));
    colorTable.append(qRgb(135, 19, 25));
    colorTable.append(qRgb(137, 20, 18));
    colorTable.append(qRgb(138, 22, 11));
    colorTable.append(qRgb(140, 25, 5));
    colorTable.append(qRgb(142, 29, 253));
    colorTable.append(qRgb(143, 33, 246));
    colorTable.append(qRgb(145, 38, 239));
    colorTable.append(qRgb(147, 44, 233));
    colorTable.append(qRgb(149, 51, 227));
    colorTable.append(qRgb(150, 58, 221));
    colorTable.append(qRgb(152, 65, 215));
    colorTable.append(qRgb(154, 73, 210));
    colorTable.append(qRgb(155, 81, 205));
    colorTable.append(qRgb(157, 89, 200));
    colorTable.append(qRgb(158, 98, 195));
    colorTable.append(qRgb(160, 107, 191));
    colorTable.append(qRgb(161, 115, 188));
    colorTable.append(qRgb(163, 124, 184));
    colorTable.append(qRgb(164, 132, 181));
    colorTable.append(qRgb(166, 140, 178));
    colorTable.append(qRgb(167, 148, 176));
    colorTable.append(qRgb(168, 155, 174));
    colorTable.append(qRgb(169, 162, 172));
    colorTable.append(qRgb(170, 169, 171));
    colorTable.append(qRgb(171, 175, 170));
    colorTable.append(qRgb(172, 180, 169));
    colorTable.append(qRgb(173, 184, 169));
    colorTable.append(qRgb(174, 188, 168));
    colorTable.append(qRgb(174, 191, 168));
    colorTable.append(qRgb(175, 193, 168));
    colorTable.append(qRgb(175, 194, 168));
    colorTable.append(qRgb(176, 195, 169));
    colorTable.append(qRgb(176, 194, 169));
    colorTable.append(qRgb(177, 193, 169));
    colorTable.append(qRgb(177, 191, 169));
    colorTable.append(qRgb(178, 188, 170));
    colorTable.append(qRgb(178, 184, 170));
    colorTable.append(qRgb(178, 180, 170));
    colorTable.append(qRgb(178, 175, 169));
    colorTable.append(qRgb(179, 169, 169));
    colorTable.append(qRgb(179, 162, 168));
    colorTable.append(qRgb(179, 155, 167));
    colorTable.append(qRgb(180, 148, 166));
    colorTable.append(qRgb(180, 140, 164));
    colorTable.append(qRgb(180, 132, 162));
    colorTable.append(qRgb(181, 124, 160));
    colorTable.append(qRgb(181, 115, 157));
    colorTable.append(qRgb(182, 106, 154));
    colorTable.append(qRgb(182, 98, 150));
    colorTable.append(qRgb(183, 89, 146));
    colorTable.append(qRgb(184, 81, 141));
    colorTable.append(qRgb(184, 73, 136));
    colorTable.append(qRgb(185, 65, 131));
    colorTable.append(qRgb(186, 58, 125));
    colorTable.append(qRgb(187, 51, 119));
    colorTable.append(qRgb(188, 44, 112));
    colorTable.append(qRgb(190, 38, 106));
    colorTable.append(qRgb(191, 33, 98));
    colorTable.append(qRgb(192, 29, 91));
    colorTable.append(qRgb(194, 25, 83));
    colorTable.append(qRgb(195, 22, 75));
    colorTable.append(qRgb(197, 20, 67));
    colorTable.append(qRgb(198, 19, 58));
    colorTable.append(qRgb(200, 19, 50));
    colorTable.append(qRgb(202, 19, 41));
    colorTable.append(qRgb(204, 20, 33));
    colorTable.append(qRgb(206, 22, 24));
    colorTable.append(qRgb(208, 25, 16));
    colorTable.append(qRgb(210, 29, 8));
    colorTable.append(qRgb(212, 33, 255));
    colorTable.append(qRgb(214, 38, 247));
    colorTable.append(qRgb(216, 44, 239));
    colorTable.append(qRgb(218, 51, 232));
    colorTable.append(qRgb(220, 58, 225));
    colorTable.append(qRgb(222, 65, 219));
    colorTable.append(qRgb(224, 73, 213));
    colorTable.append(qRgb(225, 81, 207));
    colorTable.append(qRgb(227, 89, 202));
    colorTable.append(qRgb(229, 98, 197));
    colorTable.append(qRgb(231, 107, 193));
    colorTable.append(qRgb(232, 115, 189));
    colorTable.append(qRgb(234, 124, 186));
    colorTable.append(qRgb(236, 132, 183));
    colorTable.append(qRgb(237, 140, 181));
    colorTable.append(qRgb(238, 148, 179));
    colorTable.append(qRgb(239, 155, 178));
    colorTable.append(qRgb(241, 162, 177));
    colorTable.append(qRgb(242, 169, 177));
    colorTable.append(qRgb(242, 175, 177));
    colorTable.append(qRgb(243, 180, 178));
    colorTable.append(qRgb(244, 184, 179));
    colorTable.append(qRgb(245, 188, 180));
    colorTable.append(qRgb(245, 191, 181));
    colorTable.append(qRgb(246, 193, 183));
    colorTable.append(qRgb(246, 194, 184));
    colorTable.append(qRgb(246, 195, 186));
    colorTable.append(qRgb(246, 194, 188));
    colorTable.append(qRgb(246, 193, 190));
    colorTable.append(qRgb(247, 191, 192));
    colorTable.append(qRgb(247, 188, 193));
    colorTable.append(qRgb(247, 184, 195));
    colorTable.append(qRgb(247, 180, 196));
    colorTable.append(qRgb(247, 175, 197));
    colorTable.append(qRgb(247, 169, 198));
    colorTable.append(qRgb(247, 162, 198));
    colorTable.append(qRgb(247, 155, 198));
    colorTable.append(qRgb(247, 148, 198));
    colorTable.append(qRgb(247, 140, 197));
    colorTable.append(qRgb(247, 132, 195));
    colorTable.append(qRgb(247, 124, 193));
    colorTable.append(qRgb(247, 115, 191));
    colorTable.append(qRgb(248, 106, 188));
    colorTable.append(qRgb(248, 98, 184));
    colorTable.append(qRgb(248, 89, 180));
    colorTable.append(qRgb(249, 81, 175));
    colorTable.append(qRgb(250, 73, 169));
    colorTable.append(qRgb(251, 65, 163));
    colorTable.append(qRgb(251, 58, 156));
    colorTable.append(qRgb(252, 51, 149));
    colorTable.append(qRgb(254, 44, 142));
    colorTable.append(qRgb(255, 38, 133));
    colorTable.append(qRgb(254, 33, 125));
    colorTable.append(qRgb(252, 29, 116));
    colorTable.append(qRgb(251, 25, 107));
    colorTable.append(qRgb(249, 22, 97));
    colorTable.append(qRgb(247, 20, 87));
    colorTable.append(qRgb(247, 20, 87));



    this->addToColorLookUpTables(colorTable,"Plasma");
    colorTable.clear();
}

/**
 *  @brief This method searchs in current directory after colorLookUpTables 's files, which end with .clut
 *  and loads the files. The *.clut files contains array of integers such:
 *
 *  red(0) green(0) blue(0)
 *   .....
 *  red(i) green(i) blue(i)
 *  red(i+1) green(i+1) blue(i+1)
 *  red(i+2) green(i+2) blue(i+2)
 * .....
 *  red(255) green(255) blue(255)
*/
void ImageControllerWidget::loadColorLookUpTables(QString path,bool saveGlobal){
    QDir dir(path);
    QStringList clutFiles=dir.entryList(QStringList("*.clut"));
    for (int index = 0; index < clutFiles.size(); ++index) {
        if(colorLookUpTableNames.indexOf(clutFiles.at(index))==-1){
            this->addToColorLookUpTables(readColorLookUpTable(clutFiles.at(index),path,saveGlobal),clutFiles.at(index));
        }
    }
}


QVector <QRgb> ImageControllerWidget::readColorLookUpTable(QString filename,QString path,bool saveGlobal){
    QFile file(QDir::toNativeSeparators(path+QDir::separator ()+filename));
    if(path!=QCoreApplication::applicationDirPath()&&saveGlobal)
        if(!QFile(QDir::toNativeSeparators(QCoreApplication::applicationDirPath()+QDir::separator ()+filename)).exists())
            QFile::copy(QDir::toNativeSeparators(path+QDir::separator ()+filename),QDir::toNativeSeparators(QCoreApplication::applicationDirPath()+QDir::separator ()+filename));
    QTextStream in(&file);
    file.open(QIODevice::ReadOnly);
    QVector <QRgb> rgbs;
    int values[3];
    while (!in.atEnd()) {
        in>>values[0];
        in>>values[1];
        in>>values[2];
        rgbs.append(qRgb(values[0],values[1],values[2]));
    }
    return rgbs;
}

QVector <QRgb> ImageControllerWidget::createLookUpTableLinear(int rStart, int rEnd, int gStart, int gEnd,
                                                              int bStart, int bEnd, int numberOfColors){
    QVector <QRgb> rbg;
    if((rStart>=0) && (gStart>=0) &&(bStart>=0) &&
            (rEnd<=numberOfColors) && (gEnd<=numberOfColors)&&(bEnd<=numberOfColors)){
        int divider=rEnd-rStart;
        QVector <int> red;
        int j=0;
        for(int i=0;i<numberOfColors;i++){
            if(i>rStart && i<=rEnd){
                red.append(j*255/divider);
                j++;
            }else if(i>rEnd){
                red.append(255);
            }else{
                red.append(0);
            }
        }
        divider=gEnd-gStart;
        QVector <int> green;
        j=0;
        for(int i=0;i<numberOfColors;i++){
            if(i>gStart && i<=gEnd){
                green.append(j*255/divider);
                j++;
            }else if(i>gEnd){
                green.append(255);
            }else{
                green.append(0);
            }
        }
        divider=bEnd-bStart;
        QVector <int> blue;
        j=0;
        for(int i=0;i<numberOfColors;i++){
            if(i>bStart && i<=bEnd){
                blue.append(j*255/divider);
                j++;
            }else if(i>bEnd){
                blue.append(255);
            }else{
                blue.append(0);
            }
        }

        for(int i=0;i<numberOfColors;i++){
            rbg.append(qRgb(red.at(i),green.at(i),blue.at(i)));
        }

    }
    return rbg;
}

void ImageControllerWidget::addToColorLookUpTables(QVector<QRgb> table, QString name){
    if(colorLookUpTableNames.indexOf(name)==-1){
        this->colorLookUpTables.append(table);
        this->colorLookUpTableNames.append(name);
    }
}
QString ImageControllerWidget::getColorTableSliderObjectName(){
    return this->colorTableSlider->objectName();
}

/**
 * @brief This method set stylesheet of elements in this widget adequate to given parameter
 *
 */

void ImageControllerWidget::setStyleSheets(){
    QString CheckBoxStyleSheetTemplate ="QCheckBox::indicator{width: CheckBoxIconWidthpx;height: CheckBoxIconHeightpx;}"
            "QCheckBox:indicator:checked{ image:url(:/images/%1.png);}"
            "QCheckBox:indicator:unchecked{ image:url(:/images/no%1.png);}"
            "QCheckBox::checked:hover { background-color:palette(highlight); }"
            "QCheckBox::unchecked:hover { background-color:#E0842F; }"
            // "QCheckBox{border:1px solid grey;}"
            "QCheckBox:checked{color:#%2;background-color:transparent;}"
            "QCheckBox:unchecked{color:#CheckBoxInactiveFontColor;background-color:transparent;}";
    CheckBoxStyleSheetTemplate.replace("CheckBoxIconWidth",QString::number(this->CheckBoxIconWidth));
    CheckBoxStyleSheetTemplate.replace("CheckBoxIconHeight",QString::number(this->CheckBoxIconHeight));
    CheckBoxStyleSheetTemplate.replace("CheckBoxInactiveFontColor",this->CheckBoxInactiveFontColor);
    this->starCheckBox->setStyleSheet(CheckBoxStyleSheetTemplate.arg(this->starCheckBox->objectName(),this->labelColorMap->value(this->starCheckBox->objectName())));
    this->contourCheckBox->setStyleSheet(CheckBoxStyleSheetTemplate.arg(this->contourCheckBox->objectName(),this->labelColorMap->value(this->contourCheckBox->objectName())));
    this->previewCheckBox->setStyleSheet(CheckBoxStyleSheetTemplate.arg(this->previewCheckBox->objectName(),this->labelColorMap->value(this->previewCheckBox->objectName())));
    this->localCheckBox->setStyleSheet(CheckBoxStyleSheetTemplate.arg(this->localCheckBox->objectName(),this->labelColorMap->value(this->pipeCheckBox->objectName())));
    this->linearCheckBox->setStyleSheet(CheckBoxStyleSheetTemplate.arg(this->linearCheckBox->objectName(),this->labelColorMap->value(this->linearCheckBox->objectName())));
    this->secondOrderCheckBox->setStyleSheet(CheckBoxStyleSheetTemplate.arg(this->secondOrderCheckBox->objectName(),this->labelColorMap->value(this->secondOrderCheckBox->objectName())));
    this->tauCheckBox->setStyleSheet(CheckBoxStyleSheetTemplate.arg(this->tauCheckBox->objectName(),this->labelColorMap->value(this->tauCheckBox->objectName())));
    this->relobCheckbox->setStyleSheet(CheckBoxStyleSheetTemplate.arg(this->relobCheckbox->objectName(),this->labelColorMap->value(this->tauCheckBox->objectName())));
    this->colorbarCheckBox->setStyleSheet(CheckBoxStyleSheetTemplate.arg(this->colorbarCheckBox->objectName(),this->labelColorMap->value(this->colorbarCheckBox->objectName())));
    this->invertColorsCheckBox->setStyleSheet(CheckBoxStyleSheetTemplate.arg(this->invertColorsCheckBox->objectName(),this->labelColorMap->value(this->invertColorsCheckBox->objectName())));
    this->pipeCheckBox->setStyleSheet(CheckBoxStyleSheetTemplate.arg(this->pipeCheckBox->objectName(),this->labelColorMap->value(this->pipeCheckBox->objectName())));
    this->antiAliasingCheckBox->setStyleSheet(CheckBoxStyleSheetTemplate.arg(this->antiAliasingCheckBox->objectName(),this->labelColorMap->value(this->antiAliasingCheckBox->objectName())));
    this->dopplerCatchCheckBox->setStyleSheet(CheckBoxStyleSheetTemplate.arg(this->dopplerCatchCheckBox->objectName(),this->labelColorMap->value(this->dopplerCatchCheckBox->objectName())));
    this->absoluteScaleCheckBox->setStyleSheet(CheckBoxStyleSheetTemplate.arg(this->absoluteScaleCheckBox->objectName(),this->labelColorMap->value(this->absoluteScaleCheckBox->objectName())));

    QString groupboxStyleSheetTemplate="QGroupBox{border:1px solid #%1;}";
    this->checkBoxGroupBox->setStyleSheet(groupboxStyleSheetTemplate.arg(this->northWestGroupBoxBorderColor));
    this->spinBoxGroupBox->setStyleSheet(groupboxStyleSheetTemplate.arg(this->northWestGroupBoxBorderColor));
    this->pushButtonGroupBox->setStyleSheet(groupboxStyleSheetTemplate.arg(this->northWestGroupBoxBorderColor));
    this->specpushButtonGroupBox->setStyleSheet(groupboxStyleSheetTemplate.arg(this->northWestGroupBoxBorderColor));
    groupboxStyleSheetTemplate="QGroupBox{border:0px ;}";
    this->pipeAntiGroupBox->setStyleSheet(groupboxStyleSheetTemplate);

    QString labelStyleSheetTemplate="QLabel{color:#%1;}";
    this->colorTableLabel->setStyleSheet(labelStyleSheetTemplate.arg(this->labelFontColor));
    //    this->colorTableMaxLabel->setStyleSheet(labelStyleSheetTemplate.arg(this->labelFontColor));
    //    this->colorTableMinLabel->setStyleSheet(labelStyleSheetTemplate.arg(this->labelFontColor));
    this->positionLabel->setStyleSheet(labelStyleSheetTemplate.arg(this->labelFontColor));
    this->contrastLabel->setStyleSheet(labelStyleSheetTemplate.arg(this->labelFontColor));
    this->positionAngleLabel->setStyleSheet(labelStyleSheetTemplate.arg(this->labelFontColor));
    this->numberOfPixelsLabel->setStyleSheet(labelStyleSheetTemplate.arg(this->labelFontColor));
    this->inclinationLabel->setStyleSheet(labelStyleSheetTemplate.arg(this->labelFontColor));
    this->redLambdaLabel->setStyleSheet(labelStyleSheetTemplate.arg(this->labelFontColor));
    this->blueLambdaLabel->setStyleSheet(labelStyleSheetTemplate.arg(this->labelFontColor));
    this->greenLambdaLabel->setStyleSheet(labelStyleSheetTemplate.arg(this->labelFontColor));
    this->velocityLabel->setStyleSheet(labelStyleSheetTemplate.arg(this->labelFontColor));
    this->moleculeLabel->setStyleSheet(labelStyleSheetTemplate.arg(this->labelFontColor));
    this->lineLabel->setStyleSheet(labelStyleSheetTemplate.arg(this->labelFontColor));
    this->viewangleLabel->setStyleSheet(labelStyleSheetTemplate.arg(this->labelFontColor));
    this->pointingLabel->setStyleSheet(labelStyleSheetTemplate.arg(this->labelFontColor));
    this->observerposLabel->setStyleSheet(labelStyleSheetTemplate.arg(this->labelFontColor));
    labelStyleSheetTemplate="QLabel{color:#%1;qproperty-alignment: 'Aligncenter'}";
    this->numberOfContoursLabel->setStyleSheet(labelStyleSheetTemplate.arg(this->labelFontColor));
    this->sizeLabel->setStyleSheet(labelStyleSheetTemplate.arg(this->labelFontColor));
    this->phiLabel->setStyleSheet(labelStyleSheetTemplate.arg(this->labelFontColor));
    this->saturationLabel->setStyleSheet(labelStyleSheetTemplate.arg(this->labelFontColor));


    QString spinBoxSytleSheetTemplate="QSpinBox{color:#%1;}";
    this->numberOfContoursSpinBox->setStyleSheet(spinBoxSytleSheetTemplate.arg(this->labelFontColor));
    this->numberOfPixelsSpinBox->setStyleSheet(spinBoxSytleSheetTemplate.arg(this->labelFontColor));
    this->lineSpinBox->setStyleSheet(spinBoxSytleSheetTemplate.arg(this->labelFontColor));

    QString doubleSpinBoxSytleSheetTemplate="QDoubleSpinBox{color:#%1;}";
    this->contrastSpinBox->setStyleSheet(doubleSpinBoxSytleSheetTemplate.arg(this->labelFontColor));
    this->saturationSpinBox->setStyleSheet(doubleSpinBoxSytleSheetTemplate.arg(this->labelFontColor));
    this->positionAngleSpinBox->setStyleSheet(doubleSpinBoxSytleSheetTemplate.arg(this->labelFontColor));
    this->inclinationSpinBox->setStyleSheet(doubleSpinBoxSytleSheetTemplate.arg(this->labelFontColor));
    this->phiSpinBox->setStyleSheet(doubleSpinBoxSytleSheetTemplate.arg(this->labelFontColor));
    this->velocitySpinBox->setStyleSheet(doubleSpinBoxSytleSheetTemplate.arg(this->labelFontColor));
    this->viewangleSpinBox->setStyleSheet(doubleSpinBoxSytleSheetTemplate.arg(this->labelFontColor));
    this->scaleSpinBox->setStyleSheet(doubleSpinBoxSytleSheetTemplate.arg(this->labelFontColor));

    QString comboBoxSytleSheetTemplate="QComboxBox{Color:#%1}";
    this->colorTableComboxBox->setStyleSheet(comboBoxSytleSheetTemplate.arg(this->labelFontColor));
    this->positionComboxBox->setStyleSheet(comboBoxSytleSheetTemplate.arg(this->labelFontColor));
    this->moleculeComboBox->setStyleSheet(comboBoxSytleSheetTemplate.arg(this->labelFontColor));

    QString lineEditSytleSheetTemplate="QLineEdit{Color:#%1}";
    this->sizeLineEdit->setStyleSheet(lineEditSytleSheetTemplate.arg(this->labelFontColor));
    this->redLambdaValueLineEdit->setStyleSheet(lineEditSytleSheetTemplate.arg(this->labelFontColor));
    this->blueLambdaValueLineEdit->setStyleSheet(lineEditSytleSheetTemplate.arg(this->labelFontColor));
    this->greenLambdaValueLineEdit->setStyleSheet(lineEditSytleSheetTemplate.arg(this->labelFontColor));

    this->redColorTuneValueLineEdit->setStyleSheet(lineEditSytleSheetTemplate.arg(this->labelFontColor));
    this->blueColorTuneValueLineEdit->setStyleSheet(lineEditSytleSheetTemplate.arg(this->labelFontColor));
    this->greenColorTuneValueLineEdit->setStyleSheet(lineEditSytleSheetTemplate.arg(this->labelFontColor));

    this->observerposXLineEdit->setStyleSheet(lineEditSytleSheetTemplate.arg(this->labelFontColor));
    this->observerposYLineEdit->setStyleSheet(lineEditSytleSheetTemplate.arg(this->labelFontColor));
    this->observerposZLineEdit->setStyleSheet(lineEditSytleSheetTemplate.arg(this->labelFontColor));


    this->pointXLineEdit->setStyleSheet(lineEditSytleSheetTemplate.arg(this->labelFontColor));
    this->pointYLineEdit->setStyleSheet(lineEditSytleSheetTemplate.arg(this->labelFontColor));
    this->pointZLineEdit->setStyleSheet(lineEditSytleSheetTemplate.arg(this->labelFontColor));





}


/**
 * @brief This method set invertPixels's tooltip adequate to given parameter
 *
 * @param active, this method gets the active-parameter from it's slot.
 *
 * @see activateInvertPixels()
 */

void ImageControllerWidget::setInvertColorsToolTip(bool active){
    if(active){
        this->invertColorsCheckBox->setToolTip(trUtf8("Invert colors is on"));
    }else{
        this->invertColorsCheckBox->setToolTip(trUtf8("Invert colors is off"));
    }
}
/**
 * @brief This method set pipe's tooltip adequate to given parameter
 *
 * @param active, this method gets the active-parameter from it's slot.
 *
 * @see activatePipe()
 */

void ImageControllerWidget::setPipeToolTip(bool active){
    if(active){
        this->pipeCheckBox->setToolTip(trUtf8("Communication over Pipe"));
    }else{
        this->pipeCheckBox->setToolTip(trUtf8("Communication over commandline"));
    }
}
/**
 * @brief This method set preview's tooltip adequate to given parameter
 *
 * @param active, this method gets the active-parameter from it's slot.
 *
 * @see activatePreview()
 */

void ImageControllerWidget::setPreviewToolTip(bool active){
    if(active){
        this->previewCheckBox->setToolTip(trUtf8("Preview is on"));
    }else{
        this->previewCheckBox->setToolTip(trUtf8("Preview is off"));
    }
}
/**
 * @brief This method set preview's tooltip adequate to given parameter
 *
 * @param active, this method gets the active-parameter from it's slot.
 *
 * @see activateLocal()
 */

void ImageControllerWidget::setLocalToolTip(bool active){
    if(active){
        this->localCheckBox->setToolTip(trUtf8("local observer is on"));
    }else{
        this->localCheckBox->setToolTip(trUtf8("local observer is off"));
    }
}
/**
 * @brief This method set second-order's tooltip adequate to given parameter
 *
 * @param active, this method gets the active-parameter from it's slot.
 *
 * @see activateSecondOrder()
 */

void ImageControllerWidget::setSecondOrderToolTip(bool active){
    if(active){
        this->secondOrderCheckBox->setToolTip(trUtf8("Second-order integration"));
    }else{
        this->secondOrderCheckBox->setToolTip(trUtf8("First-order integration"));
    }
}
/**
 * @brief This method set tau's tooltip adequate to given parameter
 *
 * @param active, this method gets the active-parameter from it's slot.
 *
 * @see activateTau()
 */
void ImageControllerWidget::setTauToolTip(bool active){
    if(active){
        this->tauCheckBox->setToolTip(trUtf8("computing the optical depth"));
    }else{
        this->tauCheckBox->setToolTip(trUtf8("ray-tracing a true image"));
    }
}

/**
 * @brief This method set Colorbar's tooltip adequate to given parameter
 *
 * @param active, this method gets the active-parameter from it's slot.
 *
 * @see activateColorbar()
 */
void ImageControllerWidget::setColorbarToolTip(bool active){
    if(active){
        this->colorbarCheckBox->setToolTip(trUtf8("Colorbar is on"));
    }else{
        this->colorbarCheckBox->setToolTip(trUtf8("colorbar is off"));
    }
}
/**
 * @brief This method set linear's tooltip adequate to given parameter
 *
 * @param active, this method gets the active-parameter from it's slot.
 *
 * @see activatelinear()
 *
 */
void ImageControllerWidget::setLinearToolTip(bool active){
    if(active){
        this->linearCheckBox->setToolTip(trUtf8("Linear scale"));
    }else{
        this->linearCheckBox->setToolTip(trUtf8("Logarithmic scale"));
    }
}
/**
 * @brief This method set star's tooltip adequate to given parameter
 *
 * @param active, this method gets the active-parameter from it's slot.
 *
 * @see activateStar()
 *
 */
void ImageControllerWidget::setStarToolTip(bool active){
    if(active)
        this->starCheckBox->setToolTip(trUtf8("Star is included"));
    else
        this->starCheckBox->setToolTip(trUtf8("Star is excluded"));
}
/**
 * @brief This method set relob's tooltip adequate to given parameter
 *
 * @param active, this method gets the active-parameter from it's slot.
 *
 * @see activateStar()
 *
 */
void ImageControllerWidget::setRelobToolTip(bool active){
    if(active)
        this->relobCheckbox->setToolTip(trUtf8("Observer position relative wrt pointing"));
    else
        this->relobCheckbox->setToolTip(trUtf8("Observer position absolute wrt pointing"));
}

/**
 * @brief This method set star's tooltip adequate to given parameter
 *
 * @param active, this method gets the active-parameter from it's slot.
 *
 * @see activateAntiAliasing();
 *
 */
void ImageControllerWidget::setAntiAliasingToolTip(bool active){
    if(active)
        this->antiAliasingCheckBox->setToolTip(trUtf8("Anti-aliasing is on"));
    else
        this->antiAliasingCheckBox->setToolTip(trUtf8("Anti-aliasing is off"));
}
/**
 * @brief This method set star's tooltip adequate to given parameter
 *
 * @param active, this method gets the active-parameter from it's slot.
 *
 * @see activateDopplerCatching();
 *
 */

void ImageControllerWidget::setDopplerCatchingToolTip(bool active){
    if(active)
        this->dopplerCatchCheckBox->setToolTip(trUtf8("Doppler cacthing method is on"));
    else
        this->dopplerCatchCheckBox->setToolTip(trUtf8("Doppler cacthing method is off"));
}
/**
 * @brief This method set contour's tooltip adequate to given parameter
 *
 * @param active, this method gets the active-parameter from it's slot.
 *
 * @see activateContour()
 */
void ImageControllerWidget::setContourToolTip(bool active){
    if(active)
        this->contourCheckBox->setToolTip(trUtf8("contour is on"));
    else
        this->contourCheckBox->setToolTip(trUtf8("contour is off"));
}


/**
 * @brief This method set numberOfContours's tooltip
 *
 *
 */
void ImageControllerWidget::setNumberOfContoursToolTip(int value){
    QString toolTipString(trUtf8("current number of contours: %1").arg(QString::number(value)));
    this->numberOfContoursLabel->setToolTip(toolTipString);
    this->numberOfContoursSpinBox->setToolTip(toolTipString);
}

/**
 * @brief This method set numberOfPixels's tooltip
 *
 *
 */
void ImageControllerWidget::setNumberOfPixelsToolTip(int value){
    QString toolTipString(trUtf8("current number of pixels: %1").arg(QString::number(value)));
    this->numberOfPixelsLabel->setToolTip(toolTipString);
    this->numberOfPixelsSpinBox->setToolTip(toolTipString);
}
/**
 * @brief This method set spositionAngleT's tooltip
 *
 *
 */
void ImageControllerWidget::setPositionAngleToolTip(double value){
    QString toolTipString(trUtf8("Camera Position'angle: %1 \u00BA").arg(QString::number(value)));
    this->positionAngleLabel->setToolTip(toolTipString);
    this->positionAngleSpinBox->setToolTip(toolTipString);
}
/**
 * @brief This method set phi's tooltip
 *
 *
 */
void ImageControllerWidget::setPhiToolTip(double value){
    QString toolTipString(trUtf8("\u03A6 : %1 \u00BA").arg(QString::number(value)));
    this->phiLabel->setToolTip(toolTipString);
    this->phiSpinBox->setDecimals(2);
    this->phiSpinBox->setToolTip(toolTipString);
}
/**
 * @brief This method set inclination's tooltip
 *
 *
 */
void ImageControllerWidget::setInclinationToolTip(double value){
    QString toolTipString(trUtf8("\u2220: %1 \u00BA").arg(QString::number(value)));
    this->inclinationLabel->setToolTip(toolTipString);
    this->inclinationSpinBox->setToolTip(toolTipString);
}
/**
 * @brief This method set saturation's tooltip
 *
 *
 */
void ImageControllerWidget::setSaturationToolTip(double value){
    QString toolTipString(trUtf8("Saturation'rate %1%").arg(QString::number(value)));
    this->saturationSpinBox->setToolTip(toolTipString);
    this->saturationLabel->setToolTip(toolTipString);
}
/**
 * @brief This method set sizeLineEdit's tooltip
 *
 *
 */
void ImageControllerWidget::setSizeToolTip(QString value){
    double dValue=value.toDouble();
    QString toolTipString(trUtf8("current Size: %1 %2").arg(QString::number(dValue,'e'),this->getSizeUnit()));
    this->sizeLineEdit->setToolTip(toolTipString);
    this->sizeLabel->setToolTip(toolTipString);
}
/**
 * @brief This method set redLambdaLineEdit's tooltip
 *
 *
 */
void ImageControllerWidget::setRedLambdaToolTip(QString value){
    double dValue=value.toDouble();
    QString toolTipString(trUtf8("current lambda: %1 %2").arg(QString::number(dValue),this->getLambdaUnit()));
    this->redLambdaIndexSlider->setToolTip(toolTipString);
    this->redLambdaValueLineEdit->setToolTip(toolTipString);
    this->redLambdaLabel->setToolTip(toolTipString);
}

/**
 * @brief This method set redLambdaLineEdit's tooltip
 *
 *
 */
void ImageControllerWidget::setBlueLambdaToolTip(QString value){
    double dValue=value.toDouble();
    QString toolTipString(trUtf8("current lambda: %1 %2").arg(QString::number(dValue),this->getLambdaUnit()));
    this->blueLambdaIndexSlider->setToolTip(toolTipString);
    this->blueLambdaValueLineEdit->setToolTip(toolTipString);
    this->blueLambdaLabel->setToolTip(toolTipString);
}
/**
 * @brief This method set redLambdaLineEdit's tooltip
 *
 *
 */
void ImageControllerWidget::setGreenLambdaToolTip(QString value){
    double dValue=value.toDouble();
    QString toolTipString(trUtf8("current lambda: %1 %2").arg(QString::number(dValue),this->getLambdaUnit()));
    this->greenLambdaIndexSlider->setToolTip(toolTipString);
    this->greenLambdaValueLineEdit->setToolTip(toolTipString);
    this->greenLambdaLabel->setToolTip(toolTipString);
}


/**
 * @brief This method set contrast's tooltip
 *
 *
 */
void ImageControllerWidget::setContrastToolTip(double value){
    if(!this->linearCheckBox->isChecked()){
        QString toolTipString(trUtf8("MaxLog's currrent Value:%1").arg(QString::number(value)));
        this->contrastSpinBox->setToolTip(toolTipString);
        this->contrastLabel->setToolTip(toolTipString);
    }
    else{
        QString toolTipString(trUtf8("because Linear scale is on"));
        this->contrastSpinBox->setToolTip(toolTipString);
        this->contrastLabel->setToolTip(toolTipString);
    }
}

void ImageControllerWidget::setRedLambdaLineEdit(double value){
    this->redLambdaValueLineEdit->setText(QString::number(value)+" "+getLambdaUnit());
}

/**
 * @brief This Slot change  colorTableSlider value;
 *
 * @param active, this method gets the active-parameter from signal valueChanged(int) of colortableSlider.
 * This method emit the sendCommand signal.
 *
 * @sendCommand
 */

void ImageControllerWidget::colorTableSliderValueChanged(int value){
    this->colorTableComboxBox->setCurrentIndex(colorTableComboxBox->findText(QString::number(value),Qt::MatchStartsWith));
    emit sendCommand(this->colorTableSlider->objectName());
}
/**
 * @brief This Slot change  redlambdaIndexSlider value;
 *
 * @param active, this method gets the active-parameter from signal valueChanged(int) of lambdaIndexSlider.
 * This method emit the sendCommand signal.
 *
 * @sendCommand
 */

void ImageControllerWidget::redLambdaIndexSliderValueChanged(int value){
    if(lambdaArray.size()>(value-1)){
        this->redLambdaValueLineEdit->setText(QString::number(this->lambdaArray.at(value-1))+" "+getLambdaUnit());
        this->redLambdaValueLineEdit->setCursorPosition(0);
        this->setRedLambdaToolTip(QString::number(this->lambdaArray.at(value-1)));
    }
    this->setLambdaLabelText();
    //  emit sendCommand(this->lambdaValueLineEdit->objectName());
}
/**
 * @brief This Slot change  greenlambdaIndexSlider value;
 *
 * @param active, this method gets the active-parameter from signal valueChanged(int) of lambdaIndexSlider.
 * This method emit the sendCommand signal.
 *
 * @sendCommand
 */

void ImageControllerWidget::greenLambdaIndexSliderValueChanged(int value){
    if(lambdaArray.size()>(value-1)){
        this->greenLambdaValueLineEdit->setText(QString::number(this->lambdaArray.at(value-1))+" "+getLambdaUnit());
        this->setGreenLambdaToolTip(QString::number(this->lambdaArray.at(value-1)));
        this->greenLambdaValueLineEdit->setCursorPosition(0);
        this->setLambdaLabelText();
    }
    //  emit sendCommand(this->lambdaValueLineEdit->objectName());
}
/**
 * @brief This Slot change  bluelambdaIndexSlider value;
 *
 * @param active, this method gets the active-parameter from signal valueChanged(int) of lambdaIndexSlider.
 * This method emit the sendCommand signal.
 *
 * @sendCommand
 */

void ImageControllerWidget::blueLambdaIndexSliderValueChanged(int value){
    if(lambdaArray.size()>(value-1)){
        this->blueLambdaValueLineEdit->setText(QString::number(this->lambdaArray.at(value-1))+" "+getLambdaUnit());
        this->setBlueLambdaToolTip(QString::number(this->lambdaArray.at(value-1)));
        this->blueLambdaValueLineEdit->setCursorPosition(0);
        this->setLambdaLabelText();
    }
    //  emit sendCommand(this->lambdaValueLineEdit->objectName());
}

/**
 * @brief This Slot change  colorTableComboxBox value;
 *
 * @param index, this method gets the index-parameter from signal valueChanged(index) of colorTableComboxBox.
 * This method emit the sendCommand signal.
 *
 * @sendCommand
 */
void ImageControllerWidget::colorTableComboxBoxValueChanged(int index){
    this->colorTableSlider->setValue(index+1);
    QString toolTipString(trUtf8("Color table %1 is chosen").arg(this->colorTableComboxBox->itemText(index)));
    this->colorTableSlider->setToolTip(toolTipString);
    this->colorTableComboxBox->setToolTip(toolTipString);
}
/**
 * @brief This Slot change directionComboxBox value;
 *
 * @param index, this method gets the index-parameter from signal valueChanged(index) of directionComboxBox.
 * This method emit the sendCommand signal.
 *
 * @sendCommand
 */
void ImageControllerWidget::positionComboxBoxValueChanged(int index){
    bool vertical=false;
    QString toolTipString(trUtf8("Direction %1 is chosen").arg(this->positionComboxBox->itemText(index)));
    if((index==1) |(index==0))
        vertical=true;
    if(this->vertical!=vertical)addWidgetsToMainLayout(vertical);
    this->positionComboxBox->setToolTip(toolTipString);

    emit sendCommand(this->positionComboxBox->objectName());
}
/**
 * @brief This Slot changes saturationSpinBox's value;
 *
 * @param value, this method gets the value-parameter from signal editingFinished() of saturationSpinBox.
 * This method emit the sendCommand signal.
 *
 * @sendCommand
 */


void ImageControllerWidget::saturationSpinBoxValueChanged(){
    double value=this->saturationSpinBox->value();
    if(this->saturation!=value){
        this->setSaturationToolTip(value);
        this->saturation=value;
        emit sendCommand(this->saturationSpinBox->objectName());

    }
}
/**
 * @brief This Slot changes positionAngleSpinBox's value;
 *
 * @param value, this method gets the value-parameter from signal editingFinished() of positionAngleSpinBox.
 * This method emit the sendCommand signal.
 *
 * @sendCommand
 */
void ImageControllerWidget::positionAngleSpinBoxValueChanged(){
    double value=this->positionAngleSpinBox->value();
    if(this->posAngle!=value){
        this->posAngle=value;
        this->setPositionAngleToolTip(value);
        emit sendCommand(this->positionAngleSpinBox->objectName());
    }
}
/**
 * @brief This Slot changes inclinationSpinBox's value;
 *
 * @param value, this method gets the value-parameter from signal editingFinished() of inclinationSpinBox.
 * This method emit the sendCommand signal.
 *
 * @sendCommand
 */
void ImageControllerWidget::inclinationSpinBoxValueChanged(){
    double value=this->inclinationSpinBox->value();
    if(this->inclination!=value){
        this->inclination=value;
        this->setInclinationToolTip(value);
        emit sendCommand(this->inclinationSpinBox->objectName());
    }
}
/**
 * @brief This Slot changes phiSpinBox's value;
 *
 * @param value, this method gets the value-parameter from signal editingFinished() of phiSpinBox.
 * This method emit the sendCommand signal.
 *
 * @sendCommand
 */
void ImageControllerWidget::phiSpinBoxValueChanged(){
    double value=this->phiSpinBox->value();
    if(this->phi!=value){
        this->phi=value;
        this->setPhiToolTip(value);
        emit sendCommand(this->phiSpinBox->objectName());
    }
}
/**
 * @brief This Slot is connected to signal clicked of renderImagePushButton;
 *
 * This method emit the sendCommand signal.
 *
 * @sendCommand
 */
void ImageControllerWidget::renderImagePushButtonClicked(){
    emit sendCommand(this->renderImagePushButton->objectName());
}
/**
 * @brief This Slot is connected to signal clicked of renderSpectrumPushButton;
 *
 * This method emit the sendCommand signal.
 *
 * @sendCommand
 */
void ImageControllerWidget::renderSpectrumPushButtonClicked(){
    this->setRenderSpectrumActive(true);
    emit sendCommand(this->renderSpectrumPushButton->objectName());
}


/**
 * @brief This Slot is connected to signal clicked of rightPushButton;
 *
 * This method emit the sendCommand signal.
 *
 * @sendCommand
 */
void ImageControllerWidget::rightPushButtonClicked(){
    emit removeRect();
    if(!this->isRenderSpectrumActive()){
        if(this->isZoomActive()){
            double ddd=0.1;
            double dx=zoomBox.at(1)-zoomBox.at(0);
            zoomBox[0]=zoomBox.at(0)+dx*ddd;
            zoomBox[1]=zoomBox.at(1)+dx*ddd;
            emit sendCommand(this->rightPushButton->objectName());

        }
    }else{
        emit sendCommand(this->rightPushButton->objectName());

    }

}
/**
 * @brief This Slot is connected to signal clicked of resultPushButton;
 *
 * This method emit the sendCommand signal.
 *
 * @sendCommand
 */
void ImageControllerWidget::leftPushButtonClicked(){
    emit removeRect();
    if(!this->isRenderSpectrumActive()){
        if(this->isZoomActive()){
            double ddd=0.1;
            double dx=zoomBox.at(1)-zoomBox.at(0);
            zoomBox[1]=zoomBox.at(1)-dx*ddd;
            zoomBox[0]=zoomBox.at(0)-dx*ddd;
            emit sendCommand(this->leftPushButton->objectName());
        }
    }else{
        emit sendCommand(this->leftPushButton->objectName());

    }
}
/**
 * @brief This Slot is connected to signal clicked of downPushButton;
 *
 * This method emit the sendCommand signal.
 *
 * @sendCommand
 */
void ImageControllerWidget::downPushButtonClicked(){
    emit removeRect();
    if(!this->isRenderSpectrumActive()){
        if(this->isZoomActive()){
            double ddd=0.1;
            double dx=zoomBox.at(3)-zoomBox.at(2);
            zoomBox[2]=zoomBox.at(2)-dx*ddd;
            zoomBox[3]=zoomBox.at(3)-dx*ddd;
            emit sendCommand(this->downPushButton->objectName());
        }
    }else{
        emit sendCommand(this->downPushButton->objectName());

    }
}
/**
 * @brief This Slot is connected to signal clicked of upPushButton;
 *
 * This method emit the sendCommand signal.
 *
 * @sendCommand
 */
void ImageControllerWidget::upPushButtonClicked(){
    emit removeRect();
    if(!this->isRenderSpectrumActive()){
        if(this->isZoomActive()){
            double ddd=0.1;
            double dx=zoomBox.at(3)-zoomBox.at(2);
            zoomBox[2]=zoomBox.at(2)+dx*ddd;
            zoomBox[3]=zoomBox.at(3)+dx*ddd;
            emit sendCommand(this->upPushButton->objectName());
        }
    }else{
        emit sendCommand(this->upPushButton->objectName());

    }
}
/**
 * @brief This Slot is connected to signal clicked of unzoomPushButton;
 *
 * This method emit the sendCommand signal.
 *
 * @sendCommand
 */
void ImageControllerWidget::unzoomImagePushButtonClicked(){
    emit sendCommand(this->unzoomImagePushButton->objectName());
}

/**
 * @brief This Slot changes sizeLineEdit's value und is connected to signal editingFinished of sizeLineEdit;
 *
 * This method emit the sendCommand signal.
 *
 * @sendCommand
 */
void ImageControllerWidget::sizeLineEditValueChanged(){
    bool ok;
    QString value=this->sizeLineEdit->text();
    value.remove(QRegExp(sizeRegExp));
    double dValue=value.toDouble(&ok);
    if(ok ){
        this->sizeLineEdit->setText(QString::number(dValue,'e')+"  "+this->getSizeUnit());
        this->setSizeToolTip(value);
        if( (dValue!=this->lastSizeValue)){
            this->lastSizeValue=dValue;
            emit sendCommand(this->sizeLineEdit->objectName());
        }
    }else{
        this->sizeLineEdit->setText(QString::number(this->lastSizeValue,'e')+"  "+this->getSizeUnit());
        this->setSizeToolTip(sizeLineEdit->text());
    }
    this->sizeLineEdit->setCursorPosition(0);
}
/**
 * @brief This Slot changes redlambdaValueLineEdit's value  und is connected to signal editingFinished of lambdaValueLineEdit;
 *
 * This method emit the sendCommand signal.
 *
 * @sendCommand
 */
void ImageControllerWidget::redLambdaValueLineEditValueChanged(){
    bool ok;
    QString value=this->redLambdaValueLineEdit->text();
    value.remove(QRegExp(lambdaRegExp));
    double dValue=value.toDouble(&ok);
    double lambda=0;
    if(ok){
        int index=lambdaSize;
        for (int var = 0; var < this->lambdaSize; ++var) {
            lambda=QString::number(this->lambdaArray.at(var)).toDouble();
            if(lambda>=dValue){
                index=var;
                break;
            }
        }
        this->redLambdaValueLineEdit->setCursorPosition(0);
        if((index+1)!=this->redLambdaIndexSlider->value())
            this->redLambdaIndexSlider->setValue(index+1);
        else
            this->redLambdaValueLineEdit->setText(QString::number(lambda)+" "+getLambdaUnit());
    }else{
        this->redLambdaValueLineEdit->setText(QString::number(this->lambdaArray.at(this->redLambdaIndexSlider->value()-1))+" "+getLambdaUnit());
    }
}


void ImageControllerWidget::blueColorTuneValueLineEditValueChanged(){
    bool ok=false;
    QString value=this->blueColorTuneValueLineEdit->text();
    double dValue=value.toDouble(&ok);
    if(ok){
        if(dValue!=this->lastBlueColorTuneValue){
            this->lastBlueColorTuneValue=dValue;
            emit sendCommand(this->blueColorTuneValueLineEdit->objectName());
        }
    }else
        this->blueColorTuneValueLineEdit->setText(QString::number(this->lastBlueColorTuneValue));
}


void ImageControllerWidget::greenColorTuneValueLineEditValueChanged(){
    bool ok=false;
    QString value=this->greenColorTuneValueLineEdit->text();
    double dValue= value.toDouble(&ok);
    if(ok){
        if(dValue!=this->lastGreenColorTuneValue){
            this->lastGreenColorTuneValue=dValue;
            emit sendCommand(this->greenColorTuneValueLineEdit->objectName());
        }
    }else
        this->greenColorTuneValueLineEdit->setText(QString::number(this->lastGreenColorTuneValue));
}
void ImageControllerWidget::setLineModus(bool active){
    this->lineModus=active;
    this->lineLabel->setVisible(active);
    this->lineSpinBox->setVisible(active);
    this->velocityLabel->setVisible(active);
    this->velocitySpinBox->setVisible(active);
    this->moleculeLabel->setVisible(active);
    this->moleculeComboBox->setVisible(active);
    if(active){
        this->setSendCommandLineList(true);
        QFile file(this->setup->getLinesFilname());
        QString line;
        if(file.exists()){
            file.open(QIODevice::ReadOnly);
            QTextStream in(&file);
            line=in.readLine();
            line=in.readLine();
            int total=line.toInt();
            this->blockMoleculeSpinBoxSignal(true);
            for (int var = 0; var < total; ++var) {
                line=in.readLine();
                QStringList list=line.split(QRegExp("\\s+"),QString::SkipEmptyParts);
                if(list.size()>0){
                    if(list.size()>1)
                        moleculeComboBox->addItem(list.at(0)+" ( "+list.at(1)+" ) ");
                    else
                        moleculeComboBox->addItem(list.at(0));
                }
            }
            this->blockMoleculeSpinBoxSignal(false);
            file.close();
        }
    }

}

void ImageControllerWidget::redColorTuneValueLineEditValueChanged(){
    bool ok=false;
    QString value=this->redColorTuneValueLineEdit->text();
    double dValue=value.toDouble(&ok);
    if(ok){
        if(dValue!=this->lastRedColorTuneValue){
            this->lastRedColorTuneValue=dValue;
            emit sendCommand(this->redColorTuneValueLineEdit->objectName());
        }
    }else
        this->redColorTuneValueLineEdit->setText(QString::number(this->lastRedColorTuneValue));
}

void ImageControllerWidget::velocitySpinBoxValueChanged(){
    double value=this->velocitySpinBox->value();
    this->lastVelocityValue=value;
    this->setDoLine(true);
    this->velocitySpinBox->blockSignals(true);
    if(velocitySpinBox->isEnabled()){
        this->velocitySpinBox->setEnabled(false);
        emit sendCommand(this->velocitySpinBox->objectName());
    }
}
void ImageControllerWidget::lineSpinBoxValueChanged(){
    int value=this->lineSpinBox->value();
    this->setDoLine(true);
    this->lastLineValue=value;
    this->lineSpinBox->blockSignals(true);
    if(lineSpinBox->isEnabled()){
        this->lineSpinBox->setEnabled(false);
        emit sendCommand(this->lineSpinBox->objectName());
    }

}
void ImageControllerWidget::moleculeComboBoxValueChanged(int value){
    this->setDoLine(true);
    lastMoleculeValue=value;
    this->moleculeComboBox->blockSignals(true);
    emit sendCommand(this->moleculeComboBox->objectName());

}

void ImageControllerWidget::viewangleSpinBoxValueChanged(){
    emit sendCommand(this->viewangleSpinBox->objectName());

}

void ImageControllerWidget::scaleSpinBoxValueChanged(){
    double scale=scaleSpinBox->value();
    if(scale!=0.0){
        double x=this->observerposXLineEdit->text().remove(QRegExp(sizeRegExp)).toDouble()*scale;
        double y=this->observerposYLineEdit->text().remove(QRegExp(sizeRegExp)).toDouble()*scale;
        double z=this->observerposZLineEdit->text().remove(QRegExp(sizeRegExp)).toDouble()*scale;
        this->observerposXLineEdit->blockSignals(true);
        this->observerposYLineEdit->blockSignals(true);
        this->observerposZLineEdit->blockSignals(true);
        this->observerposXLineEdit->setText(QString::number(x,'e')+" "+getSizeUnit());
        this->observerposYLineEdit->setText(QString::number(y,'e')+" "+getSizeUnit());
        this->observerposZLineEdit->setText(QString::number(z,'e')+" "+getSizeUnit());
        this->observerposXLineEdit->blockSignals(false);
        this->observerposYLineEdit->blockSignals(false);
        this->observerposZLineEdit->blockSignals(false);
        emit sendCommand(this->relobCheckbox->objectName());
    }else{
        QMessageBox::information(this,trUtf8("Invalid scale"),trUtf8("Zero is invalid scale"));
    }
}

/**
 * @brief This Slot changes bluelambdaValueLineEdit's value  und is connected to signal editingFinished of lambdaValueLineEdit;
 *
 * This method emit the sendCommand signal.
 *
 * @sendCommand
 */
void ImageControllerWidget::blueLambdaValueLineEditValueChanged(){
    bool ok;
    QString value=this->blueLambdaValueLineEdit->text();
    value.remove(QRegExp(lambdaRegExp));
    double dValue=value.toDouble(&ok);
    double lambda=0;
    if(ok){
        int index=lambdaSize;
        for (int var = 0; var < this->lambdaSize; ++var) {
            lambda=QString::number(this->lambdaArray.at(var)).toDouble();
            if(lambda>=dValue){
                index=var;
                break;
            }
        }
        this->blueLambdaValueLineEdit->setCursorPosition(0);
        if((index+1)!=this->blueLambdaIndexSlider->value())
            this->blueLambdaIndexSlider->setValue(index+1);
        else
            this->blueLambdaValueLineEdit->setText(QString::number(lambda)+" "+getLambdaUnit());
    }else{
        this->blueLambdaValueLineEdit->setText(QString::number(this->lambdaArray.at(this->blueLambdaIndexSlider->value()-1))+" "+getLambdaUnit());
    }
}
/**
 * @brief This Slot changes greenlambdaValueLineEdit's value  und is connected to signal editingFinished of lambdaValueLineEdit;
 *
 * This method emit the sendCommand signal.
 *
 * @sendCommand
 */
void ImageControllerWidget::greenLambdaValueLineEditValueChanged(){
    bool ok;
    QString value=this->greenLambdaValueLineEdit->text();
    value.remove(QRegExp(lambdaRegExp));
    double dValue=value.toDouble(&ok);
    double lambda=0;
    if(ok){
        int index=lambdaSize;
        for (int var = 0; var < this->lambdaSize; ++var) {
            lambda=QString::number(this->lambdaArray.at(var)).toDouble();
            if(lambda>=dValue){
                index=var;
                break;
            }
        }
        this->greenLambdaValueLineEdit->setCursorPosition(0);
        if((index+1)!=this->greenLambdaIndexSlider->value())
            this->greenLambdaIndexSlider->setValue(index+1);
        else
            this->greenLambdaValueLineEdit->setText(QString::number(lambda)+" "+getLambdaUnit());
    }else{
        this->greenLambdaValueLineEdit->setText(QString::number(this->lambdaArray.at(this->greenLambdaIndexSlider->value()-1))+" "+getLambdaUnit());
    }
}

/**
 * @brief This slot  will set the scroll of debugger's to Bottom, if the debugger's tab is selected.
 */
void ImageControllerWidget::middleTabWidgetCurrentIndex(int index){
    if(this->middleTabWidget->indexOf(this->debuggerGroupBox)==index)this->debugger->scrollToBottom();
}

/**
 * @brief This slot sets the scroll of debugger's to Bottom.
 */
void ImageControllerWidget::scrollDebuggerToBottom(){
    this->debugger->scrollToBottom();
}

/**
 * @brief This slot sets the scroll of debugger's to top.
 */
void ImageControllerWidget::scrollDebuggerToTop(){
    this->debugger->scrollToTop();
}
/**
 * @brief This slot clears debugger's messages.
 */

void ImageControllerWidget::clearDebugger(){
    this->debugger->model()->removeRows(0,this->debugger->model()->rowCount());
}

/**
 * @brief This Slot changes contrastSpinBox's value;
 *
 * @param value, this method gets the value-parameter from signal editingFinished() of constrastSpinBox.
 * This method emit the sendCommand signal.
 *
 * @sendCommand
 */
void ImageControllerWidget::contrastSpinBoxValueChanged(){
    double value=this->contrastSpinBox->value();
    if(this->maxlog!=value){
        this->maxlog=value;
        this->setContrastToolTip(value);
        emit sendCommand(this->contrastSpinBox->objectName());
    }
}
/**
 * @brief This Slot changes numbrOfContours's value;
 *
 * @param value, this method gets the value-parameter from signal editingFinished() of numberOfContourSpinBox.
 * This method emit the sendCommand signal.
 *
 * @sendCommand
 */
void ImageControllerWidget::numberOfContoursSpinBoxValueChanged(){
    int value=this->numberOfContoursSpinBox->value();
    if(numberOfContours!=value){
        numberOfContours=value;
        this->setNumberOfContoursToolTip(value);
        this->contourCheckBox->blockSignals(true);
        if(!contourCheckBox->isChecked())this->contourCheckBox->setChecked(true);
        this->contourCheckBox->blockSignals(false);
        emit sendCommand(this->numberOfContoursSpinBox->objectName());
    }
}
/**
 * @brief This Slot changes numbrOfPixels's value;
 *
 * @param value, this method gets the value-parameter from signal  editingFinished() of numberOfPixels.
 * This method emit the sendCommand signal.
 *
 * @sendCommand
 */
void ImageControllerWidget::numberOfPixelsSpinBoxValueChanged(){
    int value=this->numberOfPixelsSpinBox->value();
    if(numberOfPixels!=value){
        numberOfPixels=value;
        this->setNumberOfPixelsToolTip(value);
        emit sendCommand(this->numberOfPixelsSpinBox->objectName());
    }
}
/**
 * @brief This Slot changes pointXLineEdit's value;
 *
 * @param value, this method gets the value-parameter from signal  editingFinished() of pointXLineEdit.
 * This method emit the sendCommand signal.
 *
 * @sendCommand
 */
void ImageControllerWidget::pointXLineEditValueChanged(){
    bool ok;
    double value=this->pointXLineEdit->text().remove(QRegExp(this->sizeRegExp)).trimmed().toDouble(&ok);
    if(value!=lastPointX && ok ){
        this->lastPointX=value;
        this->pointXLineEdit->setText(QString::number(value,'e')+"  "+this->getSizeUnit());
        this->pointXLineEdit->setCursorPosition(0);
        emit this->sendCommand(this->pointXLineEdit->objectName());
    }

}
/**
 * @brief This Slot changes pointYLineEdit's value;
 *
 * @param value, this method gets the value-parameter from signal  editingFinished() of pointYLineEdit.
 * This method emit the sendCommand signal.
 *
 * @sendCommand
 */
void ImageControllerWidget::pointYLineEditValueChanged(){
    bool ok;
    double value=this->pointYLineEdit->text().remove(QRegExp(this->sizeRegExp)).trimmed().toDouble(&ok);;
    if(value!=lastPointY && ok ){
        this->lastPointY=value;
        this->pointYLineEdit->setText(QString::number(value,'e')+"  "+this->getSizeUnit());
        this->pointYLineEdit->setCursorPosition(0);
        emit this->sendCommand(this->pointYLineEdit->objectName());
    }

}
/**
 * @brief This Slot changes pointZLineEdit's value;
 *
 * @param value, this method gets the value-parameter from signal  editingFinished() of pointZLineEdit.
 * This method emit the sendCommand signal.
 *
 * @sendCommand
 */
void ImageControllerWidget::pointZLineEditValueChanged(){
    bool ok;
    double value=this->pointZLineEdit->text().remove(QRegExp(this->sizeRegExp)).trimmed().toDouble(&ok);;
    if(value!=lastPointZ && ok ){
        this->lastPointZ=value;
        this->pointZLineEdit->setText(QString::number(value,'e')+"  "+this->getSizeUnit());
        this->pointZLineEdit->setCursorPosition(0);
        emit this->sendCommand(this->pointZLineEdit->objectName());
    }

}


/**
 * @brief This Slot changes observerposXLineEdit's value;
 *
 * @param value, this method gets the value-parameter from signal  editingFinished() of observerposXLineEdit.
 * This method emit the sendCommand signal.
 *
 * @sendCommand
 */
void ImageControllerWidget::observerposXLineEditValueChanged(){
    bool ok;
    double value=this->observerposXLineEdit->text().remove(QRegExp(this->sizeRegExp)).trimmed().toDouble(&ok);;
    if(value!=lastobsX && ok ){
        this->lastobsX=value;
        this->observerposXLineEdit->setText(QString::number(value,'e')+"  "+this->getSizeUnit());
        this->observerposXLineEdit->setCursorPosition(0);
        emit this->sendCommand(this->observerposXLineEdit->objectName());;
    }

}

/**
 * @brief This Slot changes observerposYLineEdit's value;
 *
 * @param value, this method gets the value-parameter from signal  editingFinished() of observerposYLineEdit.
 * This method emit the sendCommand signal.
 *
 * @sendCommand
 */
void ImageControllerWidget::observerposYLineEditValueChanged(){
    bool ok;
    double value=this->observerposYLineEdit->text().remove(QRegExp(this->sizeRegExp)).trimmed().toDouble(&ok);;
    if(value!=lastobsY && ok ){
        this->lastobsY=value;
        this->observerposYLineEdit->setText(QString::number(value,'e')+"  "+this->getSizeUnit());
        this->observerposYLineEdit->setCursorPosition(0);
        emit this->sendCommand(this->observerposYLineEdit->objectName());
    }

}

/**
 * @brief This Slot changes observerposZLineEdit's value;
 *
 * @param value, this method gets the value-parameter from signal  editingFinished() of observerposZLineEdit.
 * This method emit the sendCommand signal.
 *
 * @sendCommand
 */
void ImageControllerWidget::observerposZLineEditValueChanged(){
    bool ok;
    double value=this->observerposZLineEdit->text().remove(QRegExp(this->sizeRegExp)).trimmed().toDouble(&ok);;
    if(value!=lastobsZ && ok ){
        this->lastobsZ=value;
        this->observerposZLineEdit->setText(QString::number(value,'e')+"  "+this->getSizeUnit());
        this->observerposZLineEdit->setCursorPosition(0);
        emit this->sendCommand(this->observerposZLineEdit->objectName());
    }

}

/**
 * @brief This method sets constrastSpinBox enable/disable
 *
*/
void ImageControllerWidget::setConstrastSpinBoxEnabled(bool active){
    if(active)
        this->setContrastToolTip(this->contrastSpinBox->value());
    //this->contrastSpinBox->setEnabled(active);


}




/**
 * @brief This Slot in/excludes the star.
 *
 * @param active, this method gets the active-parameter from signal clicked(bool) of starCheckBox.
 * This method emit the sendCommand signal.
 *
 * @sendCommand
 */
void ImageControllerWidget::activateColorbar(bool active){
    this->setColorbarToolTip(active);
    emit sendCommand(this->colorbarCheckBox->objectName());
    emit checktStateChanged(this->colorbarCheckBox->objectName());
}
/**
 * @brief This Slot turns communication over the pipe on/off.
 *
 * @param active, this method gets the active-parameter from signal clicked(bool) of pipeCheckBox.
 * This method emit the sendCommand signal.
 *
 * @sendCommand
 */
void ImageControllerWidget::activatePipe(bool active){
    if(active)
        this->logger->writeToLogFile(trUtf8("The communication is now over the pipe"));
    else
        this->logger->writeToLogFile(trUtf8("The communication is now over filesystem"));
    this->setPipeToolTip(active);
    if(lineModus)this->setSendCommandLineList(true);
    emit checktStateChanged(this->pipeCheckBox->objectName());
    emit sendCommand(this->pipeCheckBox->objectName());
}
void ImageControllerWidget::setColorTuneEnabled(bool active){
    if(!active){
        this->redColorTuneValueLineEdit->setToolTip(trUtf8("Active just by logarithmical modus"));
        this->blueColorTuneValueLineEdit->setToolTip(trUtf8("Active just by logarithmical modus"));
        this->greenColorTuneValueLineEdit->setToolTip(trUtf8("Active just by logarithmical modus"));
    }else{
        this->redColorTuneValueLineEdit->setToolTip(this->redColorTuneValueLineEdit->text());
        this->blueColorTuneValueLineEdit->setToolTip(this->blueColorTuneValueLineEdit->text());
        this->greenColorTuneValueLineEdit->setToolTip(this->greenColorTuneValueLineEdit->text());
    }
    this->redColorTuneValueLineEdit->setEnabled(active);
    this->blueColorTuneValueLineEdit->setEnabled(active);
    this->greenColorTuneValueLineEdit->setEnabled(active);
}

void ImageControllerWidget::activateAbsoluteScale(bool active){
    if(active){
        this->lastColorTune.resize(3);
        this->lastColorTune[0]=this->lastRedColorTuneValue;
        this->lastColorTune[1]=this->lastGreenColorTuneValue;
        this->lastColorTune[2]=this->lastBlueColorTuneValue;
        this->lastRedColorTuneValue=lastColorMaximumArray[0];
        this->lastGreenColorTuneValue=lastColorMaximumArray[1];
        this->lastBlueColorTuneValue=lastColorMaximumArray[2];
        this->blockSignals(true);
        this->redColorTuneValueLineEdit->setText(QString::number(this->lastRedColorTuneValue));
        this->greenColorTuneValueLineEdit->setText(QString::number(this->lastGreenColorTuneValue));
        this->blueColorTuneValueLineEdit->setText(QString::number(this->lastBlueColorTuneValue));
        this->blockSignals(false);
        this->logger->writeToLogFile("Absolute scale is now on");
    }else{
        if(lastColorTune.size()==3){
            this->blockSignals(true);
            this->lastRedColorTuneValue=lastColorTune[0];
            this->lastGreenColorTuneValue=lastColorTune[1];
            this->lastBlueColorTuneValue=lastColorTune[2];
            this->redColorTuneValueLineEdit->setText(QString::number(this->lastRedColorTuneValue));
            this->greenColorTuneValueLineEdit->setText(QString::number(this->lastGreenColorTuneValue));
            this->blueColorTuneValueLineEdit->setText(QString::number(this->lastBlueColorTuneValue));
            this->blockSignals(false);
        }
        this->logger->writeToLogFile("Absolute scale is now off");
    }
    emit checktStateChanged(this->absoluteScaleCheckBox->objectName());
    emit sendCommand(this->absoluteScaleCheckBox->objectName());
}

void ImageControllerWidget::activateLocalObserver(bool active){
    localModus=active;
    this->pointingLabel->setVisible(active);
    this->pointXLineEdit->setVisible(active);
    this->pointYLineEdit->setVisible(active);
    this->pointZLineEdit->setVisible(active);

    this->observerposLabel->setVisible(active);
    this->observerposXLineEdit->setVisible(active);
    this->observerposYLineEdit->setVisible(active);
    this->observerposZLineEdit->setVisible(active);

    this->inclinationLabel->setVisible(!active);
    this->inclinationSpinBox->setVisible(!active);
    this->phiLabel->setVisible(!active);
    this->phiSpinBox->setVisible(!active);


    this->relobCheckbox->setVisible(active);
    this->scaleSpinBox->setVisible(active);
    this->viewangleLabel->setVisible(active);
    this->viewangleSpinBox->setVisible(active);
    if(viewangleSpinBox->value()==0)this->viewangleSpinBox->setValue(0.8);
    if(observerposYLineEdit->text().remove(QRegExp(sizeRegExp)).toDouble()==0)setSizeLineEdit();
    if(active){
        this->previewCheckBox->setVisible(false);
        this->secondOrderCheckBox->setVisible(false);
        this->localCheckBox->setVisible(true);
        this->localCheckBox->blockSignals(true);
        this->localCheckBox->setChecked(true);
        this->localCheckBox->blockSignals(false);
    }else{
        this->localCheckBox->setVisible(false);
        this->previewCheckBox->setVisible(true);
        this->secondOrderCheckBox->setVisible(true);
    }
    updateImageMenuActionVisibilties();


}

/**
 * @brief This Slot turn invert colors on/off.
 *
 * @param active, this method gets the active-parameter from signal clicked(bool) of invertColorsCheckBox.
 * This method emit the sendCommand signal.
 *
 * @sendCommand
 */
void ImageControllerWidget::activateInvertColors(bool active){
    this->setInvertColorsToolTip(active);
    emit checktStateChanged(this->invertColorsCheckBox->objectName());
    emit sendCommand(this->invertColorsCheckBox->objectName());
}
/**
 * @brief This Slot in/excludes the star.
 *
 * @param active, this method gets the active-parameter from signal clicked(bool) of starCheckBox.
 * This method emit the sendCommand signal.
 *
 * @sendCommand
 */
void ImageControllerWidget::activateStar(bool active){
    this->setStarToolTip(active);
    emit checktStateChanged(this->starCheckBox->objectName());
    emit sendCommand(this->starCheckBox->objectName());
}
/**
 * @brief This Slot does observer position relative or absolute wrt pointing.
 *
 * @param active, this method gets the active-parameter from signal clicked(bool) of relobCheckBox.
 * This method emit the sendCommand signal.
 *
 * @sendCommand
 */
void ImageControllerWidget::activateRelob(bool active){
    this->setRelobToolTip(active);
    emit checktStateChanged(this->relobCheckbox->objectName());
    emit sendCommand(this->relobCheckbox->objectName());
}
/**
 * @brief This Slot set antialiasing.
 *
 * @param active, this method gets the active-parameter from signal clicked(bool) of antiAliasingCheckBox.
 * This method emit the sendCommand signal.
 *
 * @sendCommand
 */
void ImageControllerWidget::activateAntiAliasing(bool active){
    setAntiAliasingToolTip(active);
    emit checktStateChanged(this->antiAliasingCheckBox->objectName());
    emit sendCommand(this->antiAliasingCheckBox->objectName());
}
/**
 * @brief This Slot set doppler Cacthing.
 *
 * @param active, this method gets the active-parameter from signal clicked(bool) of dopplerCatchCheckBox.
 * This method emit the sendCommand signal.
 *
 * @sendCommand
 */
void ImageControllerWidget::activateDopplerCatching(bool active){
    setDopplerCatchingToolTip(active);
    emit checktStateChanged(this->dopplerCatchCheckBox->objectName());
    emit sendCommand(this->dopplerCatchCheckBox->objectName());
}
/**
 * @brief This Slot set linear/logarithmic scale
 *
 * @param active, this method gets the active-parameter from signal clicked(bool) of linearCheckBox.
 * This method emit the sendCommand signal.
 *
 * @sendCommand
 */
void ImageControllerWidget::activateLinear(bool active){

    this->setLinearToolTip(active);
    this->setConstrastSpinBoxEnabled(!active);
    emit checktStateChanged(this->linearCheckBox->objectName());
    emit sendCommand(this->linearCheckBox->objectName());
}
/**
 * @brief This Slot sets Second/ First Order Integration
 *
 * @param active, this method gets the active-parameter from signal clicked(bool) of secondOrderCheckBox.
 * This method emit the sendCommand signal.
 *
 * @sendCommand
 */

void ImageControllerWidget::activateSecondOrder(bool active){


    this->setSecondOrderToolTip(active);
    emit checktStateChanged(this->secondOrderCheckBox->objectName());
    emit this->sendCommand(this->secondOrderCheckBox->objectName());

}
/**
 * @brief This Slot sets Second/ First Order Integration
 *
 * @param active, this method gets the active-parameter from signal clicked(bool) of secondOrderCheckBox.
 * This method emit the sendCommand signal.
 *
 * @sendCommand
 */

void ImageControllerWidget::activatePreview(bool active){
    this->setPreviewToolTip(active);
    emit checktStateChanged(this->previewCheckBox->objectName());
    emit this->sendCommand(this->previewCheckBox->objectName());

}
/**
 * @brief This Slot sets localObserver
 *
 * @param active, this method gets the active-parameter from signal clicked(bool) of localCheckBox.
 * This method emit the sendCommand signal.
 *
 * @sendCommand
 */

void ImageControllerWidget::activateLocal(bool active){
    this->setLocalToolTip(active);
    this->setLocalModus(active);
    this->previewCheckBox->setVisible(!active);
    this->secondOrderCheckBox->setVisible(!active);
    this->inclinationLabel->setVisible(!active);
    this->inclinationSpinBox->setVisible(!active);
    this->phiLabel->setVisible(!active);
    this->phiSpinBox->setVisible(!active);
    updateImageMenuActionVisibilties();
    emit checktStateChanged(this->localCheckBox->objectName());
    emit this->sendCommand(this->localCheckBox->objectName());

}

/**
 * @brief This Slot (de)activates computing the coptical depth at the wavelength
 *
 * @param active, this method getd the active-parameter from signal clicked(bool) of tauCheckBox.
 * This method emit the sendCommand signal.
 *
 * @sendCommand
 */

void ImageControllerWidget::activateTau(bool active){

    this->setTauToolTip(active);
    emit checktStateChanged(this->tauCheckBox->objectName());
    emit this->sendCommand(this->tauCheckBox->objectName());

}

/**
 * @brief This Slot sets contour
 *
 * @param active, this method gets the active-parameter from signal clicked(bool) of contourCheckBox.
 * This method emit the sendCommand signal.
 *
 * @sendCommand
 */

void ImageControllerWidget::activateContour(bool active){

    this->setContourToolTip(active);
    emit checktStateChanged(this->contourCheckBox->objectName());
    emit this->sendCommand(this->contourCheckBox->objectName());

}
/**
 * @brief This method set the shortcuts.
 *
 */
void ImageControllerWidget::setShortcuts(){
    this->starCheckBox->setShortcut(QKeySequence("Ctrl+s"));
    this->contourCheckBox->setShortcut(QKeySequence("Ctrl+c"));
    this->previewCheckBox->setShortcut(QKeySequence("Ctrl+w"));
    this->tauCheckBox->setShortcut(QKeySequence("Ctrl+t"));
    this->secondOrderCheckBox->setShortcut(QKeySequence("Ctrl+n"));
    this->linearCheckBox->setShortcut(QKeySequence("Ctrl+l"));
    this->localCheckBox->setShortcut(QKeySequence("alt+l"));
    this->colorbarCheckBox->setShortcut(QKeySequence("Ctrl+b"));
    this->invertColorsCheckBox->setShortcut(QKeySequence("Ctrl+i"));
    this->phiShortCut=new QShortcut(QKeySequence("ALT+p"),this);
    this->connect(this->phiShortCut,SIGNAL(activated()),this->phiSpinBox,SLOT(setFocus()));
    this->inclinationShortCut=new QShortcut(QKeySequence("ALT+a"),this);
    this->connect(this->inclinationShortCut,SIGNAL(activated()),this->inclinationSpinBox,SLOT(setFocus()));
    this->redLambdaShortCut =new QShortcut(QKeySequence("ALT+r"),this);
    this->connect(this->redLambdaShortCut,SIGNAL(activated()),this->redLambdaValueLineEdit,SLOT(setFocus()));
    this->blueLambdaShortCut =new QShortcut(QKeySequence("ALT+b"),this);
    this->connect(this->blueLambdaShortCut,SIGNAL(activated()),this->blueLambdaValueLineEdit,SLOT(setFocus()));
    this->greenLambdaShortCut =new QShortcut(QKeySequence("ALT+e"),this);
    this->connect(this->greenLambdaShortCut,SIGNAL(activated()),this->greenLambdaValueLineEdit,SLOT(setFocus()));
    this->pipeCheckBox->setShortcut(QKeySequence("Ctrl+p"));
    this->antiAliasingCheckBox->setShortcut(QKeySequence("Ctrl+g"));
    this->dopplerCatchCheckBox->setShortcut(QKeySequence("Ctrl+e"));
    this->upPushButton->setShortcut(QKeySequence(Qt::Key_Up));
    this->downPushButton->setShortcut(QKeySequence(Qt::Key_Down));
    this->leftPushButton->setShortcut(QKeySequence(Qt::Key_Left));
    this->rightPushButton->setShortcut(QKeySequence(Qt::Key_Right));
}
void ImageControllerWidget::sizeUnitChanged(QAction* action){
    SetUp::unit unit;
    unit=SetUp::CM;
    if(action->text()=="pc")unit=SetUp::PC;
    else if(action->text()=="au")unit=SetUp::AU;
    this->setSizeUnit(unit);
    this->setup->setImageUnit(unit);
    this->setSizeLineEdit();
}
void ImageControllerWidget::selectDebugger(){
    this->middleTabWidget->setCurrentIndex(this->middleTabWidget->indexOf(this->debuggerGroupBox));
    emit this->scollDebuggerToBottomPushButton->animateClick();
}

/**
 *  @brief This is the setup method.
 *
 *  @details This method sets the GUI-elements's attribute icons,
 *  labels, shortcuts and slots.
 *
 *  @see setIcons()
 *  @see translate()
 *  @see setSlots()
 *  @see setShortcuts()
 *  @see setUpColorTableComboBoxValues()
 *  @see setUpDirectionComboBoxValues()
 *
 */
void ImageControllerWidget::setUp(){
    this->setUpLocalTab();
    this->setIcons();
    this->setStyleSheets();
    this->setToolTips();
    this->setUpPositionComboBoxValues();
    this->setUpContrastSpinBox();
    this->setUpSaturationSpinBbox();
    this->setUpNumberOfContoursSpinBox();
    this->setUpNumberOfPixelsSpinBox();
    this->setUpLineSpinBox();
    this->setUpVelocityBoxSpinBbox();
    this->setUpMoleculeSpinBox();
    this->setUpScaleSpinBox();
    this->setUpViewangleSpinBox();
    this->setUpPositionAngleSpinBoxSpinBbox();
    this->setSizeUnit(this->setup->getImageUnit());
    this->sizeRegExp="[^0-9.e+-]+";
    this->lambdaRegExp="[^0-9.-]+";
    this->lastSizeValue=this->setup->calculateImageSize();
    this->sizeLineEdit->setText(QString::number( lastSizeValue,'e')+"  "+getSizeUnit());
    this->sizeLineEdit->setCursorPosition(0);
    this->lastobsY=-2*lastSizeValue;
    this->observerposYLineEdit->setText(QString::number(lastobsY,'e')+"  "+getSizeUnit());

    this->lastPointY=0;
    this->lastPointX=0;
    this->lastobsZ=0;
    this->lastPointZ=0;
    this->lastobsX=0;
    this->pointXLineEdit->setText(QString::number(lastPointX,'e')+"  "+getSizeUnit());
    this->pointZLineEdit->setText(QString::number(lastPointZ,'e')+"  "+getSizeUnit());
    this->pointYLineEdit->setText(QString::number(lastPointY,'e')+"  "+getSizeUnit());

    this->observerposXLineEdit->setText(QString::number(lastobsX,'e')+"  "+getSizeUnit());
    this->observerposZLineEdit->setText(QString::number(lastobsZ,'e')+"  "+getSizeUnit());
    this->unzoomImagePushButton->setEnabled(true);
    this->scaleSpinBox->setValue(1.0);


    this->setUpPhiSpinBoxSpinBbox();
    this->setUpInclinationSpinBoxSpinBbox();

    this->setLambdaUnit(trUtf8("\u03BCm"));


    this->setSlots();
    this->lastBlueColorTuneValue=1.0;
    this->lastRedColorTuneValue=1.0;
    this->lastGreenColorTuneValue=1.0;
    this->redColorTuneValueLineEdit->setText(QString::number(this->lastRedColorTuneValue));
    this->blueColorTuneValueLineEdit->setText(QString::number(this->lastBlueColorTuneValue));
    this->greenColorTuneValueLineEdit->setText(QString::number(this->lastGreenColorTuneValue));

    if(this->setup->IsLambaRead()){

        this->setUpLambdaIndexSlider();
    }
    this->translate();

}
void ImageControllerWidget::setSizeLineEdit(){
    this->lastSizeValue=this->setup->calculateImageSize();
    this->sizeLineEdit->setText(QString::number( this->lastSizeValue,'e')+"  "+getSizeUnit());
    this->sizeLineEdit->setCursorPosition(0);
    this->lastobsY=-2*lastSizeValue;
    this->observerposYLineEdit->setText(QString::number(lastobsY,'e')+"  "+getSizeUnit());
    QString txt=this->observerposXLineEdit->text();
    PhysicalConstants *constants=PhysicalConstants::getInstance();
    if(txt.contains("CM")){
        if(getSizeUnit()=="PC"){
            double parsec=constants->getParsec();
            this->lastPointX=lastPointX/parsec;
            this->lastPointY=lastPointY/parsec;
            this->lastPointZ=lastPointZ/parsec;
            this->lastobsX=this->lastobsX/parsec;
            this->lastobsY=this->lastobsY/parsec;
            this->lastobsZ=this->lastobsZ/parsec;

            this->pointXLineEdit->setText(QString::number(lastPointX,'e')+"  "+getSizeUnit());
            this->pointYLineEdit->setText(QString::number(lastPointY,'e')+"  "+getSizeUnit());
            this->pointZLineEdit->setText(QString::number(lastPointZ,'e')+"  "+getSizeUnit());

            this->observerposXLineEdit->setText(QString::number(lastobsX,'e')+"  "+getSizeUnit());
            this->observerposZLineEdit->setText(QString::number(lastobsZ,'e')+"  "+getSizeUnit());
        }else if(getSizeUnit()=="AU"){
            double au=constants->getAstronomicalUnit();
            this->lastPointX=lastPointX/au;
            this->lastPointY=lastPointY/au;
            this->lastPointZ=lastPointZ/au;
            this->lastobsX=this->lastobsX/au;
            this->lastobsY=this->lastobsY/au;
            this->lastobsZ=this->lastobsZ/au;

            this->pointXLineEdit->setText(QString::number(lastPointX,'e')+"  "+getSizeUnit());
            this->pointYLineEdit->setText(QString::number(lastPointY,'e')+"  "+getSizeUnit());
            this->pointZLineEdit->setText(QString::number(lastPointZ,'e')+"  "+getSizeUnit());

            this->observerposXLineEdit->setText(QString::number(lastobsX,'e')+"  "+getSizeUnit());
            this->observerposZLineEdit->setText(QString::number(lastobsZ,'e')+"  "+getSizeUnit());

        }
    }else if(txt.contains("AU")){
        if(getSizeUnit()=="PC"){
            double parsec=constants->getParsec();
            double au=constants->getAstronomicalUnit();
            double erg=au/parsec;
            this->lastPointX=lastPointX*erg;
            this->lastPointY=lastPointY*erg;
            this->lastPointZ=lastPointZ*erg;
            this->lastobsX=this->lastobsX*erg;
            this->lastobsY=this->lastobsY*erg;
            this->lastobsZ=this->lastobsZ*erg;

            this->pointXLineEdit->setText(QString::number(lastPointX,'e')+"  "+getSizeUnit());
            this->pointYLineEdit->setText(QString::number(lastPointY,'e')+"  "+getSizeUnit());
            this->pointZLineEdit->setText(QString::number(lastPointZ,'e')+"  "+getSizeUnit());

            this->observerposXLineEdit->setText(QString::number(lastobsX,'e')+"  "+getSizeUnit());
            this->observerposZLineEdit->setText(QString::number(lastobsZ,'e')+"  "+getSizeUnit());
        }else if(getSizeUnit()=="CM"){
            double au=constants->getAstronomicalUnit();
            this->lastPointX=lastPointX*au;
            this->lastPointY=lastPointY*au;
            this->lastPointZ=lastPointZ*au;
            this->lastobsX=this->lastobsX*au;
            this->lastobsY=this->lastobsY*au;
            this->lastobsZ=this->lastobsZ*au;

            this->pointXLineEdit->setText(QString::number(lastPointX,'e')+"  "+getSizeUnit());
            this->pointYLineEdit->setText(QString::number(lastPointY,'e')+"  "+getSizeUnit());
            this->pointZLineEdit->setText(QString::number(lastPointZ,'e')+"  "+getSizeUnit());

            this->observerposXLineEdit->setText(QString::number(lastobsX,'e')+"  "+getSizeUnit());
            this->observerposZLineEdit->setText(QString::number(lastobsZ,'e')+"  "+getSizeUnit());

        }
    }else if(txt.contains("PC")){
        if(getSizeUnit()=="AU"){
            double parsec=constants->getParsec();
            double au=constants->getAstronomicalUnit();
            double erg=parsec/au;
            this->lastPointX=lastPointX*erg;
            this->lastPointY=lastPointY*erg;
            this->lastPointZ=lastPointZ*erg;
            this->lastobsX=this->lastobsX*erg;
            this->lastobsY=this->lastobsY*erg;
            this->lastobsZ=this->lastobsZ*erg;

            this->pointXLineEdit->setText(QString::number(lastPointX,'e')+"  "+getSizeUnit());
            this->pointYLineEdit->setText(QString::number(lastPointY,'e')+"  "+getSizeUnit());
            this->pointZLineEdit->setText(QString::number(lastPointZ,'e')+"  "+getSizeUnit());

            this->observerposXLineEdit->setText(QString::number(lastobsX,'e')+"  "+getSizeUnit());
            this->observerposZLineEdit->setText(QString::number(lastobsZ,'e')+"  "+getSizeUnit());
        }else if(getSizeUnit()=="CM"){
            double parsec=constants->getParsec();
            this->lastPointX=lastPointX*parsec;
            this->lastPointY=lastPointY*parsec;
            this->lastPointZ=lastPointZ*parsec;
            this->lastobsX=this->lastobsX*parsec;
            this->lastobsY=this->lastobsY*parsec;
            this->lastobsZ=this->lastobsZ*parsec;

            this->pointXLineEdit->setText(QString::number(lastPointX,'e')+"  "+getSizeUnit());
            this->pointYLineEdit->setText(QString::number(lastPointY,'e')+"  "+getSizeUnit());
            this->pointZLineEdit->setText(QString::number(lastPointZ,'e')+"  "+getSizeUnit());

            this->observerposXLineEdit->setText(QString::number(lastobsX,'e')+"  "+getSizeUnit());
            this->observerposZLineEdit->setText(QString::number(lastobsZ,'e')+"  "+getSizeUnit());

        }
    }
    this->pointXLineEdit->setCursorPosition(0);
    this->pointYLineEdit->setCursorPosition(0);
    this->pointZLineEdit->setCursorPosition(0);
    this->observerposXLineEdit->setCursorPosition(0);
    this->observerposYLineEdit->setCursorPosition(0);
    this->observerposZLineEdit->setCursorPosition(0);

}


/**
 *  @brief This method sets lambda's size.
 *
 */
void ImageControllerWidget::setLambdaSize(int size){
    this->lambdaSize=size;
}

/**
 *  @brief This method returns lambda's size.
 *  @return lambdaSize
 */
int ImageControllerWidget::getLambdaSize(){
    return this->lambdaSize;
}

/**
 *  @brief This method returns lambda's array.
 *  @return lambdaArray
 */
QVector<double> ImageControllerWidget::getLambdaArray(){
    return this->lambdaArray;
}
/**
 *  @brief This method sets lambda's array
 *  @param lambda's array.
 */
void ImageControllerWidget::setLambdaArray(QVector<double> array){
    this->lambdaArray=array;
}




/**
 *  @brief This method sets the size-unit's string
 *
 *  @param sizeUnit is a QString-parameter.
 */
void ImageControllerWidget::setSizeUnit(SetUp::unit iunit){
    QString unit;
    if(iunit==SetUp::CM){
        unit="cm";
    }else if(iunit==SetUp::AU){
        unit="au";
    }else if(iunit==SetUp::PC){
        unit="pc";
    }
    this->sizeUnit=unit;
}

/**
 *  @brief This method sets the lambda-unit's string
 *
 *  @param sizeUnit is a QString-parameter.
 */
void ImageControllerWidget::setLambdaUnit(QString unit){
    this->lambdaUnit=unit;
}

/**
 * @brief This method returns the sizeunit's string
 *
*/
QString ImageControllerWidget::getSizeUnit(){
    return this->sizeUnit;
}
/**
 * @brief This method returns the lambdaunit's string
 *
*/
QString ImageControllerWidget::getLambdaUnit(){
    return this->lambdaUnit;
}

/**
 *  This method sets the setting of saturationSpinBox.
 *
 *
 *
 */
void ImageControllerWidget::setUpSaturationSpinBbox(){
    this->saturationSpinBox->setDecimals(2);
    this->saturationSpinBox->setSuffix("%");
    this->saturationSpinBox->setMinimum(this->saturationSpinBoxMinimum);
    this->saturationSpinBox->setMaximum(this->saturationSpinBoxMaximum);
    this->saturationSpinBox->setSingleStep(this->saturationSpinBoxStep);
    this->saturationSpinBox->setValue(saturationSpinBoxMaximum);
}
/**
 *  This method sets the setting of PositionAngleSpinBox.
 *
 *
 *
 */
void ImageControllerWidget::setUpPositionAngleSpinBoxSpinBbox(){
    this->positionAngleSpinBox->setDecimals(2);
    this->positionAngleSpinBox->setSuffix(trUtf8("\u00BA"));
    this->positionAngleSpinBox->setMinimum(this->positionAngleSpinBoxMinimum);
    this->positionAngleSpinBox->setMaximum(this->positionAngleSpinBoxMaximum);
    this->positionAngleSpinBox->setSingleStep(this->positionAngleSpinBoxStep);
}
/**
 *  This method sets the setting of inclinationSpinBox.
 *
 *
 *
 */
void ImageControllerWidget::setUpInclinationSpinBoxSpinBbox(){
    this->inclinationSpinBox->setDecimals(2);
    this->inclinationSpinBox->setSuffix(trUtf8("\u00BA"));
    this->inclinationSpinBox->setMinimum(this->inclinationSpinBoxMinimum);
    this->inclinationSpinBox->setMaximum(this->inclinationSpinBoxMaximum);
    this->inclinationSpinBox->setSingleStep(this->inclinationSpinBoxStep);
    this->inclinationSpinBox->setValue(this->inclinationSpinBoxDefaultValue);
}
/**
 *  This method sets the setting of phiSpinBox.
 *
 *
 *
 */
void ImageControllerWidget::setUpPhiSpinBoxSpinBbox(){
    this->phiSpinBox->setDecimals(2);
    this->phiSpinBox->setSuffix(trUtf8("\u00BA"));
    this->phiSpinBox->setMinimum(this->phiSpinBoxMinimum);
    this->phiSpinBox->setMaximum(this->phiSpinBoxMaximum);
    this->phiSpinBox->setSingleStep(this->phiSpinBoxStep);
    this->phiSpinBox->setValue(this->phiSpinBoxDefaultValue);
}
/**
 *  This method sets the setting of VelocitySpinBox.
 *
 *
 *
 */
void ImageControllerWidget::setUpVelocityBoxSpinBbox(){
    this->velocitySpinBox->setDecimals(6);
    this->velocitySpinBox->setSuffix(" km/s");
    this->velocitySpinBox->setMinimum(this->velocitySpinBoxMinimum);
    this->velocitySpinBox->setMaximum(this->velocitySpinBoxMaximum);
    this->velocitySpinBox->setValue(this->velocitySpinBoxDefaultValue);
    this->velocitySpinBox->setSingleStep(this->velocitySpinBoxStep);
}
/**
 *  This method sets the setting of scaleSpinBox.
 *
 *
 *
 */
void ImageControllerWidget::setUpScaleSpinBox(){
    this->scaleSpinBox->setDecimals(6);
    this->scaleSpinBox->setMinimum(this->scaleSpinBoxMinimum);
    this->scaleSpinBox->setMaximum(this->scaleSpinBoxMaximum);
    this->scaleSpinBox->setValue(this->scaleSpinBoxDefaultValue);
    this->scaleSpinBox->setSingleStep(this->scaleSpinBoxStep);
}
/**
 *  This method sets the setting of viewangleSpinBox.
 *
 *
 *
 */
void ImageControllerWidget::setUpViewangleSpinBox(){
    this->viewangleSpinBox->setDecimals(6);
    this->viewangleSpinBox->setSuffix(trUtf8(" \u33AD"));
    this->viewangleSpinBox->setMinimum(this->viewangleSpinBoxMinimum);
    this->viewangleSpinBox->setMaximum(this->viewangleSpinBoxMaximum);
    this->viewangleSpinBox->setValue(this->viewangleSpinBoxDefaultValue);
    this->viewangleSpinBox->setSingleStep(this->viewangleSpinBoxStep);
}

/**
 *  This method sets the setting of contrastSpinBox.
 *
 *
 *
 */
void ImageControllerWidget::setUpContrastSpinBox(){
    this->contrastSpinBox->setMinimum(this->contrastSpinBoxMinimum);
    this->contrastSpinBox->setMaximum(this->contrastSpinBoxMaximum);
    this->contrastSpinBox->setSingleStep(this->contrastSpinBoxStep);
    this->contrastSpinBox->setValue(6.0);
    this->setConstrastSpinBoxEnabled(!this->linearCheckBox->isChecked());
}
/**
 *  This method sets the setting of numberOfContoursSpinBox.
 *
 *
 *
 */
void ImageControllerWidget::setUpNumberOfContoursSpinBox(){
    this->numberOfContoursSpinBox->setMinimum(this->numberOfContourSpinBoxMinimum);
    this->numberOfContoursSpinBox->setMaximum(this->numberOfContourSpinBoxMaximum);
    this->numberOfContoursSpinBox->setValue(20.);
    this->numberOfContoursSpinBox->setSingleStep(this->numberOfContourSpinBoxStep);
}
/**
 *  This method sets the setting of lineSpinBox.
 *
 *
 *
 */
void ImageControllerWidget::setUpLineSpinBox(){
    this->lineSpinBox->setMinimum(1);
    this->lineSpinBox->setMaximum(10000);
    this->lineSpinBox->setSingleStep(1);
}
/**
 *  This method sets the setting of moleculeSpinBox.
 *
 *
 *
 */
void ImageControllerWidget::setUpMoleculeSpinBox(){

}


/**
 *  This method sets the setting of numberOfPixelsSpinBox.
 *
 *
 *
 */
void ImageControllerWidget::setUpNumberOfPixelsSpinBox(){
    this->numberOfPixelsSpinBox->setMinimum(this->numberOfPixelsSpinBoxMinimum);
    this->numberOfPixelsSpinBox->setMaximum(this->numberOfPixelsSpinBoxMaximum);
    this->numberOfPixelsSpinBox->setSingleStep(this->numberOfPixelsSpinBoxStep);
    this->numberOfPixelsSpinBox->setValue(this->numberOfPixelsSpinBoxDefaultValue);
}

/**
 * This method sets the values of colorTableComboxBoxValues.
 *
 */
void ImageControllerWidget::setUpColorTableComboBoxValues(){
    this->colorTableComboxBox->clear();
    this->colorTableComboxBox->setIconSize(QSize(comboBoxIconWidth,comboBoxIconWidth));
    for (int var = 1; var <= this->colorLookUpTables.size(); ++var) {
        this->colorTableComboxBox->addItem(QIcon(QString(":/images/%1.png").arg(
                                                     this->colorTableSlider->objectName())),QString::number(var)+" "
                                           +this->colorLookUpTableNames.at(var-1));
    }
    this->colorTableComboxBox->setToolTip(colorTableComboxBox->itemText(this->colorTableComboxBox->currentIndex()));
    this->colorTableSlider->setToolTip(colorTableComboxBox->itemText(this->colorTableComboxBox->currentIndex()));
}
/**
 * This method sets the values of directionComboxBoxValues.
 *
 */
void ImageControllerWidget::setUpPositionComboBoxValues(){
    this->positionComboxBox->setIconSize(QSize(comboBoxIconWidth,comboBoxIconWidth));
    this->positionComboxBox->addItem(qApp->style()->standardIcon(QStyle::SP_ArrowRight),"Right");
    this->positionComboxBox->addItem(qApp->style()->standardIcon(QStyle::SP_ArrowLeft),"Left");
    this->positionComboxBox->addItem(qApp->style()->standardIcon(QStyle::SP_ArrowDown),"Bottom");
    this->positionComboxBox->addItem(qApp->style()->standardIcon(QStyle::SP_ArrowUp),"Top");
    this->positionComboxBox->setToolTip(this->positionComboxBox->itemText(this->positionComboxBox->currentIndex()));
}



/**
 * @brief This method set the slots.
 *
 */
void ImageControllerWidget::setSlots(){
    this->connect(this->starCheckBox,SIGNAL(toggled(bool)),this,SLOT(activateStar(bool)));
    this->connect(this->secondOrderCheckBox,SIGNAL(toggled(bool)),this,SLOT(activateSecondOrder(bool)));
    this->connect(this->absoluteScaleCheckBox,SIGNAL(toggled(bool)),this,SLOT(activateAbsoluteScale(bool)));
    this->connect(this->tauCheckBox,SIGNAL(toggled(bool)),this,SLOT(activateTau(bool)));
    this->connect(this->previewCheckBox,SIGNAL(toggled(bool)),this,SLOT(activatePreview(bool)));
    this->connect(this->localCheckBox,SIGNAL(toggled(bool)),this,SLOT(activateLocal(bool)));
    this->connect(this->invertColorsCheckBox,SIGNAL(toggled(bool)),this,SLOT(activateInvertColors(bool)));
    this->connect(this->contourCheckBox,SIGNAL(toggled(bool)),this,SLOT(activateContour(bool)));
    this->connect(this->linearCheckBox,SIGNAL(toggled(bool)),this,SLOT(activateLinear(bool)));
    this->connect(this->colorbarCheckBox,SIGNAL(toggled(bool)),this,SLOT(activateColorbar(bool)));
    this->connect(this->pipeCheckBox,SIGNAL(toggled(bool)),this,SLOT(activatePipe(bool)));
    this->connect(this->colorTableSlider,SIGNAL(valueChanged(int)),this,SLOT(colorTableSliderValueChanged(int)));
    this->connect(this->redLambdaIndexSlider,SIGNAL(valueChanged(int)),this,SLOT(redLambdaIndexSliderValueChanged(int)));
    this->connect(this->greenLambdaIndexSlider,SIGNAL(valueChanged(int)),this,SLOT(greenLambdaIndexSliderValueChanged(int)));
    this->connect(this->blueLambdaIndexSlider,SIGNAL(valueChanged(int)),this,SLOT(blueLambdaIndexSliderValueChanged(int)));
    this->connect(this->colorTableComboxBox,SIGNAL(currentIndexChanged(int)),this,SLOT(colorTableComboxBoxValueChanged(int)));
    this->connect(this->positionComboxBox,SIGNAL(currentIndexChanged(int)),this,SLOT(positionComboxBoxValueChanged(int)));
    this->connect(this->contrastSpinBox,SIGNAL(editingFinished()),this,SLOT(contrastSpinBoxValueChanged()));
    this->connect(this->pointXLineEdit,SIGNAL(editingFinished()),this,SLOT(pointXLineEditValueChanged()));
    this->connect(this->pointYLineEdit,SIGNAL(editingFinished()),this,SLOT(pointYLineEditValueChanged()));
    this->connect(this->pointZLineEdit,SIGNAL(editingFinished()),this,SLOT(pointZLineEditValueChanged()));
    this->connect(this->observerposXLineEdit,SIGNAL(editingFinished()),this,SLOT(observerposXLineEditValueChanged()));
    this->connect(this->observerposYLineEdit,SIGNAL(editingFinished()),this,SLOT(observerposYLineEditValueChanged()));
    this->connect(this->observerposZLineEdit,SIGNAL(editingFinished()),this,SLOT(observerposZLineEditValueChanged()));
    this->connect(this->saturationSpinBox,SIGNAL(editingFinished()),this,SLOT(saturationSpinBoxValueChanged()));
    this->connect(this->numberOfContoursSpinBox,SIGNAL(editingFinished()),this,SLOT(numberOfContoursSpinBoxValueChanged()));
    this->connect(this->positionAngleSpinBox,SIGNAL(editingFinished()),this,SLOT(positionAngleSpinBoxValueChanged()));
    this->connect(this->numberOfPixelsSpinBox,SIGNAL(editingFinished()),this,SLOT(numberOfPixelsSpinBoxValueChanged()));
    this->connect(this->sizeLineEdit,SIGNAL(editingFinished()),this,SLOT(sizeLineEditValueChanged()));
    this->connect(this->redLambdaValueLineEdit,SIGNAL(editingFinished()),this,SLOT(redLambdaValueLineEditValueChanged()));
    this->connect(this->greenLambdaValueLineEdit,SIGNAL(editingFinished()),this,SLOT(greenLambdaValueLineEditValueChanged()));
    this->connect(this->blueLambdaValueLineEdit,SIGNAL(editingFinished()),this,SLOT(blueLambdaValueLineEditValueChanged()));
    this->connect(this->phiSpinBox,SIGNAL(editingFinished()),this,SLOT(phiSpinBoxValueChanged()));
    this->connect(this->scaleSpinBox,SIGNAL(editingFinished()),this,SLOT(scaleSpinBoxValueChanged()));
    this->connect(this->viewangleSpinBox,SIGNAL(editingFinished()),this,SLOT(viewangleSpinBoxValueChanged()));
    this->connect(this->inclinationSpinBox,SIGNAL(editingFinished()),this,SLOT(inclinationSpinBoxValueChanged()));
    this->connect(this->unzoomImagePushButton,SIGNAL(clicked()),this,SLOT(unzoomImagePushButtonClicked()));
    this->connect(this->renderImagePushButton,SIGNAL(clicked()),this,SLOT(renderImagePushButtonClicked()));
    this->connect(this->renderSpectrumPushButton,SIGNAL(clicked()),this,SLOT(renderSpectrumPushButtonClicked()));
    this->connect(this->rightPushButton,SIGNAL(clicked()),this,SLOT(rightPushButtonClicked()));
    this->connect(this->leftPushButton,SIGNAL(clicked()),this,SLOT(leftPushButtonClicked()));
    this->connect(this->upPushButton,SIGNAL(clicked()),this,SLOT(upPushButtonClicked()));
    this->connect(this->downPushButton,SIGNAL(clicked()),this,SLOT(downPushButtonClicked()));
    this->connect(this->clearDebuggerPushButton,SIGNAL(clicked()),this,SLOT(clearDebugger()));
    this->connect(this->scollDebuggerToBottomPushButton,SIGNAL(clicked()),this,SLOT(scrollDebuggerToBottom()));
    this->connect(this->scollDebuggerToTopPushButton,SIGNAL(clicked()),this,SLOT(scrollDebuggerToTop()));
    this->connect(this->middleTabWidget,SIGNAL(currentChanged(int)),this,SLOT(middleTabWidgetCurrentIndex(int)));
    this->connect(this->antiAliasingCheckBox,SIGNAL(toggled(bool)),this,SLOT(activateAntiAliasing(bool)));
    this->connect(this->relobCheckbox,SIGNAL(toggled(bool)),this,SLOT(activateRelob(bool)));
    this->connect(this->dopplerCatchCheckBox,SIGNAL(toggled(bool)),this,SLOT(activateDopplerCatching(bool)));
    this->connect(this->showRadmc3dPushButton,SIGNAL(clicked()),this,SLOT(showRadmc3dPushButtonClicked()));
    this->connect(this->redColorTuneValueLineEdit,SIGNAL(editingFinished()),this,SLOT(redColorTuneValueLineEditValueChanged()));
    this->connect(this->greenColorTuneValueLineEdit,SIGNAL(editingFinished()),this,SLOT(greenColorTuneValueLineEditValueChanged()));
    this->connect(this->blueColorTuneValueLineEdit,SIGNAL(editingFinished()),this,SLOT(blueColorTuneValueLineEditValueChanged()));
    this->connect(this->moleculeComboBox,SIGNAL(currentIndexChanged(int)),this,SLOT(moleculeComboBoxValueChanged(int)));
    this->connect(this->lineSpinBox,SIGNAL(editingFinished()),this,SLOT(lineSpinBoxValueChanged()));
    this->connect(this->velocitySpinBox,SIGNAL(editingFinished()),this,SLOT(velocitySpinBoxValueChanged()));


}

void ImageControllerWidget::showRadmc3dPushButtonClicked(){
    emit sendCommand(this->showRadmc3dPushButton->objectName());
}

/**
 * @brief This method set the GUI-elements's labels.
 *
 * @details This method should be called, if the GUI-language is changed, so that the translations will be loaded.
 */
void ImageControllerWidget::translate(){
    this->setToolTips();
    this->tauCheckBox->setText(trUtf8("&\u03C4"));
    this->previewCheckBox->setToolTip(trUtf8("preview"));
    this->localCheckBox->setToolTip(trUtf8("Local Observer"));
    this->relobCheckbox->setText(trUtf8("Relative &Observer"));
    this->previewCheckBox->setText(trUtf8("Previe&w"));
    this->antiAliasingCheckBox->setText(trUtf8("antiAliasing"));
    this->antiAliasingCheckBox->setToolTip(trUtf8("antiAliasing"));
    this->dopplerCatchCheckBox->setText(trUtf8("Doppler Catching"));
    this->dopplerCatchCheckBox->setToolTip(trUtf8("Doppler Catching method"));
    this->showRadmc3dPushButton->setText(trUtf8("Show RADMC3D's logfile"));
    this->contourCheckBox->setText(trUtf8("&Contour"));
    this->secondOrderCheckBox->setText(trUtf8("&ND"));
    this->starCheckBox->setText(trUtf8("&Star"));
    this->linearCheckBox->setText(trUtf8("&Linear"));
    this->colorbarCheckBox->setText(trUtf8("Color&bar"));
    this->colorTableLabel->setText(trUtf8("&Colortable"));
    this->invertColorsCheckBox->setText(trUtf8("&Invert colors"));
    this->pipeCheckBox->setText(trUtf8("&Pipe"));
    this->positionLabel->setText(trUtf8("Posi&tion"));
    this->contrastLabel->setTextFormat(Qt::RichText);
    this->contrastLabel->setText(trUtf8("<html><img src=\":/images/%1.png\" width=\"%2\" height=\"%3\" >Ma&xLog</html>")
                                 .arg(this->contrastSpinBox->objectName(),QString::number(this->labelIconWidth)
                                      ,QString::number(this->labelIconHeight)));
    this->saturationLabel->setTextFormat(Qt::RichText);
    this->saturationLabel->setText(trUtf8("<html><img src=\":/images/%1.png\" width=\"%2\" height=\"%3\" >Sat&uration</html>")
                                   .arg(this->saturationSpinBox->objectName(),QString::number(this->labelIconWidth)
                                        ,QString::number(this->labelIconHeight)));
    this->numberOfContoursLabel->setTextFormat(Qt::RichText);
    this->numberOfContoursLabel->setText(trUtf8("<html><img src=\":/images/%1.png\" width=\"%2\" height=\"%3\" >Number &of contours</html>")
                                         .arg(this->numberOfContoursSpinBox->objectName(),QString::number(this->labelIconWidth)
                                              ,QString::number(this->labelIconHeight)));

    this->positionAngleLabel->setTextFormat(Qt::RichText);
    this->positionAngleLabel->setText(trUtf8("<html><img src=\":/images/%1.png\" width=\"%2\" height=\"%3\"><font size=3>Position's an&gle</html>")
                                      .arg(this->positionAngleSpinBox->objectName(),QString::number(this->labelIconWidth*1.25)
                                           ,QString::number(this->labelIconHeight*1.25)));

    this->numberOfPixelsLabel->setTextFormat(Qt::RichText);
    this->numberOfPixelsLabel->setText(trUtf8("<html><img src=\":/images/%1.png\" width=\"%2\" height=\"%3\" >&Number of pixels</html>")
                                       .arg(this->numberOfPixelsSpinBox->objectName(),QString::number(this->labelIconWidth)
                                            ,QString::number(this->labelIconHeight)));

    this->sizeLabel->setTextFormat(Qt::RichText);
    this->sizeLabel->setText(trUtf8("<html><img src=\":/images/%1.png\" width=\"%2\" height=\"%3\" >Si&ze</html>")
                             .arg(this->sizeLineEdit->objectName(),QString::number(this->labelIconWidth)
                                  ,QString::number(this->labelIconHeight)));

    this->inclinationLabel->setTextFormat(Qt::RichText);
    this->inclinationLabel->setText(trUtf8("<html><img src=\":/images/%1.png\" width=\"%2\" height=\"%3\" >&ang;</html>")
                                    .arg(this->inclinationSpinBox->objectName(),QString::number(this->labelIconWidth)
                                         ,QString::number(this->labelIconHeight)));

    this->phiLabel->setTextFormat(Qt::RichText);
    this->phiLabel->setText(trUtf8("<html><img src=\":/images/%1.png\" width=\"%2\" height=\"%3\" >%4</html>")
                            .arg(this->phiSpinBox->objectName(),QString::number(this->labelIconWidth)
                                 ,QString::number(this->labelIconHeight),trUtf8("\u03A6")));

    this->moleculeLabel->setTextFormat(Qt::RichText);
    this->moleculeLabel->setText(trUtf8("<html><img src=\":/images/%1.png\" width=\"%2\" height=\"%3\" >%4</html>")
                                 .arg(this->moleculeComboBox->objectName(),QString::number(this->labelIconWidth*2)
                                      ,QString::number(this->labelIconHeight*2),trUtf8("iMolecule")));
    this->lineLabel->setTextFormat(Qt::RichText);
    this->lineLabel->setText(trUtf8("<html><img src=\":/images/%1.png\" width=\"%2\" height=\"%3\" >%4</html>")
                             .arg(this->lineSpinBox->objectName(),QString::number(this->labelIconWidth*2)
                                  ,QString::number(this->labelIconHeight*2),trUtf8("iLine")));

    this->velocityLabel->setTextFormat(Qt::RichText);
    this->velocityLabel->setText(trUtf8("<html><img src=\":/images/%1.png\" width=\"%2\" height=\"%3\" >%4</html>")
                                 .arg(this->velocitySpinBox->objectName(),QString::number(this->labelIconWidth*2)
                                      ,QString::number(this->labelIconHeight*2),trUtf8("Velocity")));

    this->pointingLabel->setTextFormat(Qt::RichText);
    this->pointingLabel->setText(trUtf8("<html><img src=\":/images/%1.png\" width=\"%2\" height=\"%3\" >%4</html>")
                                 .arg(this->pointingLabel->objectName(),QString::number(this->labelIconWidth*2)
                                      ,QString::number(this->labelIconHeight*2),trUtf8("Pointing")));


    this->observerposLabel->setTextFormat(Qt::RichText);
    this->observerposLabel->setText(trUtf8("<html><img src=\":/images/%1.png\" width=\"%2\" height=\"%3\" >%4</html>")
                                    .arg(this->observerposLabel->objectName(),QString::number(this->labelIconWidth*2)
                                         ,QString::number(this->labelIconHeight*2),trUtf8("Observer Position")));
    this->viewangleLabel->setTextFormat(Qt::RichText);
    this->viewangleLabel->setText(trUtf8("<html><img src=\":/images/%1.png\" width=\"%2\" height=\"%3\"><font size=3>View's angle</html>")
                                  .arg(this->positionAngleSpinBox->objectName(),QString::number(this->labelIconWidth*1.25)
                                       ,QString::number(this->labelIconHeight*1.25)));

    this->localCheckBox->setText(trUtf8("Local Observer"));

    this->setLambdaLabelText();
    this->absoluteScaleCheckBox->setText(trUtf8("Absol&ute fixed brightness scale"));

    this->setShortcuts();
    this->setUpWhatsThis();
}

void ImageControllerWidget::setLambdaLabelText(){
    this->redLambdaLabel->setTextFormat(Qt::RichText);
    this->redLambdaLabel->setText(trUtf8("<html><img src=\":/images/%1.png\" width=\"%2\" height=\"%3\" >Red &lambda;(%4)</html>")
                                  .arg("lambda",QString::number(this->labelIconWidth)
                                       ,QString::number(this->labelIconHeight),QString::number(this->redLambdaIndexSlider->value())));
    this->blueLambdaLabel->setTextFormat(Qt::RichText);
    this->blueLambdaLabel->setText(trUtf8("<html><img src=\":/images/%1.png\" width=\"%2\" height=\"%3\" >Blue &lambda;(%4)</html>")
                                   .arg("lambda",QString::number(this->labelIconWidth)
                                        ,QString::number(this->labelIconHeight),QString::number(this->blueLambdaIndexSlider->value())));

    this->greenLambdaLabel->setTextFormat(Qt::RichText);
    this->greenLambdaLabel->setText(trUtf8("<html><img src=\":/images/%1.png\" width=\"%2\" height=\"%3\" >Green &lambda;(%4)</html>")
                                    .arg("lambda",QString::number(this->labelIconWidth)
                                         ,QString::number(this->labelIconHeight),QString::number(this->greenLambdaIndexSlider->value())));

}

void ImageControllerWidget::setUpWhatsThis(){

    this->starCheckBox->setWhatsThis(trUtf8("<html><font color=blue><img src=\":/images/%1.png\" width=\"%2\" height=\"%3\">"
                                            "Include stars in images.<br/>"
                                            "<img src=\":/images/no%1.png\" width=\"%2\" height=\"%3\">Do not include stars in images. Only the circumstellar /"
                                            "interstellar material is imaged as if a perfect coronograph is used.<br/>shortcut: %4"
                                            "</html>").arg(this->starCheckBox->objectName(),QString::number(this->whatsThisIconWidth),
                                                           QString::number(this->whatsThisIconHeight),this->starCheckBox->shortcut().toString()));

    this->secondOrderCheckBox->setWhatsThis(trUtf8("<html><font color=blue><img src=\":/images/%1.png\"  width=\"%2\" height=\"%3\">"
                                                   " Use second order integration ray tracing.<br/>"
                                                   "<img src=\":/images/no%1.png\"  width=\"%2\" height=\"%3\">Use first order integration ray tracing (See Section 9.8)."
                                                   "<br/>shortcut: %4"
                                                   "</html>").arg(this->secondOrderCheckBox->objectName(),QString::number(this->whatsThisIconWidth),
                                                                  QString::number(this->whatsThisIconHeight),this->secondOrderCheckBox->shortcut().toString()));


    this->tauCheckBox->setWhatsThis(trUtf8("<html><font color=blue><img src=\":/images/%1.png\"  width=\"%2\" height=\"%3\">"
                                           " If this option is set, then instead of ray-tracing a true image, the camera will "
                                           "compute the optical depth at the wavelength given by e.g. inu and puts this into an image "
                                           "output as if it were a true image. Can be useful for analysis of models.<br/>"
                                           "<img src=\":/images/no%1.png\"  width=\"%2\" height=\"%3\">reset to normal imaging mod."
                                           "<br/>shortcut: %4"
                                           "</html>").arg(this->tauCheckBox->objectName(),QString::number(this->whatsThisIconWidth),
                                                          QString::number(this->whatsThisIconHeight),this->tauCheckBox->shortcut().toString()));

    this->contourCheckBox->setWhatsThis(trUtf8("<html><font color=blue><img src=\":/images/%1.png\"  width=\"%2\" height=\"%3\">"
                                               "Overplot %4 contours over the image.<br/>"
                                               "<img src=\":/images/no%1.png\"  width=\"%2\" height=\"%3\">it does not plot contours ."
                                               "<br/>shortcut: %5"
                                               "</html>").arg(this->contourCheckBox->objectName(),QString::number(this->whatsThisIconWidth),
                                                              QString::number(this->whatsThisIconHeight),this->numberOfContoursSpinBox->text()
                                                              ,this->contourCheckBox->shortcut().toString()));

    this->previewCheckBox->setWhatsThis(trUtf8("<html><font color=blue><img src=\":/images/%1.png\"  width=\"%2\" height=\"%3\">"
                                               "Puts nrrefine (Specifies a maximum depth of refinement of the pixels) to a large value to assue flux conservation.<br/>"
                                               "<img src=\":/images/no%1.png\"  width=\"%2\" height=\"%3\"> Puts nrrefine to 0 so that each pixel of the image corre-"
                                               "sponds only to 1 ray. This is fast but not reliable and therefore not recommended(See section 9.6)."
                                               "<br/>shortcut: %4"
                                               "</html>").arg(this->previewCheckBox->objectName(),QString::number(this->whatsThisIconWidth),
                                                              QString::number(this->whatsThisIconHeight),this->previewCheckBox->shortcut().toString()));

    this->linearCheckBox->setWhatsThis(trUtf8("<html><font color=blue><img src=\":/images/%1.png\"  width=\"%2\" height=\"%3\">"
                                              "Linear scale(you can not use absoulte scale in this modus).<br/><br/><br/>"
                                              "<img src=\":/images/no%1.png\"  width=\"%2\" height=\"%3\">Logarithmic scale."
                                              "<br/>shortcut: %4"
                                              "</html>").arg(this->linearCheckBox->objectName(),QString::number(this->whatsThisIconWidth),
                                                             QString::number(this->whatsThisIconHeight),this->linearCheckBox->shortcut().toString()));

    this->localCheckBox->setWhatsThis(trUtf8("<html><font color=blue><img src=\":/images/%1.png\"  width=\"%2\" height=\"%3\">"
                                             "Local Observer is on.<br/><br/><br/>"
                                             "<img src=\":/images/no%1.png\"  width=\"%2\" height=\"%3\">Local Observer is off."
                                             "<br/>shortcut: %4"
                                             "</html>").arg(this->localCheckBox->objectName(),QString::number(this->whatsThisIconWidth),
                                                            QString::number(this->whatsThisIconHeight),this->localCheckBox->shortcut().toString()));

    this->absoluteScaleCheckBox->setWhatsThis(trUtf8("<html><font color=blue><img src=\":/images/%1.png\"  width=\"%2\" height=\"%3\">"
                                                     "Enable to set specifically the range of the image brightness in terms of the log of the image values by each color tune.<br/><br/><br/>"
                                                     "<img src=\":/images/no%1.png\"  width=\"%2\" height=\"%3\">inactive."
                                                     "<br/>shortcut: %4"
                                                     "</html>").arg(this->absoluteScaleCheckBox->objectName(),QString::number(this->whatsThisIconWidth),
                                                                    QString::number(this->whatsThisIconHeight),this->absoluteScaleCheckBox->shortcut().toString()));


    this->colorbarCheckBox->setWhatsThis(trUtf8("<html><font color=blue><img src=\":/images/%1.png\"  width=\"%2\" height=\"%3\">"
                                                "Shows colorbar.<br/><br/><br/>"
                                                "<img src=\":/images/no%1.png\"  width=\"%2\" height=\"%3\">Hides colorbar."
                                                "<br/>shortcut: %4"
                                                "</html>").arg(this->colorbarCheckBox->objectName(),QString::number(this->whatsThisIconWidth),
                                                               QString::number(this->whatsThisIconHeight),this->colorbarCheckBox->shortcut().toString()));

    this->invertColorsCheckBox->setWhatsThis(trUtf8("<html><font color=blue><img src=\":/images/%1.png\"  width=\"%2\" height=\"%3\">"
                                                    "Invert colors.<br/><br/><br/>"
                                                    "<img src=\":/images/no%1.png\"  width=\"%2\" height=\"%3\">Do not invert colors."
                                                    "<br/>shortcut: %4"
                                                    "</html>").arg(this->invertColorsCheckBox->objectName(),QString::number(this->whatsThisIconWidth),
                                                                   QString::number(this->whatsThisIconHeight),this->invertColorsCheckBox->shortcut().toString()));

    this->pipeCheckBox->setWhatsThis(trUtf8("<html><font color=blue><img src=\":/images/%1.png\"  width=\"%2\" height=\"%3\">"
                                            "Child modus is on.This prevents RADMC-3D from exiting after each main command is done. Instead, RADMC-"
                                            "3D will wait for further commands being given on the RADMC-3D internal command line."
                                            " This can be useful if multiple actions are to be taken, and the user does not want to wait for"
                                            " long file input reading.<br/><br/><br/>"
                                            "<img src=\":/images/no%1.png\"  width=\"%2\" height=\"%3\">Child modus is off."
                                            "<br/>shortcut: %4"
                                            "</html>").arg(this->pipeCheckBox->objectName(),QString::number(this->whatsThisIconWidth),
                                                           QString::number(this->whatsThisIconHeight),this->pipeCheckBox->shortcut().toString()));
    this->antiAliasingCheckBox->setWhatsThis(trUtf8("<html><font color=blue><img src=\":/images/%1.png\"  width=\"%2\" height=\"%3\">"
                                                    "Anti aliasing is on.<br/><br/><br/>"
                                                    "<img src=\":/images/no%1.png\"  width=\"%2\" height=\"%3\">anti aliasing is off."
                                                    "<br/>shortcut: %4"
                                                    "</html>").arg(this->antiAliasingCheckBox->objectName(),QString::number(this->whatsThisIconWidth),
                                                                   QString::number(this->whatsThisIconHeight),this->antiAliasingCheckBox->shortcut().toString()));
    this->dopplerCatchCheckBox->setWhatsThis(trUtf8("<html><font color=blue><img src=\":/images/%1.png\"  width=\"%2\" height=\"%3\">"
                                                    "Doppler catch method is on(See 7.12 section).<br/><br/><br/>"
                                                    "<img src=\":/images/no%1.png\"  width=\"%2\" height=\"%3\">doppler catch method is off."
                                                    "<br/>shortcut: %4"
                                                    "</html>").arg(this->dopplerCatchCheckBox->objectName(),QString::number(this->whatsThisIconWidth),
                                                                   QString::number(this->whatsThisIconHeight),this->dopplerCatchCheckBox->shortcut().toString()));

    this->contrastSpinBox->setWhatsThis(trUtf8("<html><font color=blue><img src=\":/images/%1.png\"  width=\"%2\" height=\"%3\">"
                                               "Set the maximum number of factors of 10 the log brightness color coding will span.(max=%4% min=%5%)."
                                               " Active just by logarithmic scale"
                                               "<br/>shortcut: %6"
                                               "</html>").arg(this->contrastSpinBox->objectName(),QString::number(this->whatsThisIconWidth),
                                                              QString::number(this->whatsThisIconHeight),QString::number(contrastSpinBoxMaximum)
                                                              ,QString::number(contrastSpinBoxMinimum),QKeySequence::mnemonic(this->contrastLabel->text()).toString()));

    this->saturationSpinBox->setWhatsThis(trUtf8("<html><font color=blue><img src=\":/images/%1.png\"  width=\"%2\" height=\"%3\">"
                                                 "Allows you to enhance the contrast of very weak emission regions by saturating bright regionsn.(default=%4% ,max=%4% min=%5%)."
                                                 "<br/>shortcut: %6"
                                                 "</html>").arg(this->saturationSpinBox->objectName(),QString::number(this->whatsThisIconWidth),
                                                                QString::number(this->whatsThisIconHeight),QString::number(saturationSpinBoxMaximum)
                                                                ,QString::number(saturationSpinBoxMinimum),QKeySequence::mnemonic(this->saturationLabel->text()).toString()));


    this->numberOfContoursSpinBox->setWhatsThis(trUtf8("<html><font color=blue><img src=\":/images/%1.png\"  width=\"%2\" height=\"%3\">"
                                                       "This spinbox sets the number of contours, that should be ploted.(max=%4 min=%5)."
                                                       "<br/>shortcut: %6"
                                                       "</html>").arg(this->numberOfContoursSpinBox->objectName(),QString::number(this->whatsThisIconWidth),
                                                                      QString::number(this->whatsThisIconHeight),QString::number(numberOfContourSpinBoxMaximum)
                                                                      ,QString::number(numberOfContourSpinBoxMinimum),
                                                                      QKeySequence::mnemonic(this->numberOfContoursLabel->text()).toString()));

    this->positionAngleSpinBox->setWhatsThis(trUtf8("<html><font color=blue><img src=\":/images/%1.png\"  width=\"%2\" height=\"%3\">"
                                                    "This rotates the camera itself around the (0,0) point in the image plane(max=%4\u00BA min=%5\u00BA)."
                                                    "<br/>shortcut: %6"
                                                    "</html>").arg(this->positionAngleSpinBox->objectName(),QString::number(this->whatsThisIconWidth),
                                                                   QString::number(this->whatsThisIconHeight),
                                                                   QString::number(positionAngleSpinBoxMaximum),
                                                                   QString::number(positionAngleSpinBoxMinimum),
                                                                   QKeySequence::mnemonic(this->positionAngleLabel->text()).toString()));

    this->numberOfPixelsSpinBox->setWhatsThis(trUtf8("<html><font color=blue><img src=\":/images/%1.png\"  width=\"%2\" height=\"%3\">"
                                                     "This spinBox sets the number of pixels (max=%4 min=%5)."
                                                     "<br/>shortcut: %6"
                                                     "</html>").arg(this->numberOfPixelsSpinBox->objectName(),QString::number(this->whatsThisIconWidth),
                                                                    QString::number(this->whatsThisIconHeight),
                                                                    QString::number(numberOfPixelsSpinBoxMaximum),QString::number(numberOfPixelsSpinBoxMinimum),
                                                                    QKeySequence::mnemonic(this->numberOfPixelsLabel->text()).toString()));

    this->sizeLineEdit->setWhatsThis(trUtf8("<html><font color=blue><img src=\":/images/%1.png\"  width=\"%2\" height=\"%3\">"
                                            "This spinBox sets the size(%4) of image.This image size is measured from the image center to the left or right"
                                            " or top or bottom. This gives always square images."
                                            "<br/>shortcut: %5"
                                            "</html>").arg(this->sizeLineEdit->objectName(),QString::number(this->whatsThisIconWidth),
                                                           QString::number(this->whatsThisIconHeight),this->getSizeUnit()
                                                           ,QKeySequence::mnemonic(this->sizeLabel->text()).toString()));


    this->inclinationSpinBox->setWhatsThis(trUtf8("<html><font color=blue><img src=\":/images/%1.png\"  width=\"%2\" height=\"%3\">"
                                                  "For the case when the camera is at infinity (i.e. at a large distance so"
                                                  " that no local perspective has to be taken into account) this inclination specifies the direction"
                                                  "toward which the camera for images and spectra is positioned. Incl = 0 means toward the"
                                                  " positive z-axis (in cartesian space), incl=90 means toward a position in the x-y-plane and"
                                                  "incl=180 means toward the negative z-axis. The angle is given in degrees."
                                                  "<br/>shortcut: %4"
                                                  "</html>").arg(this->inclinationSpinBox->objectName(),QString::number(this->whatsThisIconWidth),
                                                                 QString::number(this->whatsThisIconHeight),this->inclinationShortCut->key().toString()));

    this->phiSpinBox->setWhatsThis(trUtf8("<html><font color=blue><img src=\":/images/%1.png\"  width=\"%2\" height=\"%3\">"
                                          "Like incl, but now the remaining angle, also given in degrees."
                                          "Examples: incl=90 and phi=0 means that the observer is located at infinity toward the"
                                          "negative y axis; incl=90 and phi=90 means that the observer is located at infinity toward"
                                          "the negative x axis; incl=90 and phi=180 means that the observer is located at infinity"
                                          "toward the positive y axis (looking back in negative y direction). Rotation of the observer"
                                          "around the object around the z-axis goes in clockwise direction. The starting point of this"
                                          "rotation is such that for incl=0 and phi=0 the (x, y) in the image plane correspond to the"
                                          "(x, y) in the 3-D space, with x pointing toward the right and y pointing upward. Examples: ]"
                                          "if we fix the position of the observer at for instance incl=0 (i.e. we look at the object from"
                                          "the top from the positive z-axis at infinity downward), then increasing phi means rotating"
                                          "the object counter-clockwise in the image plane."
                                          "<br/>shortcut: %4"
                                          "</html>").arg(this->phiSpinBox->objectName(),QString::number(this->whatsThisIconWidth),
                                                         QString::number(this->whatsThisIconHeight),this->phiShortCut->key().toString()));

    this->redLambdaValueLineEdit->setWhatsThis(trUtf8("<html><font color=blue><img src=\":/images/%1.png\"  width=\"%2\" height=\"%3\">"
                                                      "Specify the the wavelength(red) from the wavelength micron.inp file"
                                                      " for which a ray-trace image should be made."
                                                      "<br/>shortcut: %4"
                                                      "</html>").arg("lambda",QString::number(this->whatsThisIconWidth),
                                                                     QString::number(this->whatsThisIconHeight),this->redLambdaShortCut->key().toString()));
    this->blueLambdaValueLineEdit->setWhatsThis(trUtf8("<html><font color=blue><img src=\":/images/%1.png\"  width=\"%2\" height=\"%3\">"
                                                       "Specify the the wavelength(blue) from the wavelength micron.inp file"
                                                       " for which a ray-trace image should be made."
                                                       "<br/>shortcut: %4"
                                                       "</html>").arg("lambda",QString::number(this->whatsThisIconWidth),
                                                                      QString::number(this->whatsThisIconHeight),this->blueLambdaShortCut->key().toString()));
    this->greenLambdaValueLineEdit->setWhatsThis(trUtf8("<html><font color=blue><img src=\":/images/%1.png\"  width=\"%2\" height=\"%3\">"
                                                        "Specify the the wavelength(green) from the wavelength micron.inp file"
                                                        " for which a ray-trace image should be made."
                                                        "<br/>shortcut: %4"
                                                        "</html>").arg("lambda",QString::number(this->whatsThisIconWidth),
                                                                       QString::number(this->whatsThisIconHeight),this->greenLambdaShortCut->key().toString()));
    this->unzoomImagePushButton->setWhatsThis(trUtf8("<html><font color=blue><img src=\":/images/%1.png\" width=\"%2\"height=\"%3\">"
                                                     "This pushbutton zooms out. <br/>"
                                                     "</html>").arg(this->unzoomImagePushButton->objectName(),
                                                                    QString::number(this->pushButtonIconWidth),QString::number(this->pushButtonIconheight)));
    this->renderImagePushButton->setWhatsThis(trUtf8("<html><font color=blue><img src=\":/images/%1.png\" width=\"%2\"height=\"%3\">"
                                                     "This pushbutton sends rendering image command to radmc3d. <br/>"
                                                     "</html>").arg(this->renderImagePushButton->objectName(),
                                                                    QString::number(this->pushButtonIconWidth*2),QString::number(this->pushButtonIconheight*2)));
    this->renderSpectrumPushButton->setWhatsThis(trUtf8("<html><font color=blue><img src=\":/images/%1.png\" width=\"%2\"height=\"%3\">"
                                                        "This pushbutton sends rendering spectrum command to radmc3d. <br/>"
                                                        "</html>").arg(this->renderSpectrumPushButton->objectName(),
                                                                       QString::number(this->pushButtonIconWidth*2),QString::number(this->pushButtonIconheight*2)));

    this->upPushButton->setWhatsThis(trUtf8("<html><font color=blue><img src=\":/images/%1.png\" width=\"%2\"height=\"%3\">"
                                            "This pushbutton moves image up<br/>"
                                            "</br>ShortCut:%4"
                                            "</html>").arg(this->upPushButton->objectName(),
                                                           QString::number(this->pushButtonIconWidth)
                                                           ,QString::number(this->pushButtonIconheight),this->upPushButton->shortcut().toString()));

    this->downPushButton->setWhatsThis(trUtf8("<html><font color=blue><img src=\":/images/%1.png\" width=\"%2\"height=\"%3\">"
                                              "This pushbutton moves image dwon. <br/>"
                                              "</br>ShortCut:%4"
                                              "</html>").arg(this->downPushButton->objectName(),
                                                             QString::number(this->pushButtonIconWidth)
                                                             ,QString::number(this->pushButtonIconheight),this->downPushButton->shortcut().toString()));

    this->leftPushButton->setWhatsThis(trUtf8("<html><font color=blue><img src=\":/images/%1.png\" width=\"%2\"height=\"%3\">"
                                              "This pushbutton moves image left. <br/>"
                                              "</br>ShortCut:%4"
                                              "</html>").arg(this->leftPushButton->objectName(),
                                                             QString::number(this->pushButtonIconWidth)
                                                             ,QString::number(this->pushButtonIconheight),this->leftPushButton->shortcut().toString()));

    this->rightPushButton->setWhatsThis(trUtf8("<html><font color=blue><img src=\":/images/%1.png\" width=\"%2\"height=\"%3\">"
                                               "This pushbutton moves image right. <br/>"
                                               "</br>ShortCut:%4"
                                               "</html>").arg(this->rightPushButton->objectName(),
                                                              QString::number(this->pushButtonIconWidth)
                                                              ,QString::number(this->pushButtonIconheight),this->rightPushButton->shortcut().toString()));

    this->clearDebuggerPushButton->setWhatsThis(trUtf8("<html><font color=blue><img src=\":/images/%1.png\" width=\"%2\"height=\"%3\">"
                                                       "This pushbutton clear debugger. <br/>"
                                                       "</html>").arg(this->clearDebuggerPushButton->objectName(),
                                                                      QString::number(this->pushButtonIconWidth),QString::number(this->pushButtonIconheight)));

    this->scollDebuggerToBottomPushButton->setWhatsThis(trUtf8("<html><font color=blue><img src=\":/images/%1.png\" width=\"%2\"height=\"%3\">"
                                                               "This pushbutton scroll debugger to bottom. <br/>"
                                                               "</html>").arg(this->scollDebuggerToBottomPushButton->objectName(),
                                                                              QString::number(this->pushButtonIconWidth),QString::number(this->pushButtonIconheight)));
    this->scollDebuggerToTopPushButton->setWhatsThis(trUtf8("<html><font color=blue><img src=\":/images/%1.png\" width=\"%2\"height=\"%3\">"
                                                            "This pushbutton  scroll debugger to top. <br/>"
                                                            "</html>").arg(this->scollDebuggerToTopPushButton->objectName(),
                                                                           QString::number(this->pushButtonIconWidth),QString::number(this->pushButtonIconheight)));
    this->colorTableComboxBox->setWhatsThis(trUtf8("<html><font color=blue><img src=\":/images/%1.png\" width=\"%2\"height=\"%3\">"
                                                   "Loading colortables with ending *.clut in current directory, which contains<br/>"
                                                   "red(i) green(i) blue(i) <br/> red(i+1) green(i+1) blue(i+1) <br/> as integer array<br/>"
                                                   "</html>").arg(this->colorTableSlider->objectName(),
                                                                  QString::number(this->pushButtonIconWidth),QString::number(this->pushButtonIconheight)));
    this->colorTableSlider->setWhatsThis(trUtf8("<html><font color=blue><img src=\":/images/%1.png\" width=\"%2\"height=\"%3\">"
                                                "Loading colortables with ending *.clut in current directory, which contains<br/>"
                                                "red(i) green(i) blue(i)<br/>red(i+1) green(i+1) blue(i+1)<br/> as integer array<br/>"
                                                "</html>").arg(this->colorTableSlider->objectName(),
                                                               QString::number(this->pushButtonIconWidth),QString::number(this->pushButtonIconheight)));
    this->showRadmc3dPushButton->setWhatsThis(trUtf8("<html><font color=blue><img src=\":/images/%1.png\" width=\"%2\"height=\"%3\">"
                                                     "This pushbutton  show radmc3d's log. <br/>"
                                                     "</html>").arg("radmc3d",
                                                                    QString::number(this->pushButtonIconWidth),QString::number(this->pushButtonIconheight)));

    this->redColorTuneValueLineEdit->setWhatsThis(trUtf8("<html><font color=blue>greenColorTune:</br><br/>"
                                                         "<img src=\":/images/no%1.png\" width=\"%2\"height=\"%3\">(Absolute scale is off) If set to 1, then rescale the brightness of"
                                                         " all channels the same value, to get the best"
                                                         " color depth.you can directly specify the weight of red color."
                                                         "<br/><br/><img src=\":/images/%1.png\" width=\"%2\"height=\"%3\">(Absolute scale is on)Set specifically the range of the image brightness "
                                                         "in terms of the log of the image values.</font><html>").arg(
                                                      this->absoluteScaleCheckBox->objectName(),
                                                      QString::number(this->pushButtonIconWidth),QString::number(this->pushButtonIconheight)));

    this->blueColorTuneValueLineEdit->setWhatsThis(trUtf8("<html><font color=blue>greenColorTune:</br><br/>"
                                                          "<img src=\":/images/no%1.png\" width=\"%2\"height=\"%3\">(Absolute scale is off) If set to 1, then rescale the brightness of"
                                                          " all channels the same value, to get the best"
                                                          " color depth.you can directly specify the weight of blue color."
                                                          "<br/><br/><img src=\":/images/%1.png\" width=\"%2\"height=\"%3\">(Absolute scale is on)Set specifically the range of the image brightness "
                                                          "in terms of the log of the image values.</font><html>").arg(
                                                       this->absoluteScaleCheckBox->objectName(),
                                                       QString::number(this->pushButtonIconWidth),QString::number(this->pushButtonIconheight)));

    this->greenColorTuneValueLineEdit->setWhatsThis(trUtf8("<html><font color=blue>greenColorTune:</br><br/>"
                                                           "<img src=\":/images/no%1.png\" width=\"%2\"height=\"%3\">(Absolute scale is off) If set to 1, then rescale the brightness of"
                                                           " all channels the same value, to get the best"
                                                           " color depth.you can directly specify the weight of green color ."
                                                           "<br/><br/><img src=\":/images/%1.png\" width=\"%2\"height=\"%3\">(Absolute scale is on)Set specifically the range of the image brightness "
                                                           "in terms of the log of the image values.</font><html>").arg(
                                                        this->absoluteScaleCheckBox->objectName(),
                                                        QString::number(this->pushButtonIconWidth),QString::number(this->pushButtonIconheight)));

    this->moleculeComboBox->setWhatsThis(trUtf8("<html><font color=blue>"
                                                "<img src=\":/images/%1.png\" width=\"%2\"height=\"%3\"> specify the wavelength by specifying which line of which <font color=red>molecule<font color=blue>, and at which "
                                                "velocity you want to render.</font><html>").arg(
                                             this->moleculeComboBox->objectName(),
                                             QString::number(this->pushButtonIconWidth),QString::number(this->pushButtonIconheight)));
    this->lineSpinBox->setWhatsThis(trUtf8("<html><font color=blue>"
                                           "<img src=\":/images/%1.png\" width=\"%2\"height=\"%3\"> specify the wavelength by specifying which<font color=red> line <font color=blue> of which molecule, and at which "
                                           "velocity you want to render.</font><html>").arg(
                                        this->lineSpinBox->objectName(),
                                        QString::number(this->pushButtonIconWidth),QString::number(this->pushButtonIconheight)));
    this->velocitySpinBox->setWhatsThis(trUtf8("<html><font color=blue>"
                                               "<img src=\":/images/%1.png\" width=\"%2\"height=\"%3\"> specify the wavelength by specifying which line of which molecule, and at which "
                                               "<font color=red>velocity<font color=blue> you want to render.</font><html>").arg(
                                            this->velocitySpinBox->objectName(),
                                            QString::number(this->pushButtonIconWidth),QString::number(this->pushButtonIconheight)));

    this->relobCheckbox->setWhatsThis(trUtf8("<html><font color=blue>"
                                             "<img src=\":/images/%1.png\" width=\"%2\"height=\"%3\">"
                                             "Do observer position relative <br/><img src=\":/images/no%1.png\" width=\"%2\"height=\"%3\"> Do observer position absolute  with respect to pointing.</font><html>").arg(
                                          this->relobCheckbox->objectName(),
                                          QString::number(this->pushButtonIconWidth),QString::number(this->pushButtonIconheight)));
    this->viewangleSpinBox->setWhatsThis(tr("<html><font color=blue>Viewing width.</font><html>"));
    this->scaleSpinBox->setWhatsThis(trUtf8("<html><font color=blue>Optional scaling factor for obs pos.</font><html>"));

}

/**
 * @brief This method set the icons.
 *
 */
void ImageControllerWidget::setIcons(){


    this->unzoomImagePushButton->setIcon(QIcon(QString(":/images/%1.png").arg(this->unzoomImagePushButton->objectName())));
    this->unzoomImagePushButton->setIconSize(QSize(this->pushButtonIconWidth,this->pushButtonIconheight));
    this->unzoomImagePushButton->setFixedSize(this->pushButtonIconWidth,this->pushButtonIconheight);

    this->renderImagePushButton->setIcon(QIcon(QString(":/images/%1.png").arg(this->renderImagePushButton->objectName())));
    this->renderImagePushButton->setIconSize(QSize(this->pushButtonIconWidth-10,this->pushButtonIconheight-10));
    this->renderImagePushButton->setFixedSize(this->pushButtonIconWidth,this->pushButtonIconheight);

    this->renderSpectrumPushButton->setIcon(QIcon(QString(":/images/%1.png").arg(this->renderSpectrumPushButton->objectName())));
    this->renderSpectrumPushButton->setIconSize(QSize(this->pushButtonIconWidth*2-10,this->pushButtonIconheight*2-10));
    this->renderSpectrumPushButton->setFixedSize(this->pushButtonIconWidth*2,this->pushButtonIconheight*2);

    this->upPushButton->setIcon(QIcon(QString(":/images/%1.png").arg(this->upPushButton->objectName())));
    this->upPushButton->setIconSize(QSize(this->pushButtonIconWidth/2,this->pushButtonIconheight/2));
    this->upPushButton->setFixedSize(this->pushButtonIconWidth,this->pushButtonIconheight);

    this->downPushButton->setIcon(QIcon(QString(":/images/%1.png").arg(this->downPushButton->objectName())));
    this->downPushButton->setIconSize(QSize(this->pushButtonIconWidth/2,this->pushButtonIconheight/2));
    this->downPushButton->setFixedSize(this->pushButtonIconWidth,this->pushButtonIconheight);

    this->leftPushButton->setIcon(QIcon(QString(":/images/%1.png").arg(this->leftPushButton->objectName())));
    this->leftPushButton->setIconSize(QSize(this->pushButtonIconWidth/2,this->pushButtonIconheight/2));
    this->leftPushButton->setFixedSize(this->pushButtonIconWidth,this->pushButtonIconheight);

    this->rightPushButton->setIcon(QIcon(QString(":/images/%1.png").arg(this->rightPushButton->objectName())));
    this->rightPushButton->setIconSize(QSize(this->pushButtonIconWidth/2,this->pushButtonIconheight/2));
    this->rightPushButton->setFixedSize(this->pushButtonIconWidth,this->pushButtonIconheight);

    this->clearDebuggerPushButton->setIcon(QIcon(QString(":/images/%1.png").arg(this->clearDebuggerPushButton->objectName())));
    this->clearDebuggerPushButton->setIconSize(QSize(this->comboBoxIconWidth,this->comboBoxIconHeight));
    this->clearDebuggerPushButton->setFixedSize(this->comboBoxIconWidth+2,this->comboBoxIconHeight+2);

    this->scollDebuggerToBottomPushButton->setIcon(QIcon(QString(":/images/%1.png").arg(this->scollDebuggerToBottomPushButton->objectName())));
    this->scollDebuggerToBottomPushButton->setIconSize(QSize(this->comboBoxIconWidth,this->comboBoxIconHeight));
    this->scollDebuggerToBottomPushButton->setFixedSize(this->comboBoxIconWidth+2,this->comboBoxIconHeight+2);

    this->scollDebuggerToTopPushButton->setIcon(QIcon(QString(":/images/%1.png").arg(this->scollDebuggerToTopPushButton->objectName())));
    this->scollDebuggerToTopPushButton->setIconSize(QSize(this->comboBoxIconWidth,this->comboBoxIconHeight));
    this->scollDebuggerToTopPushButton->setFixedSize(this->comboBoxIconWidth+2,this->comboBoxIconHeight+2);
    table=NULL;

}

void ImageControllerWidget::changePhiAndInclination(double phi,double inclination){
    this->blockSignals(true);
    this->phiSpinBox->setValue(phi);
    this->inclinationSpinBox->setValue(inclination);
    this->blockSignals(false);
    emit sendCommand(this->renderImagePushButton->objectName());
}

/**
 * @brief This method set the tooltips.
 *
 */
void ImageControllerWidget::setToolTips(){
    this->setStarToolTip(this->starCheckBox->isChecked());
    this->setSecondOrderToolTip(this->secondOrderCheckBox->isChecked());
    this->setContourToolTip(this->contourCheckBox->isChecked());
    this->setLinearToolTip(this->linearCheckBox->isChecked());
    this->setColorbarToolTip(this->colorbarCheckBox->isChecked());
    this->setTauToolTip(this->tauCheckBox->isChecked());
    this->setPreviewToolTip(this->previewCheckBox->isChecked());
    this->setInvertColorsToolTip(this->invertColorsCheckBox->isChecked());
    this->setPipeToolTip(this->pipeCheckBox->isChecked());
    this->setContrastToolTip(this->contrastSpinBox->value());
    this->setSaturationToolTip(this->saturationSpinBox->value());
    this->setNumberOfContoursToolTip(this->numberOfContoursSpinBox->value());
    this->setPositionAngleToolTip(this->positionAngleSpinBox->value());
    this->setNumberOfPixelsToolTip(this->numberOfPixelsSpinBox->value());
    this->setSizeToolTip(this->sizeLineEdit->text().remove(QRegExp(sizeRegExp)));
    this->setRedLambdaToolTip(this->redLambdaValueLineEdit->text().remove(QRegExp(lambdaRegExp)));
    this->setBlueLambdaToolTip(this->blueLambdaValueLineEdit->text().remove(QRegExp(lambdaRegExp)));
    this->setGreenLambdaToolTip(this->greenLambdaValueLineEdit->text().remove(QRegExp(lambdaRegExp)));
    this->setPhiToolTip(this->phiSpinBox->value());
    this->setInclinationToolTip(this->inclinationSpinBox->value());
    this->unzoomImagePushButton->setToolTip(trUtf8("zoom out"));
    this->renderImagePushButton->setToolTip(trUtf8("Render image"));
    this->upPushButton->setToolTip(trUtf8("Move up"));
    this->downPushButton->setToolTip(trUtf8("Move down"));
    this->rightPushButton->setToolTip(trUtf8("Move right"));
    this->leftPushButton->setToolTip(trUtf8("Move left"));
    this->clearDebuggerPushButton->setToolTip(trUtf8("Clear"));
    this->scollDebuggerToBottomPushButton->setToolTip(trUtf8("Scroll to bottom"));
    this->scollDebuggerToTopPushButton->setToolTip(trUtf8("scroll to top"));
    this->showRadmc3dPushButton->setToolTip(trUtf8("Show RADMC3D's Logfile"));
    this->viewangleSpinBox->setToolTip(tr("Viewing width"));
    this->scaleSpinBox->setToolTip(trUtf8("Optional scaling factor for obs pos"));

}

void ImageControllerWidget::setUpLocalTab(){
    QHBoxLayout *pointingHBox=new QHBoxLayout();
    this->pointingLabel=new QLabel();
    this->pointingLabel->setObjectName(trUtf8("pointing"));
    this->pointXLineEdit=new QLineEdit();
    this->pointXLineEdit->setObjectName("pointx");
    this->pointYLineEdit=new QLineEdit();
    this->pointYLineEdit->setObjectName("pointy");
    this->pointZLineEdit=new QLineEdit();
    this->pointZLineEdit->setObjectName("pointz");
    pointingHBox->addWidget(this->pointingLabel);
    pointingHBox->addWidget(this->pointXLineEdit);
    pointingHBox->addWidget(this->pointYLineEdit);
    pointingHBox->addWidget(this->pointZLineEdit);
    QHBoxLayout *observerposHBox=new QHBoxLayout();
    this->observerposLabel=new QLabel();
    this->observerposLabel->setObjectName(trUtf8("localAction"));
    this->observerposXLineEdit=new QLineEdit();
    this->observerposXLineEdit->setObjectName("obsx");
    this->observerposYLineEdit=new QLineEdit();
    this->observerposYLineEdit->setObjectName("obsy");
    this->observerposZLineEdit=new QLineEdit();
    this->observerposZLineEdit->setObjectName("obsz");
    observerposHBox->addWidget(this->observerposLabel);
    observerposHBox->addWidget(this->observerposXLineEdit);
    observerposHBox->addWidget(this->observerposYLineEdit);
    observerposHBox->addWidget(this->observerposZLineEdit);
    QHBoxLayout *viewHbox=new QHBoxLayout();
    this->viewangleLabel=new QLabel();
    this->viewangleSpinBox=new QDoubleSpinBox();
    this->viewangleSpinBox->setObjectName("viewangle");
    this->viewangleSpinBox->setKeyboardTracking(false);
    this->relobCheckbox=new QCheckBox();
    this->relobCheckbox->setObjectName("relob");
    this->scaleSpinBox=new QDoubleSpinBox();
    this->scaleSpinBox->setObjectName("scaleSpinBox");
    this->scaleSpinBox->setKeyboardTracking(false);
    viewHbox->addWidget(this->viewangleLabel);
    viewHbox->addWidget(this->viewangleSpinBox);
    viewHbox->addWidget(this->relobCheckbox);
    viewHbox->addWidget(this->scaleSpinBox);
    this->spinBoxLayout->addLayout(pointingHBox,9,0,1,5);
    this->spinBoxLayout->addLayout(observerposHBox,10,0,1,5);
    this->spinBoxLayout->addLayout(viewHbox,11,0,1,5);


}


void ImageControllerWidget::readTransferFile(bool load){
    QFile file(this->setup->getTransferDefaultFilename());
    map.clear();
    appendCommandVector.clear();
    if(file.exists()){
        file.open(QIODevice::ReadOnly|QIODevice::Text);
        QTextStream in(&file);
        while (!in.atEnd()) {
            QStringList strlist=in.readLine().split("=");
            if(strlist.size()==2){
                map.insert(strlist.at(0).trimmed(),strlist.at(1).trimmed());
                appendCommandVector.append(false);

            }
        }
        file.close();
    }
    if(map.size()>0){
        this->logger->writeToLogFile(trUtf8("userTransfer option is on"));
        this->usertransferOn=true;
        if(table==NULL){
            table=new QTableWidget();
            QStringList headers;
            table->setColumnCount(4);
            headers.append("Key");
            headers.append("Value");
            headers.append("Default");
            headers.append("Add to command");
            table->setHorizontalHeaderLabels(headers);
            QMapIterator<QString, QString> mapIterator(map);
            int var=0;
            QTableWidgetItem *item;
            table->setRowCount(map.size());
            while(mapIterator.hasNext()){
                mapIterator.next();
                item = new QTableWidgetItem(mapIterator.key());
                item->setFlags(item->flags() & ~Qt::ItemIsEditable);
                table->setItem(var,0,item);
                item->setFlags(item->flags() & Qt::ItemIsEditable);
                item = new QTableWidgetItem(mapIterator.value());
                table->setItem(var,1,item);
                item = new QTableWidgetItem(mapIterator.value());
                item->setFlags(Qt::NoItemFlags);
                table->setItem(var,2,item);
                item = new QTableWidgetItem("");
                item->setCheckState(Qt::Unchecked);
                table->setItem(var,3,item);
                var++;
            }
            if(load){
                QFile file(this->setup->getTransferFilename());
                file.open(QIODevice::ReadOnly|QIODevice::Text);
                QTextStream in(&file);
                while (!in.atEnd()) {
                    QStringList strlist=in.readLine().split("=");
                    if(strlist.size()==2){
                        if(map.contains(strlist.at(0))){
                            QList<QTableWidgetItem *> items=table->findItems(strlist.at(0),Qt::MatchCaseSensitive);
                            foreach (item, items) {
                                if(item->column()==0){
                                    QTableWidgetItem *valueItem=table->item(item->row(),1);
                                    valueItem->setText(strlist.at(1).trimmed());
                                }
                            }
                        }
                    }
                }
                file.close();
            }
            QHeaderView *headerView = table->horizontalHeader();
#if QT_VERSION < 0x050000
            headerView->setResizeMode(QHeaderView::Stretch);
            headerView->setResizeMode(1, QHeaderView::Stretch);
#else
            headerView->setSectionResizeMode(QHeaderView::Stretch);
            headerView->setSectionResizeMode(1,QHeaderView::Stretch);
#endif
            int index=middleTabWidget->indexOf(this->debuggerGroupBox);
            this->usertransferLayout=new QVBoxLayout();
            this->usertransferLayout->addWidget(table);
            this->usertransferGroupBox=new QGroupBox();
            this->usertransferGroupBox->setLayout(usertransferLayout);
            this->resetUserTransferPushButton=new QPushButton(trUtf8("Reset to default values"));
            this->resetUserTransferPushButton->setObjectName(trUtf8("reset"));
            this->resetUserTransferPushButton->setIcon(QIcon(trUtf8(":/images/%1.png").arg(this->resetUserTransferPushButton->objectName())));
            this->usertransferLayout->addWidget(resetUserTransferPushButton);
            connect(table,SIGNAL(itemChanged(QTableWidgetItem*)),this,SLOT(itemChanged(QTableWidgetItem*)));
            connect(resetUserTransferPushButton,SIGNAL(clicked()),this,SLOT(resetTransferFile()));
            this->middleTabWidget->insertTab(index,usertransferGroupBox,QIcon(":/images/transfer.png"),trUtf8("usertransfer"));
        }
    }
}

 QMap<QString,QString> ImageControllerWidget::getAppendCommand(){
     QMap<QString,QString> mapCommand;
     mapCommand.clear();
      QMapIterator<QString, QString> mapIterator(map);
     for (int var = 0; var < appendCommandVector.size(); ++var) {
           mapIterator.next();
         if(appendCommandVector.at(var)){
             mapCommand.insert(mapIterator.key(),mapIterator.value());
         }
     }
     return mapCommand;
}

void ImageControllerWidget::resetTransferFile(){
    int rowCount=this->table->rowCount();
    this->table->blockSignals(true);
    for (int var = 0; var < rowCount; ++var) {
        this->table->item(var,1)->setText(this->table->item(var,2)->text());
    }
    this->table->blockSignals(false);
    QFile file(this->setup->getTransferFilename());
    file.open(QIODevice::WriteOnly|QIODevice::Text);
    QTextStream in(&file);
    QMapIterator<QString, QString> mapIterator(map);
    while(mapIterator.hasNext()){
        mapIterator.next();
        in<<mapIterator.key()<<"="<<mapIterator.value()<<endl;
    }
    file.close();
    this->logger->writeToLogFile(trUtf8("transfer.inp is reseted to default"));
    emit readMyAction();
}

void ImageControllerWidget::removeTransferTab(){
    if(this->table!=NULL){
        int index=this->middleTabWidget->indexOf(this->usertransferGroupBox);
        this->logger->writeToLogFile(trUtf8("userTransfer option is off"));
        this->middleTabWidget->removeTab(index);
        this->usertransferOn=false;
        this->table=NULL;
    }
}

void ImageControllerWidget::itemChanged(QTableWidgetItem *item){
    if( item->text().trimmed()!=""){
        map.insert(table->item(item->row(),0)->text(),item->text());
        QFile file(this->setup->getTransferFilename());
        file.open(QIODevice::WriteOnly|QIODevice::Text);
        QTextStream in(&file);
        QMapIterator<QString, QString> mapIterator(map);
        while(mapIterator.hasNext()){
            mapIterator.next();
            in<<mapIterator.key()<<"="<<mapIterator.value()<<endl;
        }
        file.close();
        emit readMyAction();
    }else{
        if(item->column()!=3){
            QString key=table->item(item->row(),0)->text();
            item->setText(map.value(key));
            QMessageBox::information(this,trUtf8("Invalid value"),trUtf8("Value of %1 was  empty, value will reset").arg(key));
        }else{
            if(item->row()<appendCommandVector.size()){
                appendCommandVector[item->row()]=item->checkState();
            }
        }
    }

}
void ImageControllerWidget::reset(){
    this->blockSignals(true);
    this->starCheckBox->setChecked(false);
    this->secondOrderCheckBox->setChecked(false);
    this->tauCheckBox->setChecked(false);
    this->contourCheckBox->setChecked(false);
    this->previewCheckBox->setChecked(false);
    this->colorbarCheckBox->setChecked(false);
    this->linearCheckBox->setChecked(false);
    this->invertColorsCheckBox->setChecked(false);
    this->pipeCheckBox->setChecked(false);
    this->antiAliasingCheckBox->setChecked(false);
    this->dopplerCatchCheckBox->setChecked(false);
    this->absoluteScaleCheckBox->blockSignals(true);
    this->absoluteScaleCheckBox->setChecked(false);
    this->absoluteScaleCheckBox->blockSignals(false);
    this->colorTableSlider->setValue(0);

    this->contrastSpinBox->setValue(6.0);
    this->scaleSpinBox->setValue(1.0);
    this->numberOfContoursSpinBox->setValue(20.);

    this->numberOfPixelsSpinBox->setValue(100);

    this->phiSpinBox->setValue(phiSpinBoxDefaultValue);

    this->saturationSpinBox->setValue(saturationSpinBoxMaximum);

    this->positionAngleSpinBox->setValue(0);
    this->viewangleSpinBox->setValue(0.8);
    this->relobCheckbox->setChecked(false);

    this->inclinationSpinBox->setValue(inclinationSpinBoxDefaultValue);

    this->redLambdaIndexSlider->setValue(this->setup->getRedLambdaIndex());
    this->blueLambdaIndexSlider->setValue(this->setup->getBlueLambdaIndex());
    this->greenLambdaIndexSlider->setValue(this->setup->getGreenLambdaIndex());


    this->lastRedColorTuneValue=1.0;
    this->lastBlueColorTuneValue=1.0;
    this->lastGreenColorTuneValue=1.0;

    this->redColorTuneValueLineEdit->setText(QString::number(this->lastRedColorTuneValue));
    this->blueColorTuneValueLineEdit->setText(QString::number(this->lastBlueColorTuneValue));
    this->greenColorTuneValueLineEdit->setText(QString::number(this->lastGreenColorTuneValue));

    this->lastMoleculeValue=0;
    this->lastLineValue=1.0;
    this->lastVelocityValue=0.0;

    this->velocitySpinBox->setValue(this->lastVelocityValue);
    this->lineSpinBox->setValue(this->lastLineValue);
    this->moleculeComboBox->setCurrentIndex(this->lastMoleculeValue);

    this->lastSizeValue=this->setup->calculateImageSize();
    this->sizeLineEdit->setText(QString::number( lastSizeValue,'e')+"  "+getSizeUnit());
    this->sizeLineEdit->setCursorPosition(0);
    this->lastobsY=-2*lastSizeValue;
    this->observerposYLineEdit->setText(QString::number(lastobsY,'e')+"  "+getSizeUnit());

    this->lastPointY=0;
    this->lastPointX=0;
    this->lastobsZ=0;
    this->lastPointZ=0;
    this->lastobsX=0;
    this->pointXLineEdit->setText(QString::number(lastPointX,'e')+"  "+getSizeUnit());
    this->pointZLineEdit->setText(QString::number(lastPointZ,'e')+"  "+getSizeUnit());
    this->pointYLineEdit->setText(QString::number(lastPointY,'e')+"  "+getSizeUnit());

    this->observerposXLineEdit->setText(QString::number(lastobsX,'e')+"  "+getSizeUnit());
    this->observerposZLineEdit->setText(QString::number(lastobsZ,'e')+"  "+getSizeUnit());

    loadAsDefaults();
    this->blockSignals(false);
}

/**
 * @brief This method loads the last state of GUI-Elements, if the stated are saved.
 * @param loading, if true, the sendCommand is interrupted.
 * @see save
*/
void ImageControllerWidget::load(){

    QSettings *settings=setup->getLocalSettings();
    settings->beginGroup("ImageControllerWidget");
    this->blockSignals(true);
    this->starCheckBox->setChecked(settings->value(this->starCheckBox->objectName()).toBool());
    this->secondOrderCheckBox->setChecked(settings->value(this->secondOrderCheckBox->objectName()).toBool());
    this->tauCheckBox->setChecked(settings->value(this->tauCheckBox->objectName()).toBool());
    this->contourCheckBox->setChecked(settings->value(this->contourCheckBox->objectName()).toBool());
    this->previewCheckBox->setChecked(settings->value(this->previewCheckBox->objectName()).toBool());
    this->colorbarCheckBox->setChecked(settings->value(this->colorbarCheckBox->objectName()).toBool());
    this->linearCheckBox->setChecked(settings->value(this->linearCheckBox->objectName()).toBool());
    this->invertColorsCheckBox->setChecked(settings->value(this->invertColorsCheckBox->objectName()).toBool());
    this->pipeCheckBox->setChecked(settings->value(this->pipeCheckBox->objectName()).toBool());
    this->antiAliasingCheckBox->setChecked(settings->value(this->antiAliasingCheckBox->objectName()).toBool());
    this->renderSpectrumActive=(settings->value("spectromIsOn").toBool());
    this->dopplerCatchCheckBox->setChecked(settings->value(this->dopplerCatchCheckBox->objectName()).toBool());
    this->absoluteScaleCheckBox->blockSignals(true);
    this->absoluteScaleCheckBox->setChecked(settings->value(this->absoluteScaleCheckBox->objectName()).toBool());
    this->relobCheckbox->setChecked(settings->value(this->relobCheckbox->objectName()).toBool());
    this->absoluteScaleCheckBox->blockSignals(false);
    this->colorTableSlider->setValue(settings->value(this->colorTableSlider->objectName()).toInt());
    this->blockSignals(false);
    this->positionComboxBox->setCurrentIndex(settings->value(this->positionComboxBox->objectName()).toInt());
    this->blockSignals(true);

    this->contrastSpinBox->setValue(settings->value(this->contrastSpinBox->objectName()).toDouble());

    this->numberOfContoursSpinBox->setValue(settings->value(this->numberOfContoursSpinBox->objectName()).toDouble());

    this->numberOfPixelsSpinBox->setValue(settings->value(this->numberOfPixelsSpinBox->objectName()).toDouble());

    this->phiSpinBox->setValue(settings->value(this->phiSpinBox->objectName()).toDouble());

    this->saturationSpinBox->setValue(settings->value(this->saturationSpinBox->objectName()).toDouble());
    this->scaleSpinBox->setValue(settings->value(this->scaleSpinBox->objectName()).toDouble());

    this->positionAngleSpinBox->setValue(settings->value(this->positionAngleSpinBox->objectName()).toDouble());
    this->viewangleSpinBox->setValue(settings->value(this->viewangleSpinBox->objectName()).toDouble());


    this->inclinationSpinBox->setValue(settings->value(this->inclinationSpinBox->objectName()).toDouble());

    this->redLambdaIndexSlider->setValue(settings->value(this->redLambdaValueLineEdit->objectName()).toInt());
    this->blueLambdaIndexSlider->setValue(settings->value(this->blueLambdaValueLineEdit->objectName()).toInt());
    this->greenLambdaIndexSlider->setValue(settings->value(this->greenLambdaValueLineEdit->objectName()).toInt());

    this->lastSizeValue=settings->value(this->sizeLineEdit->objectName()).toDouble();
    this->sizeLineEdit->setText(QString::number(lastSizeValue,'e')+" "+getSizeUnit());
    this->sizeLineEdit->setCursorPosition(0);


    this->lastPointX=settings->value(this->pointXLineEdit->objectName()).toDouble();
    this->pointXLineEdit->setText(QString::number(lastPointX,'e')+" "+getSizeUnit());

    this->lastPointY=settings->value(this->pointYLineEdit->objectName()).toDouble();
    this->pointYLineEdit->setText(QString::number(lastPointY,'e')+" "+getSizeUnit());

    this->lastPointZ=settings->value(this->pointZLineEdit->objectName()).toDouble();
    this->pointZLineEdit->setText(QString::number(lastPointZ,'e')+" "+getSizeUnit());

    this->lastobsX=settings->value(this->observerposXLineEdit->objectName()).toDouble();
    this->observerposXLineEdit->setText(QString::number(lastobsX,'e')+" "+getSizeUnit());

    this->lastobsY=settings->value(this->observerposYLineEdit->objectName()).toDouble();
    this->observerposYLineEdit->setText(QString::number(lastobsY,'e')+" "+getSizeUnit());

    this->lastobsZ=settings->value(this->observerposZLineEdit->objectName()).toDouble();
    this->observerposZLineEdit->setText(QString::number(lastobsZ,'e')+" "+getSizeUnit());




    this->lastRedColorTuneValue=settings->value(this->redColorTuneValueLineEdit->objectName()).toDouble();
    this->lastBlueColorTuneValue=settings->value(this->blueColorTuneValueLineEdit->objectName()).toDouble();
    this->lastGreenColorTuneValue=settings->value(this->greenColorTuneValueLineEdit->objectName()).toDouble();
    this->redColorTuneValueLineEdit->setText(QString::number(this->lastRedColorTuneValue));
    this->blueColorTuneValueLineEdit->setText(QString::number(this->lastBlueColorTuneValue));
    this->greenColorTuneValueLineEdit->setText(QString::number(this->lastGreenColorTuneValue));


    this->velocitySpinBox->setValue(settings->value(this->velocitySpinBox->objectName()).toDouble());
    this->lineSpinBox->setValue(settings->value(this->lineSpinBox->objectName()).toInt());
    lastMoleculeValue=settings->value(this->moleculeComboBox->objectName()).toInt();


    this->numberOfContours=this->numberOfContoursSpinBox->value();
    this->numberOfPixels=this->numberOfPixelsSpinBox->value();
    this->saturation=saturationSpinBox->value();
    this->phi=this->phiSpinBox->value();
    this->inclination=this->inclinationSpinBox->value();
    this->maxlog=this->contrastSpinBox->value();
    this->posAngle=this->positionAngleSpinBox->value();
    QString zoom=settings->value("zoomBox").toString();
    if(zoom.trimmed()!=""){
        QStringList list=zoom.split(QRegExp("\\s+"),QString::SkipEmptyParts);
        zoomBox.clear();
        zoomBox.append(list.at(0).toDouble());
        zoomBox.append(list.at(1).toDouble());
        zoomBox.append(list.at(2).toDouble());
        zoomBox.append(list.at(3).toDouble());
        zoomActive=true;
    }
    this->pointXLineEdit->setCursorPosition(0);
    this->pointYLineEdit->setCursorPosition(0);
    this->pointZLineEdit->setCursorPosition(0);
    this->observerposXLineEdit->setCursorPosition(0);
    this->observerposYLineEdit->setCursorPosition(0);
    this->observerposZLineEdit->setCursorPosition(0);
    settings->endGroup();
    this->blockSignals(false);
}
/**
 * @brief This method loads the last state of GUI-Elements, if the stated are saved.
 * @param loading, if true, the sendCommand is interrupted.
 * @see save
*/
void ImageControllerWidget::loadAsDefaults(){

    QSettings *settings=setup->getGlobalSettings();
    settings->beginGroup("ImageControllerWidget");
    if(settings->childKeys().size()>0){
        this->blockSignals(true);
        this->starCheckBox->setChecked(settings->value(this->starCheckBox->objectName()).toBool());
        this->secondOrderCheckBox->setChecked(settings->value(this->secondOrderCheckBox->objectName()).toBool());
        this->tauCheckBox->setChecked(settings->value(this->tauCheckBox->objectName()).toBool());
        this->contourCheckBox->setChecked(settings->value(this->contourCheckBox->objectName()).toBool());
        this->previewCheckBox->setChecked(settings->value(this->previewCheckBox->objectName()).toBool());
        this->colorbarCheckBox->setChecked(settings->value(this->colorbarCheckBox->objectName()).toBool());
        this->linearCheckBox->setChecked(settings->value(this->linearCheckBox->objectName()).toBool());
        this->invertColorsCheckBox->setChecked(settings->value(this->invertColorsCheckBox->objectName()).toBool());
        this->pipeCheckBox->setChecked(settings->value(this->pipeCheckBox->objectName()).toBool());
        this->antiAliasingCheckBox->setChecked(settings->value(this->antiAliasingCheckBox->objectName()).toBool());
        this->renderSpectrumActive=(settings->value("spectromIsOn").toBool());
        this->dopplerCatchCheckBox->setChecked(settings->value(this->dopplerCatchCheckBox->objectName()).toBool());
        this->absoluteScaleCheckBox->blockSignals(true);
        this->absoluteScaleCheckBox->setChecked(settings->value(this->absoluteScaleCheckBox->objectName()).toBool());
        this->relobCheckbox->setChecked(settings->value(this->relobCheckbox->objectName()).toBool());
        this->absoluteScaleCheckBox->blockSignals(false);
        this->colorTableSlider->setValue(settings->value(this->colorTableSlider->objectName()).toInt());
        this->blockSignals(false);
        this->positionComboxBox->setCurrentIndex(settings->value(this->positionComboxBox->objectName()).toInt());
        this->blockSignals(true);
        this->scaleSpinBox->setValue(settings->value(this->scaleSpinBox->objectName()).toDouble());
        this->contrastSpinBox->setValue(settings->value(this->contrastSpinBox->objectName()).toDouble());
        this->numberOfContoursSpinBox->setValue(settings->value(this->numberOfContoursSpinBox->objectName()).toDouble());
        this->numberOfPixelsSpinBox->setValue(settings->value(this->numberOfPixelsSpinBox->objectName()).toDouble());
        this->phiSpinBox->setValue(settings->value(this->phiSpinBox->objectName()).toDouble());
        this->saturationSpinBox->setValue(settings->value(this->saturationSpinBox->objectName()).toDouble());
        this->positionAngleSpinBox->setValue(settings->value(this->positionAngleSpinBox->objectName()).toDouble());
        this->viewangleSpinBox->setValue(settings->value(this->viewangleSpinBox->objectName()).toDouble());
        this->inclinationSpinBox->setValue(settings->value(this->inclinationSpinBox->objectName()).toDouble());
        this->numberOfContours=this->numberOfContoursSpinBox->value();
        this->numberOfPixels=this->numberOfPixelsSpinBox->value();
        this->saturation=saturationSpinBox->value();
        this->phi=this->phiSpinBox->value();
        this->inclination=this->inclinationSpinBox->value();
        this->maxlog=this->contrastSpinBox->value();
        this->posAngle=this->positionAngleSpinBox->value();
        this->velocitySpinBox->setValue(settings->value(this->velocitySpinBox->objectName()).toDouble());
        settings->endGroup();
    }else{
        this->setLineModus(false);
        activateLocalObserver(false);

        removeColorModusWidgets();
    }
    this->blockSignals(false);
    this->pointXLineEdit->setCursorPosition(0);
    this->pointYLineEdit->setCursorPosition(0);
    this->pointZLineEdit->setCursorPosition(0);
    this->observerposXLineEdit->setCursorPosition(0);
    this->observerposYLineEdit->setCursorPosition(0);
    this->observerposZLineEdit->setCursorPosition(0);
}
/**
 * @brief This method saves the state of GUI-Elements.
 * @see load
*/
void ImageControllerWidget::save(){
    QSettings *settings=setup->getLocalSettings();
    settings->beginGroup("ImageControllerWidget");
    settings->setValue(this->starCheckBox->objectName(),this->starCheckBox->isChecked());
    settings->setValue(this->secondOrderCheckBox->objectName(),this->secondOrderCheckBox->isChecked());
    settings->setValue(this->tauCheckBox->objectName(),this->tauCheckBox->isChecked());
    settings->setValue(this->contourCheckBox->objectName(),this->contourCheckBox->isChecked());
    settings->setValue(this->previewCheckBox->objectName(),this->previewCheckBox->isChecked());
    settings->setValue(this->colorbarCheckBox->objectName(),this->colorbarCheckBox->isChecked());
    settings->setValue(this->linearCheckBox->objectName(),this->linearCheckBox->isChecked());
    settings->setValue(this->relobCheckbox->objectName(),this->relobCheckbox->isChecked());
    settings->setValue("spectromIsOn",this->renderSpectrumActive);
    settings->setValue(this->invertColorsCheckBox->objectName(),this->invertColorsCheckBox->isChecked());
    settings->setValue(this->pipeCheckBox->objectName(),this->pipeCheckBox->isChecked());
    settings->setValue(this->antiAliasingCheckBox->objectName(),this->antiAliasingCheckBox->isChecked());
    settings->setValue(this->dopplerCatchCheckBox->objectName(),this->dopplerCatchCheckBox->isChecked());
    settings->setValue(this->absoluteScaleCheckBox->objectName(),this->absoluteScaleCheckBox->isChecked());
    settings->setValue(this->colorTableSlider->objectName(),this->colorTableSlider->value());
    settings->setValue(this->positionComboxBox->objectName(),this->positionComboxBox->currentIndex());
    settings->setValue(this->contrastSpinBox->objectName(),this->contrastSpinBox->value());
    settings->setValue(this->numberOfContoursSpinBox->objectName(),this->numberOfContoursSpinBox->value());
    settings->setValue(this->numberOfPixelsSpinBox->objectName(),this->numberOfPixelsSpinBox->value());
    settings->setValue(this->phiSpinBox->objectName(),this->phiSpinBox->value());
    settings->setValue(this->viewangleSpinBox->objectName(),this->viewangleSpinBox->value());
    settings->setValue(this->saturationSpinBox->objectName(),this->saturationSpinBox->value());
    settings->setValue(this->positionAngleSpinBox->objectName(),this->positionAngleSpinBox->value());
    settings->setValue(this->inclinationSpinBox->objectName(),this->inclinationSpinBox->value());
    settings->setValue(this->redLambdaValueLineEdit->objectName(),this->redLambdaIndexSlider->value());
    settings->setValue(this->greenLambdaValueLineEdit->objectName(),this->greenLambdaIndexSlider->value());
    settings->setValue(this->blueLambdaValueLineEdit->objectName(),this->blueLambdaIndexSlider->value());

    settings->setValue(this->sizeLineEdit->objectName(),this->sizeLineEdit->text().remove(QRegExp(sizeRegExp)));
    settings->setValue(this->pointXLineEdit->objectName(),this->pointXLineEdit->text().remove(QRegExp(sizeRegExp)));
    settings->setValue(this->pointYLineEdit->objectName(),this->pointYLineEdit->text().remove(QRegExp(sizeRegExp)));
    settings->setValue(this->pointZLineEdit->objectName(),this->pointZLineEdit->text().remove(QRegExp(sizeRegExp)));
    settings->setValue(this->observerposXLineEdit->objectName(),this->observerposXLineEdit->text().remove(QRegExp(sizeRegExp)));
    settings->setValue(this->observerposYLineEdit->objectName(),this->observerposYLineEdit->text().remove(QRegExp(sizeRegExp)));
    settings->setValue(this->observerposZLineEdit->objectName(),this->observerposZLineEdit->text().remove(QRegExp(sizeRegExp)));

    settings->setValue(this->redColorTuneValueLineEdit->objectName(),this->redColorTuneValueLineEdit->text());
    settings->setValue(this->blueColorTuneValueLineEdit->objectName(),this->blueColorTuneValueLineEdit->text());
    settings->setValue(this->greenColorTuneValueLineEdit->objectName(),this->greenColorTuneValueLineEdit->text());

    settings->setValue(this->velocitySpinBox->objectName(),this->velocitySpinBox->value());
    settings->setValue(this->lineSpinBox->objectName(),this->lineSpinBox->value());
    settings->setValue(this->moleculeComboBox->objectName(),lastMoleculeValue);
    settings->setValue(this->scaleSpinBox->objectName(),scaleSpinBox->value());
    if(this->isZoomActive()){
        QString zoom=QString::number(this->zoomBox.at(0))+" "+QString::number(this->zoomBox.at(1))
                +" "+QString::number(this->zoomBox.at(2))
                +" "+QString::number(this->zoomBox.at(3));
        settings->setValue("zoomBox",zoom);
    }
    settings->endGroup();
}

/**
 * @brief This method saves the state of GUI-Elements as defaults.
 * @see load
*/
void ImageControllerWidget::saveAsDefaults(){
    QSettings *settings=setup->getGlobalSettings();
    settings->beginGroup("ImageControllerWidget");
    settings->setValue(this->starCheckBox->objectName(),this->starCheckBox->isChecked());
    settings->setValue(this->secondOrderCheckBox->objectName(),this->secondOrderCheckBox->isChecked());
    settings->setValue(this->tauCheckBox->objectName(),this->tauCheckBox->isChecked());
    settings->setValue(this->contourCheckBox->objectName(),this->contourCheckBox->isChecked());
    settings->setValue(this->previewCheckBox->objectName(),this->previewCheckBox->isChecked());
    settings->setValue(this->colorbarCheckBox->objectName(),this->colorbarCheckBox->isChecked());
    settings->setValue(this->linearCheckBox->objectName(),this->linearCheckBox->isChecked());
    settings->setValue(this->relobCheckbox->objectName(),this->relobCheckbox->isChecked());
    settings->setValue("spectromIsOn",this->renderSpectrumActive);
    settings->setValue(this->invertColorsCheckBox->objectName(),this->invertColorsCheckBox->isChecked());
    settings->setValue(this->pipeCheckBox->objectName(),this->pipeCheckBox->isChecked());
    settings->setValue(this->antiAliasingCheckBox->objectName(),this->antiAliasingCheckBox->isChecked());
    settings->setValue(this->dopplerCatchCheckBox->objectName(),this->dopplerCatchCheckBox->isChecked());
    settings->setValue(this->absoluteScaleCheckBox->objectName(),this->absoluteScaleCheckBox->isChecked());
    settings->setValue(this->colorTableSlider->objectName(),this->colorTableSlider->value());
    settings->setValue(this->positionComboxBox->objectName(),this->positionComboxBox->currentIndex());
    settings->setValue(this->contrastSpinBox->objectName(),this->contrastSpinBox->value());
    settings->setValue(this->numberOfContoursSpinBox->objectName(),this->numberOfContoursSpinBox->value());
    settings->setValue(this->numberOfPixelsSpinBox->objectName(),this->numberOfPixelsSpinBox->value());
    settings->setValue(this->phiSpinBox->objectName(),this->phiSpinBox->value());
    settings->setValue(this->viewangleSpinBox->objectName(),this->viewangleSpinBox->value());
    settings->setValue(this->saturationSpinBox->objectName(),this->saturationSpinBox->value());
    settings->setValue(this->positionAngleSpinBox->objectName(),this->positionAngleSpinBox->value());
    settings->setValue(this->inclinationSpinBox->objectName(),this->inclinationSpinBox->value());
    settings->setValue(this->scaleSpinBox->objectName(),scaleSpinBox->value());
    settings->setValue(this->velocitySpinBox->objectName(),this->velocitySpinBox->value());
    settings->endGroup();
}
