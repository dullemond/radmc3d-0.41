#include "mainwindow.h"
#define NINT(a) ((a) >= 0.0 ? (int)((a)+0.5) : (int)((a)-0.5))

int const MainWindow::EXIT_CODE_REBOOT = -123456789;

/**
*   @file mainwindow.h
*   @brief The class Mainwindow is derived from QMainWindow. All the widget in the Gui will be placed on
*   the MainWindow
*   @author Farzin Sereshti
*   @version 1.0
*/

/**
 *  @brief This is the contsructur.
 *
 *  @param a widget parent can be given. This a is optional parameter.
 *  @details This method initializes the mainWindow's widget.
 *  @see setUp()
 *
 */
MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent)
{


    unformattedSetting=UnformattedSetting::getInstance(this);
    initialized=false;
    exitCode=0;
    waitForFinishing=false;
    this->isErrorbar=false;
    this->zoomActive=false;
    this->exit=false;
    this->loaded=false;
    this->setZoomBox=false;
    this->startnexttime=false;
    this->spectrumIsRead=false;
    this->imageIsRead=false;
    this->timer=NULL;
    this->firstTime=true;
    this->spectrum=NULL;
    this->image=NULL;
    this->showLog=true;
    this->nativeMenuBar=false;
    this->historyMax=5;
    this->setup=SetUp::getInstance();
    this->physicalConstants=PhysicalConstants::getInstance();
    this->logger=Logger::getInstance();
    QSettings *setting=this->setup->getGlobalSettings();
    startnexttime=false;
    if(setting->value("startNextTimeFromCurrentDirectory").toBool()){
        startnexttime=true;
        QString path=readLastWorkingDirectory();
        path.remove("\n");
        if((path.trimmed()!="") & (QFile(path).exists())){
            QDir::setCurrent(path);
            this->logger->createLogFile();
            this->logger->writeToLogFile(trUtf8("Current directory is %1").arg(QDir::currentPath()));
        }
    }
    this->nativeMenuBar=setting->value("NativeMenuBar").toBool();
    if(setting->contains("ByteOrder")){
        this->unformattedSetting->setByteOrder((QDataStream::ByteOrder)setting->value("ByteOrder").toInt());
        if(this->unformattedSetting->getByteOrder()==QDataStream::LittleEndian)
            this->logger->writeToLogFile("current byte order is LittleEndian");
        else
            this->logger->writeToLogFile("current byte order is BigEndian");
    }
    if(setting->contains("RecordSize")){
        this->unformattedSetting->setRecordLengthSize(setting->value("RecordSize").toInt());
        this->logger->writeToLogFile(trUtf8("current recordlength's size is %1 bytes").arg(this->unformattedSetting->getRecordLengthSize()));
    }
    mutex.lock();
    exit=true;
    mutex.unlock();
    setZoomBox=false;
    this->imageIsRead=false;
    this->imageControllerWidgetTabString=trUtf8("Render image");
    this->mainLayout=new QHBoxLayout();
    this->leftWidget=new QStackedWidget();
    this->rightWidget=new QTabWidget();
    this->imageLabel=new ImageLabelWidget(this->coordDistanceToImage);
    QHBoxLayout* imagebox=new QHBoxLayout();
    this->imageGroup=new QGroupBox();
    this->imageGroup->setLayout(imagebox);
    imagebox->setAlignment(Qt::AlignCenter);
    imagebox->addWidget(this->imageLabel);

    this->spectrumLabel=new SpectrumLabelWidget(this->coordDistanceToImage);
    QHBoxLayout *spectrumbox=new QHBoxLayout();
    this->spectrumGroup=new QGroupBox();
    this->spectrumGroup->setLayout(spectrumbox);
    spectrumbox->setAlignment(Qt::AlignCenter);
    spectrumbox->addWidget(this->spectrumLabel);

    this->mainLayout->setAlignment(Qt::AlignCenter);
    this->leftWidget->addWidget(imageGroup);
    this->leftWidget->addWidget(spectrumGroup);
    this->mainLayout->addWidget(this->leftWidget);
    this->mainLayout->addWidget(this->rightWidget);
    this->textEdit=new QTextEdit();
    this->texteditbox=new QVBoxLayout();
    this->textGroup=new QGroupBox();
    this->textGroup->setLayout(texteditbox);
    this->closeLogPushButton=new QPushButton();
    this->texteditbox->addWidget(closeLogPushButton);
    this->texteditbox->addWidget(this->textEdit);
    this->leftWidget->addWidget(this->textGroup);
    this->setCentralWidget(new QWidget());
    this->centralWidget()->setLayout(this->mainLayout);
    this->imageControllerWidget=new ImageControllerWidget();
    this->rightWidget->setDocumentMode(true);
    this->rightWidget->setMinimumSize(300,400);
    this->rightWidget->addTab(this->imageControllerWidget,QIcon(QString(":/images/%1.png").arg(this->imageControllerWidget->objectName())),this->imageControllerWidgetTabString);
    this->helpMenu=new QMenu(this);
    this->manualAction=new QAction(this);
    this->guimanualAction=new QAction(this);
    this->aboutAction=new QAction(this);
    this->whatsthisAction=QWhatsThis::createAction(this);
    this->manualAction->setIconVisibleInMenu(true);
    this->guimanualAction->setIconVisibleInMenu(true);
    this->whatsthisAction->setIconVisibleInMenu(true);
    this->aboutAction->setIconVisibleInMenu(true);
    this->helpMenu->addAction(this->manualAction);
    this->helpMenu->addAction(this->guimanualAction);
    this->helpMenu->addAction(this->aboutAction);
    this->helpMenu->addAction(this->whatsthisAction);
    this->nativeMenuBarAction=new QAction(this->helpMenu);
    this->nativeMenuBarAction->setCheckable(true);
    this->nativeMenuBarAction->setChecked(this->nativeMenuBar);
    this->helpMenu->addAction(this->nativeMenuBarAction);
    this->configMenu=new QMenu(this);
    this->escAction=new QAction(this);
    this->escAction->setShortcut(QKeySequence(Qt::Key_Escape));
    this->escAction->setShortcutContext(Qt::ApplicationShortcut);
    this->addAction(this->escAction);
    this->sizeGroup=new QActionGroup(this);
    this->sizeGroup->setExclusive(true);
    sizemenu=new QMenu(configMenu);
    sizemenu->setIcon(QIcon(":/images/size.png"));
    sizemenu->setTitle("Size unit");
    sizemenu->setObjectName("sizemenu");
    sizesList<<"cm"<<"pc"<<"au";
    foreach (QString str, sizesList) {
        QAction *action=new QAction(sizemenu);
        action->setCheckable(true);
        action->setText(str);
        action->setIconVisibleInMenu(true);
        action->setActionGroup(sizeGroup);
        sizemenu->addAction(action);
        action->setIcon(QIcon(":/images/size.png"));
    }
    sizemenu->actions().at(0)->setChecked(true);
    configMenu->addMenu(sizemenu);

    this->progressDialog=new QProgressDialog();
    this->progressDialog->setValue(0);
    this->progressDialog->setMinimum(0);
    this->progressDialog->setMaximum(0);
    this->progressDialog->setMaximumHeight(20);
    this->progressDialog->setStyleSheet("QProgressBar::chunk { background-color: #009B00; width: 3px; margin:0.7px}");


    this->statusBar()->addPermanentWidget(progressDialog,1);
    this->statusBar()->hide();

    this->userTransferAction=new QAction(configMenu);
    this->userTransferAction->setCheckable(true);
    this->userTransferAction->setObjectName("transfer");

    this->colorAction=new QAction(configMenu);
    this->colorAction->setCheckable(true);
    this->colorAction->setObjectName(trUtf8("color"));

    this->lineOptionAction=new QAction(configMenu);
    this->lineOptionAction->setCheckable(true);
    this->imageControllerWidget->setLineModus(false);
    this->lineOptionAction->setObjectName(trUtf8("line"));
    this->lineOptionAction->setIconVisibleInMenu(true);

    this->changeCurrentDirectoryAction=new QAction(configMenu);
    this->changeCurrentDirectoryAction->setIconVisibleInMenu(true);
    this->changeCurrentDirectoryAction->setObjectName("dir");

    this->localAction=new QAction(configMenu);
    this->localAction->setCheckable(true);
    this->localAction->setObjectName(trUtf8("localAction"));
    this->localAction->setIconVisibleInMenu(true);

    this->savedefaultsAction=new QAction(this->configMenu);
    this->savedefaultsAction->setObjectName("global");
    this->savedefaultsAction->setIconVisibleInMenu(true);

    this->configMenu->addAction(this->userTransferAction);
    this->configMenu->addAction(this->colorAction);
    this->configMenu->addAction(this->lineOptionAction);
    this->configMenu->addAction(this->localAction);
    this->configMenu->addAction(this->changeCurrentDirectoryAction);
    this->configMenu->addAction(this->savedefaultsAction);

    this->historyMenu=new QMenu(this);
    this->radmc3dProcess=NULL;

    this->spectrumMenu=new QMenu(this);
    this->loadPointsAction=new QAction(spectrumMenu);
    this->loadPointsAction->setCheckable(true);
    this->loadPointsAction->setObjectName("LSP");
    this->loadPointsAction->setIconVisibleInMenu(true);
    this->loadPointsAction->setChecked(false);
    this->xLinAction=new QAction(spectrumMenu);
    this->yLinAction=new QAction(spectrumMenu);
    this->xLinAction->setObjectName("xlin");
    this->xLinAction->setShortcut(QKeySequence("Ctrl+x"));
    this->xLinAction->setCheckable(true);
    this->xLinAction->setIconVisibleInMenu(true);
    this->yLinAction->setShortcut(QKeySequence("Ctrl+y"));
    this->yLinAction->setCheckable(true);
    this->yLinAction->setObjectName("ylin");
    this->yLinAction->setIconVisibleInMenu(true);
    this->setIcons();
    this->setImageMenu();
    this->imageControllerWidget->activateLocalObserver(this->localAction->isChecked());
    this->readImageAction=new QWidgetAction(imageMenu);
    this->readImageAction->setObjectName(trUtf8("readImage"));
    this->readImageAction->setShortcut(QKeySequence("Ctrl+r"));
    this->readImageAction->setIcon(QIcon(trUtf8(":/images/%1.png").arg(this->readImageAction->objectName())));
    this->readImageAction->setIconVisibleInMenu(true);
    this->readSpectrumAction=new QWidgetAction(this->spectrumMenu);
    this->readSpectrumAction->setObjectName(trUtf8("readSpectrum"));
    this->readSpectrumAction->setShortcut(QKeySequence("Ctrl+u"));
    this->readSpectrumAction->setIcon(QIcon(trUtf8(":/images/%1.png").arg(this->readImageAction->objectName())));
    this->readSpectrumAction->setIconVisibleInMenu(true);
    this->imageControllerWidget->loadColorLookUpTables(QCoreApplication::applicationDirPath(),false);
    this->imageControllerWidget->setUpColorTableComboBoxValues();
    this->imageControllerWidget->setUpColorTableSlider();
    this->setupHistory(true);


#if (QT_VERSION >= QT_VERSION_CHECK(5, 0, 0))
    QTextCodec::setCodecForLocale(QTextCodec::codecForName("UTF-8"));
#else
    QTextCodec::setCodecForCStrings(QTextCodec::codecForName("UTF-8"));
#endif
    freqUnitsList<<"\u03BB [m\u03BC]"<<"\u03BD [Hz]"<<"\u03BD [eV]"<<"\u03BD [KeV]";
    this->freqGroup=new QActionGroup(this);
    this->freqGroup->setExclusive(true);
    freqmenu=new QMenu(spectrumMenu);
    freqmenu->setIcon(QIcon(":/images/size.png"));
    freqmenu->setTitle("X-title");
    freqmenu->setObjectName("X-title");
    foreach (QString str, freqUnitsList) {
        QAction *action=new QAction(freqmenu);
        action->setCheckable(true);
        action->setText(str);
        action->setIconVisibleInMenu(true);
        action->setActionGroup(freqGroup);
        freqmenu->addAction(action);

    }
    freqmenu->actions().at(0)->setChecked(true);



    fluxUnitsList<<"\u03BDF\u1D65[erg\u00B7cm\u207B\u00B2\u00B7s\u207B\u00B9]"
                <<"\u03BDF\u1D65[JyHz]"<<"F\u1D65[erg\u00B7cm\u207B\u00B2\u00B7Hz\u207B\u00B9\u00B7s\u207B\u00B9]"
               <<"F\u1D65[HZ]"<<"\u03BDL\u1D65[erg\u00B7s\u207B\u00B9]"<<"\u03BDL\u1D65[L\u2609]"<<
                 "L\u1D65[erg\u00B7Hz\u207B\u00B9\u00B7s\u207B\u00B9]";

    this->fluxGroup=new QActionGroup(this);
    this->fluxGroup->setExclusive(true);
    fluxmenu=new QMenu(spectrumMenu);
    fluxmenu->setIcon(QIcon(":/images/size.png"));
    fluxmenu->setTitle("Y-title");
    fluxmenu->setObjectName("Y-title");
    foreach (QString str, fluxUnitsList) {
        QAction *action=new QAction(fluxmenu);
        action->setCheckable(true);
        action->setText(str);
        action->setIconVisibleInMenu(true);
        action->setActionGroup(fluxGroup);
        fluxmenu->addAction(action);
    }
    fluxmenu->actions().at(0)->setChecked(true);

    this->setUp();
    this->initialize();
    this->menuBar()->addMenu(this->historyMenu);
    this->menuBar()->addMenu(this->imageMenu);
    this->imageMenu->addAction(readImageAction);
    this->menuBar()->addMenu(this->spectrumMenu);
    this->spectrumMenu->addMenu(freqmenu);
    this->spectrumMenu->addMenu(fluxmenu);
    this->spectrumMenu->addAction(this->xLinAction);
    this->spectrumMenu->addAction(this->yLinAction);
    this->spectrumMenu->addAction(readSpectrumAction);
    this->spectrumMenu->addAction(loadPointsAction);
    this->menuBar()->addMenu(this->configMenu);
    this->menuBar()->addMenu(this->helpMenu);
}
QString MainWindow::readLastWorkingDirectory(){
    QString filename=QCoreApplication::applicationDirPath();
    if(!filename.endsWith(QDir::separator ()))filename=filename+QDir::separator ();
    filename=filename+this->setup->getHistoryFileName();
    QFile file(QDir::toNativeSeparators(filename));
    QString line;
    if(file.exists()){
        file.open(QIODevice::ReadOnly);
        line=file.readLine();
        file.close();
    }
    return line;
}

void MainWindow::setupHistory(bool read){
    QString filename=QCoreApplication::applicationDirPath();
    if(!filename.endsWith(QDir::separator ()))filename=filename+QDir::separator ();
    filename=filename+this->setup->getHistoryFileName();
    if(read){
        QFile file(filename);
        QString line;
        if(file.exists()){
            file.open(QIODevice::ReadOnly);
            while (!file.atEnd()){
                line=file.readLine();
                if(line.trimmed()!="")history.append(line.trimmed());
            }
            file.close();
        }
    }
    if(historyActions.size()>0){
        foreach(QAction *action,historyActions){
            historyMenu->removeAction(action);
        }
    }
    historyActions.clear();
    QString currentDir=QDir::currentPath().trimmed();
    if(history.size()>0){
        int index=history.indexOf(currentDir);
        if(index!=-1)history.remove(index);
    }
    history.push_front(currentDir);
    if(history.size()>historyMax)history.remove(historyMax);
    QFile file(filename);
    file.open(QIODevice::WriteOnly);
    QTextStream in(&file);
    foreach(QString path,history){
        if(QFile(path).exists()){
            QAction *action=new QAction(historyMenu);
            action->setIcon(QIcon(":/images/history.png"));
            action->setIconVisibleInMenu(true);
            in<<path<<endl;
            action->setText(path);
            historyActions.append(action);
            historyMenu->addAction(action);
        }
    }
    file.close();
}

bool MainWindow::saveDialog(){
    int result=QMessageBox::question(this,"Save config","Do you want to save configs, bevor changing directory?",QMessageBox::Yes|QMessageBox::No|QMessageBox::Cancel);
    if (result == QMessageBox::Yes) {
        QSettings *settings=setup->getLocalSettings();
        QFile file(settings->fileName());
        file.remove();
        this->imageControllerWidget->save();
        this->save();
    }else if(result == QMessageBox::Cancel){
        return false;
    }
    return true;
}

void MainWindow::initialize(){
    QFile file(this->setup->getRadmc3dLogFilename());
    file.remove();
    this->setup->setLambda();
    if(this->setup->IsLambaRead()){
        this->imageControllerWidget->setUpLambdaIndexSlider();
        radmc3dProcess=NULL;
        if(QFile(this->setup->getradmc3dFilename()).exists()){
            this->setRadmc3dProcessLocalExisting(true);
            this->logger->writeToLogFile("radmc3d is local existing");
            radmc3dWriteGridFile(this->imageControllerWidget->isPipe());
        }else{
            this->setRadmc3dProcessLocalExisting(false);
            this->setup->setAmrGrid();
        }
        this->setup->checkRadmc3dlocal();
        this->findRadmc3dStarterPath();
        if(this->setup->isAmrGridRead()){
            loaded=false;
            initialized=true;
            this->imageControllerWidget->loadColorLookUpTables(QDir::currentPath(),true);
            this->imageControllerWidget->setUpColorTableComboBoxValues();
            this->imageControllerWidget->setUpColorTableSlider();
            makeImage=new MakeImage(this->imageControllerWidget);
            makeSpectrum=new MakeSpectrum(this->imageControllerWidget);
            this->checkLaststate();
            //This method should be called after setAmrGrid, because the image'size will be calculated from amrgrid file.
            this->colorChanged(this->colorAction->isChecked(),false);
            if(!loaded) this->imageControllerWidget->setSizeLineEdit();
            this->updateImageMenu();
            this->imageControllerWidget->updateImageMenuActionVisibilties();
            emit this->imageControllerWidget->positionComboxBoxValueChanged(this->imageControllerWidget->getPosition());
            this->setUpWhatsThis();
            if(imageControllerWidget->isAntiAliasing())
                this->transformMode=Qt::SmoothTransformation;
            else
                this->transformMode=Qt::FastTransformation;

            emit sendCommand(this->imageControllerWidget->getPipeCheckBoxObjectName());
        }else{
            initialized=false;
        }
    }else{
        initialized=false;
    }
    if(!initialized)showRadmc3dLog();
}
void MainWindow::changeCurrentDirectoryByAction(QAction* action){
    if(QFile(action->text()).exists()){
        if(action->text()!=QDir::currentPath()){
            if(saveDialog()){
                changeCurrentDirectory(action->text());
            }

        }else{
            if(QMessageBox::question(this,trUtf8("Do you want to reload setting?"),trUtf8("Do you want to reload setting?"),QMessageBox::Yes|QMessageBox::No)==QMessageBox::Yes){
                changeCurrentDirectory(action->text());
            }
        }
    }else{
        QMessageBox::information(this,trUtf8("Directory is not existing"),trUtf8("Directory is not existing %1").arg(action->text()));
        setupHistory(false);
    }
}

void MainWindow::changeCurrentDirectory(QString path){
    imageIsRead=false;
    initialized=false;
    showLog=false;
    this->spectrumLabel->setZoomActive(false);
    if((radmc3dProcess!=NULL) && (radmc3dProcess->state()==QProcess::Running))killRadmc3d();
    QDir::setCurrent(path);
    this->logger->createLogFile();
    this->logger->writeToLogFile(trUtf8("Current directory is %1").arg(QDir::currentPath()));
    initialized=false;
    this->setup->setLambdaRead(false);
    this->setup->setAmrGridRead(false);
    initialize();
    setupHistory(false);
    showLog=true;
}

void MainWindow::changeCurrentDirectory(){
    QFileDialog dialog;
    dialog.setFileMode(QFileDialog::Directory);
    dialog.setDirectory(QDir::currentPath());
    // dialog.setOption(QFileDialog::ShowDirsOnly);
    int result=dialog.exec();
    if(result==dialog.Accepted){
        if(QDir::currentPath()!=dialog.directory().absolutePath()){
            if(saveDialog()){
                QMessageBox::information(this,trUtf8("Directory is choosen"),trUtf8("Directory %1 is choosen").arg(dialog.directory().absolutePath()));
                changeCurrentDirectory(dialog.directory().absolutePath());
            }
        }
    }
}

void MainWindow::radmc3dWriteGridFile(bool pipe){
    mutex.lock();
    exit=false;
    mutex.unlock();
    this->logger->redirectToLogFile(false);
    this->logger->writeToLogFile(trUtf8("radmc3d write grid file, be patient"));
    this->logger->redirectToLogFile(true);
    if(!pipe){
        QFile file;
        bool amrfileExist=false;
        QTime lastModifiedTime;
        for (int var = 0; var < this->setup->getAmrGridFilenames().size(); ++var) {
            file.setFileName(this->setup->getAmrGridFilenames().at(var));
            if(file.exists()){
                amrfileExist=true;
                QFileInfo info(file);
                lastModifiedTime=info.lastModified().time();
                break;
            }
        }
        QString command="./"+this->setup->getradmc3dFilename()+" writegridfile";
        radmc3dProcess=new QProcess();
        radmc3dProcess->start(command);
        this->logger->writeToLogFile(command);
        bool started=radmc3dProcess->waitForStarted();
        if(!started){
            this->logger->writeToLogFile(trUtf8("radmc3d could not be started"));
            this->logger->stop(EXIT_FAILURE);
        }
        radmc3dProcess->waitForFinished(-1);
        if(radmc3dProcess->exitStatus()==QProcess::NormalExit){
            QString str=radmc3dProcess->readAllStandardOutput();
            QFile file(this->setup->getRadmc3dLogFilename());
            QTextStream in(&file);
            file.open(QIODevice::WriteOnly);
            in<<str;
            file.close();
        }
        if(amrfileExist){
            QFileInfo info(file);
            if(lastModifiedTime==info.lastModified().time()){
                this->logger->writeToLogFile(trUtf8("AMR file %1 is not recreated").arg(file.fileName()));
            }
        }
        this->setup->setAmrGrid();
    }else{
        QString command="./"+this->setup->getradmc3dFilename()+" child";
        radmc3dProcess=new QProcess();
        this->radmc3dProcess->setReadChannel(QProcess::StandardOutput);
        radmc3dProcess->start(command);
        this->logger->writeToLogFile(command);
        bool started=radmc3dProcess->waitForStarted();
        if(!started){
            this->logger->writeToLogFile(trUtf8("radmc3d could not be started"));
            this->logger->stop(EXIT_FAILURE);
        }
        QTextStream in(radmc3dProcess);
        in<<"writegridfile"<<endl;
        in<<"enter"<<endl;
        in<<"respondwhenready"<<endl;
        in<<"enter"<<endl;
        in.flush();
        radmc3dProcess->waitForBytesWritten();
        bool done=false;
        while(!done){
            QCoreApplication::processEvents();
            radmc3dProcess->waitForReadyRead(200);
            QString str=QString::fromLatin1(radmc3dProcess->readAll());
            if(str.toInt()==1){
                this->logger->writeToLogFile("radmc3d is done with writing grid file");
                done=true;
            }
            if(radmc3dProcess->state()==QProcess::NotRunning){
                this->logger->writeToLogFile(trUtf8("radmc3d not running more in child modus"));
            }
        }
        this->setup->setAmrGrid();
    }
}

void MainWindow::calculateContour(){
    if(imageIsRead){
        QVector<QVector <double> > pixelValues;
        QVector <double > contourLevel;
        double min=1,max=-1;
        int total=this->imageControllerWidget->getNumberOfContours();
        if(!this->imageControllerWidget->isRGBModus()){
            if(lastImage.height()<400)
                scaleImage=lastImage.scaled(400,400,
                                            Qt::KeepAspectRatio,Qt::SmoothTransformation).
                        convertToFormat(QImage::Format_Indexed8,lastImage.colorTable());
            else
                scaleImage=lastImage.convertToFormat(QImage::Format_Indexed8,lastImage.colorTable());


            pixelValues.resize(scaleImage.height());
            for (int iy = 0; iy < scaleImage.height(); ++iy){
                pixelValues[iy].resize(scaleImage.width());
                for (int ix = 0; ix < scaleImage.width(); ++ix){
                    pixelValues[iy][ix]=scaleImage.pixelIndex(ix,iy);
                }
                double minvalue=*min_element(pixelValues[iy].begin(),pixelValues[iy].end());
                double maxvalue=*max_element(pixelValues[iy].begin(),pixelValues[iy].end());
                if(iy==0||minvalue<min)min=minvalue;
                if(iy==0||maxvalue>max)max=maxvalue;
            }
        }else{
            if(lastImage.height()<400)
                scaleImage=lastImage.scaled(400,400,
                                            Qt::KeepAspectRatio,Qt::SmoothTransformation).
                        convertToFormat(QImage::Format_RGB888);
            else
                scaleImage=lastImage.convertToFormat(QImage::Format_RGB888);
            pixelValues.resize(scaleImage.height());
            for (int iy = 0; iy < scaleImage.height(); ++iy){
                pixelValues[iy].resize(scaleImage.width());
                for (int ix = 0; ix < scaleImage.width(); ++ix){
                    QRgb rgb=scaleImage.pixel(ix,iy);
                    pixelValues[iy][ix]=qRed(rgb)*0.84375 + qGreen(rgb)*0.140625 + qBlue(rgb)*0.0234375;
                }
            }
            for (int iy = 0; iy < scaleImage.height(); ++iy){
                double minvalue=*min_element(pixelValues[iy].begin(),pixelValues[iy].end());
                double maxvalue=*max_element(pixelValues[iy].begin(),pixelValues[iy].end());
                if(iy==0||minvalue<min)min=minvalue;
                if(iy==0||maxvalue>max)max=maxvalue;
            }

        }
        if(max-min>1){
            double value;
            for (int var = 1; var <= total; ++var) {
                value=((max-min)*var*1./(total+1.))+min;
                if(contourLevel.indexOf(value)==-1)
                    contourLevel.append(value);
            }
            contourPoints.clear();
            contourPoints=this->Contour(pixelValues,0,scaleImage.height()-1,0,scaleImage.width()-1,
                                        image->createY(scaleImage.height()),
                                        image->createX(scaleImage.width()),contourLevel.size(),contourLevel);
        }
    }

}

/**
 * @ written by paul bourke, adapte to QT by Farzin Sereshti
 *  Derivation from the fortran version of CONREC by Paul Bourke
 *  d               ! matrix of data to contour
 *  ilb,iub,jlb,jub ! index bounds of data matrix
 *  x               ! data matrix column coordinates
 *  y               ! data matrix row coordinates
 *  nc              ! number of contour levels
 *  z               ! contour levels in increasing order
*/
QVector <QVector<QPoint> > MainWindow:: Contour(QVector < QVector<double> > d,int ilb,int iub,int jlb,int jub,
                                                QVector<double> x, QVector<double> y,int nc,QVector<double> z)
{
#define xsect(p1,p2) (h[p2]*xh[p1]-h[p1]*xh[p2])/(h[p2]-h[p1])
#define ysect(p1,p2) (h[p2]*yh[p1]-h[p1]*yh[p2])/(h[p2]-h[p1])
    QVector< QVector<QPoint> > points;
    int m1,m2,m3,case_value;
    double dmin,dmax,x1=0,x2=0,y1=0,y2=0;
    int i,j,k,m;
    double h[5];
    int sh[5];
    double xh[5],yh[5];
    int im[4] = {0,1,1,0},jm[4]={0,0,1,1};
    int castab[3][3][3] = {
        { {0,0,8},{0,2,5},{7,6,9} },
        { {0,3,4},{1,3,1},{4,3,0} },
        { {9,6,7},{5,2,0},{8,0,0} }
    };
    double temp1,temp2;

    for (j=(jub-1);j>=jlb;j--) {
        QCoreApplication::processEvents();
        for (i=ilb;i<=iub-1;i++) {
            temp1 = min(d[i][j],d[i][j+1]);
            temp2 = min(d[i+1][j],d[i+1][j+1]);
            dmin  = min(temp1,temp2);
            temp1 = max(d[i][j],d[i][j+1]);
            temp2 = max(d[i+1][j],d[i+1][j+1]);
            dmax  = max(temp1,temp2);
            if (dmax < z[0] || dmin > z[nc-1])
                continue;
            for (k=0;k<nc;k++) {
                QVector<QPoint> vec;
                if (z[k] < dmin || z[k] > dmax)
                    continue;
                for (m=4;m>=0;m--) {
                    if (m > 0) {
                        h[m]  = d[i+im[m-1]][j+jm[m-1]]-z[k];
                        xh[m] = x[i+im[m-1]];
                        yh[m] = y[j+jm[m-1]];
                    } else {
                        h[0]  = 0.25 * (h[1]+h[2]+h[3]+h[4]);
                        xh[0] = 0.50 * (x[i]+x[i+1]);
                        yh[0] = 0.50 * (y[j]+y[j+1]);
                    }
                    if (h[m] > 0.0)
                        sh[m] = 1;
                    else if (h[m] < 0.0)
                        sh[m] = -1;
                    else
                        sh[m] = 0;
                }

                /*
               Note: at this stage the relative heights of the corners and the
               centre are in the h array, and the corresponding coordinates are
               in the xh and yh arrays. The centre of the box is indexed by 0
               and the 4 corners by 1 to 4 as shown below.
               Each triangle is then indexed by the parameter m, and the 3
               vertices of each triangle are indexed by parameters m1,m2,and m3.
               It is assumed that the centre of the box is always vertex 2
               though this isimportant only when all 3 vertices lie exactly on
               the same contour level, in which case only the side of the box
               is drawn.
                  vertex 4 +-------------------+ vertex 3
                           | \               / |
                           |   \    m-3    /   |
                           |     \       /     |
                           |       \   /       |
                           |  m=2    X   m=2   |       the centre is vertex 0
                           |       /   \       |
                           |     /       \     |
                           |   /    m=1    \   |
                           | /               \ |
                  vertex 1 +-------------------+ vertex 2
            */
                /* Scan each triangle in the box */
                for (m=1;m<=4;m++) {
                    m1 = m;
                    m2 = 0;
                    if (m != 4)
                        m3 = m + 1;
                    else
                        m3 = 1;
                    if ((case_value = castab[sh[m1]+1][sh[m2]+1][sh[m3]+1]) == 0)
                        continue;
                    switch (case_value) {
                    case 1: /* Line between vertices 1 and 2 */
                        x1 = xh[m1];
                        y1 = yh[m1];
                        x2 = xh[m2];
                        y2 = yh[m2];
                        break;
                    case 2: /* Line between vertices 2 and 3 */
                        x1 = xh[m2];
                        y1 = yh[m2];
                        x2 = xh[m3];
                        y2 = yh[m3];
                        break;
                    case 3: /* Line between vertices 3 and 1 */
                        x1 = xh[m3];
                        y1 = yh[m3];
                        x2 = xh[m1];
                        y2 = yh[m1];
                        break;
                    case 4: /* Line between vertex 1 and side 2-3 */
                        x1 = xh[m1];
                        y1 = yh[m1];
                        x2 = xsect(m2,m3);
                        y2 = ysect(m2,m3);
                        break;
                    case 5: /* Line between vertex 2 and side 3-1 */
                        x1 = xh[m2];
                        y1 = yh[m2];
                        x2 = xsect(m3,m1);
                        y2 = ysect(m3,m1);
                        break;
                    case 6: /* Line between vertex 3 and side 1-2 */
                        x1 = xh[m3];
                        y1 = yh[m3];
                        x2 = xsect(m1,m2);
                        y2 = ysect(m1,m2);
                        break;
                    case 7: /* Line between sides 1-2 and 2-3 */
                        x1 = xsect(m1,m2);
                        y1 = ysect(m1,m2);
                        x2 = xsect(m2,m3);
                        y2 = ysect(m2,m3);
                        break;
                    case 8: /* Line between sides 2-3 and 3-1 */
                        x1 = xsect(m2,m3);
                        y1 = ysect(m2,m3);
                        x2 = xsect(m3,m1);
                        y2 = ysect(m3,m1);
                        break;
                    case 9: /* Line between sides 3-1 and 1-2 */
                        x1 = xsect(m3,m1);
                        y1 = ysect(m3,m1);
                        x2 = xsect(m1,m2);
                        y2 = ysect(m1,m2);
                        break;
                    default:
                        break;
                    }
                    int posx1=-1;
                    double nearst;
                    for (int var = 0; var < iub; ++var) {
                        if(posx1==-1|| (fabs(x.at(var)-x1)<fabs(nearst-x1))){
                            nearst=x.at(var);
                            posx1=var;
                        }
                    }
                    int posx2=-1;

                    for (int var = 0; var < iub; ++var) {
                        if(posx2==-1|| (fabs(x.at(var)-x2)<fabs(nearst-x2))){
                            nearst=x.at(var);
                            posx2=var;
                        }
                    }
                    int posy1=-1;
                    for (int var = 0; var < jub; ++var) {
                        if(posy1==-1|| (fabs(y.at(var)-y1)<fabs(nearst-y1))){
                            nearst=y.at(var);
                            posy1=var;
                        }
                    }
                    int posy2=-1;

                    for (int var = 0; var < jub; ++var) {
                        if(posy2==-1|| (fabs(y.at(var)-y2)<fabs(nearst-y2))){
                            nearst=y.at(var);
                            posy2=var;
                        }
                    }

                    vec.append(QPoint(posy1,posx1));
                    vec.append(QPoint(posy2,posx2));



                } /* m */
                points.append(vec);
            } /* k - contour */
        } /* i */
    } /* j */
    return points;
}

void MainWindow::startRadmc3dAsChild(){
    showProgressBar();
    disconnect(&futureWacther);
    if( radmc3dProcess==NULL || radmc3dProcess->state()==QProcess::NotRunning){
        mutex.lock();
        exit=false;
        mutex.unlock();
        this->radmc3dProcess=new QProcess();
        QFile radmc3dfilelocal(this->setup->getradmc3dFilename());
        QFile radmc3dfileglobal(this->radmc3dStarterPATH);
        if(!(radmc3dfilelocal.exists()||radmc3dfileglobal.exists())){
            this->logger->writeToLogFile(tr("RADMC3D executable file does not exists in current directory and in PATH"));
            this->imageIsRead=false;
        }
        this->radmc3dProcess->setReadChannel(QProcess::StandardOutput);
        if(radmc3dProcessIsLocalExisting){
            this->radmc3dProcess->start("./"+this->setup->getradmc3dFilename()+" child");
            this->logger->writeToLogFile("./"+this->setup->getradmc3dFilename()+" child");
        }
        else{
            this->radmc3dProcess->start(this->setup->getradmc3dFilename()+" child");
            this->logger->writeToLogFile(this->setup->getradmc3dFilename()+" child");
        }
        if(!radmc3dProcess->waitForStarted()){
            this->logger->writeToLogFile(tr("RADMC3D is not starting, try again"));
            this->imageIsRead=false;
        }

    }
    if(this->radmc3dProcess->state()!=QProcess::NotRunning){
        this->setCancelable(true);
        if(!this->imageControllerWidget->isRenderSpectrumActive()){
            QString command=this->setup->getradmc3dFilename()+" "+makeImage->getImageCommand();
            this->logger->writeToLogFile(command);
            this->makeImage->sendMakeImageCommandToRadmc3d(this->radmc3dProcess);
            this->makeImage->sendWriteImageCommandToRadmc3d(this->radmc3dProcess);
            imageIsRead=false;
            this->image=new Image();
            this->thread=new QThread();
            this->image->moveToThread(thread);
            this->connect(thread,SIGNAL(started()),this,SLOT(readImage()));
            this->connect(thread, SIGNAL(finished()), this->thread, SLOT(deleteLater()));
            this->connect(this->image,SIGNAL(setCancelable(bool)),this,SLOT(setCancelable(bool)));
            this->connect(image,SIGNAL(finished()),this,SLOT(hideProgressBarAndEditAndReload()));
            this->thread->start();
        }else{
            this->makeSpectrum->sendMakeSpectrumCommandToRadmc3d(this->radmc3dProcess);
            this->makeSpectrum->sendWriteSpectrumCommandToRadmc3d(this->radmc3dProcess);
            spectrumIsRead=false;
            this->spectrum=new Spectrum();
            this->thread=new QThread();
            this->spectrum->moveToThread(thread);
            this->connect(thread,SIGNAL(started()),this,SLOT(readSpectrum()));
            this->connect(thread, SIGNAL(finished()), this->thread, SLOT(deleteLater()));
            this->connect(this->spectrum,SIGNAL(setCancelable(bool)),this,SLOT(setCancelable(bool)));
            this->connect(spectrum,SIGNAL(finished()),this,SLOT(hideProgressBarAndEditAndReload()));
            this->thread->start();
        }
    }else{
        this->logger->writeToLogFile(trUtf8("RADMC3D is not running"));
        hideProgressBarAndShowRadmc3dLog();
    }
}

void MainWindow::startRadmc3d(){
    if(this->imageControllerWidget->isRenderSpectrumActive()){
        this->spectrumIsRead=false;
        this->renameSpectrumFile();
    }else{
        this->imageIsRead=false;
        this->renameImageFile();
    }
    if(exit==false){
        sendRadmc3dExitCommand();
    }
    mutex.lock();
    exit=true;
    mutex.unlock();
    this->radmc3dProcess=new QProcess();
    this->connect(radmc3dProcess,SIGNAL(started()),this,SLOT(showProgressBar()));
    this->connect(radmc3dProcess,SIGNAL(finished(int,QProcess::ExitStatus)),this,SLOT(radmc3dProcessFinished(int,QProcess::ExitStatus)));
    QString command="";
    if(!this->imageControllerWidget->isRenderSpectrumActive()){
        command=this->setup->getradmc3dFilename()+" "+makeImage->getImageCommand();
    }else{
        command=this->setup->getradmc3dFilename()+" "+makeSpectrum->getSpectrumCommand();
    }
    if(this->IsRadmc3dProcessLocalExisting())command="./"+command;
    this->logger->writeToLogFile(command);
    QFile radmc3dfilelocal(this->setup->getradmc3dFilename());
    QFile radmc3dfileglobal(this->radmc3dStarterPATH);
    if(!(radmc3dfilelocal.exists()||radmc3dfileglobal.exists())){
        this->logger->writeToLogFile(tr("RADMC3D executable file does not exists in current directory and in PATH"));
        imageIsRead=false;
        spectrumIsRead=false;
    }
    radmc3dProcess->setReadChannel(QProcess::StandardOutput);
    radmc3dProcess->start(command);
    if(!radmc3dProcess->waitForStarted()){
        this->logger->writeToLogFile(tr("RADMC3D is not starting, try again"));
        imageIsRead=false;
        spectrumIsRead=false;
        disconnect(radmc3dProcess);
        hideProgressBarAndShowRadmc3dLog();
    }
    this->setCancelable(true);
    if(radmc3dProcess->error()==QProcess::FailedToStart){
        this->logger->writeToLogFile("RADMC3DProcess could not be started");
        imageIsRead=false;
        spectrumIsRead=false;
        disconnect(radmc3dProcess);
        hideProgressBarAndShowRadmc3dLog();
    }

}


void MainWindow:: radmc3dProcessFinished(int exitCode,QProcess::ExitStatus status){
    hideProgressBar();
    if(status==QProcess::NormalExit){
        QString str=radmc3dProcess->readAllStandardOutput();
        QFile file(this->setup->getRadmc3dLogFilename());
        QTextStream in(&file);
        file.open(QIODevice::WriteOnly);
        in<<str;
        file.close();
        this->setCancelable(false);
        if(!this->imageControllerWidget->isRenderSpectrumActive()){
            readImage();
            if(this->imageIsRead){
                if(this->imageControllerWidget->isDoLine()){
                    this->imageControllerWidget->setDoLine(false);
                    double lambda=this->image->getLambdaAt(0);
                    this->imageControllerWidget->setRedLambdaLineEdit(lambda);
                }
            }
            this->imageControllerWidget->setDoLine(false);
            if(imageIsRead){
                pixelEdit();
                this->changeImageColorLookUpTable(this->imageControllerWidget->isColorInvert());
                if(this->imageControllerWidget->isContour())calculateContour();
                reloadImage();
            }else{
                showRadmc3dLog();
            }
        }else{
            if(readSpectrum()){
                this->spectrumIsRead=true;
                reloadSpectrum();
            }else{
                showRadmc3dLog();
            }
        }
    }else{
        this->logger->writeToLogFile(trUtf8("The process is crashed _CODE = %1").arg(exitCode));
        showRadmc3dLog();
    }
    waitForFinishing=false;

}
bool MainWindow::readSpectrum(){
    if(!this->imageControllerWidget->isPipe()){
        spectrum=new Spectrum();
        if(!spectrum->readSpectrumFormatted(this->setup->getSpectrumFilename(),this->setup->getSpectrumEnding())){
            spectrumIsRead=false;
            showRadmc3dLog();
        }else{
            spectrumLabel->setSpectrum(spectrum);
            this->spectrumIsRead=true;
            return true;
        }
    }else{
        QString command=this->setup->getradmc3dFilename()+" "+makeSpectrum->getSpectrumCommand();
        this->logger->writeToLogFile(command);
        emit spectrum->readSpectrumFormattedFromRadmc3d(radmc3dProcess,&exit);
    }
    return true;
}

void MainWindow::showRadmc3dLog(){
    if(showLog){
        hideProgressBar();
        QFile file(this->setup->getRadmc3dLogFilename());
        QTextStream in(&file);
        file.open(QIODevice::ReadOnly);
        QString content=in.readAll();
        textEdit->setText(content);
        textEdit->setReadOnly(true);
        this->leftWidget->setCurrentIndex(this->leftWidget->indexOf(textGroup));
        QTextCursor c = textEdit->textCursor();
        c.movePosition(QTextCursor::End);
        textEdit->setTextCursor(c);
    }
    emit this->imageControllerWidget->selectDebugger();
}
void MainWindow::renameImageFile(){
    QFile file(this->setup->getImageFilename()+this->setup->getImageFileEnding());
    QString filename=file.fileName()+".backup";
    QFile refile(filename);
    if(refile.exists())refile.remove();
    if(file.exists() && !file.rename(filename)){
        imageIsRead=false;
        this->logger->writeToLogFile(trUtf8("imagefile cannnot be renamed to  %1").arg(filename));
    }
}
void MainWindow::renameSpectrumFile(){
    QFile file(this->setup->getSpectrumFilename()+this->setup->getSpectrumEnding());
    QString filename=file.fileName()+".backup";
    QFile refile(filename);
    if(refile.exists())refile.remove();
    if(file.exists() && !file.rename(filename)){
        spectrumIsRead=false;
        this->logger->writeToLogFile(trUtf8("spectrumfile cannnot be renamed to  %1").arg(filename));
    }
}

void MainWindow::readImage(){
    imageIsRead=false;
    if(!this->imageControllerWidget->isPipe()){
        image=new Image();
        if(!image->readImageFormatted(this->setup->getImageFilename(),this->setup->getImageFileEnding())){
            imageIsRead=false;
            showRadmc3dLog();
        }else{
            this->imageIsRead=true;
            imageLabel->setImage(image);
            this->imageControllerWidget->setNumberOfPixelsSpinBox(image->getNumberOfPixelX());
        }
    }else{
        emit image->readImageFormattedFromRadmc3d(radmc3dProcess,&exit);
    }
}
void MainWindow::pixelEdit(){
    if(imageIsRead){
        hideLog();
        int nx=image->getNumberOfPixelX();
        int ny=image->getNumberOfPixelY();
        QSize size(nx,ny);
        QImage::Format format;
        if(imageControllerWidget->isRGBModus())
            format=QImage::Format_RGB32;
        else
            format=QImage::Format_Indexed8;
        QImage qImage(size,format );
        if(!imageControllerWidget->isRGBModus())
            qImage.setColorTable(this->imageControllerWidget->getCurrentColorLookUpTable());

        QVector< QVector <QVector < double> > >intensity;
        if(image->getImages().size()>0){
            intensity.append(image->getImageAt(0));
            if(imageControllerWidget->isRGBModus()){
                if(image->getNumberOfImages()==3){
                    intensity.append(image->getImageAt(1));
                    intensity.append(image->getImageAt(2));
                }
            }
        }else
            imageIsRead=false;
        if(imageIsRead){
            QVector <double> minValues;
            QVector <double> maxValues;
            double maxValue=-1;
            double minValue=-1;

            minValues.resize(intensity.size());
            maxValues.resize(intensity.size());
            QVector< QVector <double> > lgrange;
            lgrange.clear();
            if(this->imageControllerWidget->isRGBModus()&&this->imageControllerWidget->isAbsoluteScale()){

                QVector<double> vec1,vec0;
                vec1.append(log10(this->imageControllerWidget->getRedColorTuneValueLineEditText().toDouble()));
                vec1.append(log10(this->imageControllerWidget->getGreenColorTuneValueLineEditText().toDouble()));
                vec1.append(log10(this->imageControllerWidget->getBlueColorTuneValueLineEditText().toDouble()));
                vec0.append(vec1.at(0)-this->imageControllerWidget->getMaxLog());
                vec0.append(vec1.at(1)-this->imageControllerWidget->getMaxLog());
                vec0.append(vec1.at(2)-this->imageControllerWidget->getMaxLog());
                lgrange.append(vec0);
                lgrange.append(vec1);

            }

            double saturate=this->imageControllerWidget->getSaturate();
            if(saturate<1.){
                for (int at = 0; at < intensity.size(); ++at){
                    for (int index = 0; index < ny; ++index){
                        double* max=max_element(intensity[at][index].begin(),intensity[at][index].end());
                        if(index==0 && at==0)maxValue=*max;
                        if((*max)>maxValue)maxValue=*max;
                    }
                }

                saturate=saturate*maxValue;
                for (int  j= 0; j < intensity.size(); ++j){
                    for (int iy = 0; iy < ny; ++iy){
                        for (int ix = 0; ix < nx; ++ix){
                            if(intensity[j][iy][ix]>saturate)intensity[j][iy][ix]=saturate;
                        }
                    }
                }
            }

            if(this->imageControllerWidget->isLinear()){
                if(lgrange.size()==0){
                    for (int h = 0; h < intensity.size(); ++h){
                        for (int index = 0; index < ny; ++index){
                            double* max=max_element(intensity[h][index].begin(),intensity[h][index].end());
                            if(index==0 )maxValue=*max;
                            if((*max)>maxValue)maxValue=*max;
                        }
                        maxValues[h]=maxValue;
                        minValues[h]=0;
                    }
                }else{
                    minValues=lgrange.at(0);
                    maxValues=lgrange.at(1);
                    for (int v = 0; v < intensity.size(); ++v){
                        minValues[v]=pow(10.,minValues[v]);
                        maxValues[v]=pow(10.,maxValues[v]);
                    }
                }
                if(this->imageControllerWidget->isRGBModus()&& !this->imageControllerWidget->isAbsoluteScale()){
                    bool ok;
                    double redColorTune,blueColorTune,greenColorTune;
                    redColorTune=this->imageControllerWidget->getRedColorTuneValueLineEditText().toDouble(&ok);
                    if(ok)
                        blueColorTune=this->imageControllerWidget->getBlueColorTuneValueLineEditText().toDouble(&ok);
                    if(ok)
                        greenColorTune=this->imageControllerWidget->getGreenColorTuneValueLineEditText().toDouble(&ok);
                    if(ok){
                        if(minValues.size()==3){
                            minValues[0]=minValues[0]/redColorTune;
                            minValues[1]=minValues[1]/greenColorTune;
                            minValues[2]=minValues[2]/blueColorTune;
                            maxValues[0]=maxValues[0]/redColorTune;
                            maxValues[1]=maxValues[1]/greenColorTune;
                            maxValues[2]=maxValues[2]/blueColorTune;
                        }
                    }else{
                        this->logger->writeToLogFile(trUtf8("one of the colortunes value is invalid red(%1) blue(%2) green(%3)").arg(
                                                         this->imageControllerWidget->getRedColorTuneValueLineEditText(),
                                                         this->imageControllerWidget->getBlueColorTuneValueLineEditText(),
                                                         this->imageControllerWidget->getGreenColorTuneValueLineEditText()));
                    }
                }
            }else{
                for (int w = 0; w < intensity.size(); ++w){
                    for (int iy = 0; iy < ny; ++iy){
                        for (int ix = 0; ix < nx; ++ix){
                            intensity[w][iy][ix]=log10(intensity[w][iy][ix]+pow(10.,-200));
                        }
                    }
                }
                if(lgrange.size()==0){
                    for (int z = 0; z < intensity.size(); ++z){
                        for (int index = 0; index < ny; ++index){
                            double* max=max_element(intensity[z][index].begin(),intensity[z][index].end());
                            double* min=min_element(intensity[z][index].begin(),intensity[z][index].end());
                            if(index==0)maxValue=*max;
                            if(index==0)minValue=*min;
                            if((*max)>maxValue)maxValue=*max;
                            if((*min)<minValue)minValue=*min;
                        }
                        maxValues[z]=maxValue;
                        minValues[z]=minValue;
                    }
                    bool setMinValues=false;
                    if(this->imageControllerWidget->getMaxLog()>0){
                        for (int m = 0; m < intensity.size(); ++m){
                            if(fabs(maxValues[m]-minValues[m])>this->imageControllerWidget->getMaxLog()){
                                setMinValues=true;
                                break;
                            }
                        }
                    }

                    if(setMinValues){
                        for (int o = 0; o < intensity.size(); ++o){
                            minValues[o]=maxValues[o]-this->imageControllerWidget->getMaxLog();
                        }

                    }
                    if(this->imageControllerWidget->isRGBModus()&& !this->imageControllerWidget->isAbsoluteScale()){
                        bool ok;
                        double redColorTune,blueColorTune,greenColorTune;
                        redColorTune=this->imageControllerWidget->getRedColorTuneValueLineEditText().toDouble(&ok);
                        if(ok)
                            blueColorTune=this->imageControllerWidget->getBlueColorTuneValueLineEditText().toDouble(&ok);
                        if(ok)
                            greenColorTune=this->imageControllerWidget->getGreenColorTuneValueLineEditText().toDouble(&ok);
                        if(ok){
                            if(minValues.size()==3){
                                minValues[0]=minValues[0]-log10(redColorTune);
                                minValues[1]=minValues[1]-log10(greenColorTune);
                                minValues[2]=minValues[2]-log10(blueColorTune);
                                maxValues[0]=maxValues[0]-log10(redColorTune);
                                maxValues[1]=maxValues[1]-log10(greenColorTune);
                                maxValues[2]=maxValues[2]-log10(blueColorTune);
                            }
                        }else{
                            this->logger->writeToLogFile(trUtf8("one of the colortunes value is invalid red(%1) blue(%2) green(%3)").arg(
                                                             this->imageControllerWidget->getRedColorTuneValueLineEditText(),
                                                             this->imageControllerWidget->getBlueColorTuneValueLineEditText(),
                                                             this->imageControllerWidget->getGreenColorTuneValueLineEditText()));
                        }
                    }
                }else{
                    minValues=lgrange.at(0);
                    maxValues=lgrange.at(1);
                }

            }
            for (int t = 0; t < intensity.size(); ++t){
                for (int iy = 0; iy < ny; ++iy){
                    for (int ix = 0; ix < nx; ++ix){
                        if(intensity[t][iy][ix]>maxValues[t])intensity[t][iy][ix]=maxValues[t];
                        if(intensity[t][iy][ix]<minValues[t])intensity[t][iy][ix]=minValues[t];
                    }
                }
            }

            for (int q = 0; q < intensity.size(); ++q){
                for (int iy = 0; iy < ny; ++iy){
                    for (int ix = 0; ix < nx; ++ix){
                        intensity[q][iy][ix]=((intensity[q][iy][ix]-minValues[q])/(maxValues[q]-minValues[q]));
                        intensity[q][iy][ix]=intensity[q][iy][ix]*255;
                    }

                }
            }
            if(imageControllerWidget->isRGBModus()){
                if(intensity.size()==3){
                    for (int iy = 0; iy < ny; ++iy){
                        for (int ix = 0; ix < nx; ++ix){
                            QRgb  pixelvalue=qRgb(intensity[0][iy][ix],intensity[1][iy][ix],intensity[2][iy][ix]);
                            qImage.setPixel(ix,iy,pixelvalue);
                        }
                    }
                }

            }else{
                for (int iy = 0; iy < ny; ++iy){
                    for (int ix = 0; ix < nx; ++ix){
                        qImage.setPixel(ix,iy,intensity[0][iy][ix]);
                    }
                }

            }
            if(this->imageControllerWidget->isRGBModus()){
                if(this->imageControllerWidget->isLinear()){
                    this->imageControllerWidget->setLastColorMaximumArray(maxValues);
                    this->imageControllerWidget->setLastColorMinimumArray(minValues);
                }else{
                    if(maxValues.size()==3){
                        maxValues[0]=pow(10.,maxValues[0]);
                        maxValues[1]=pow(10.,maxValues[1]);
                        maxValues[2]=pow(10.,maxValues[2]);
                        minValues[0]=pow(10.,minValues[0]);
                        minValues[1]=pow(10.,minValues[1]);
                        minValues[2]=pow(10.,minValues[2]);
                        this->imageControllerWidget->setLastColorMaximumArray(maxValues);
                        this->imageControllerWidget->setLastColorMinimumArray(minValues);
                    }
                }
            }
            lastImage=qImage;
        }
    }
}
void MainWindow::changeImageColorLookUpTable(bool invert){
    if(!imageControllerWidget->isRGBModus())
        this->lastImage.setColorTable(this->imageControllerWidget->getCurrentColorLookUpTable());
    if(invert)this->lastImage.invertPixels(QImage::InvertRgb);

}
void MainWindow::reloadSpectrum(){
    if(this->spectrumIsRead){
        int index=this->leftWidget->currentIndex();
        this->imageControllerWidget->setRenderSpectrumActive(true);
        if(index!=this->leftWidget->indexOf(this->spectrumGroup))this->leftWidget->setCurrentWidget(spectrumGroup);
        int size;
        int width=this->width();
        int height=this->height();
        if(width>height){
            if(width/2>height)
                size=height;
            else
                size=width/2;
        }else{
            if(height/2>width)
                size=width;
            else
                size=height/2;
        }

        QPainter painter;
        QPixmap pixmap(size,size);
        QVector<Point> lpoints;
        lpoints.clear();
        if(loadPointsAction->isChecked())lpoints=loadedPoints;
        QVector<double> xcoord,ycoord;
        if(this->spectrumLabel->getfrequencyUnit()==Spectrum::MICRON){
            double fact=this->physicalConstants->getLightSpeed()/this->physicalConstants->getMicron();
            xcoord=this->spectrum->getFrequencies();
            int total= this->spectrum->getNumberOfFrequencies();
            for (int var = 0; var <total; ++var) {
                xcoord[var]=fact/xcoord[var];
            }
            if(loadPointsAction->isChecked()){
                for (int var = 0; var <lpoints.size(); ++var) {
                    lpoints[var].x=fact/lpoints[var].x;
                }
                if(isErrorbar){
                    for (int var = 0; var <lpoints.size(); ++var) {
                        lpoints[var].x1=fact/lpoints[var].x1;
                        lpoints[var].x2=fact/lpoints[var].x2;
                    }
                }
            }
        }else if(this->spectrumLabel->getfrequencyUnit()==Spectrum::HZ){
            xcoord=this->spectrum->getFrequencies();

        }else if(this->spectrumLabel->getfrequencyUnit()==Spectrum::EV){
            double h=this->physicalConstants->getPlanckconstant();
            double ev=this->physicalConstants->getElectronVolt();
            double h_ev=h/ev;
            xcoord=this->spectrum->getFrequencies();
            int total= this->spectrum->getNumberOfFrequencies();
            for (int var = 0; var <total; ++var) {
                xcoord[var]=h_ev*xcoord[var];
            }
            if(loadPointsAction->isChecked()){
                for (int var = 0; var <lpoints.size(); ++var) {
                    lpoints[var].x=h_ev*lpoints[var].x;
                }
                if(isErrorbar){
                    for (int var = 0; var <lpoints.size(); ++var) {
                        lpoints[var].x1=h_ev*lpoints[var].x1;
                        lpoints[var].x2=h_ev*lpoints[var].x2;
                    }
                }
            }

        }else if(this->spectrumLabel->getfrequencyUnit()==Spectrum::KEV){
            double h=this->physicalConstants->getPlanckconstant();
            double kev=this->physicalConstants->getKiloElectronVolt();
            double h_kev=h/kev;
            xcoord=this->spectrum->getFrequencies();
            int total= this->spectrum->getNumberOfFrequencies();
            for (int var = 0; var <total; ++var) {
                xcoord[var]=h_kev*xcoord[var];
            }
            if(loadPointsAction->isChecked()){
                for (int var = 0; var <lpoints.size(); ++var) {
                    lpoints[var].x=h_kev*lpoints[var].x;
                }
                if(isErrorbar){
                    for (int var = 0; var <lpoints.size(); ++var) {
                        lpoints[var].x1=h_kev*lpoints[var].x1;
                        lpoints[var].x2=h_kev*lpoints[var].x2;
                    }
                }
            }

        }
        int total= this->spectrum->getNumberOfFrequencies();
        double distanceFactor;
        if(this->spectrumLabel->getFluxUnit()>Spectrum::FNUJY){
            distanceFactor=4*physicalConstants->getPi()*pow(this->physicalConstants->getParsec(),2.);
        }else{
            distanceFactor=1.;
        }
        if(this->spectrumLabel->getFluxUnit()==Spectrum::NUFNU){
            ycoord=this->spectrum->getFrequencies();
            for (int var = 0; var <total; ++var) {
                ycoord[var]=ycoord[var]*this->spectrum->getFlux().at(var)* distanceFactor;
            }
            if(loadPointsAction->isChecked()){
                for (int var = 0; var <lpoints.size(); ++var) {
                    lpoints[var].y=loadedPoints[var].x*lpoints[var].y* distanceFactor;
                }
                if(isErrorbar){
                    for (int var = 0; var <lpoints.size(); ++var) {
                        lpoints[var].y1=loadedPoints[var].x*lpoints[var].y1* distanceFactor;
                        lpoints[var].y2=loadedPoints[var].x*lpoints[var].y2* distanceFactor;
                    }
                }
            }
        }else if(this->spectrumLabel->getFluxUnit()==Spectrum::NUFNUJY){
            ycoord=this->spectrum->getFrequencies();
            for (int var = 0; var <total; ++var) {
                ycoord[var]=ycoord[var]*this->spectrum->getFlux().at(var)*1e23* distanceFactor;
            }
            if(loadPointsAction->isChecked()){
                for (int var = 0; var <lpoints.size(); ++var) {
                    lpoints[var].y=loadedPoints[var].x*lpoints[var].y*1e23* distanceFactor;
                }
                if(isErrorbar){
                    for (int var = 0; var <lpoints.size(); ++var) {
                        lpoints[var].y1=loadedPoints[var].x*lpoints[var].y1*1e23* distanceFactor;
                        lpoints[var].y2=loadedPoints[var].x*lpoints[var].y2*1e23* distanceFactor;
                    }
                }
            }
        }else if(this->spectrumLabel->getFluxUnit()==Spectrum::FNU){

            for (int var = 0; var <total; ++var) {
                ycoord.append(this->spectrum->getFlux().at(var)* distanceFactor);
            }
            if(loadPointsAction->isChecked()){
                for (int var = 0; var <lpoints.size(); ++var) {
                    lpoints[var].y=lpoints[var].y* distanceFactor;
                }
                if(isErrorbar){
                    for (int var = 0; var <lpoints.size(); ++var) {
                        lpoints[var].y1=lpoints[var].y1* distanceFactor;
                        lpoints[var].y2=lpoints[var].y2* distanceFactor;
                    }
                }
            }
        }else if(this->spectrumLabel->getFluxUnit()==Spectrum::FNUJY){
            for (int var = 0; var <total; ++var) {
                ycoord.append(this->spectrum->getFlux().at(var)*1e23);
            }
            if(loadPointsAction->isChecked()){
                for (int var = 0; var <lpoints.size(); ++var) {
                    lpoints[var].y=lpoints[var].y *1e23;
                }
                if(isErrorbar){
                    for (int var = 0; var <lpoints.size(); ++var) {
                        lpoints[var].y1=lpoints[var].y1 *1e23;
                        lpoints[var].y2=lpoints[var].y2 *1e23;
                    }
                }
            }
        }else if(this->spectrumLabel->getFluxUnit()==Spectrum::NULNU){
            ycoord=this->spectrum->getFrequencies();
            for (int var = 0; var <total; ++var) {
                ycoord[var]=ycoord[var]*this->spectrum->getFlux().at(var)* distanceFactor;
            }
            if(loadPointsAction->isChecked()){
                for (int var = 0; var <lpoints.size(); ++var) {
                    lpoints[var].y=loadedPoints[var].x*lpoints[var].y *distanceFactor;
                }
                if(isErrorbar){
                    for (int var = 0; var <lpoints.size(); ++var) {
                        lpoints[var].y1=loadedPoints[var].x*lpoints[var].y1 *distanceFactor;
                        lpoints[var].y2=loadedPoints[var].x*lpoints[var].y2 *distanceFactor;
                    }
                }
            }
        }else if(this->spectrumLabel->getFluxUnit()==Spectrum::NULULSUN){
            distanceFactor=distanceFactor*2.5956986e-34;
            ycoord=this->spectrum->getFrequencies();
            for (int var = 0; var <total; ++var) {
                ycoord[var]=ycoord[var]*this->spectrum->getFlux().at(var)* distanceFactor;
            }
            if(loadPointsAction->isChecked()){
                for (int var = 0; var <lpoints.size(); ++var) {
                    lpoints[var].y=loadedPoints[var].x*lpoints[var].y*distanceFactor;
                }
                if(isErrorbar){
                    for (int var = 0; var <lpoints.size(); ++var) {
                        lpoints[var].y1=loadedPoints[var].x*lpoints[var].y1*distanceFactor;
                        lpoints[var].y2=loadedPoints[var].x*lpoints[var].y2*distanceFactor;
                    }
                }
            }

        }else if(this->spectrumLabel->getFluxUnit()==Spectrum::LNU){
            for (int var = 0; var <total; ++var) {
                ycoord.append(this->spectrum->getFlux().at(var)* distanceFactor);
            }
            if(loadPointsAction->isChecked()){
                for (int var = 0; var <lpoints.size(); ++var) {
                    lpoints[var].y=lpoints[var].y *distanceFactor;
                }
                if(isErrorbar){
                    for (int var = 0; var <lpoints.size(); ++var) {
                        lpoints[var].y1=lpoints[var].y1 *distanceFactor;
                        lpoints[var].y2=lpoints[var].y2 *distanceFactor;
                    }
                }
            }
        }

        if(xcoord.size()>0&& ycoord.size()>0){
            if(xcoord.at(0)>xcoord.at(xcoord.size()-1)){
                reverse(xcoord.begin(),xcoord.end());
                reverse(ycoord.begin(),ycoord.end());
            }
        }
        this->spectrum->setYcoord(ycoord);
        this->spectrum->setXcoord(xcoord);

        double xStart=0,xEnd=0;
        double yStart=0,yEnd=0;
        if(!this->spectrumLabel->isZoomActive()){
            xStart=*min_element(xcoord.begin(),xcoord.end());
            xEnd=*max_element(xcoord.begin(),xcoord.end());

            yEnd=*max_element(ycoord.begin(),ycoord.end());
        }


        pixmap.fill(Qt::white);

        painter.begin(&pixmap);
        painter.setPen(Qt::black);
        painter.setRenderHint(QPainter::Antialiasing);
        QVector<QPoint> points;
        int x,y;
        if(!this->spectrumLabel->isZoomActive()){
            if(this->spectrumLabel->isXLogOn()){
                xStart=log10(xStart);
                xEnd=log10(xEnd);
            }
            if(this->spectrumLabel->isYLogOn()){
                yEnd=log10(yEnd);
                yStart=yEnd-this->imageControllerWidget->getMaxLog();
            }else{
                yStart=*min_element(ycoord.begin(),ycoord.end());
            }
        }
        if(this->spectrumLabel->isZoomActive()){
            spectrumSetting *spectrumsetting=this->spectrumSettingVector.at(this->spectrumSettingVector.size()-1);
            xStart=spectrumsetting->getXStart();
            yStart=spectrumsetting->getYStart();
            xEnd=spectrumsetting->getXEnd();
            yEnd=spectrumsetting->getYEnd();
        }
        this->spectrum->setXStart(xStart);
        this->spectrum->setYStart(yStart);
        this->spectrum->setXEnd(xEnd);
        this->spectrum->setYEnd(yEnd);
        QVector<double>xVector,yVector;
        if(painter.isActive()){
            int total= this->spectrum->getNumberOfFrequencies();
            for (int var = 0; var <total; ++var) {
                if(this->spectrumLabel->isXLogOn())
                    x=((log10(xcoord[var])-xStart)/(xEnd-xStart))*(size-2*coordDistanceToImage);
                else
                    x=((xcoord[var]-xStart)/(xEnd-xStart))*(size-2*coordDistanceToImage);
                if(this->spectrumLabel->isYLogOn())
                    y=((log10(ycoord[var])-yStart)/(yEnd-yStart))*(size-2*coordDistanceToImage);
                else
                    y=((ycoord[var]-yStart)/(yEnd-yStart))*(size-2*coordDistanceToImage);
                xVector.append(x);
                yVector.append(y);
            }

            for (int var = 0; var <yVector.size(); ++var) {
                if(yVector.at(var)>=0){
                    points.append(QPoint(xVector.at(var)+coordDistanceToImage
                                         ,yVector.at(var)+coordDistanceToImage));
                }else{
                    points.append(QPoint(xVector.at(var)+coordDistanceToImage
                                         ,coordDistanceToImage));
                }
            }

            painter.setWindow(0,size,size,-size);
            painter.drawPolyline(points);
            QVector<QLine> xerrorLines,yerrorLines;
            if(loadPointsAction->isChecked()){
                painter.setPen(Qt::red);
                QPen pen = painter.pen() ;
                int width=pen.width();
                pen.setWidth(5);
                pen.setCapStyle(Qt::RoundCap);
                painter.setPen(pen);
                points.clear();
                double x1,x2,y1,y2;
                for (int var = 0; var < lpoints.size(); ++var) {
                    if(this->spectrumLabel->isXLogOn())
                        x=((log10(lpoints[var].x)-xStart)/(xEnd-xStart))*(size-2*coordDistanceToImage);
                    else
                        x=((lpoints[var].x-xStart)/(xEnd-xStart))*(size-2*coordDistanceToImage);
                    if(this->spectrumLabel->isYLogOn())
                        y=((log10(lpoints[var].y)-yStart)/(yEnd-yStart))*(size-2*coordDistanceToImage);
                    else
                        y=((lpoints[var].y-yStart)/(yEnd-yStart))*(size-2*coordDistanceToImage);
                    points.append(QPoint(x+coordDistanceToImage,y+coordDistanceToImage));
                    if(isErrorbar){
                        if(this->spectrumLabel->isXLogOn()){
                            x1=((log10(lpoints[var].x1)-xStart)/(xEnd-xStart))*(size-2*coordDistanceToImage);
                            x2=((log10(lpoints[var].x2)-xStart)/(xEnd-xStart))*(size-2*coordDistanceToImage);
                        }else{
                            x1=((lpoints[var].x1-xStart)/(xEnd-xStart))*(size-2*coordDistanceToImage);
                            x2=((lpoints[var].x2-xStart)/(xEnd-xStart))*(size-2*coordDistanceToImage);
                        }
                        QLine xline,yline;
                        if(x2-x1!=0){
                            xline.setP1(QPoint(x1+coordDistanceToImage,y+coordDistanceToImage));
                            xline.setP2(QPoint(x2+coordDistanceToImage,y+coordDistanceToImage));
                            xerrorLines.append(xline);
                        }
                        if(this->spectrumLabel->isYLogOn()){
                            y1=((log10(lpoints[var].y1)-yStart)/(yEnd-yStart))*(size-2*coordDistanceToImage);
                            y2=((log10(lpoints[var].y2)-yStart)/(yEnd-yStart))*(size-2*coordDistanceToImage);
                        }else{
                            y1=((lpoints[var].y1-yStart)/(yEnd-yStart))*(size-2*coordDistanceToImage);
                            y2=((lpoints[var].y2-yStart)/(yEnd-yStart))*(size-2*coordDistanceToImage);
                        }

                        if(y2-y1!=0){
                            yline.setP1(QPoint(x+coordDistanceToImage,y1+coordDistanceToImage));
                            yline.setP2(QPoint(x+coordDistanceToImage,y2+coordDistanceToImage));
                            yerrorLines.append(yline);
                        }
                    }
                }
                painter.drawPoints(points);
                pen.setCapStyle(Qt::SquareCap);
                pen.setWidth(width);;
                painter.setPen(pen);
                painter.setPen(Qt::black);
                int errorbarDelta=3;
                if(yerrorLines.size()>0){
                    painter.drawLines(yerrorLines);
                    foreach(QLine line,yerrorLines){
                        if(line.y1()>0){
                            painter.drawLine(line.x1()-errorbarDelta,line.y1(),line.x1()+errorbarDelta,line.y1());
                            painter.drawLine(line.x1()-errorbarDelta,line.y2(),line.x1()+errorbarDelta,line.y2());
                        }
                    }
                }
                if(xerrorLines.size()>0){
                    painter.drawLines(xerrorLines);
                    foreach(QLine line,xerrorLines){
                        if(line.y1()>0){
                            painter.drawLine(line.x1(),line.y1()-errorbarDelta,line.x1(),line.y1()+errorbarDelta);
                            painter.drawLine(line.x2(),line.y1()-errorbarDelta,line.x2(),line.y1()+errorbarDelta);
                        }
                    }
                }
            }
            painter.setWindow(0,0,size,size);
            painter.fillRect(0,0,size,this->coordDistanceToImage,Qt::white);
            painter.fillRect(0,0,this->coordDistanceToImage,size,Qt::white);
            painter.fillRect(size-this->coordDistanceToImage, 0,this->coordDistanceToImage,size,Qt::white);
            painter.fillRect(0,size-this->coordDistanceToImage,size,this->coordDistanceToImage,Qt::white);

            painter.drawRect(this->coordDistanceToImage,this->coordDistanceToImage,
                             size-this->coordDistanceToImage*2,size-this->coordDistanceToImage*2);
            int divider=8;
            int x1,y1,x2,y2;
            double xvalue=xStart;
            double yvalue=yEnd;
            QFont font = painter.font() ;
            font.setPointSizeF(font.pointSizeF() * 1.25);
            painter.setFont(font);
            QTransform transform;
            transform.rotate(+270);
            painter.setTransform(transform);
            int move=20;

            int var;
            for ( var= 0; var < fluxUnitsList.size(); ++var) {
                if(fluxGroup->actions().at(var)->isChecked())break;
            }
            if(var==0||var==2||var==6){
                move=100;
            }
            QString fluxText=fluxUnitsList.at(var);

            painter.drawText(-size/2-move,coordDistanceToImage-50,fluxText);
            transform.rotate(-270);
            painter.setTransform(transform);

            for ( var= 0; var < freqUnitsList.size(); ++var) {
                if(freqGroup->actions().at(var)->isChecked())break;
            }
            QString freqText=freqUnitsList.at(var);
            painter.drawText(size/2-50,size-coordDistanceToImage+40,freqText);
            font = painter.font() ;
            font.setPointSizeF(font.pointSizeF() / 1.25);
            painter.setFont(font);
            for (int index = 0; index <= divider; ++index) {
                x1=this->coordDistanceToImage+(size-this->coordDistanceToImage*2.)/8.*index;
                x2=x1;
                y1=this->coordDistanceToImage;
                y2=y1-5;
                painter.drawLine(x1,y1,x2,y2);
                painter.drawLine(y1,x1,y2,x2);
                if(index%2==0){
                    if(this->spectrumLabel->isYLogOn())
                        painter.drawText(y1-45,x2,trUtf8("%1").arg(QString::number(pow(10.,yvalue),'e',0)));
                    else
                        painter.drawText(y1-45,x2,trUtf8("%1").arg(QString::number((yvalue),'e',0)));
                    yvalue=yvalue-(yEnd-yStart)/(divider/2.);

                }
                y1=size-this->coordDistanceToImage;
                y2=y1+5;
                painter.drawLine(x1,y1,x2,y2);
                painter.drawLine(y1,x1,y2,x2);
                if(index%2==0){
                    if(this->spectrumLabel->isXLogOn())
                        painter.drawText(x1-10,y2+15,trUtf8("%1").arg(QString::number(pow(10.,xvalue),'e',0)));
                    else
                        painter.drawText(x1-10,y2+15,trUtf8("%1").arg(QString::number((xvalue),'e',0)));
                    xvalue=xvalue+(xEnd-xStart)/(divider/2.);
                }


            }

        }

        painter.end();
        this->spectrum->setPixmap(pixmap);
        this->spectrumLabel->setPixmap(pixmap);
    }
}

void MainWindow::reloadImage(){
    if(imageIsRead){
        int index=this->leftWidget->currentIndex();
        this->imageControllerWidget->setRenderSpectrumActive(false);
        if(index!=this->leftWidget->indexOf(this->imageGroup))this->leftWidget->setCurrentWidget(imageGroup);
        int size;
        int width=this->width();
        int height=this->height();
        if(width>height){
            if(width/2>height)
                size=height;
            else
                size=width/2;
        }else{
            if(height/2>width)
                size=width;
            else
                size=height/2;
        }

        QPainter painter;
        QPixmap pixmap(size,size);

        QPixmap imagePixmap(lastImage.width(),lastImage.height());


        painter.begin(&imagePixmap);
        if(painter.isActive()){
            painter.drawImage(0,0,this->lastImage);

            painter.end();
        }

        pixmap.fill(this->leftWidget->foregroundRole());

        painter.begin(&pixmap);
        painter.setPen(Qt::white);
        if(painter.isActive()){
            double maxValue[3];
            if(this->imageControllerWidget->isTau()||this->image->isLocalObserver()){
                QVector< QVector< QVector< double> > > intensity;
                intensity.append(image->getImageAt(0));
                if(this->imageControllerWidget->isRGBModus()){
                    intensity.append(image->getImageAt(1));
                    intensity.append(image->getImageAt(2));
                }
                for(int i=0;i<intensity.size();i++){
                    for (int index = 0; index < intensity.at(0).size(); ++index){
                        double  max=*max_element(intensity[i][index].begin(),intensity[i][index].end());
                        if(maxValue[i]<max||(i==0&&index==0))maxValue[i]=max;
                    }
                }
            }
            if(this->imageLabel->isPlot()){
                imagePixmap=imagePixmap.scaled(size-this->coordDistanceToImage*2,size-this->coordDistanceToImage*2,
                                               Qt::KeepAspectRatio,this->transformMode);
                painter.drawPixmap(this->coordDistanceToImage,this->coordDistanceToImage,
                                   imagePixmap);
            }else{
                imagePixmap=imagePixmap.scaled(size,size,
                                               Qt::KeepAspectRatio,this->transformMode);
                painter.drawPixmap(0,0,imagePixmap);
            }
            if(this->imageControllerWidget->isTau()||this->image->isLocalObserver())
                if(!this->imageControllerWidget->isRGBModus())
                    painter.drawText(pixmap.width()-170,pixmap.height()-5,trUtf8("Max  \u03C4au = %1 ")
                                     .arg(QString::number(maxValue[0],'e')));
                else{
                    painter.setPen(Qt::red);
                    painter.drawText(pixmap.width()-170,pixmap.height()-35,trUtf8("Max  \u03C4au = %1 ")
                                     .arg(QString::number(maxValue[0],'e')));

                    painter.setPen(Qt::green);
                    painter.drawText(pixmap.width()-170,pixmap.height()-20,trUtf8("Max  \u03C4au = %1 ")
                                     .arg(QString::number(maxValue[1],'e')));

                    painter.setPen(QColor(30,140,255));
                    painter.drawText(pixmap.width()-170,pixmap.height()-5,trUtf8("Max  \u03C4au = %1 ")
                                     .arg(QString::number(maxValue[2],'e')));
                    painter.setPen(Qt::white);
                }
            else
                if(!this->imageControllerWidget->isRGBModus())
                    painter.drawText(pixmap.width()-200,pixmap.height()-5,trUtf8("Flux = %1 Jy @ 1pc")
                                     .arg(QString::number(image->flux.at(0)*pow(10.,23),'e')));
                else{
                    painter.setPen(Qt::red);
                    painter.drawText(pixmap.width()-200,pixmap.height()-35,trUtf8("Flux = %1 Jy @ 1pc")
                                     .arg(QString::number(image->flux.at(0)*pow(10.,23),'e')));

                    painter.setPen(Qt::green);
                    painter.drawText(pixmap.width()-200,pixmap.height()-20,trUtf8("Flux = %1 Jy @ 1pc")
                                     .arg(QString::number(image->flux.at(1)*pow(10.,23),'e')));

                    painter.setPen(QColor(30,140,255));
                    painter.drawText(pixmap.width()-200,pixmap.height()-5,trUtf8("Flux = %1 Jy @ 1pc")
                                     .arg(QString::number(image->flux.at(2)*pow(10.,23),'e')));
                    painter.setPen(Qt::white);
                }

            QPen pen=painter.pen();
            pen.setWidth(2);
            painter.setPen(pen);
            if(this->imageLabel->isPlot()){
                painter.drawRect(this->coordDistanceToImage,this->coordDistanceToImage,
                                 size-this->coordDistanceToImage*2,size-this->coordDistanceToImage*2);
            }

            if(this->imageControllerWidget->isColorBar()){
                int colorBarSchrinkSize=5;
                QFont font = painter.font() ;
                font.setPointSize(font.pointSize() / 1.2);
                painter.setFont(font);
                painter.drawRect(size-this->coordDistanceToImage+10,this->coordDistanceToImage+colorBarSchrinkSize,
                                 this->coordDistanceToImage/8,size-this->coordDistanceToImage*2-colorBarSchrinkSize*2);
                if(!this->imageControllerWidget->isRGBModus()){
                    QVector<QRgb> tble=this->imageControllerWidget->getCurrentColorLookUpTable();
                    if(this->imageControllerWidget->isColorInvert())reverse(tble.begin(),tble.end());
                    int total=(size-this->coordDistanceToImage*2-colorBarSchrinkSize*2);
                    for (int var = 0; var <total; ++var) {
                        painter.fillRect(size-this->coordDistanceToImage+10,this->coordDistanceToImage+(total-1-var)+colorBarSchrinkSize, this->coordDistanceToImage/8,
                                         1,QColor::fromRgb(tble.at(tble.size()*var/total)));
                    }
                }else{
                    double max,min;
                    max=image->getMax(0);
                    min=image->getMin(0);
                    int total=(size-this->coordDistanceToImage*2-colorBarSchrinkSize*2);
                    int maxColorNumber=0;
                    int mult=-1;
                    if(this->imageControllerWidget->isColorInvert()){
                        mult=1;
                        maxColorNumber=255;
                    }
                    for (int var = 0; var <total; ++var) {
                        painter.fillRect(size-this->coordDistanceToImage+10,this->coordDistanceToImage+(total-1-var)+colorBarSchrinkSize, this->coordDistanceToImage/8,
                                         1,QColor::fromRgb(qRgb(maxColorNumber-((max-min)/max*var*255./total)*mult,0,0)));
                    }

                    painter.drawRect(size-this->coordDistanceToImage+20+this->coordDistanceToImage/8,this->coordDistanceToImage+colorBarSchrinkSize,
                                     this->coordDistanceToImage/8,size-this->coordDistanceToImage*2-colorBarSchrinkSize*2);
                    max=image->getMax(1);
                    min=image->getMin(1);
                    for (int var = 0; var <total; ++var) {
                        painter.fillRect(size-this->coordDistanceToImage+20+this->coordDistanceToImage/8,this->coordDistanceToImage+(total-1-var)+colorBarSchrinkSize, this->coordDistanceToImage/8,
                                         1,QColor::fromRgb(qRgb(0,maxColorNumber-((max-min)/max*var*255./total)*mult,0)));
                    }
                    painter.drawText(size-this->coordDistanceToImage+10+this->coordDistanceToImage/8,this->coordDistanceToImage-12,
                                     QString::number(max,'e',0));
                    painter.drawText(size-this->coordDistanceToImage+15+this->coordDistanceToImage/8,size-this->coordDistanceToImage+15,
                                     QString::number(min,'e',0));


                    painter.drawRect(size-this->coordDistanceToImage+30+this->coordDistanceToImage/4,this->coordDistanceToImage+colorBarSchrinkSize,
                                     this->coordDistanceToImage/8,size-this->coordDistanceToImage*2-colorBarSchrinkSize*2);
                    max=image->getMax(2);
                    min=image->getMin(2);
                    for (int var = 0; var <total; ++var) {
                        painter.fillRect(size-this->coordDistanceToImage+30+this->coordDistanceToImage/4,this->coordDistanceToImage+(total-1-var)+colorBarSchrinkSize, this->coordDistanceToImage/8,
                                         1,QColor::fromRgb(qRgb(0,0,maxColorNumber-((max-min)/max*var*255./total)*mult)));
                    }
                    painter.drawText(size-this->coordDistanceToImage+30+this->coordDistanceToImage/8,this->coordDistanceToImage,
                                     QString::number(max,'e',0));
                    painter.drawText(size-this->coordDistanceToImage+35+this->coordDistanceToImage/8,size-this->coordDistanceToImage+5,
                                     QString::number(min,'e',0));


                }
                if(this->imageControllerWidget->isLinear()){
                painter.drawText(size-this->coordDistanceToImage+5,this->coordDistanceToImage,
                                 QString::number(image->getMax(0),'e',0));
                painter.drawText(size-this->coordDistanceToImage+5,size-this->coordDistanceToImage+5,
                                 QString::number(image->getMin(0),'e',0));
                }else{
                    painter.drawText(size-this->coordDistanceToImage+5,this->coordDistanceToImage,
                                     QString::number(log10(image->getMax(0)),'g',3));
                    painter.drawText(size-this->coordDistanceToImage+5,size-this->coordDistanceToImage+5,
                                     QString::number(log10(image->getMin(0)),'g',3));
                }
                font = painter.font() ;
                font.setPointSize(font.pointSize() * 1.2);
                painter.setFont(font);


            }
            int divider=8;
            int x1,y1,x2,y2;
            double xend,yend,xstart,ystart;
            if(!this->imageControllerWidget->isLocalModus()){
                if(!this->imageControllerWidget->isZoomActive()){
                    xend=image->getNumberOfPixelX()*image->getSizeOfPixelX()/2.;
                    yend=image->getNumberOfPixelY()*image->getSizeOfPixelY()/2;
                    xstart=xend*-1;
                    ystart=yend*-1;
                }else{
                    QVector<double> zoombox=this->imageControllerWidget->getZoomBox();
                    xstart=zoombox.at(0);
                    xend=zoombox.at(1);
                    ystart=zoombox.at(2);
                    yend=zoombox.at(3);
                }
                if(this->setup->getImageUnit()==this->setup->PC){
                    xend=xend/this->physicalConstants->getParsec();
                    yend=yend/this->physicalConstants->getParsec();
                    ystart=ystart/this->physicalConstants->getParsec();
                    xstart=xstart/this->physicalConstants->getParsec();
                }else if(this->setup->getImageUnit()==this->setup->AU){
                    xend=xend/this->physicalConstants->getAstronomicalUnit();
                    yend=yend/this->physicalConstants->getAstronomicalUnit();
                    ystart=ystart/this->physicalConstants->getAstronomicalUnit();
                    xstart=xstart/this->physicalConstants->getAstronomicalUnit();
                }
            }else{
                xend=this->imageControllerWidget->getViewangleSpinBoxValue()/2;
                yend=this->imageControllerWidget->getViewangleSpinBoxValue()/2;
                xstart=xend*-1;
                ystart=yend*-1;
            }
            if(this->imageLabel->isPlot()){
                double xvalue=xstart;
                double yvalue=yend;
                QTransform transform;
                transform.rotate(90);
                painter.setTransform(transform);
                painter.drawText(size/2,40-coordDistanceToImage,this->imageControllerWidget->getSizeUnit());
                transform.rotate(-90);
                painter.setTransform(transform);
                painter.drawText(size/2-20,size-coordDistanceToImage+40,this->imageControllerWidget->getSizeUnit());
                for (int index = 0; index <= divider; ++index) {
                    x1=this->coordDistanceToImage+(size-this->coordDistanceToImage*2.)/8.*index;
                    x2=x1;
                    y1=this->coordDistanceToImage;
                    y2=y1-5;
                    painter.drawLine(x1,y1,x2,y2);
                    painter.drawLine(y1,x1,y2,x2);
                    if(index%2==0){
                        painter.drawText(y1-60,x2,trUtf8("%1").arg(QString::number(yvalue,'e',0)));
                        yvalue=yvalue-(yend-ystart)/(divider/2.);
                    }
                    y1=size-this->coordDistanceToImage;
                    y2=y1+5;
                    painter.drawLine(x1,y1,x2,y2);
                    painter.drawLine(y1,x1,y2,x2);
                    if(index%2==0){
                        painter.drawText(x1-20,y2+15,trUtf8("%1").arg(QString::number(xvalue,'e',0)));
                        xvalue=xvalue+(xend-xstart)/(divider/2.);
                    }


                }
            }

            if(this->imageControllerWidget->isContour()){
                QPen pen=painter.pen();
                pen.setWidthF(0.8);
                pen.setColor(palette().color(QPalette::Highlight));
                painter.setRenderHint(QPainter::Antialiasing);
                painter.setPen(pen);
                for (int var = 0; var < contourPoints.size(); ++var) {
                    QVector < QPoint > poly;
                    poly=contourPoints[var];
                    QPolygon polygen(poly);
                    QTransform transpoly;
                    transpoly=transpoly.scale(imagePixmap.width()*1./scaleImage.width(),imagePixmap.height()*1./scaleImage.height());
                    poly=transpoly.map(poly);
                    int dis=coordDistanceToImage;
                    if(!this->imageLabel->isPlot())dis=0;
                    for (int index = 0; index < poly.size(); ++index) {
                        poly[index].setX(poly[index].x()+dis);
                        poly[index].setY(poly[index].y()+dis);
                    }
                    painter.drawPolygon(poly,Qt::WindingFill);
                }
            }
            painter.end();
        }
        this->image->setPixmap(pixmap);
        this->imageLabel->setPixmap(pixmap);
        hideLog();
    }
}

QVector<QPoint> MainWindow::PointInterpolation(const QVector<QPoint> &points){
    QVector<QPoint> points2;
    for(int i=0; i<points.size()-1; i++){
        QPoint p1(points.at(i)), p2(points.at(i+1));

        int dx = p2.x() - p1.x();
        int dy = p2.y() - p1.y();

        points2 << p1;

        if(abs(dx) > 1 || abs(dy) > 1){
            // linear interpolation
            if(dx == 0){
                // vertical
                int increment = dy > 0 ? 1: -1;
                for(int y=p1.y()+increment; y!= p2.y(); y+= increment){
                    points2 << QPoint(p1.x(),y);
                }
            }else{
                // calculate the gradient
                double k = double(dy)/double(dx);

                if(abs(dx) >= abs(dy)){
                    int increment = dx > 0 ? 1 : -1;
                    for(int x=p1.x() + increment; x!=p2.x(); x+= increment){
                        double y = p1.y() + (x-p1.x()) * k;
                        points2 << QPoint(x, round(y));
                    }
                }else{
                    int increment = dy > 0 ? 1 : -1;
                    for(int y=p1.y()+increment; y!=p2.y(); y+= increment){
                        // (y - p1.y)/(x - p1.x) = k
                        double x = (y-p1.y())/k + p1.x();
                        points2 << QPoint(round(x), y);
                    }
                }
            }
        }
    }

    return points2;
}

void   MainWindow:: resizeEvent(QResizeEvent* event){
    if(!this->imageControllerWidget->isRenderSpectrumActive()){
        if(imageIsRead ){
            reloadImage();
        }
    }else{
        if(spectrumIsRead){
            reloadSpectrum();
        }
    }
    QMainWindow::resizeEvent(event);
}
void MainWindow::calculateZoom(QPoint first,QPoint last){
    double yEnd,xEnd,xStart,yStart;
    QVector<double> zoombox;
    if(!this->imageControllerWidget->isRenderSpectrumActive()){
        if(!this->imageControllerWidget->isZoomActive()){
            yEnd=this->image->getSizeOfPixelY()*this->image->getNumberOfPixelY()/2;
            xEnd=this->image->getSizeOfPixelX()*this->image->getNumberOfPixelX()/2;
            xStart=-xEnd;
            yStart=-yEnd;
        }else{
            zoombox=this->imageControllerWidget->getZoomBox();
            xStart=zoombox.at(0);
            xEnd=zoombox.at(1);
            yStart=zoombox.at(2);
            yEnd=zoombox.at(3);
        }
        int distance=this->coordDistanceToImage;
        if(!this->imageLabel->isPlot())distance=0;
        int ix=first.x()-distance;
        int iy=first.y()-distance;
        int lx=last.x()-distance;
        int ly=last.y()-distance;
        int width,height;
        width=this->image->getPixmapWidth()-2*distance;
        height=this->image->getPixmapHeight()-2*distance;
        zoombox.clear();
        double buffer;
        zoombox.append(xStart+((xEnd-xStart)*(ix*1./width)));
        zoombox.append(xStart+((xEnd-xStart)*(lx*1./width)));
        if(zoombox.at(0)>zoombox.at(1)){
            buffer=zoombox.at(0);
            zoombox[0]=zoombox.at(1);
            zoombox[1]=buffer;
        }
        zoombox.append(yStart+((yEnd-yStart)*(1.-(ly*1./height))));
        zoombox.append(yStart+((yEnd-yStart)*(1.-(iy*1./height))));
        if(zoombox.at(2)>zoombox.at(3)){
            buffer=zoombox.at(2);
            zoombox[2]=zoombox.at(3);
            zoombox[3]=buffer;
        }
        this->zoomBox=zoombox;
        setZoomBox=true;
    }else{
        spectrumSetting* spectrumsetting=new spectrumSetting();
        if(!this->spectrumLabel->isZoomActive()){
            xStart=this->spectrum->getXStart();
            yStart=this->spectrum->getYStart();
            xEnd=this->spectrum->getXEnd();
            yEnd=this->spectrum->getYEnd();
            this->spectrumLabel->setZoomActive(true);
            spectrumSettingVector.clear();
            emit setZoomOutVisible(true);
        }else{
            xStart=this->spectrum->getXStart();
            yStart=this->spectrum->getYStart();
            xEnd=this->spectrum->getXEnd();
            yEnd=this->spectrum->getYEnd();
        }
        int distance=this->coordDistanceToImage;
        int ix=first.x()-distance;
        int iy=first.y()-distance;
        int lx=last.x()-distance;
        int ly=last.y()-distance;
        int width,height;
        width=this->spectrum->getPixmapWidth()-2*distance;
        height=this->spectrum->getPixmapHeight()-2*distance;
        zoombox.clear();
        double buffer;
        zoombox.append(xStart+((xEnd-xStart)*(ix*1./width)));
        zoombox.append(xStart+((xEnd-xStart)*(lx*1./width)));
        if(zoombox.at(0)>zoombox.at(1)){
            buffer=zoombox.at(0);
            zoombox[0]=zoombox.at(1);
            zoombox[1]=buffer;
        }
        zoombox.append(yStart+((yEnd-yStart)*(1.-(ly*1./height))));
        zoombox.append(yStart+((yEnd-yStart)*(1.-(iy*1./height))));
        if(zoombox.at(2)>zoombox.at(3)){
            buffer=zoombox.at(2);
            zoombox[2]=zoombox.at(3);
            zoombox[3]=buffer;
        }
        spectrumsetting->setXStart(zoombox[0]);
        spectrumsetting->setYStart(zoombox[2]);
        spectrumsetting->setXEnd(zoombox[1]);
        spectrumsetting->setYEnd(zoombox[3]);
        spectrumSettingVector.append(spectrumsetting);
        reloadSpectrum();

    }
}




void MainWindow::setUpWhatsThis(){
    this->sizemenu->setWhatsThis(trUtf8("<html><font color=blue><img src=\":/images/size.png\" width=\"%1\" height=\"%2\"/>"
                                        "Scale can given as parallax second or astronomical unit or cm </html>").arg(
                                     QString::number(this->imageControllerWidget->whatsThisIconWidth), QString::number(this->imageControllerWidget->whatsThisIconHeight)));
    this->colorAction->setWhatsThis(trUtf8("<html><font color=blue>Switch color Modus between<br/><img src=\":/images/%3.png\" width=\"%1\" height=\"%2\"/>"
                                           "RGB and"
                                           "<font color=blue><img src=\":/images/no%3.png\" width=\"%1\" height=\"%2\"/>Indexed8</html>").arg(
                                        QString::number(this->imageControllerWidget->whatsThisIconWidth),
                                        QString::number(this->imageControllerWidget->whatsThisIconHeight),this->colorAction->objectName()));
    this->lineOptionAction->setWhatsThis(trUtf8("<html><font color=blue><img src=\":/images/%3.png\" width=\"%1\" height=\"%2\"/>"
                                                "Include Line options( This includes entry fields such as imol,"
                                                "iline and vkms to make it easier to specify the precise "
                                                "wavelength of the image in case of lines)  </br>"
                                                "<img src=\":/images/no%3.png\" width=\"%1\" height=\"%2\"/>exclude Line options</html>").arg(
                                             QString::number(this->imageControllerWidget->whatsThisIconWidth),
                                             QString::number(this->imageControllerWidget->whatsThisIconHeight),this->lineOptionAction->objectName()));
    this->userTransferAction->setWhatsThis(trUtf8("<html><font color=blue><img src=\":/images/%3.png\" width=\"%1\" height=\"%2\"/>"
                                                  "Load/unload %4 and %5 (if available) in the current directory</html>").arg(
                                               QString::number(this->imageControllerWidget->whatsThisIconWidth),
                                               QString::number(this->imageControllerWidget->whatsThisIconHeight)
                                               ,this->userTransferAction->objectName(),
                                               this->setup->getTransferDefaultFilename(),this->setup->getTransferFilename()));
    this->changeCurrentDirectoryAction->setWhatsThis(trUtf8("<html><font color=blue><img src=\":/images/%3.png\" width=\"%1\" height=\"%2\"/>"
                                                            "Change current working directory</html>").arg(
                                                         QString::number(this->imageControllerWidget->whatsThisIconWidth),
                                                         QString::number(this->imageControllerWidget->whatsThisIconHeight),
                                                         this->changeCurrentDirectoryAction->objectName()));
    this->localAction->setWhatsThis(trUtf8("<html><font color=blue><img src=\":/images/%3.png\" width=\"%1\" height=\"%2\"/>"
                                           "Local observer is on. <img src=\":/images/no%3.png\" width=\"%1\" height=\"%2\"/>"
                                           "Local observer is off.(see 9.11 section)</html>").arg(
                                        QString::number(this->imageControllerWidget->whatsThisIconWidth),
                                        QString::number(this->imageControllerWidget->whatsThisIconHeight),
                                        this->localAction->objectName()));
    this->savedefaultsAction->setWhatsThis(trUtf8("<html><font color=blue><img src=\":/images/%3.png\" width=\"%1\" height=\"%2\"/>"
                                                  "Save the Gui state as application's defaults.</html>").arg(
                                               QString::number(this->imageControllerWidget->whatsThisIconWidth),
                                               QString::number(this->imageControllerWidget->whatsThisIconHeight),
                                               this->savedefaultsAction->objectName()));
    this->readImageAction->setWhatsThis(trUtf8("<html><font color=blue><img src=\":/images/%3.png\" width=\"%1\" height=\"%2\"/>"
                                               "Load image file(formtted or g77 unformatted).<br/><br/>  Attention:<font color=red>if you choose wrong format, GUI will be crashed.</html>").arg(
                                            QString::number(this->imageControllerWidget->whatsThisIconWidth),
                                            QString::number(this->imageControllerWidget->whatsThisIconHeight)
                                            ,this->readImageAction->objectName()));
    this->readSpectrumAction->setWhatsThis(trUtf8("<html><font color=blue><img src=\":/images/%3.png\" width=\"%1\" height=\"%2\"/>"
                                                  "Load spectrum file(formatted).</html>").arg(
                                               QString::number(this->imageControllerWidget->whatsThisIconWidth),
                                               QString::number(this->imageControllerWidget->whatsThisIconHeight)
                                               ,this->readImageAction->objectName()));

    QList<QPushButton *> list=this->progressDialog->findChildren<QPushButton *>();
    list.at(0)->setWhatsThis(trUtf8("<html><font color=blue><img src=\":/images/kill.png\" width=\"%1\" height=\"%2\"/>"
                                    "This Button kill RADMC3D Process. It means, if you are in child(pipe)modus, RADMC3D will be killed "
                                    "and by your next command restarted again.</html>")
                             .arg(QString::number(this->imageControllerWidget->whatsThisIconWidth),
                                  QString::number(this->imageControllerWidget->whatsThisIconHeight)));

}

/**
 * @brief This method loads the last saved state of GUI, If the state of GUI was saved before.
 * The user has the possibilty to delete the last saved state of GUI.
 *
 */
void MainWindow::checkLaststate(){
    SetUp *setup=SetUp::getInstance();
    QFile file(setup->getLocalSettingsFilename());
    if(file.exists()){
        msgBox=new QMessageBox(this);
        connect(msgBox,SIGNAL(finished(int)),this,SLOT(closeMsg(int)));
        msgBox->setText(trUtf8("Do you want to load the last saved state of GUI?"));
        msgBox->setWindowTitle(trUtf8("Load last state"));
        msgBox->addButton(trUtf8("Yes"), QMessageBox::ActionRole);
        msgBox->addButton(trUtf8("No"), QMessageBox::ActionRole);
        msgBox->addButton(trUtf8("No, delete saved state"), QMessageBox::ActionRole);
        msgBox->addButton(trUtf8("No,reset"), QMessageBox::ActionRole);
        msgBox->setDefaultButton(QMessageBox::No);
        msgBox->setIcon(QMessageBox::Question);
        msgBox->exec();
    }else{
        if(!firstTime){
            msgBox=new QMessageBox(this);
            connect(msgBox,SIGNAL(finished(int)),this,SLOT(closeMessage(int)));
            msgBox->setText(trUtf8("Do you want to reset state of GUI?"));
            msgBox->addButton(trUtf8("Yes,reset"), QMessageBox::ActionRole);
            msgBox->addButton(trUtf8("No"), QMessageBox::ActionRole);
            msgBox->setIcon(QMessageBox::Question);
            msgBox->exec();
        }else{
            this->imageControllerWidget->loadAsDefaults();
            loadAsdefaults();
        }
    }
    if(firstTime)firstTime=false;
}

void MainWindow::closeMessage(int result){
    if(result==0){
        this->colorAction->blockSignals(true);
        this->colorAction->setChecked(false);
        this->colorAction->blockSignals(false);
        this->localAction->blockSignals(true);
        this->localAction->setIcon(QIcon(trUtf8(":/images/no%1.png").arg(this->localAction->objectName())));
        this->localAction->setChecked(false);
        this->localAction->blockSignals(false);
        this->lineOptionAction->blockSignals(true);
        this->lineOptionAction->setChecked(false);
        this->lineOptionAction->blockSignals(false);
        this->userTransferAction->blockSignals(true);
        this->userTransferAction->setChecked(false);
        this->userTransferAction->blockSignals(false);
        this->imageControllerWidget->setZoomActivity(false);
        this->setZoomBox=false;
        this->sizemenu->actions().at(0)->blockSignals(true);
        this->sizemenu->actions().at(0)->setCheckable(true);
        this->sizemenu->actions().at(0)->blockSignals(false);
        loadAsdefaults();
        this->imageControllerWidget->reset();
    }else if(result==1){
        this->imageControllerWidget->setZoomActivity(false);
        this->setZoomBox=false;
        if(this->userTransferAction->isChecked()){
            if(!QFile(this->setup->getTransferDefaultFilename()).exists()){
                this->userTransferAction->blockSignals(true);
                this->userTransferAction->setChecked(false);
                this->userTransferAction->setIcon(QIcon(trUtf8(":/images/no%1.png").arg(this->userTransferAction->objectName())));
                emit removeTransferTab();
                this->userTransferAction->blockSignals(false);
            }
        }
        if(this->imageControllerWidget->isLineModus()){
            if(!QFile(this->setup->getLinesFilname()).exists()){
                this->lineOptionAction->blockSignals(true);
                this->lineOptionAction->setChecked(false);
                this->lineOptionAction->setIcon(QIcon(trUtf8(":/images/no%1.png").arg(this->lineOptionAction->objectName())));
                this->imageControllerWidget->setLineModus(false);
                this->lineOptionAction->blockSignals(false);
            }
        }
    }
}

void MainWindow::closeMsg(int result){
    if(result==0){
        this->imageControllerWidget->setZoomActivity(false);
        this->setZoomBox=false;
        this->imageControllerWidget->reset();
        this->localAction->blockSignals(true);
        this->localAction->setIcon(QIcon(trUtf8(":/images/no%1.png").arg(this->localAction->objectName())));
        this->localAction->setChecked(false);
        this->localAction->blockSignals(false);
        this->load();
        this->imageControllerWidget->load();
        loaded=true;
    }else if(result==2){
        QFile file(setup->getLocalSettingsFilename());
        file.remove();
    }else if(result==3){
        this->localAction->blockSignals(true);
        this->localAction->setIcon(QIcon(trUtf8(":/images/no%1.png").arg(this->localAction->objectName())));
        this->localAction->setChecked(false);
        this->localAction->blockSignals(false);
        this->colorAction->blockSignals(true);
        this->colorAction->setChecked(false);
        this->colorAction->blockSignals(false);
        this->lineOptionAction->blockSignals(true);
        this->lineOptionAction->setChecked(false);
        this->lineOptionAction->blockSignals(false);
        this->userTransferAction->blockSignals(true);
        this->userTransferAction->setChecked(false);
        this->userTransferAction->blockSignals(false);
        this->imageControllerWidget->setZoomActivity(false);
        this->setZoomBox=false;
        this->sizemenu->actions().at(0)->blockSignals(true);
        this->sizemenu->actions().at(0)->setCheckable(true);
        this->sizemenu->actions().at(0)->blockSignals(false);
        loadAsdefaults();
        this->imageControllerWidget->reset();
    }else if(result==1){
        this->imageControllerWidget->setZoomActivity(false);
        this->setZoomBox=false;
        if(this->userTransferAction->isChecked()){
            if(!QFile(this->setup->getTransferDefaultFilename()).exists()){
                this->userTransferAction->blockSignals(true);
                this->userTransferAction->setChecked(false);
                this->userTransferAction->setIcon(QIcon(trUtf8(":/images/no%1.png").arg(this->userTransferAction->objectName())));
                emit removeTransferTab();
                this->userTransferAction->blockSignals(false);
            }
        }
        if(this->imageControllerWidget->isLineModus()){
            if(!QFile(this->setup->getLinesFilname()).exists()){
                this->lineOptionAction->setChecked(false);
                this->lineOptionAction->setIcon(QIcon(trUtf8(":/images/no%1.png").arg(this->lineOptionAction->objectName())));
                this->imageControllerWidget->setLineModus(false);
            }
        }
    }
    this->msgBox->close();
}
void MainWindow::updateImageMenu(){
    QList<QCheckBox *> list=this->imageControllerWidget->findChildren<QCheckBox *>();
    int size=list.size();
    QCheckBox *checkBox;
    QString no="";
    for(int i = 0; i <size ; ++i) {
        checkBox=list.at(i);
        QString objectName=checkBox->objectName()+"Action";
        QWidgetAction *action=this->findChild<QWidgetAction *>(objectName);

        if(!checkBox->isChecked()){
            no="no";
        }else{
            no="";
        }
        action->setIcon(QIcon(QString(":/images/%1%2.png").arg(no,checkBox->objectName())));
    }
}

void MainWindow::setImageMenu(){
    QList<QCheckBox *> list=this->imageControllerWidget->findChildren<QCheckBox *>();
    this->imageMenu=new QMenu();
    int size=list.size();
    QWidgetAction *action;
    QString no="";
    QCheckBox *checkBox;
    for (int i = 0; i <size ; ++i) {
        checkBox=list.at(i);
        action=new QWidgetAction(this);
        action->setToolTip(checkBox->shortcut().toString().trimmed());
        action->setShortcut(checkBox->shortcut());
        action->setText(checkBox->text());
        action->setWhatsThis(checkBox->whatsThis());
        action->setVisible(true);
        action->setShortcutContext(Qt::WidgetShortcut);
        if(!checkBox->isChecked()){
            no="no";
        }else{
            no="";
        }
        action->setObjectName(checkBox->objectName()+"Action");
        action->setIcon(QIcon(QString(":/images/%1%2.png").arg(no,checkBox->objectName())));
        action->setIconVisibleInMenu(true);
        this->imageMenu->addAction(action);
        this->connect(action,SIGNAL(triggered()),checkBox,SLOT(animateClick()));


    }

    this->menuBar()->setNativeMenuBar(this->nativeMenuBar);
}


void MainWindow::changeIcon(QString checkBoxobjectName){
    QString  objectName=checkBoxobjectName;
    QWidgetAction *action;
    objectName.append("Action");
    action=this->findChild<QWidgetAction*>(objectName);
    QCheckBox *checkBox=this->imageControllerWidget->findChild<QCheckBox*>(checkBoxobjectName);
    QString no="";
    if((checkBox!=NULL) && (action!=NULL)){
        if(!checkBox->isChecked()){
            this->imageMenu->setStyleSheet(trUtf8("QMenu::item#%1{color :%2;}").arg(checkBox->objectName(),this->imageControllerWidget->CheckBoxInactiveFontColor));
            no="no";
        }else{
            this->imageMenu->setStyleSheet(trUtf8("QMenu::item#%1{color :%2;}").arg(objectName,this->imageControllerWidget->labelColorMap->value(checkBoxobjectName)));
            no="";
        }
        action->setIcon(QIcon(QString(":/images/%1%2.png").arg(no,checkBoxobjectName)));
    }
}
void MainWindow::setWindowSize(bool checked){
    if(checked){
        if(this->mainLayout->direction()==QBoxLayout::LeftToRight||this->mainLayout->direction()==QBoxLayout::RightToLeft)
            this->setFixedSize(800,600);

        else
            this->setFixedSize(600,800);
    }else{
        this->layout()->setSizeConstraint(QLayout::SetMinAndMaxSize);
    }
}

/**
 *  @brief This method sets the widget direction
 *
 *  @details This method set mainLayout's direction adquat to
 *  parameter direction
 *
 */
void MainWindow::setMainLayoutDirection(QBoxLayout::Direction direction){
    if(this->mainLayout->direction()==QBoxLayout::LeftToRight||this->mainLayout->direction()==QBoxLayout::RightToLeft){
        if(direction==QBoxLayout::TopToBottom||direction==QBoxLayout::BottomToTop){
            this->setFixedSize(600,800);
        }
    }else{
        if(direction==QBoxLayout::LeftToRight||direction==QBoxLayout::RightToLeft){
            this->setFixedSize(800,600);
        }

    }
    this->mainLayout->setDirection(direction);
    this->mainLayout->removeWidget(this->leftWidget);
    this->mainLayout->removeWidget(this->rightWidget);
    this->mainLayout->addWidget(this->leftWidget);
    this->mainLayout->addWidget(this->rightWidget);
    this->layout()->setSizeConstraint(QLayout::SetMinAndMaxSize);

}
void   MainWindow:: userTransferChanged(bool checked){
    if(checked){
        if(!QFile(this->setup->getTransferDefaultFilename()).exists()){
            QMessageBox::information(this,trUtf8("Transfer file is not found"),trUtf8("there is no %1 in current directory").arg(this->setup->getTransferDefaultFilename()));
            this->userTransferAction->setChecked(false);
        }else{
            this->userTransferAction->setIcon(QIcon(trUtf8(":/images/%1.png").arg(this->userTransferAction->objectName())));
            emit readTransferFile();
        }
    }else{
        this->userTransferAction->setIcon(QIcon(trUtf8(":/images/no%1.png").arg(this->userTransferAction->objectName())));
        emit removeTransferTab();
    }
}
void MainWindow::localObserverChanged(bool checked,bool reload){
    if(checked){
        this->imageControllerWidget->activateLocalObserver(true);
        this->localAction->setIcon(QIcon(trUtf8(":/images/%1.png").arg(this->localAction->objectName())));
        if(reload)sendCommand(this->imageControllerWidget->getRenderImagePushButtonObjectName());
    }else{
        this->imageControllerWidget->activateLocalObserver(false);
        this->localAction->setIcon(QIcon(trUtf8(":/images/no%1.png").arg(this->localAction->objectName())));
        if(reload)sendCommand(this->imageControllerWidget->getRenderImagePushButtonObjectName());
    }
}

void MainWindow::colorChanged(bool checked,bool reload){
    if(checked){
        if(this->imageControllerWidget->isLineModus()){
            QMessageBox::information(this,trUtf8("incompatibility"),trUtf8("The line modus is incompatible with color modus. Line modus is now deactivated"));
            this->lineOptionAction->blockSignals(true);
            this->lineOptionAction->setChecked(false);
            this->lineOptionAction->blockSignals(false);
            this->lineOptionAction->setIcon(QIcon(trUtf8(":/images/no%1.png").arg(this->lineOptionAction->objectName())));
            this->imageControllerWidget->setLineModus(false);
        }
        this->colorAction->setIcon(QIcon(trUtf8(":/images/%1.png").arg(this->colorAction->objectName())));
        QAction *action=this->findChild<QWidgetAction*>(imageControllerWidget->getAbsoluteScaleCheckBoxObjectName()+"Action");
        this->colorAction->blockSignals(true);
        this->colorAction->setChecked(true);
        this->colorAction->blockSignals(false);
        if(action!=NULL)action->setVisible(true);
        this->imageControllerWidget->addColorModusWidgets();
        if(reload)sendCommand(this->imageControllerWidget->getRenderImagePushButtonObjectName());
    }else{
        this->colorAction->setIcon(QIcon(trUtf8(":/images/no%1.png").arg(this->colorAction->objectName())));
        QAction *action=this->findChild<QWidgetAction*>(imageControllerWidget->getAbsoluteScaleCheckBoxObjectName()+"Action");
        this->colorAction->blockSignals(true);
        this->colorAction->setChecked(false);
        this->colorAction->blockSignals(false);
        if(action!=NULL)action->setVisible(false);
        this->imageControllerWidget->removeColorModusWidgets();
        if(reload)sendCommand(this->imageControllerWidget->getRenderImagePushButtonObjectName());
    }
}
void MainWindow::lineOptionChanged(bool checked,bool reload){
    if(checked){
        if(!this->imageControllerWidget->isRGBModus() && QFile(this->setup->getLinesFilname()).exists()){
            this->lineOptionAction->setIcon(QIcon(trUtf8(":/images/%1.png").arg(this->lineOptionAction->objectName())));
            this->imageControllerWidget->setLineModus(true);
            if(reload)sendCommand(this->imageControllerWidget->getRenderImagePushButtonObjectName());
        }else{
            if(this->imageControllerWidget->isRGBModus())
                if( !QFile(this->setup->getLinesFilname()).exists())
                    QMessageBox::information(this,trUtf8("Line modus can not be activated"),trUtf8("You can not activate Line modus with color modus and there is no %1").arg(this->setup->getLinesFilname()));
                else
                    QMessageBox::information(this,trUtf8("Line modus can not be activated"),trUtf8("You can not activate Line modus with color modus"));
            else
                QMessageBox::information(this,trUtf8("Line modus can not be activated"),trUtf8("Line modus can not be activated, because there is no %1").arg(this->setup->getLinesFilname()));
            this->lineOptionAction->setChecked(false);
        }
    }else{
        this->lineOptionAction->setIcon(QIcon(trUtf8(":/images/no%1.png").arg(this->lineOptionAction->objectName())));
        this->imageControllerWidget->setLineModus(false);
        if(reload)sendCommand(this->imageControllerWidget->getRenderImagePushButtonObjectName());
    }
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
 *
 */
void MainWindow::setUp(){

    this->translate();
    this->setSlots();
    this->setShortcuts();
}
/**
 * @brief This method overrides the closeEvent. This method shows a messagebox to give a possibility
 * to save the state of GUI / abort closing / close the program without saving state of GUI.
 *
 */
void MainWindow::closeEvent(QCloseEvent *event){
    QMessageBox msgBox;
    msgBox.setText(trUtf8("<html>Are you sure you want to exit?<br/><br/><font color=blue>Current directory : %1</html>").arg(QDir::currentPath()));
    msgBox.setWindowTitle(trUtf8("Are you sure?"));
    QCheckBox *startNextTimeCheckbox=new QCheckBox(&msgBox);
    startNextTimeCheckbox->setChecked(startnexttime);
    startNextTimeCheckbox->setText(trUtf8("Start next time from current directory"));
    QAbstractButton *yesButton = msgBox.addButton(trUtf8("Yes"), QMessageBox::ActionRole);
    QAbstractButton *saveButton = msgBox.addButton(trUtf8("Yes, save GUI"), QMessageBox::ActionRole);
    QAbstractButton *noButton  =msgBox.addButton(trUtf8("No"), QMessageBox::ActionRole);
    startNextTimeCheckbox->blockSignals(true);
    msgBox.addButton(startNextTimeCheckbox, QMessageBox::NoRole);
    msgBox.setDefaultButton(QMessageBox::No);
    msgBox.setIcon(QMessageBox::Question);
    msgBox.exec();
    if(msgBox.clickedButton()==yesButton){
        if(radmc3dProcess!=NULL){
            if(radmc3dProcess->state()==QProcess::Running){
                this->statusBar()->show();
                this->statusBar()->showMessage(trUtf8("Be patient, radmc3d closing"));
                this->logger->redirectToLogFile(false);
                this->logger->writeToLogFile(trUtf8("Be patient, radmc3d closing"));
                this->statusBar()->show();
                this->killRadmc3d(event);
            }
        }
    }else if(msgBox.clickedButton()==saveButton){
        QSettings *settings=setup->getLocalSettings();
        QFile file(settings->fileName());
        file.remove();
        this->imageControllerWidget->save();
        this->save();
        if(radmc3dProcess!=NULL){
            if(radmc3dProcess->state()==QProcess::Running){
                this->statusBar()->show();
                this->statusBar()->showMessage(trUtf8("Be patient, radmc3d closing"));
                this->logger->redirectToLogFile(false);
                this->logger->writeToLogFile(trUtf8("Be patient, radmc3d closing"));
                this->statusBar()->show();
                this->killRadmc3d(event);
            }
        }
    }else  if(msgBox.clickedButton()==noButton){
        QCoreApplication::processEvents();
        if(event!=NULL)event->ignore();
        exitCode=0;
    }
    QSettings *setting=this->setup->getGlobalSettings();
    setting->setValue("startNextTimeFromCurrentDirectory",startNextTimeCheckbox->isChecked());
    setting->setValue("NativeMenuBar",nativeMenuBar);
    setting->setValue("ByteOrder",unformattedSetting->getByteOrder());
    setting->setValue("RecordSize",unformattedSetting->getRecordLengthSize());
    msgBox.close();
}
void MainWindow::sendRadmc3dExitCommand(){
    if(radmc3dProcess!=NULL &&radmc3dProcess->state()==QProcess::Running){
        QTextStream in(radmc3dProcess);
        in<<"exit"<<endl;
        in.flush();
        radmc3dProcess->waitForBytesWritten();
    }
    if(radmc3dProcess!=NULL)radmc3dProcess->waitForFinished(200);
    if(radmc3dProcess!=NULL && radmc3dProcess->state()==QProcess::Running){
        this->logger->writeToLogFile(trUtf8("The RADMC3D was terminated with \"exit\" command"));
        this->logger->writeToLogFile(trUtf8("The RADMC3D will killed now by force"));
        radmc3dProcess->kill();
    }
}
void MainWindow::killRadmc3d(){
    this->logger->writeToLogFile(trUtf8("You killed RADMC3D"));
    this->showLog=false;
    killRadmc3d(NULL);
}

void MainWindow::killRadmc3d(QCloseEvent *event){
    mutex.lock();
    exit=true;
    mutex.unlock();
    if(radmc3dProcess!=NULL){
        radmc3dProcess->close();
        radmc3dProcess->closeWriteChannel();
        int tryNr=1;
        while(radmc3dProcess->state()==QProcess::Running){
            radmc3dProcess->kill();
            radmc3dProcess->waitForFinished(300);
            this->logger->writeToLogFile(trUtf8("Try %1 to kill RADMC3D").arg(QString::number(tryNr)));
            tryNr++;
        }if(event!=NULL)event->accept();

    }
}

void MainWindow::freqUnitChanged(QAction *action){
    int var;
    unzoomSpectrum();
    for ( var= 0; var < freqUnitsList.size(); ++var) {
        if(freqUnitsList.at(var)==action->text())break;
    }
    spectrumLabel->setFrequencyUnit(static_cast<Spectrum::frequencyUnit>(var));
    if(spectrum!=NULL){
        reloadSpectrum();
    }
}
void MainWindow::fluxUnitChanged(QAction *action){
    int var;
    unzoomSpectrum();
    for ( var= 0; var < fluxUnitsList.size(); ++var) {
        if(fluxUnitsList.at(var)==action->text())break;
    }
    spectrumLabel->setFluxUnit(static_cast<Spectrum::fluxUnit>(var));
    if(spectrum!=NULL){
        reloadSpectrum();
    }
}
void MainWindow::saveAsDefaultsSlot(){
    this->saveAsDefaults();
    this->imageControllerWidget->saveAsDefaults();
    QMessageBox::information(this,trUtf8("Saved as defaults"),trUtf8("Saving defaults was successfully"));
}

/**
 * @brief This method saves the state of GUI-Elements.
 * @see load
*/
void MainWindow::saveAsDefaults(){
    QSettings *settings=setup->getGlobalSettings();
    settings->beginGroup("MainWindow");
    int i=0;
    for(i=0;i<sizemenu->actions().size();i++){
        if(sizemenu->actions().at(i)->isChecked())break;
    }
    if(i>=sizemenu->actions().size())i=0;
    settings->setValue(this->sizemenu->objectName(),QString::number(i));

    for(i=0;i<freqmenu->actions().size();i++){
        if(freqmenu->actions().at(i)->isChecked())break;
    }
    if(i>=freqmenu->actions().size())i=0;
    settings->setValue(this->freqmenu->objectName(),QString::number(i));

    for(i=0;i<fluxmenu->actions().size();i++){
        if(fluxmenu->actions().at(i)->isChecked())break;
    }
    if(i>=fluxmenu->actions().size())i=0;
    settings->setValue(this->fluxmenu->objectName(),QString::number(i));

    settings->setValue(this->userTransferAction->objectName(),this->userTransferAction->isChecked());
    settings->setValue(this->localAction->objectName(),this->localAction->isChecked());
    settings->setValue(this->colorAction->objectName(),this->colorAction->isChecked());
    settings->setValue(this->lineOptionAction->objectName(),this->lineOptionAction->isChecked());
    settings->setValue(this->xLinAction->objectName(),this->xLinAction->isChecked());
    settings->setValue(this->yLinAction->objectName(),this->yLinAction->isChecked());
    settings->endGroup();
}

/**
 * @brief This method saves the state of GUI-Elements.
 * @see load
*/
void MainWindow::save(){
    QSettings *settings=setup->getLocalSettings();
    settings->beginGroup("MainWindow");
    int i=0;
    for(i=0;i<sizemenu->actions().size();i++){
        if(sizemenu->actions().at(i)->isChecked())break;
    }
    if(i>=sizemenu->actions().size())i=0;
    settings->setValue(this->sizemenu->objectName(),QString::number(i));

    for(i=0;i<freqmenu->actions().size();i++){
        if(freqmenu->actions().at(i)->isChecked())break;
    }
    if(i>=freqmenu->actions().size())i=0;
    settings->setValue(this->freqmenu->objectName(),QString::number(i));

    for(i=0;i<fluxmenu->actions().size();i++){
        if(fluxmenu->actions().at(i)->isChecked())break;
    }
    if(i>=fluxmenu->actions().size())i=0;
    settings->setValue(this->fluxmenu->objectName(),QString::number(i));

    settings->setValue(this->userTransferAction->objectName(),this->userTransferAction->isChecked());
    settings->setValue(this->localAction->objectName(),this->localAction->isChecked());
    settings->setValue(this->colorAction->objectName(),this->colorAction->isChecked());
    settings->setValue(this->lineOptionAction->objectName(),this->lineOptionAction->isChecked());
    settings->setValue(this->xLinAction->objectName(),this->xLinAction->isChecked());
    settings->setValue(this->yLinAction->objectName(),this->yLinAction->isChecked());
    settings->endGroup();
}
/**
 * @brief This method loads the state of GUI-Elements.
 * @see save
*/
void MainWindow::load(){
    QSettings *settings=setup->getLocalSettings();
    settings->beginGroup("MainWindow");


    int value= settings->value(this->sizemenu->objectName()).toInt();
    this->sizemenu->actions().at(value)->setChecked(true);

    value= settings->value(this->freqmenu->objectName()).toInt();
    this->freqmenu->actions().at(value)->setChecked(true);

    value= settings->value(this->fluxmenu->objectName()).toInt();
    this->fluxmenu->actions().at(value)->setChecked(true);

    this->setup->setImageUnit(static_cast<SetUp::unit>(value));
    this->imageControllerWidget->setSizeUnit(static_cast<SetUp::unit>(value));


    value= settings->value(this->sizemenu->objectName()).toInt();
    this->sizemenu->actions().at(value)->setChecked(true);



    this->userTransferAction->setChecked(settings->value(this->userTransferAction->objectName()).toBool());
    this->localAction->setChecked(settings->value(this->localAction->objectName()).toBool());
    this->userTransferChanged(this->userTransferAction->isChecked());
    this->imageControllerWidget->activateLocalObserver(this->localAction->isChecked());
    if(userTransferAction->isChecked())this->imageControllerWidget->readTransferFile(true);
    this->localObserverChanged(this->localAction->isChecked(),false);
    this->colorAction->setChecked(settings->value(this->colorAction->objectName()).toBool());
    this->lineOptionAction->setChecked(settings->value(this->lineOptionAction->objectName()).toBool());
    this->lineOptionChanged(this->lineOptionAction->isChecked(),false);
    this->xLinAction->setChecked(settings->value(this->xLinAction->objectName()).toBool());
    this->yLinAction->setChecked(settings->value(this->yLinAction->objectName()).toBool());
    if(this->xLinAction->isChecked())
        this->xLinAction->setIcon(QIcon(trUtf8(":/images/%1.png").arg("linear")));
    else
        this->xLinAction->setIcon(QIcon(trUtf8(":/images/no%1.png").arg("linear")));
    if(this->yLinAction->isChecked())
        this->yLinAction->setIcon(QIcon(trUtf8(":/images/%1.png").arg("linear")));
    else
        this->yLinAction->setIcon(QIcon(trUtf8(":/images/no%1.png").arg("linear")));
    this->spectrumLabel->setXLogOn(!xLinAction->isChecked());
    this->spectrumLabel->setYLogOn(!yLinAction->isChecked());

    settings->endGroup();
}
/**
 * @brief This method loads the state of GUI-Elements as defaults.
 * @see save
*/
void MainWindow::loadAsdefaults(){
    QSettings *settings=setup->getGlobalSettings();
    settings->beginGroup("MainWindow");
    if(settings->childKeys().size()>0){

        int value= settings->value(this->sizemenu->objectName()).toInt();
        this->sizemenu->actions().at(value)->setChecked(true);

        value= settings->value(this->freqmenu->objectName()).toInt();
        this->freqmenu->actions().at(value)->setChecked(true);

        value= settings->value(this->fluxmenu->objectName()).toInt();
        this->fluxmenu->actions().at(value)->setChecked(true);

        this->setup->setImageUnit(static_cast<SetUp::unit>(value));
        this->imageControllerWidget->setSizeUnit(static_cast<SetUp::unit>(value));


        value= settings->value(this->sizemenu->objectName()).toInt();
        this->sizemenu->actions().at(value)->setChecked(true);



        this->userTransferAction->setChecked(settings->value(this->userTransferAction->objectName()).toBool());
        this->localAction->setChecked(settings->value(this->localAction->objectName()).toBool());
        this->userTransferChanged(this->userTransferAction->isChecked());
        this->imageControllerWidget->activateLocalObserver(this->localAction->isChecked());
        if(userTransferAction->isChecked())this->imageControllerWidget->readTransferFile(true);
        this->localObserverChanged(this->localAction->isChecked(),false);
        this->colorAction->setChecked(settings->value(this->colorAction->objectName()).toBool());
        this->lineOptionAction->setChecked(settings->value(this->lineOptionAction->objectName()).toBool());
        this->lineOptionChanged(this->lineOptionAction->isChecked(),false);
        this->xLinAction->setChecked(settings->value(this->xLinAction->objectName()).toBool());
        this->yLinAction->setChecked(settings->value(this->yLinAction->objectName()).toBool());
        if(this->xLinAction->isChecked())
            this->xLinAction->setIcon(QIcon(trUtf8(":/images/%1.png").arg("linear")));
        else
            this->xLinAction->setIcon(QIcon(trUtf8(":/images/no%1.png").arg("linear")));
        if(this->yLinAction->isChecked())
            this->yLinAction->setIcon(QIcon(trUtf8(":/images/%1.png").arg("linear")));
        else
            this->yLinAction->setIcon(QIcon(trUtf8(":/images/no%1.png").arg("linear")));
        this->spectrumLabel->setXLogOn(!xLinAction->isChecked());
        this->spectrumLabel->setYLogOn(!yLinAction->isChecked());
    }
    settings->endGroup();
}
/**
 * @brief This method set the shortcuts.
 *
 */
void MainWindow::setShortcuts(){
    this->manualAction->setShortcut(QKeySequence("CTRL+M"));
    this->guimanualAction->setShortcut(QKeySequence("SHIFT+M"));
    this->aboutAction->setShortcut(QKeySequence("CTRL+A"));
    this->whatsthisAction->setShortcut(QKeySequence(Qt::Key_F1));
    this->loadPointsAction->setShortcut(QKeySequence("SHIFT+P"));
}
void MainWindow::sizeUnitChanged( ){
    reloadImage();
}
void MainWindow::hideLog(){
    if(!this->imageControllerWidget->isRenderSpectrumActive())
        this->leftWidget->setCurrentIndex(this->leftWidget->indexOf(this->imageGroup));
    else
        this->leftWidget->setCurrentIndex(this->leftWidget->indexOf(this->spectrumGroup));

}
void MainWindow::readMyAction(){
    if(this->imageControllerWidget->isPipe()){
        if((radmc3dProcess!=NULL) & (radmc3dProcess->state()==QProcess::Running)){
            QTextStream in(radmc3dProcess);
            in<<"myaction"<<endl;
            in<<"enter"<<endl;
            in.flush();
            radmc3dProcess->waitForBytesWritten();
        }
        this->logger->writeToLogFile(trUtf8("Send command myAction to radmc3d"));
        emit sendCommand(this->imageControllerWidget->getRenderImagePushButtonObjectName());
    }
}

void MainWindow::setCancelable(bool cancelable){

    QList<QPushButton *> list=this->progressDialog->findChildren<QPushButton *>();
    if(cancelable)
        list.at(0)->show();
    else
        list.at(0)->hide();
}
void MainWindow::nativeMenuBarChanged(bool active){
    int result=QMessageBox::question(this,trUtf8("Do you want to restart the GUI now?"),trUtf8("Changes will take effect when you restart,Do you want to do that now?"),QMessageBox::Yes|QMessageBox::No);
    nativeMenuBar=active;
    if(result==QMessageBox::Yes){
        exitCode=EXIT_CODE_REBOOT;
        close();

    }
}

/**
 * @brief This method set the slots.
 *
 */
void MainWindow::setSlots(){
    this->connect(this->progressDialog,SIGNAL(canceled()),this,SLOT(killRadmc3d()));
    this->connect(this->escAction,SIGNAL(triggered()),this,SLOT(removeRect()));
    this->connect(this->fluxGroup,SIGNAL(triggered(QAction*)),this,SLOT(fluxUnitChanged(QAction*)));
    this->connect(this->nativeMenuBarAction,SIGNAL(triggered(bool)),this,SLOT(nativeMenuBarChanged(bool)));
    this->connect(this->aboutAction,SIGNAL(triggered()),this,SLOT(showAbout()));
    this->connect(this->manualAction,SIGNAL(triggered()),this,SLOT(showManual()));
    this->connect(this->guimanualAction,SIGNAL(triggered()),this,SLOT(showGUIManual()));
    this->connect(this->imageControllerWidget,SIGNAL(sendCommand(QString)),this,SLOT(sendCommand(QString)));
    this->connect(this->sizeGroup,SIGNAL(triggered(QAction*)),this->imageControllerWidget,SLOT(sizeUnitChanged(QAction*)));
    this->connect(this->sizeGroup,SIGNAL(triggered(QAction*)),this,SLOT(sizeUnitChanged()));
    this->connect(this->freqGroup,SIGNAL(triggered(QAction*)),this,SLOT(freqUnitChanged(QAction*)));
    this->connect(this->closeLogPushButton,SIGNAL(clicked()),this,SLOT(hideLog()));
    this->connect(this->imageLabel,SIGNAL(pixelValueChange(QString)),this,SLOT(showPixelValue(QString)));
    this->connect(this->spectrumLabel,SIGNAL(pixelValueChange(QString)),this,SLOT(showPixelValue(QString)));
    this->connect(this->spectrumLabel,SIGNAL(drawCircle(int,int)),this,SLOT(drawCircle(int,int)));
    this->connect(this->imageLabel,SIGNAL(drawRect(QPoint,QPoint)),this,SLOT(drawRect(QPoint,QPoint)));
    this->connect(this->spectrumLabel,SIGNAL(drawRect(QPoint,QPoint)),this,SLOT(drawRect(QPoint,QPoint)));
    this->connect(this->spectrumLabel,SIGNAL(removeZoom()),this,SLOT(removeZoom()));
    this->connect(this->imageLabel,SIGNAL(removeRect()),this,SLOT(removeRect()));
    this->connect(this->spectrumLabel,SIGNAL(removeRect()),this,SLOT(removeRect()));
    this->connect(this->imageLabel,SIGNAL(reloadImage()),this,SLOT(reloadImage()));
    this->connect(this->imageControllerWidget,SIGNAL(removeRect()),this,SLOT(removeRect()));
    this->connect(this->imageLabel,SIGNAL(calculateZoom(QPoint,QPoint)),this,SLOT(calculateZoom(QPoint,QPoint)));
    this->connect(this->spectrumLabel,SIGNAL(calculateZoom(QPoint,QPoint)),this,SLOT(calculateZoom(QPoint,QPoint)));
    this->connect(this->imageLabel,SIGNAL(phiAndinclinactionValueChanged(double,double)),this->imageControllerWidget,SLOT(changePhiAndInclination(double,double)));
    this->connect(this->userTransferAction,SIGNAL(triggered(bool)),this,SLOT(userTransferChanged(bool)));
    this->connect(this->localAction,SIGNAL(triggered(bool)),this,SLOT(localObserverChanged(bool)));
    this->connect(this->changeCurrentDirectoryAction,SIGNAL(triggered()),this,SLOT(changeCurrentDirectory()));
    this->connect(this->colorAction,SIGNAL(triggered(bool)),this,SLOT(colorChanged(bool)));
    this->connect(this->lineOptionAction,SIGNAL(triggered(bool)),this,SLOT(lineOptionChanged(bool)));
    this->connect(this,SIGNAL(readTransferFile()),this->imageControllerWidget,SLOT(readTransferFile()));
    this->connect(this,SIGNAL(removeTransferTab()),this->imageControllerWidget,SLOT(removeTransferTab()));
    this->connect(imageControllerWidget,SIGNAL(readMyAction()),this,SLOT(readMyAction()));
    this->connect(this->imageControllerWidget,SIGNAL(setActionVisible(QString,bool)),this,SLOT(setActionVisible(QString,bool)));
    this->connect(historyMenu,SIGNAL(triggered(QAction*)),this,SLOT(changeCurrentDirectoryByAction(QAction*)));
    this->connect(this->readImageAction,SIGNAL(triggered()),this,SLOT(readAndReloadImage()));
    this->connect(this->readSpectrumAction,SIGNAL(triggered()),this,SLOT(readAndReloadSpectrum()));
    this->connect(this,SIGNAL(activateLocalObserver(bool)),this->imageControllerWidget,SLOT(activateLocalObserver(bool)));
    this->connect(this->imageControllerWidget,SIGNAL(checktStateChanged(QString)),this,SLOT(changeIcon(QString)));
    this->connect(this->xLinAction,SIGNAL(triggered(bool)),this,SLOT(xLinearActive(bool)));
    this->connect(this->yLinAction,SIGNAL(triggered(bool)),this,SLOT(yLinearActive(bool)));
    this->connect(this,SIGNAL(setZoomOutVisible(bool)),this->spectrumLabel,SLOT(setZoomOutVisible(bool)));
    this->connect(this->savedefaultsAction,SIGNAL(triggered()),this,SLOT(saveAsDefaultsSlot()));
    this->connect(this->loadPointsAction,SIGNAL(triggered(bool)),this,SLOT(loadPointsSlot(bool)));

}

void MainWindow::unzoomSpectrum(){
    spectrumSettingVector.clear();
    this->spectrumLabel->setZoomActive(false);
    emit setZoomOutVisible(false);
}

void MainWindow::removeZoom(){
    if(this->spectrumSettingVector.size()>0){
        this->spectrumSettingVector.remove(this->spectrumSettingVector.size()-1);
        if(this->spectrumSettingVector.size()==0){
            this->spectrumLabel->setZoomActive(false);
            emit setZoomOutVisible(false);
        }
    }
    reloadSpectrum();
}

void MainWindow::xLinearActive(bool active){
    unzoomSpectrum();
    if(this->xLinAction->isChecked())
        this->xLinAction->setIcon(QIcon(trUtf8(":/images/%1.png").arg("linear")));
    else
        this->xLinAction->setIcon(QIcon(trUtf8(":/images/no%1.png").arg("linear")));
    this->spectrumLabel->setXLogOn(!active);
    if(spectrum!=NULL)
        reloadSpectrum();
}

void MainWindow::yLinearActive(bool active){
    unzoomSpectrum();
    if(this->yLinAction->isChecked())
        this->yLinAction->setIcon(QIcon(trUtf8(":/images/%1.png").arg("linear")));
    else
        this->yLinAction->setIcon(QIcon(trUtf8(":/images/no%1.png").arg("linear")));
    this->spectrumLabel->setYLogOn(!active);
    if(spectrum!=NULL)
        reloadSpectrum();
}

void MainWindow::readAndReloadSpectrum(){
    QString filename;
    QTextStream in;
    QFile file;
    QStringList list;
    if(spectrum==NULL)spectrum=new Spectrum();
    this->spectrum->findSpectrumFile(&filename,this->setup->getSpectrumFilename(),this->setup->getSpectrumEnding(),Qt::CaseSensitive,&list);
    if(list.size()>0){
        QDialog dialog;
        QVBoxLayout  box;
        label=new QLabel();
        label->setTextFormat(Qt::RichText);
        label->setText(trUtf8("Choose the spectrum file, you want to read?"));
        box.addWidget(label);
        QComboBox combo;
        combo.addItems(list);
        box.addWidget(&combo);
        QHBoxLayout hbox;
        QPushButton okButton;
        okButton.setText("Ok");
        QPushButton cancelButton;
        cancelButton.setText("Cancel");
        hbox.addWidget(&okButton);
        hbox.addWidget(&cancelButton);
        box.addLayout(&hbox);
        dialog.setLayout(&box);
        this->connect(&okButton,SIGNAL(clicked()),&dialog,SLOT(accept()));
        this->connect(&cancelButton,SIGNAL(clicked()),&dialog,SLOT(reject()));
        dialog.exec();
        if(dialog.result()==QDialog::Accepted){
            showProgressBar();
            spectrum=new Spectrum();
            spectrumIsRead=spectrum->readSpectrumFormatted(combo.itemText(combo.currentIndex()),"");
            if(!spectrumIsRead){
                showRadmc3dLog();
            }else{
                this->logger->writeToLogFile(trUtf8("The Spectrum file %1 is read").arg(combo.itemText(combo.currentIndex())));
                spectrumLabel->setSpectrum(spectrum);
                reloadSpectrum();
                hideProgressBar();
            }
        }else{
            this->disconnect(&okButton);
            this->disconnect(&cancelButton);
            this->disconnect(&combo);
            this->logger->writeToLogFile("Reading spectrum is Canceled");
        }
    }else{
        QMessageBox::information(this,trUtf8("Cannot load spectrum"),trUtf8("Cannot load spectrum, see logfile"));
        emit this->imageControllerWidget->selectDebugger();
    }
}
void MainWindow::loadPointsSlot(bool load){
    if(this->imageControllerWidget->isRenderSpectrumActive()){
        QString fileName;
        bool loadSuccessfully=load;
        if(load){
            fileName = QFileDialog::getOpenFileName(this,
                                                    tr("Load points from a file"), QDir::currentPath(), tr("*"));
            if(!fileName.trimmed().isEmpty()){
                this->loadPointsAction->setIcon(QIcon(trUtf8(":/images/%1.png").arg(this->loadPointsAction->objectName())));
                QFile file(fileName);
                if(file.exists()){
                    QString line;
                    file.open(QIODevice::ReadOnly);
                    QTextStream in(&file);
                    loadedPoints.clear();
                    bool check=false;
                    double velocityDivideByMicron=this->physicalConstants->getLightSpeed()/this->physicalConstants->getMicron();
                    QStringList list;
                    int size=0;
                    while(!in.atEnd()){
                        line=in.readLine();
                        list=line.split(QRegExp("\\s+"),QString::SkipEmptyParts);
                        size=list.size();
                        if(!check){
                            check=true;
                            if(size!=2 && size!=6){
                                loadSuccessfully=false;
                                QMessageBox::information(this,trUtf8("False format"),trUtf8("<html>File %1 should contains 2 columns(x y)<br/>or"
                                                                                            "should contains 6 columns(x y &Delta;x<sub>1</sub> &Delta;x<sub>2</sub> &Delta;y<sub>1</sub> &Delta;y<sub>2</sub>)<html>").arg(fileName));
                                break;
                            }else if(size==2)isErrorbar=false;
                            else isErrorbar=true;
                        }
                        if(!isErrorbar){
                            Point p;
                            p.x=velocityDivideByMicron/list.at(0).toDouble();
                            p.y=list.at(1).toDouble();
                            loadedPoints.append(p);
                        }else{
                            Point p;
                            p.x=velocityDivideByMicron/list.at(0).toDouble();
                            p.y=list.at(1).toDouble();
                            p.x1=velocityDivideByMicron/(list.at(0).toDouble()-list.at(2).toDouble());
                            p.x2=velocityDivideByMicron/(list.at(0).toDouble()+list.at(3).toDouble());
                            p.y1=p.y-list.at(4).toDouble();
                            p.y2=p.y+list.at(5).toDouble();
                            loadedPoints.append(p);
                        }
                    }
                }else{
                    loadSuccessfully=false;
                }
            }else{
                loadSuccessfully=false;
            }
        }
        if(!loadSuccessfully){
            this->loadPointsAction->setIcon(QIcon(trUtf8(":/images/no%1.png").arg(this->loadPointsAction->objectName())));
            this->loadPointsAction->blockSignals(true);
            this->loadPointsAction->setChecked(false);
            this->loadPointsAction->blockSignals(false);
        }
        reloadSpectrum();
    }else{
        this->loadPointsAction->setIcon(QIcon(trUtf8(":/images/no%1.png").arg(this->loadPointsAction->objectName())));
        this->loadPointsAction->blockSignals(true);
        this->loadPointsAction->setChecked(false);
        this->loadPointsAction->blockSignals(false);
        QMessageBox::information(this,trUtf8("Load a spectrum first"),trUtf8("You should load a spectrum first,bevor loading points"));
    }
}

void MainWindow::readAndReloadImage(){
    QString filename;
    QTextStream in;
    QFile file;
    QStringList list;
    if(image==NULL)image=new Image();
    this->image->findImageFile(&filename,this->setup->getImageFilename(),this->setup->getImageFileEnding(),Qt::CaseSensitive,&list);
    if(list.size()>0){
        QDialog dialog;
        formatted=new QCheckBox();
        formatted->setText(trUtf8("Formatted"));
        QVBoxLayout  box;
        label=new QLabel();
        label->setTextFormat(Qt::RichText);
        box.addWidget(label);
        box.addWidget(formatted);
        QComboBox combo;
        combo.addItems(list);
        box.addWidget(&combo);
        QHBoxLayout hbox;
        QPushButton okButton;
        okButton.setText("Ok");
        QPushButton cancelButton;
        cancelButton.setText("Cancel");
        hbox.addWidget(&okButton);
        hbox.addWidget(&cancelButton);
        box.addLayout(&hbox);
        dialog.setLayout(&box);
        checkFormatted(list.at(0));
        this->connect(&okButton,SIGNAL(clicked()),&dialog,SLOT(accept()));
        this->connect(&combo,SIGNAL(currentIndexChanged(QString)),this,SLOT(checkFormatted(QString)));
        this->connect(&cancelButton,SIGNAL(clicked()),&dialog,SLOT(reject()));
        formatted->blockSignals(true);
        dialog.exec();
        if(dialog.result()==QDialog::Accepted){
            showProgressBar();
            image=new Image();
            if(formatted->isChecked()){
                imageIsRead=image->readImageFormatted(combo.itemText(combo.currentIndex()),"");
            }else{
                imageIsRead=image->readImageUnformatted(combo.itemText(combo.currentIndex()),"");
            }
            if(!imageIsRead){
                showRadmc3dLog();
            }else{
                this->logger->writeToLogFile(trUtf8("The image file %1 is read").arg(combo.itemText(combo.currentIndex())));
                this->setZoomBox=false;
                this->imageControllerWidget->setZoomActivity(false);
                if(image->getNumberOfImages()==3){
                    colorChanged(true,false);
                }else{
                    colorChanged(false,false);
                }
                pixelEdit();
                this->changeImageColorLookUpTable(this->imageControllerWidget->isColorInvert());
                if(this->imageControllerWidget->isContour())calculateContour();
                reloadImage();
                imageLabel->setImage(image);
                this->imageControllerWidget->setNumberOfPixelsSpinBox(image->getNumberOfPixelX());
                reloadImage();
                hideProgressBar();
            }
        }else{
            this->disconnect(&okButton);
            this->disconnect(&cancelButton);
            this->disconnect(&combo);
            this->logger->writeToLogFile("Reading image is Canceled");
        }
    }else{
        QMessageBox::information(this,trUtf8("Cannot load image"),trUtf8("Cannot load image, see logfile"));
        emit this->imageControllerWidget->selectDebugger();
    }
}
void MainWindow::checkFormatted(QString filename){
    QFile file(filename);
    if(file.exists()){
        file.open(QFile::ReadOnly);
        QString line=file.readLine();
        formatted->setChecked(line.contains("\n"));
        this->label->setTextFormat(Qt::RichText);
        if(formatted->isChecked())
            this->label->setText(trUtf8("<html>Choose the image, you want to read(<font color=blue>Recommanded formatted</font>)?<br/><br/><font color=red>"
                                        "Attention: wrong format will crash GUI</html>"));
        else
            this->label->setText(trUtf8("<html>Choose the image, you want to read(<font color=blue>Recommanded unformatted</font>)?<br/><br/><font color=red>"
                                        "Attention: wrong format will crash GUI</html>"));
        file.close();
    }
}

void MainWindow::setActionVisible(const QString objectName, bool visibility){
    QAction *action=this->findChild<QAction*>(objectName+"Action");
    if(action!=NULL){
        action->setVisible(visibility);
    }
    this->imageMenu->repaint();
}

void MainWindow::drawRect(QPoint first, QPoint last){
    if(!this->imageControllerWidget->isRenderSpectrumActive()){
        QPainter painter;
        QPixmap  pixmap=this->image->getPixmap();
        painter.begin(&pixmap);
        if(painter.isActive()){
            QPen pen=painter.pen();
            pen.setColor(palette().color(QPalette::Highlight));
            painter.setPen(pen);
            zoomActive=true;
            setZoomBox=true;
            QRect zoomRect=QRect(first.x(),first.y(),last.x()-first.x(),last.y()-first.y());
            painter.drawRect(zoomRect);
            painter.end();
        }
        this->imageLabel->setPixmap(pixmap);
    }else{
        QPainter painter;
        QPixmap  pixmap=this->spectrum->getPixmap();
        painter.begin(&pixmap);
        if(painter.isActive()){
            QPen pen=painter.pen();
            pen.setColor(palette().color(QPalette::Highlight));
            painter.setPen(pen);
            setZoomBox=true;
            QRect zoomRect=QRect(first.x(),first.y(),last.x()-first.x(),last.y()-first.y());
            painter.drawRect(zoomRect);
            painter.end();
        }
        this->spectrumLabel->setPixmap(pixmap);
    }
}
void MainWindow::removeRect(){
    if(!this->imageControllerWidget->isRenderSpectrumActive()){
        if(setZoomBox){
            setZoomBox=false;
            this->imageLabel->setPixmap(this->image->getPixmap());
        }
    }else{
        if(setZoomBox){
            setZoomBox=false;
            this->spectrumLabel->setPixmap(this->spectrum->getPixmap());
        }
    }
}

void MainWindow::showPixelValue(QString value){
    this->statusBar()->show();
    this->statusBar()->showMessage(value);
}
void MainWindow::drawCircle(int x,int y){
    QPainter painter;
    QPixmap pixmap=spectrum->getPixmap();
    if(!pixmap.isNull()){
        int size=pixmap.width();
        painter.begin(&pixmap);
        if(painter.isActive()){
            painter.setPen(Qt::blue);
            painter.setWindow(0,size,size,-size);
            painter.drawEllipse(x-5,y-5,10,10);
            painter.setWindow(0,0,size,size);
            painter.end();
        }
        spectrumLabel->setPixmap(pixmap);
        spectrumLabel->repaint();
    }
}

void MainWindow::showProgressBar(){
    timer=new QTimer();
    this->changeCurrentDirectoryAction->setEnabled(false);
    QList<QAction*> historyActions=this->historyMenu->findChildren<QAction*>();
    foreach(QAction *action,historyActions){
        action->setEnabled(false);
    }
    QList<QAction*> imageActions=this->imageMenu->findChildren<QAction*>();
    foreach(QAction *action,imageActions){
        action->setEnabled(false);
    }

    this->statusBar()->show();
    this->localAction->setEnabled(false);
    this->userTransferAction->setEnabled(false);
    this->progressDialog->setVisible(true);
    this->readImageAction->setEnabled(false);
    this->sizemenu->setEnabled(false);
    this->freqmenu->setEnabled(false);
    this->fluxmenu->setEnabled(false);
    this->imageControllerWidget->setEnabled(false);
    this->imageControllerWidget->blockSignals(true);
    this->imageLabel->blockSignals(true);
    this->readSpectrumAction->setEnabled(false);
    this->colorAction->setEnabled(false);
    this->lineOptionAction->setEnabled(false);
    this->colorAction->blockSignals(true);
    this->lineOptionAction->blockSignals(true);
    this->loadPointsAction->setEnabled(false);
    this->imageLabel->setMouseTracking(false);
    this->statusBar()->clearMessage();
    connect(timer,SIGNAL(timeout()),this->statusBar(),SLOT(repaint()));
    timer->start(200);

}

void MainWindow::hideProgressBar(){
    this->imageControllerWidget->blockSignals(false);
    this->changeCurrentDirectoryAction->setEnabled(true);
    QList<QAction*> historyActions=this->historyMenu->findChildren<QAction*>();
    foreach(QAction *action,historyActions){
        action->setEnabled(true);
    }
    QList<QAction*> imageActions=this->imageMenu->findChildren<QAction*>();
    foreach(QAction *action,imageActions){
        action->setEnabled(true);
    }

    this->readImageAction->setEnabled(true);
    this->readSpectrumAction->setEnabled(true);
    this->sizemenu->setEnabled(true);
    this->freqmenu->setEnabled(true);
    this->fluxmenu->setEnabled(true);
    this->userTransferAction->setEnabled(true);
    this->loadPointsAction->setEnabled(true);
    this->localAction->setEnabled(true);
    this->statusBar()->hide();
    if(timer!=NULL)timer->stop();
    this->progressDialog->setVisible(false);
    this->colorAction->setEnabled(true);
    this->lineOptionAction->setEnabled(true);
    this->imageControllerWidget->blockMoleculeSpinBoxSignal(false);
    this->imageControllerWidget->blockLineSpinBoxSignal(false);
    this->imageControllerWidget->blockVelocitySpinBoxSignal(false);
    this->imageLabel->blockSignals(false);
    this->colorAction->blockSignals(false);
    this->lineOptionAction->blockSignals(false);
    this->imageLabel->setMouseTracking(true);
    this->imageControllerWidget->setEnabled(true);
}
void MainWindow::hideProgressBarAndShowRadmc3dLog(){
    hideProgressBar();
    showRadmc3dLog();
}
void MainWindow::hideProgressBarAndReload(){
    hideProgressBar();
    reloadImage();
    this->imageControllerWidget->setEnabled(true);
}
void MainWindow::hideProgressBarAndEditAndReload(){
    if(!this->imageControllerWidget->isRenderSpectrumActive()){
        imageIsRead=false;
        if(!image->isError())imageIsRead=true;
        if(imageIsRead){
            if(this->imageControllerWidget->isDoLine()){
                this->imageControllerWidget->setDoLine(false);
                double lambda=this->image->getLambdaAt(0);
                this->imageControllerWidget->setRedLambdaLineEdit(lambda);
            }
            imageLabel->setImage(image);
            hideProgressBar();
            pixelEdit();
            this->changeImageColorLookUpTable(this->imageControllerWidget->isColorInvert());
            if(this->imageControllerWidget->isContour())calculateContour();
            reloadImage();
            this->imageControllerWidget->setEnabled(true);
        }else{
            this->logger->writeToLogFile(trUtf8("RADMC3D is crashed, bevor writting  to the PIPE was completed"));
            hideProgressBarAndShowRadmc3dLog();
        }
        waitForFinishing=false;
    }else{
        spectrumIsRead=false;
        if(!spectrum->isError())spectrumIsRead=true;
        if(spectrumIsRead){
            spectrumLabel->setSpectrum(spectrum);
            hideProgressBar();
            this->spectrumIsRead=true;
            this->reloadSpectrum();
            this->imageControllerWidget->setEnabled(true);
        }else{
            this->logger->writeToLogFile(trUtf8("RADMC3D is crashed, bevor writting  to the PIPE was completed"));
            hideProgressBarAndShowRadmc3dLog();
        }
        waitForFinishing=false;
    }
}
void MainWindow::imageEdit(bool changeColorTable,bool  invertPixel,bool calculateContour,bool pixelEdit){
    if(pixelEdit)this->pixelEdit();
    if(changeColorTable){
        this->changeImageColorLookUpTable(invertPixel);
    }
    if(calculateContour)this->calculateContour();
}

void MainWindow::sendCommand(const QString objectName){
    if(!this->imageControllerWidget->isObjectNameBelongsToturnToRenderImageList(objectName)&& image!=NULL){
        this->imageControllerWidget->setRenderSpectrumActive(false);
    }

    this->hideLog();
    disconnect(&futureWacther);
    if(objectName==this->imageControllerWidget->getInvertColorsCheckBoxObjectName()){
        if(image!=NULL){
            this->showProgressBar();
            future=QtConcurrent::run(this,&MainWindow::imageEdit,true,true,false,false);
            connect(&futureWacther,SIGNAL(finished()),SLOT(hideProgressBarAndReload()));
            futureWacther.setFuture(future);
        }

    }else if(objectName==this->imageControllerWidget->getColorTableSliderObjectName()){
        if(image!=NULL){
            this->showProgressBar();
            bool invertpixel=this->imageControllerWidget->isColorInvert();
            future=QtConcurrent::run(this,&MainWindow::imageEdit,true,invertpixel,this->imageControllerWidget->isContour(),true);
            connect(&futureWacther,SIGNAL(finished()),SLOT(hideProgressBarAndReload()));
            futureWacther.setFuture(future);
        }

    }else if(objectName==this->imageControllerWidget->getAntiAliasingCheckBoxObjectName()){
        if(image!=NULL){
            if(this->imageControllerWidget->isAntiAliasing()){
                this->transformMode=Qt::SmoothTransformation;
            }else{
                this->transformMode=Qt::FastTransformation;
            }
            this->reloadImage();
        }

    }else if(objectName==this->imageControllerWidget->getLinearObjectName()||
             this->imageControllerWidget->getContrastObjectName()==objectName||
             objectName==this->imageControllerWidget->getSaturateObjectName()){
        if(((this->imageControllerWidget->isRenderSpectrumActive())&& (this->imageControllerWidget->getContrastObjectName()==objectName))){
            if(spectrum!=NULL){
                reloadSpectrum();
            }
        }else{
            if(image!=NULL){
                this->showProgressBar();
                future=QtConcurrent::run(this,&MainWindow::imageEdit,false,false,true,true);
                connect(&futureWacther,SIGNAL(finished()),SLOT(hideProgressBarAndReload()));
                futureWacther.setFuture(future);
            }
        }

    }else if(objectName==this->imageControllerWidget->getColorBarObjectName()){
        if(image!=NULL){
            this->reloadImage();
        }

    }else if(objectName==this->imageControllerWidget->getNumberOfContoursObjectName()||
             this->imageControllerWidget->getContourCheckBoxObjectName()==objectName){
        if(image!=NULL){
            this->showProgressBar();
            future=QtConcurrent::run(this,&MainWindow::imageEdit,false,false,true,false);
            connect(&futureWacther,SIGNAL(finished()),SLOT(hideProgressBarAndReload()));
            futureWacther.setFuture(future);
        }

    }else if(objectName==this->imageControllerWidget->getRedColorTuneValueLineEditObjectName()||
             objectName==this->imageControllerWidget->getBlueColorTuneValueLineEditObjectName()||
             objectName==this->imageControllerWidget->getGreenColorTuneValueLineEditObjectName()||
             objectName==this->imageControllerWidget->getAbsoluteScaleCheckBoxObjectName()){
        if(image!=NULL){
            this->showProgressBar();
            bool invertpixel=this->imageControllerWidget->isColorInvert();
            future=QtConcurrent::run(this,&MainWindow::imageEdit,true,invertpixel,false,true);
            connect(&futureWacther,SIGNAL(finished()),SLOT(hideProgressBarAndReload()));
            futureWacther.setFuture(future);
        }
    } else if(objectName==this->imageControllerWidget->getPositionComboxBoxObjectName()){
        int pos=this->imageControllerWidget->getPosition();
        if(pos==0)
            this->setMainLayoutDirection(QBoxLayout::LeftToRight);
        else if(pos==1)
            this->setMainLayoutDirection(QBoxLayout::RightToLeft);
        else if(pos==2)
            this->setMainLayoutDirection(QBoxLayout::TopToBottom);
        else
            this->setMainLayoutDirection(QBoxLayout::BottomToTop);

    }else if(this->imageControllerWidget->getShowRadmc3dPushButtonObjectName()==objectName){
        showRadmc3dLog();

    }else if(this->imageControllerWidget->getPipeCheckBoxObjectName()==objectName){
        this->imageLabel->setPhiAndInclination(this->imageControllerWidget->getPhi(),this->imageControllerWidget->getInclination());
        if(setZoomBox && zoomBox.size()==4){
            this->imageControllerWidget->setZoomBox(zoomBox);
            zoomBox.clear();
        }
        if(initialized){
            if(!waitForFinishing){
                if(!this->imageControllerWidget->isObjectNameBelongsToturnToRenderImageList(objectName)){
                    this->imageControllerWidget->setRenderSpectrumActive(false);
                }
                waitForFinishing=true;
                if(this->imageControllerWidget->isPipe())
                    startRadmc3dAsChild();
                else
                    startRadmc3d();
            }

        }else{
            initialize();
            this->progressDialog->setVisible(false);
        }
    }else {
        if(setZoomBox && zoomBox.size()==4){
            this->imageControllerWidget->setZoomBox(zoomBox);
            zoomBox.clear();
        }
        if(objectName==this->imageControllerWidget->getUnzoomImagePushButtonObjectName()||objectName==this->imageControllerWidget->getSizeLineEditObjectName())this->imageControllerWidget->setZoomActivity(false);
        this->imageLabel->setPhiAndInclination(this->imageControllerWidget->getPhi(),this->imageControllerWidget->getInclination());
        if(initialized){
            if(!waitForFinishing){
                waitForFinishing=true;
                if(!this->imageControllerWidget->isObjectNameBelongsToturnToRenderImageList(objectName)){
                    this->imageControllerWidget->setRenderSpectrumActive(false);
                }
                if(this->imageControllerWidget->isRenderSpectrumActive()&& (objectName==this->imageControllerWidget->getUpPushButtonObjectName()||
                                                                            objectName==this->imageControllerWidget->getDownPushButtonObjectName()||
                                                                            objectName==this->imageControllerWidget->getLeftPushButtonObjectName()||
                                                                            objectName==this->imageControllerWidget->getRightPushButtonObjectName())){
                    if(this->spectrumLabel->isZoomActive()){
                        if(objectName==this->imageControllerWidget->getUpPushButtonObjectName()){
                            double ddd=0.1;
                            spectrumSetting *spectrumsetting=this->spectrumSettingVector.at(this->spectrumSettingVector.size()-1);
                            double dy=spectrumsetting->getYEnd()-spectrumsetting->getYStart();
                            spectrumsetting->setYStart(spectrumsetting->getYStart()+dy*ddd);
                            spectrumsetting->setYEnd(spectrumsetting->getYEnd()+dy*ddd);
                            this->spectrumSettingVector.replace(this->spectrumSettingVector.size()-1,spectrumsetting);
                        }else if(objectName==this->imageControllerWidget->getDownPushButtonObjectName()){
                            double ddd=0.1;
                            spectrumSetting *spectrumsetting=this->spectrumSettingVector.at(this->spectrumSettingVector.size()-1);
                            double dy=spectrumsetting->getYEnd()-spectrumsetting->getYStart();
                            spectrumsetting->setYStart(spectrumsetting->getYStart()-dy*ddd);
                            spectrumsetting->setYEnd(spectrumsetting->getYEnd()-dy*ddd);
                            this->spectrumSettingVector.replace(this->spectrumSettingVector.size()-1,spectrumsetting);

                        }else if(objectName==this->imageControllerWidget->getLeftPushButtonObjectName()){
                            double ddd=0.1;
                            spectrumSetting *spectrumsetting=this->spectrumSettingVector.at(this->spectrumSettingVector.size()-1);
                            double dx=spectrumsetting->getXEnd()-spectrumsetting->getXStart();
                            spectrumsetting->setXStart(spectrumsetting->getXStart()-dx*ddd);
                            spectrumsetting->setXEnd(spectrumsetting->getXEnd()-dx*ddd);
                            this->spectrumSettingVector.replace(this->spectrumSettingVector.size()-1,spectrumsetting);

                        }else if(objectName==this->imageControllerWidget->getRightPushButtonObjectName()){
                            double ddd=0.1;
                            spectrumSetting *spectrumsetting=this->spectrumSettingVector.at(this->spectrumSettingVector.size()-1);
                            double dx=spectrumsetting->getXEnd()-spectrumsetting->getXStart();
                            spectrumsetting->setXStart(spectrumsetting->getXStart()+dx*ddd);
                            spectrumsetting->setXEnd(spectrumsetting->getXEnd()+dx*ddd);
                            this->spectrumSettingVector.replace(this->spectrumSettingVector.size()-1,spectrumsetting);
                        }
                    }
                    waitForFinishing=false;
                    reloadSpectrum();
                }else{
                    if(this->imageControllerWidget->isPipe())
                        startRadmc3dAsChild();
                    else
                        startRadmc3d();
                }
            }

        }else{
            this->progressDialog->setVisible(false);
            initialize();
        }
    }

}


/**
 * @brief This method set the GUI-elements's labels.
 *
 * @details This method should be called, if the GUI-language is changed, so that the translations will be loaded.
 */
void MainWindow::translate(){
    this->imageMenu->setTitle(trUtf8("&Image"));
    this->spectrumMenu->setTitle(trUtf8("&Spectrum"));
    this->configMenu->setTitle(trUtf8("&Config"));
    this->helpMenu->setTitle(trUtf8("&Help"));
    this->whatsthisAction->setText(trUtf8("What's This?"));
    this->setWindowTitle(trUtf8("RADMC3DGUI"));
    this->manualAction->setText(trUtf8("RADMC3D's& manual"));
    this->guimanualAction->setText(trUtf8("GUI's man&ual"));
    this->aboutAction->setText("&About");
    this->userTransferAction->setText("Trans&fer");
    this->localAction->setText("&LocalObserver");
    this->colorAction->setText(trUtf8("RGB mo&dus"));
    this->lineOptionAction->setText(trUtf8("Li&ne modus"));
    this->changeCurrentDirectoryAction->setText(trUtf8("Change direc&tory"));
    this->rightWidget->setTabText(0,this->imageControllerWidgetTabString);
    this->historyMenu->setTitle(trUtf8("Hist&ory"));
    this->readImageAction->setText(trUtf8("&Read Image"));
    this->readSpectrumAction->setText(trUtf8("Read Spectr&um"));
    this->xLinAction->setText(trUtf8("&XLinear"));
    this->yLinAction->setText(trUtf8("&YLinear"));
    this->nativeMenuBarAction->setText(trUtf8("Native menuBar"));
    this->savedefaultsAction->setText(trUtf8("Save as defaults"));
    this->loadPointsAction->setText(trUtf8("Load Points"));
}
/**
 * @brief This method set the icons.
 *
 */
void MainWindow::setIcons(){
    QIcon icon(":/images/radmc3d.png");
    this->setMainWindowIcon(icon);
    this->manualAction->setIcon(icon);
    this->guimanualAction->setIcon(icon);
    this->aboutAction->setIcon(QApplication::style()->standardIcon(QStyle::SP_MessageBoxInformation));
    this->whatsthisAction->setIcon(QApplication::style()->standardIcon(QStyle::SP_DialogHelpButton));
    this->closeLogPushButton->setFixedSize(24,24);
    this->closeLogPushButton->setIcon(QApplication::style()->standardIcon(QStyle::SP_DialogCloseButton));
    this->userTransferAction->setIconVisibleInMenu(true);
    this->colorAction->setIconVisibleInMenu(true);
    this->savedefaultsAction->setIcon(QIcon(trUtf8(":/images/%1.png").arg(this->savedefaultsAction->objectName())));
    if(userTransferAction->isChecked())
        this->userTransferAction->setIcon(QIcon(trUtf8(":/images/%1.png").arg(this->userTransferAction->objectName())));
    else
        this->userTransferAction->setIcon(QIcon(trUtf8(":/images/no%1.png").arg(this->userTransferAction->objectName())));
    if(colorAction->isChecked())
        this->colorAction->setIcon(QIcon(trUtf8(":/images/%1.png").arg(this->colorAction->objectName())));
    else
        this->colorAction->setIcon(QIcon(trUtf8(":/images/no%1.png").arg(this->colorAction->objectName())));
    if(lineOptionAction->isChecked())
        this->lineOptionAction->setIcon(QIcon(trUtf8(":/images/%1.png").arg(this->lineOptionAction->objectName())));
    else
        this->lineOptionAction->setIcon(QIcon(trUtf8(":/images/no%1.png").arg(this->lineOptionAction->objectName())));

    if(this->localAction->isChecked())
        this->localAction->setIcon(QIcon(trUtf8(":/images/%1.png").arg(this->localAction->objectName())));
    else
        this->localAction->setIcon(QIcon(trUtf8(":/images/no%1.png").arg(this->localAction->objectName())));
    if(this->xLinAction->isChecked())
        this->xLinAction->setIcon(QIcon(trUtf8(":/images/%1.png").arg("linear")));
    else
        this->xLinAction->setIcon(QIcon(trUtf8(":/images/no%1.png").arg("linear")));
    if(this->yLinAction->isChecked())
        this->yLinAction->setIcon(QIcon(trUtf8(":/images/%1.png").arg("linear")));
    else
        this->yLinAction->setIcon(QIcon(trUtf8(":/images/no%1.png").arg("linear")));


    this->loadPointsAction->setIcon(QIcon(trUtf8(":/images/no%1.png").arg(loadPointsAction->objectName())));
    this->changeCurrentDirectoryAction->setIcon(QIcon(trUtf8(":/images/%1.png").arg(this->changeCurrentDirectoryAction->objectName())));
    QList<QPushButton *> list=this->progressDialog->findChildren<QPushButton *>();
    list.at(0)->setIconSize(QSize(25,25));
    list.at(0)->setIcon(QIcon(trUtf8(":/images/kill.png")));


}
/**
 * @brief This method set the mainWindow icon.
 * @param the mainWindow will have the given icon.
 */
void MainWindow::setMainWindowIcon(const QIcon& icon){
    this->setWindowIcon(icon);
}

/**
 * @brief This slot show the about's messagebox
 *
 */
void MainWindow::showAbout(){
    QMessageBox::about(this,trUtf8("RADMC3D"),
                       trUtf8("<html><img src=\":/images/copyleft.png\"> RADMC3D written by <a href=\"http://www.ita.uni-heidelberg.de/~dullemond/\">Prof. Dr. Cornelis P. Dullemond</a>"
                              "<br/><br/><img src=\":/images/copyleft.png\"> GUI written by <a href=\"mailto:sereshti@uni-heidelberg.de\">Farzin Sereshti</html>"));
}

/**
 * @brief This method looking after radmc3dstarterscrip.
 * @details This method searching radmc3dstarterscript, which is usually
 * created in ~/bin/ folder. This method looks after  radmc3dstarterscrip in all folder, which are defiened
 * in environment variable PATH.
 *
 * @return  true, if the radmc3dstarterscript is found. Otherwise returns false
*/
bool MainWindow::findRadmc3dStarterPath(){
    QStringList environmentVariables=QProcess::systemEnvironment();
    QStringListIterator iterator(environmentVariables);
    QString systemEnvironment;
    QDir radmc3dDirectory;
    QString  binDirectoryString= QDir::homePath();
    if(!binDirectoryString.endsWith(QDir::separator())){
        binDirectoryString.append(QDir::separator());
    }
    binDirectoryString.append("bin");
    radmc3dDirectory.setPath(binDirectoryString);
    QString filePath=radmc3dDirectory.filePath("radmc3d");
    if(QFile(filePath).exists()){
        this->radmc3dStarterPATH=filePath;
        return true;
    }
     while(iterator.hasNext()){
        systemEnvironment= iterator.next();
        if(systemEnvironment.startsWith("PATH",Qt::CaseSensitive)){
            systemEnvironment.remove(0,systemEnvironment.indexOf("=")+1);
            this->PATHDirectories= systemEnvironment.split(":");
            for(int i=0;i<PATHDirectories.size();i++){
                radmc3dDirectory.setPath(this->PATHDirectories.at(i));
                QString filePath=radmc3dDirectory.filePath("radmc3d");
                if(QFile(filePath).exists()){
                    this->radmc3dStarterPATH=filePath;
                    return true;
                }
            }
        }
    }
    return false;
}
void MainWindow::showGUIManual(){
    QString path=QCoreApplication::applicationDirPath();
    int index=path.indexOf("GUI")+3;
    path=path.remove(index,path.length()-index);
    path.append("/manual/viewimage_manual.pdf");
    path=QDir::toNativeSeparators(path);
    if(!QDesktopServices::openUrl(QUrl::fromLocalFile(QDir::toNativeSeparators(path))))
        QMessageBox::critical(this,trUtf8("RADMC3D's manual could not be opened"),trUtf8("QT could not find a pdf reader to open %1").arg(path));
    else{		
        this->statusBar()->show();
        this->statusBar()->showMessage(trUtf8("GUI-Manual file %1 is loaded").arg(path));
        this->progressDialog->hide();
    }
}

/**
 *  @brief This method show the current version radmc3d's manual.
 *
 *  @details This method can just show radmc3d's manual, if the radmc3dstarterscript is found and the manual's name of
 *  RADMC3D-pdffile is defined as radmc-3d*.pdf in radmc3d's manual folder.
 *
 *  @see readRadmc3dPath()
 *  @see findRadmc3dStarterPath()
 */
void MainWindow::showManual(){
    if(findRadmc3dStarterPath()){
        if(readRadmc3dPath()){
            QString radmc3dManualDir=this->radmc3dPATH+"manual"+QDir::separator();
            QDir manualDir(radmc3dManualDir);
            if(manualDir.exists()){
                QDirIterator iterator(manualDir);
                bool found=false;
                QString radmc3dManualFile;
                while(!found && iterator.hasNext()){
                    iterator.next();
                    radmc3dManualFile=iterator.fileName();
                    if(radmc3dManualFile.startsWith("radmc-3d") && radmc3dManualFile.endsWith("pdf")){
                        found=true;
                        radmc3dManualFile=iterator.filePath();
                    }
                }
                if(found){
                    if(!QDesktopServices::openUrl(QUrl::fromLocalFile(radmc3dManualFile))){
                        QMessageBox::critical(this,trUtf8("RADMC3D's manual could not be opened"),trUtf8("QT could not find a pdf reader to open %1").arg(radmc3dManualFile));
                    }else{
                        this->statusBar()->show();
                        this->statusBar()->showMessage(trUtf8("RADMC3D's Manual file %1 is loaded").arg(radmc3dManualFile));
                        this->progressDialog->hide();
                    }
                }else
                    QMessageBox::critical(this,trUtf8("RADMC3D's manual is missing"),trUtf8("There is no manual PDF-file in the following path:<br/> %1").arg(radmc3dManualDir));
            }
        }
    }else{

        QMessageBox::critical(this,trUtf8("RADMC3D's manual is missing"),trUtf8("RADMC3D's starter script could not be found in the following directories(system PATH):<br/> %1").arg(PATHDirectories.join(" ")));
    }

}
/**
 * @brief This method finds out the current Path of installed RADMC3D.
 *
 * @details This method reads the the current Path of RADMC3D in radmc3dstarterscript, if the script is found.
 * This script is important for showing  PDF-file later. The radmc3d's pathstring in radmc3dstarterscript should start with " and contains substring src.
 *
 * @return true, if the radmc3d's path is found and the radmc3dstartpath is readable. Otherwise returns false.
 *
 * @see findRadmc3dStarterPath()
 */
bool MainWindow::readRadmc3dPath(){
    QFile file(this->radmc3dStarterPATH);
    if(file.open(QFile::ReadOnly)){
        QByteArray byte=file.readAll();
        QString string=byte.constData();
        int indexofPathBegin=string.indexOf('"'+QDir::separator())+1;
        int indexofPathEnd=string.indexOf("src");
        this->radmc3dPATH=string.mid(indexofPathBegin,indexofPathEnd-indexofPathBegin).trimmed();
        file.close();
    }else{
        QMessageBox::critical(this,trUtf8("RADMC3D's manual is missing"),trUtf8("RADMC3D's starterscript %1 could not be opened").arg(this->radmc3dStarterPATH));
        return false;
    }
    return true;
}
