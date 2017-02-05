#include "makeimage.h"

MakeImage::MakeImage(ImageControllerWidget *imageControllerWidget,QObject *parent) :
    QObject(parent)
{
    this->imageControllerWidget=imageControllerWidget;
    this->setup=SetUp::getInstance();
    this->physicalConstants=PhysicalConstants::getInstance();
    this->logger=Logger::getInstance();
}


void MakeImage::sendWriteImageCommandToRadmc3d(QProcess *process){

    QTextStream in(process);
    in<<"writeimage"<<endl;
    in.flush();

    process->waitForBytesWritten();

}


void MakeImage::sendMakeImageCommandToRadmc3d(QProcess *process){
    bool localObserver=this->imageControllerWidget->isLocalModus();
    QString numberOfPixels=QString::number(this->imageControllerWidget->getNumberOfPixels());
    QString lambda=QString::number(this->imageControllerWidget->getLambdaRed());
    QString posang=QString::number(this->imageControllerWidget->getPositionAngle());
    QString incl=QString::number(this->imageControllerWidget->getInclination());
    QString phi=QString::number(this->imageControllerWidget->getPhi());
    if(imageControllerWidget->isRGBModus()) imageControllerWidget->writeLambdasToCamerWavelengthFilename();
    QVector<double> sizes;
    double size;
    if(!this->imageControllerWidget->isZoomActive()){
        size=this->imageControllerWidget->getSize()/2.;
        sizes.append(size*-1.);
        sizes.append(size);
        sizes.append(size*-1.);
        sizes.append(size);
    }else{
        sizes=this->imageControllerWidget->getZoomBox();
    }
    QVector<double> points;
    QVector<double> obs;
    QTextStream in(process);
    in<<"image"<<endl;
    in<<"npix"<<endl;
    in<<numberOfPixels.toStdString().c_str()<<endl;
    if(!imageControllerWidget->isRGBModus()){
        if(!this->imageControllerWidget->isDoLine()){
            in<<"lambda"<<endl;
            in<<lambda.toStdString().c_str()<<endl;
        }
    }else{
        in<<"loadlambda"<<endl;
    }
    in<<"posang"<<endl;
    in<<posang.toStdString().c_str()<<endl;
    if(!localObserver){
        in<<"incl"<<endl;
        in<<incl.toStdString().c_str()<<endl;
        in<<"phi"<<endl;
        in<<phi.toStdString().c_str()<<endl;
    }else{
        points.resize(3);
        obs.resize(3);
        points[0]=this->imageControllerWidget->getPointXLineEditValue();
        points[1]=this->imageControllerWidget->getPointYLineEditValue();
        points[2]=this->imageControllerWidget->getPointZLineEditValue();
        obs[0]=this->imageControllerWidget->getObserverposXLineEditValue();
        obs[1]=this->imageControllerWidget->getObserverposYLineEditValue();
        obs[2]=this->imageControllerWidget->getObserverposZLineEditValue();
        if(this->imageControllerWidget->isRelativeObserver()){
            obs[0]+=points[0];
            obs[1]+=points[1];
            obs[2]+=points[2];
        }
        in<<"sizeradian"<<endl;
        in<<QString::number(this->imageControllerWidget->getViewangleSpinBoxValue()).toStdString().c_str()<<endl;
    }
    if(setup->getImageUnit()==SetUp::CM){
        sizes[0]=sizes.at(0)/this->physicalConstants->getAstronomicalUnit();
        sizes[1]=sizes.at(1)/this->physicalConstants->getAstronomicalUnit();
        sizes[2]=sizes.at(2)/this->physicalConstants->getAstronomicalUnit();
        sizes[3]=sizes.at(3)/this->physicalConstants->getAstronomicalUnit();
        if(!localObserver){
            in<<"zoomau"<<endl;
            in<<QString::number(sizes.at(0))<<endl;
            in<<QString::number(sizes.at(1))<<endl;
            in<<QString::number(sizes.at(2))<<endl;
            in<<QString::number(sizes.at(3))<<endl;
            in<<"pointau"<<endl;
            in<<"0.0"<<endl;
            in<<"0.0"<<endl;
            in<<"0.0"<<endl;
        }else{
            points[0]=points.at(0)/this->physicalConstants->getAstronomicalUnit();
            points[1]=points.at(1)/this->physicalConstants->getAstronomicalUnit();
            points[2]=points.at(2)/this->physicalConstants->getAstronomicalUnit();

            obs[0]=obs.at(0)/this->physicalConstants->getAstronomicalUnit();
            obs[1]=obs.at(1)/this->physicalConstants->getAstronomicalUnit();
            obs[2]=obs.at(2)/this->physicalConstants->getAstronomicalUnit();
            in<<"locobsau"<<endl;
            in<<QString::number(obs.at(0))<<endl;
            in<<QString::number(obs.at(1))<<endl;
            in<<QString::number(obs.at(2))<<endl;
            in<<"pointau"<<endl;
            in<<QString::number(points.at(0))<<endl;
            in<<QString::number(points.at(1))<<endl;
            in<<QString::number(points.at(2))<<endl;
        }
    }else if(setup->getImageUnit()==SetUp::AU){
        if(this->imageControllerWidget->isZoomActive()){
            sizes[0]=sizes.at(0)/this->physicalConstants->getAstronomicalUnit();
            sizes[1]=sizes.at(1)/this->physicalConstants->getAstronomicalUnit();
            sizes[2]=sizes.at(2)/this->physicalConstants->getAstronomicalUnit();
            sizes[3]=sizes.at(3)/this->physicalConstants->getAstronomicalUnit();
        }
        if(!localObserver){
            in<<"zoomau"<<endl;
            in<<QString::number(sizes.at(0))<<endl;
            in<<QString::number(sizes.at(1))<<endl;
            in<<QString::number(sizes.at(2))<<endl;
            in<<QString::number(sizes.at(3))<<endl;
            in<<"pointau"<<endl;
            in<<"0.0"<<endl;
            in<<"0.0"<<endl;
            in<<"0.0"<<endl;
        }else{
            in<<"locobsau"<<endl;
            in<<QString::number(obs.at(0))<<endl;
            in<<QString::number(obs.at(1))<<endl;
            in<<QString::number(obs.at(2))<<endl;
            in<<"pointau"<<endl;
            in<<QString::number(points.at(0))<<endl;
            in<<QString::number(points.at(1))<<endl;
            in<<QString::number(points.at(2))<<endl;
        }
    }else if(setup->getImageUnit()==SetUp::PC){
        if(this->imageControllerWidget->isZoomActive()){
            sizes[0]=sizes.at(0)/this->physicalConstants->getParsec();
            sizes[1]=sizes.at(1)/this->physicalConstants->getParsec();
            sizes[2]=sizes.at(2)/this->physicalConstants->getParsec();
            sizes[3]=sizes.at(3)/this->physicalConstants->getParsec();
        }
        if(!localObserver){
            in<<"zoompc"<<endl;
            in<<QString::number(sizes.at(0))<<endl;
            in<<QString::number(sizes.at(1))<<endl;
            in<<QString::number(sizes.at(2))<<endl;
            in<<QString::number(sizes.at(3))<<endl;
            in<<"pointpc"<<endl;
            in<<"0.0"<<endl;
            in<<"0.0"<<endl;
            in<<"0.0"<<endl;
        }else{
            in<<"locobspc"<<endl;
            in<<QString::number(obs.at(0))<<endl;
            in<<QString::number(obs.at(1))<<endl;
            in<<QString::number(obs.at(2))<<endl;
            in<<"pointpc"<<endl;
            in<<QString::number(points.at(0))<<endl;
            in<<QString::number(points.at(1))<<endl;
            in<<QString::number(points.at(2))<<endl;
        }
    }
    if(this->imageControllerWidget->isPreview() ||localObserver){
        in<<"nofluxcons"<<endl;
    }else{
        in<<"fluxcons"<<endl;
    }
    if(this->imageControllerWidget->isStar()){
        in<<"inclstar"<<endl;
    }else{
        in<<"nostar"<<endl;
    }
    if(this->imageControllerWidget->isSecondOrder()&&!localObserver){
        in<<"secondorder"<<endl;
    }
    if(this->imageControllerWidget->isTau()){
        in<<"tracetau"<<endl;
    }else{
        in<<"tracenormal"<<endl;
    }
    if(this->imageControllerWidget->isDopplerCatchingAvailable()){
        if(this->imageControllerWidget->isDopplerCatchingActive())
            in<<"doppcatch"<<endl;
    }
    if(this->imageControllerWidget->isLineModus() && !this->imageControllerWidget->isRGBModus()){
        in<<"inclline"<<endl;
        if(this->imageControllerWidget->isDoLine()){
            in<<"iline"<<endl;
            QString iline=QString::number(this->imageControllerWidget->getILine());
            in<<iline.toStdString().c_str()<<endl;

            in<<"imolspec"<<endl;
            QString imol=QString::number(this->imageControllerWidget->getIMolecule());
            in<<imol.toStdString().c_str()<<endl;


            QString velocity=QString::number(this->imageControllerWidget->getVelocity());
            in<<"vkms"<<endl;
            in<<velocity.toStdString().c_str()<<endl;
        }else{
            if(this->imageControllerWidget->getsendCommandLineList()){
                in<<"linelist"<<endl;
                this->imageControllerWidget->setSendCommandLineList(false);
            }
        }
    }
    if(this->imageControllerWidget->isUsertransferOn()){
        QMapIterator<QString, QString> mapIterator( this->imageControllerWidget->getAppendCommand());
        while(mapIterator.hasNext()){
            mapIterator.next();
            in<<mapIterator.key()<<endl;
            in<<mapIterator.value()<<endl;
        }

    }

    in<<"enter"<<endl;
    in.flush();
    process->waitForBytesWritten();

}


QString MakeImage::getImageCommand(){
    bool localObserver=this->imageControllerWidget->isLocalModus();
    if(imageControllerWidget->isRGBModus()) imageControllerWidget->writeLambdasToCamerWavelengthFilename();
    QString command="image";
    QString numberOfPixels=QString::number(this->imageControllerWidget->getNumberOfPixels());
    command=command+" npix "+numberOfPixels;
    QString lambda=QString::number(this->imageControllerWidget->getLambdaRed());
    if(!imageControllerWidget->isRGBModus()&& !this->imageControllerWidget->isDoLine())
        command=command+" lambda "+lambda;
    QString posang=QString::number(this->imageControllerWidget->getPositionAngle());
    command=command+" posang "+posang;
    QVector<double> points;
    QVector<double> obs;
    if(!localObserver){
        QString incl=QString::number(this->imageControllerWidget->getInclination());
        command=command+" incl "+incl;
        QString phi=QString::number(this->imageControllerWidget->getPhi());
        command=command+" phi "+phi;
    }else{
        points.resize(3);
        obs.resize(3);
        points[0]=this->imageControllerWidget->getPointXLineEditValue();
        points[1]=this->imageControllerWidget->getPointYLineEditValue();
        points[2]=this->imageControllerWidget->getPointZLineEditValue();
        obs[0]=this->imageControllerWidget->getObserverposXLineEditValue();
        obs[1]=this->imageControllerWidget->getObserverposYLineEditValue();
        obs[2]=this->imageControllerWidget->getObserverposZLineEditValue();
        if(this->imageControllerWidget->isRelativeObserver()){
            obs[0]+=points[0];
            obs[1]+=points[1];
            obs[2]+=points[2];
        }
        command=command+" sizeradian "+QString::number(this->imageControllerWidget->getViewangleSpinBoxValue());
    }
    QVector<double> sizes;
    double size;
    if(!this->imageControllerWidget->isZoomActive()){
        size=this->imageControllerWidget->getSize()/2.;
        sizes.append(size*-1.);
        sizes.append(size);
        sizes.append(size*-1.);
        sizes.append(size);
    }else{
        sizes=this->imageControllerWidget->getZoomBox();
    }
    if(setup->getImageUnit()==SetUp::CM){
        sizes[0]=sizes.at(0)/this->physicalConstants->getAstronomicalUnit();
        sizes[1]=sizes.at(1)/this->physicalConstants->getAstronomicalUnit();
        sizes[2]=sizes.at(2)/this->physicalConstants->getAstronomicalUnit();
        sizes[3]=sizes.at(3)/this->physicalConstants->getAstronomicalUnit();
        if(!localObserver){
            command=command+" zoomau "+QString::number(sizes.at(0))+" "+
                    QString::number(sizes.at(1))+" "+QString::number(sizes.at(2))+" "+QString::number(sizes.at(3));
            command=command+" pointau 0.0  0.0  0.0";
        }else{
            points[0]=points.at(0)/this->physicalConstants->getAstronomicalUnit();
            points[1]=points.at(1)/this->physicalConstants->getAstronomicalUnit();
            points[2]=points.at(2)/this->physicalConstants->getAstronomicalUnit();

            obs[0]=obs.at(0)/this->physicalConstants->getAstronomicalUnit();
            obs[1]=obs.at(1)/this->physicalConstants->getAstronomicalUnit();
            obs[2]=obs.at(2)/this->physicalConstants->getAstronomicalUnit();
            command=command+" pointau " +QString::number(points.at(0))+" "+QString::number(points.at(1))+" "+QString::number(points.at(2));
            command=command+" locobsau " +QString::number(obs.at(0))+" "+QString::number(obs.at(1))+" "+QString::number(obs.at(2));
        }
    }else if(setup->getImageUnit()==SetUp::AU){
        if(this->imageControllerWidget->isZoomActive()){
            sizes[0]=sizes.at(0)/this->physicalConstants->getAstronomicalUnit();
            sizes[1]=sizes.at(1)/this->physicalConstants->getAstronomicalUnit();
            sizes[2]=sizes.at(2)/this->physicalConstants->getAstronomicalUnit();
            sizes[3]=sizes.at(3)/this->physicalConstants->getAstronomicalUnit();
        }
        if(!localObserver){
            command=command+" zoomau "+QString::number(sizes.at(0))+" "+
                    QString::number(sizes.at(1))+" "+QString::number(sizes.at(2))+" "+QString::number(sizes.at(3));
            command=command+" pointau 0.0  0.0  0.0";
        }else{
            command=command+" pointau " +QString::number(points.at(0))+" "+QString::number(points.at(1))+" "+QString::number(points.at(2));
            command=command+" locobsau " +QString::number(obs.at(0))+" "+QString::number(obs.at(1))+" "+QString::number(obs.at(2));
        }
    }else if(setup->getImageUnit()==SetUp::PC){
        if(this->imageControllerWidget->isZoomActive()){
            sizes[0]=sizes.at(0)/this->physicalConstants->getParsec();
            sizes[1]=sizes.at(1)/this->physicalConstants->getParsec();
            sizes[2]=sizes.at(2)/this->physicalConstants->getParsec();
            sizes[3]=sizes.at(3)/this->physicalConstants->getParsec();
        }

        if(!localObserver){
            command=command+" zoompc "+QString::number(sizes.at(0))+" "+
                    QString::number(sizes.at(1))+" "+QString::number(sizes.at(2))+" "+QString::number(sizes.at(3));
            command=command+" pointpc 0.0  0.0  0.0";
        }else{
            command=command+" pointpc " +QString::number(points.at(0))+" "+QString::number(points.at(1))+" "+QString::number(points.at(2));
            command=command+" locobspc " +QString::number(obs.at(0))+" "+QString::number(obs.at(1))+" "+QString::number(obs.at(2));
        }
    }
    if(this->imageControllerWidget->isPreview()||localObserver){
        command=command+" nofluxcons";
    }else{
        command=command+" fluxcons";
    }
    if(this->imageControllerWidget->isStar()){
        command=command+" inclstar";
    }else{
        command=command+" nostar";
    }
    if(this->imageControllerWidget->isSecondOrder()&&!localObserver)command=command+" secondorder ";
    if(this->imageControllerWidget->isTau()){
        command=command+" tracetau";
    }else{
        command=command+" tracenormal";
    }
    if(imageControllerWidget->isRGBModus()){
        command=command+" loadlambda";
    }
    if(this->imageControllerWidget->isDopplerCatchingAvailable()){
        if(this->imageControllerWidget->isDopplerCatchingActive())
            command=command+" doppcatch";
    }

    if(this->imageControllerWidget->isLineModus()&& !this->imageControllerWidget->isRGBModus()){
        command=command+" inclline";
        if(this->imageControllerWidget->isDoLine()){
            QString iline=QString::number(this->imageControllerWidget->getILine());
            command=command+" iline "+iline;

            QString imol=QString::number(this->imageControllerWidget->getIMolecule());
            command=command+" imolspec "+imol;

            QString velocity=QString::number(this->imageControllerWidget->getVelocity());
            command=command+" vkms "+velocity;

        }else{
            if(this->imageControllerWidget->getsendCommandLineList()){
                command=command+" linelist";
                if(!this->imageControllerWidget->isPipe())this->imageControllerWidget->setSendCommandLineList(false);
            }
        }
    }
    if(this->imageControllerWidget->isUsertransferOn()){
        QMapIterator<QString, QString> mapIterator( this->imageControllerWidget->getAppendCommand());
        while(mapIterator.hasNext()){
            mapIterator.next();
            command=command+" "+mapIterator.key()+" "+mapIterator.value();
        }

    }
    return command;
}
