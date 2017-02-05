#include "makespectrum.h"

MakeSpectrum::MakeSpectrum(ImageControllerWidget *imageControllerWidget,QObject *parent) :
    QObject(parent)
{
    this->imageControllerWidget=imageControllerWidget;
    this->setup=SetUp::getInstance();
    this->physicalConstants=PhysicalConstants::getInstance();
    this->logger=Logger::getInstance();
}

void MakeSpectrum::sendMakeSpectrumCommandToRadmc3d(QProcess *process){
    QVector<double> sizes;
    double size;
    QString incl=QString::number(this->imageControllerWidget->getInclination());
    QString phi=QString::number(this->imageControllerWidget->getPhi());
    if(!this->imageControllerWidget->isZoomActive()){
        size=this->imageControllerWidget->getSize()/2.;
        sizes.append(size*-1.);
        sizes.append(size);
        sizes.append(size*-1.);
        sizes.append(size);
    }else{
        sizes=this->imageControllerWidget->getZoomBox();
    }
    QTextStream in(process);
    in<<"sed"<<endl;
    in<<"incl"<<endl;
    in<<incl.toStdString().c_str()<<endl;
    in<<"phi"<<endl;
    in<<phi.toStdString().c_str()<<endl;
    if(setup->getImageUnit()==SetUp::CM){
        sizes[0]=sizes.at(0)/this->physicalConstants->getAstronomicalUnit();
        sizes[1]=sizes.at(1)/this->physicalConstants->getAstronomicalUnit();
        sizes[2]=sizes.at(2)/this->physicalConstants->getAstronomicalUnit();
        sizes[3]=sizes.at(3)/this->physicalConstants->getAstronomicalUnit();
        in<<"zoomau"<<endl;
        in<<QString::number(sizes.at(0))<<endl;
        in<<QString::number(sizes.at(1))<<endl;
        in<<QString::number(sizes.at(2))<<endl;
        in<<QString::number(sizes.at(3))<<endl;
        in<<"pointau"<<endl;
        in<<"0.0"<<endl;
        in<<"0.0"<<endl;
        in<<"0.0"<<endl;

    }else if(setup->getImageUnit()==SetUp::AU){
        if(this->imageControllerWidget->isZoomActive()){
            sizes[0]=sizes.at(0)/this->physicalConstants->getAstronomicalUnit();
            sizes[1]=sizes.at(1)/this->physicalConstants->getAstronomicalUnit();
            sizes[2]=sizes.at(2)/this->physicalConstants->getAstronomicalUnit();
            sizes[3]=sizes.at(3)/this->physicalConstants->getAstronomicalUnit();
        }
        in<<"zoomau"<<endl;
        in<<QString::number(sizes.at(0))<<endl;
        in<<QString::number(sizes.at(1))<<endl;
        in<<QString::number(sizes.at(2))<<endl;
        in<<QString::number(sizes.at(3))<<endl;
        in<<"pointau"<<endl;
        in<<"0.0"<<endl;
        in<<"0.0"<<endl;
        in<<"0.0"<<endl;

    }else if(setup->getImageUnit()==SetUp::PC){
        if(this->imageControllerWidget->isZoomActive()){
            sizes[0]=sizes.at(0)/this->physicalConstants->getParsec();
            sizes[1]=sizes.at(1)/this->physicalConstants->getParsec();
            sizes[2]=sizes.at(2)/this->physicalConstants->getParsec();
            sizes[3]=sizes.at(3)/this->physicalConstants->getParsec();
        }
        in<<"zoompc"<<endl;
        in<<QString::number(sizes.at(0))<<endl;
        in<<QString::number(sizes.at(1))<<endl;
        in<<QString::number(sizes.at(2))<<endl;
        in<<QString::number(sizes.at(3))<<endl;
        in<<"pointpc"<<endl;
        in<<"0.0"<<endl;
        in<<"0.0"<<endl;
        in<<"0.0"<<endl;

    }

    if(this->imageControllerWidget->isStar()){
        in<<"inclstar"<<endl;
    }else{
        in<<"nostar"<<endl;
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

void MakeSpectrum::sendWriteSpectrumCommandToRadmc3d(QProcess *process){

    QTextStream in(process);
    in<<"writespec"<<endl;
    in.flush();

    process->waitForBytesWritten();

}
QString MakeSpectrum::getSpectrumCommand(){
    QString command="sed";
    QString incl=QString::number(this->imageControllerWidget->getInclination());
    command=command+" incl "+incl;
    QString phi=QString::number(this->imageControllerWidget->getPhi());
    command=command+" phi "+phi;
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
        command=command+" zoomau "+QString::number(sizes.at(0))+" "+
                QString::number(sizes.at(1))+" "+QString::number(sizes.at(2))+" "+QString::number(sizes.at(3));
        command=command+" pointau 0.0  0.0  0.0";

    }else if(setup->getImageUnit()==SetUp::AU){
        if(this->imageControllerWidget->isZoomActive()){
            sizes[0]=sizes.at(0)/this->physicalConstants->getAstronomicalUnit();
            sizes[1]=sizes.at(1)/this->physicalConstants->getAstronomicalUnit();
            sizes[2]=sizes.at(2)/this->physicalConstants->getAstronomicalUnit();
            sizes[3]=sizes.at(3)/this->physicalConstants->getAstronomicalUnit();
        }

        command=command+" zoomau "+QString::number(sizes.at(0))+" "+
                QString::number(sizes.at(1))+" "+QString::number(sizes.at(2))+" "+QString::number(sizes.at(3));
        command=command+" pointau 0.0  0.0  0.0";

    }else if(setup->getImageUnit()==SetUp::PC){
        if(this->imageControllerWidget->isZoomActive()){
            sizes[0]=sizes.at(0)/this->physicalConstants->getParsec();
            sizes[1]=sizes.at(1)/this->physicalConstants->getParsec();
            sizes[2]=sizes.at(2)/this->physicalConstants->getParsec();
            sizes[3]=sizes.at(3)/this->physicalConstants->getParsec();
        }
        command=command+" zoompc "+QString::number(sizes.at(0))+" "+
                QString::number(sizes.at(1))+" "+QString::number(sizes.at(2))+" "+QString::number(sizes.at(3));
        command=command+" pointpc 0.0  0.0  0.0";

    }
    if(this->imageControllerWidget->isStar()){
        command=command+" inclstar";
    }else{
        command=command+" nostar";
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
