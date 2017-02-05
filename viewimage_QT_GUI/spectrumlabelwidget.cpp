#include "spectrumlabelwidget.h"

SpectrumLabelWidget::SpectrumLabelWidget(int coordDistanceToImage,QWidget *parent) :
    QLabel(parent)
{
    this->mouseMoved=false;
    this->mousePressed=false;
    this->zoomActive=false;
    this->zoomOutButton=new QToolButton(this);
    this->zoomOutButton->setVisible(false);
    this->zoomOutButton->setIcon(QIcon(trUtf8(":/images/zoomOut.ico")));
    this->connect(this->zoomOutButton,SIGNAL(clicked()),this,SIGNAL(removeZoom()));
    this->zoomOutButton->setShortcut(QKeySequence(Qt::Key_Minus));
    this->zoomOutButton->setWhatsThis(trUtf8("<html><font color=\"blue\">This button zooms out the spectrum(shortcut: %1 ).</html>").arg(this->zoomOutButton->shortcut().toString()));
    this->zoomOutButton->adjustSize();
    this->zoomOutButton->move(0,0);
    this->coordDistanceToImage=coordDistanceToImage;
    this->unit=Spectrum::MICRON;
    this->fluxunit=Spectrum::NUFNU;
    this->xLog=true;
    this->yLog=true;
    setMouseTracking(true);
}
void SpectrumLabelWidget::setZoomOutVisible(bool visible){
    this->zoomOutButton->setVisible(visible);
}

void SpectrumLabelWidget::mouseMoveEvent(QMouseEvent *event){
    if(event->type()==QMouseEvent::MouseMove){
        if(!mousePressed){
            if(this->pixmap()!=NULL){
                if(this->spectrum!=NULL){
                    QPoint point=event->pos();
                    QSize size=this->pixmap()->size();
                    if(point.x()>this->coordDistanceToImage && point.y()>this->coordDistanceToImage&&
                            point.x()<size.width()-this->coordDistanceToImage && point.y()<size.height()-this->coordDistanceToImage){
                        point.setX(point.x()-this->coordDistanceToImage);
                        point.setY(point.y()-this->coordDistanceToImage);
                        if(point.x()>0 && point.y()>0){
                            size.setHeight(size.height()-2*coordDistanceToImage);
                            size.setWidth(size.width()-2*coordDistanceToImage);
                            double x=point.x()*((this->spectrum->getXEnd()-this->spectrum->getXStart()*1.)/(size.width()))+this->spectrum->getXStart();
                            if(this->isXLogOn())
                                x=pow(10.,x);
                            for (int var = 0; var < this->spectrum->getNumberOfFrequencies(); ++var) {
                                double xat=this->spectrum->getXcoord().at(var);
                                if(xat==x){
                                    int pY=(this->spectrum->getYcoord().at(x)-this->spectrum->getYStart())
                                            /(this->spectrum->getYEnd()-this->spectrum->getYStart())*size.height();
                                    if((pY>0)&&(pY<(pixmap()->height()-2*coordDistanceToImage)))
                                        emit drawCircle(point.x()+this->coordDistanceToImage,pY+this->coordDistanceToImage);
                                    emit  pixelValueChange(trUtf8("X=%1 Y=%2").arg(QString::number(xat),QString::number(this->spectrum->getYcoord().at(x))));
                                    break;
                                }else if(xat>x){
                                    if(var>0){
                                        double xDiff=xat-this->spectrum->getXcoord().at(var-1);
                                        double yDiff=this->spectrum->getYcoord().at(var)-this->spectrum->getYcoord().at(var-1);
                                        double m=yDiff/xDiff;
                                        double b=this->spectrum->getYcoord().at(var)-m*xat;
                                        double y=m*x+b;
                                        int pY;
                                        if(this->yLog)
                                            pY=(log10(y)-this->spectrum->getYStart())
                                                    /(this->spectrum->getYEnd()-this->spectrum->getYStart())*size.height();
                                        else
                                            pY=((y)-this->spectrum->getYStart())
                                                    /(this->spectrum->getYEnd()-this->spectrum->getYStart())*size.height();
                                        if((pY>0)&&(pY<(pixmap()->height()-2*coordDistanceToImage)))
                                            emit drawCircle(point.x()+this->coordDistanceToImage,pY+this->coordDistanceToImage);
                                        emit  pixelValueChange(trUtf8("X=%1 Y=%2").arg(QString::number(x),QString::number(y)));
                                        break;
                                    }
                                }
                            }


                        }
                    }
                }
            }
        }else{
            this->last=event->pos();
            mouseMoved=true;
            QSize size=this->pixmap()->size();
            int distance=last.y()-first.y();
            int yLength;
            if(distance>0){
                yLength=last.x()-first.x();
                if(yLength<0)yLength=yLength*-1;
            }else{
                yLength=last.x()-first.x();
                if(yLength>0)yLength=yLength*-1;
            }
            last.setY(first.y()+yLength);
            distance=this->coordDistanceToImage;
            if(last.x()>distance && last.y()>distance&&
                    last.x()<size.width()-distance && last.y()<size.height()-distance){
                emit drawRect(first,last);
                emit pixelValueChange("");
            }
        }
    }
}
void SpectrumLabelWidget::mousePressEvent(QMouseEvent *event){
    if(event->button()==Qt::RightButton && event->type()==QMouseEvent::MouseButtonPress){
        QString no="";
        if(!this->signalsBlocked()){
            QMenu *menu=new QMenu(this);


            QAction *saveAction=new QAction(menu);
            saveAction->setObjectName("save");
            saveAction->setIcon(QIcon(trUtf8(":/images/%1.png").arg(saveAction->objectName())));
            saveAction->setIconVisibleInMenu(true);
            saveAction->setText(saveAction->objectName());


            QAction *printAction=new QAction(menu);
            printAction->setObjectName("print");
            printAction->setIcon(QIcon(trUtf8(":/images/%1.png").arg(printAction->objectName())));
            printAction->setIconVisibleInMenu(true);
            printAction->setText(printAction->objectName());

            menu->addAction(saveAction);
            menu->addAction(printAction);


            connect(saveAction,SIGNAL(triggered()),this,SLOT(savePixmap()));
            connect(printAction,SIGNAL(triggered()),this,SLOT(printPixmap()));
            menu->exec(event->globalPos());
        }
    }else if(event->button()==Qt::LeftButton && event->type()==QMouseEvent::MouseButtonPress){
        first=event->pos();
        QSize size=this->pixmap()->size();
        int distance=this->coordDistanceToImage;
        if(first.x()>distance && first.y()>distance&&
                first.x()<size.width()-distance && first.y()<size.height()-distance){
            this->mousePressed=true;
            this->mouseMoved=false;
        }
    }
}
void SpectrumLabelWidget::mouseReleaseEvent(QMouseEvent *event){
    if(event!=NULL){
        if(this->mousePressed){
            this->mousePressed=false;
            if(mouseMoved){
                this->mouseMoved=false;
                emit calculateZoom(first,last);
            }
        }
    }
}

void SpectrumLabelWidget::printPixmap(){
    QPrinter printer;
    QPrintDialog dlg( &printer );
    if( dlg.exec() == QDialog::Accepted )
    {
        QPainter painter( &printer );
        QPixmap pixmap=this->spectrum->getPixmap();
        pixmap=pixmap.scaled(printer.pageRect().width(), printer.pageRect().height(), Qt::KeepAspectRatio);
        painter.drawPixmap(0,0,pixmap);
        painter.end();

    }
    dlg.close();
}

void SpectrumLabelWidget::savePixmap(){
    emit pixelValueChange("*.BMP *.png *.jpg *.jpeg *.xbm *.xpm are supported");
    QString filename = QFileDialog::getSaveFileName(this,"save file dialog",QDir::currentPath(),
                                                    "Images (*.BMP *.png *.jpg *.jpeg *.xbm *.xpm)");
    this->spectrum->getPixmap().save(filename);
}
