#include "imagelabelwidget.h"

ImageLabelWidget::ImageLabelWidget(int coordDistanceToImag,QWidget *parent) :
    QLabel(parent)
{

    this->mouseMoved=false;
    this->plot=true;
    this->setCoordDistanceToImage(coordDistanceToImag);
    setMouseTracking(true);

    this->zoomActive=true;

}
void ImageLabelWidget::setPhiAndInclination(double phi,double inclination){
    this->inclination=inclination;
    this->phi=phi;
}

void ImageLabelWidget::mouseMoveEvent(QMouseEvent *event){
    if(event->type()==QMouseEvent::MouseMove){
        if(!mousePressed){
            if(this->pixmap()!=NULL){
                if(this->image!=NULL){
                    QPoint point=event->pos();
                    QSize size=this->pixmap()->size();
                    if(point.x()>this->coordDistanceToImage && point.y()>this->coordDistanceToImage&&
                            point.x()<size.width()-this->coordDistanceToImage && point.y()<size.height()-this->coordDistanceToImage){
                        point.setX(point.x()-this->coordDistanceToImage);
                        point.setY(point.y()-this->coordDistanceToImage);
                        if(point.x()>0 && point.y()>0){
                            size.setHeight(size.height()-2*coordDistanceToImage);
                            size.setWidth(size.width()-2*coordDistanceToImage);
                            int y=point.y()*(this->image->getNumberOfPixelY()*1./size.height());
                            int x=point.x()*(this->image->getNumberOfPixelX()*1./size.width());
                            if(this->image->getNumberOfImages()==1)
                                emit  pixelValueChange(QString::number(this->image->getImageAt(0).at(y)
                                                                       .at(x)));
                            else if(this->image->getNumberOfImages()==3){
                                QString str="Red = "+QString::number(this->image->getImageAt(0).at(y)
                                                                     .at(x));
                                str=str+"   Green = "+QString::number(this->image->getImageAt(1).at(y)
                                                                      .at(x));
                                str=str+"   Blue = "+QString::number(this->image->getImageAt(2).at(y)
                                                                     .at(x));
                                emit pixelValueChange(str);

                            }
                        }
                    }
                }
            }
        }else{
            if(!this->image->isLocalObserver()){
                this->last=event->pos();
                mouseMoved=true;
                QSize size=this->pixmap()->size();
                if(zoomActive){
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
                    if(!this->isPlot())distance=0;
                    if(last.x()>distance && last.y()>distance&&
                            last.x()<size.width()-distance && last.y()<size.height()-distance){
                        emit drawRect(first,last);
                        emit pixelValueChange("");
                    }
                }
            }
        }
    }
}
void ImageLabelWidget::mousePressEvent(QMouseEvent *event){
    first=event->pos();
    if(event->button()==Qt::LeftButton && event->type()==QMouseEvent::MouseButtonPress){
        if(!this->image->isLocalObserver()){
            QSize size=this->pixmap()->size();
            int distance=this->coordDistanceToImage;
            if(!this->isPlot())distance=0;
            if(first.x()>distance && first.y()>distance&&
                    first.x()<size.width()-distance && first.y()<size.height()-distance){
                this->mousePressed=true;
                emit removeRect();
                mouseMoved=false;
            }
        }
    }else if(event->button()==Qt::RightButton && event->type()==QMouseEvent::MouseButtonPress){
        QString no="";
        if(!zoomActive)no="no";
        if(!this->signalsBlocked()){
            QMenu *menu=new QMenu(this);
            if(!this->image->isLocalObserver()){
                zoomAction=new QAction(menu);
                zoomAction->setObjectName("zoomIn");
                zoomAction->setIcon(QIcon(QString(":/images/%2%1.png").arg(zoomAction->objectName(),no)));
                zoomAction->setIconVisibleInMenu(true);
                zoomAction->setText(zoomAction->objectName());
                zoomAction->setCheckable(true);
                zoomAction->setChecked(zoomActive);
                menu->addAction(zoomAction);
                no="";
                if(zoomActive)no="no";
                rotateAction=new QAction(menu);
                rotateAction->setObjectName("rotate");
                rotateAction->setIcon(QIcon(QString(":/images/%2%1.png").arg(rotateAction->objectName(),no)));
                rotateAction->setText(rotateAction->objectName());
                rotateAction->setCheckable(true);
                rotateAction->setChecked(!zoomActive);
                rotateAction->setIconVisibleInMenu(true);
                menu->addAction(rotateAction);
            }
            QAction *plotAction=new QAction(menu);
            plotAction->setObjectName("plot");
            plotAction->setCheckable(true);
            plotAction->setChecked(this->plot);
            no="";
            if(!plot)no="no";
            plotAction->setIcon(QIcon(trUtf8(":/images/%2%1.png").arg(plotAction->objectName(),no)));
            plotAction->setIconVisibleInMenu(true);
            plotAction->setText(plotAction->objectName());

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

            menu->addAction(plotAction);
            menu->addAction(saveAction);
            menu->addAction(printAction);

            if(!this->image->isLocalObserver())connect(rotateAction,SIGNAL(triggered(bool)),this,SLOT(rotateClicked(bool)));
            connect(plotAction,SIGNAL(triggered(bool)),this,SLOT(plotClicked(bool)));
            connect(saveAction,SIGNAL(triggered()),this,SLOT(savePixmap()));
            connect(printAction,SIGNAL(triggered()),this,SLOT(printPixmap()));
            if(!this->image->isLocalObserver()) connect(zoomAction,SIGNAL(triggered(bool)),this,SLOT(zoomClicked(bool)));
            menu->exec(event->globalPos());
        }
    }
}
void ImageLabelWidget::printPixmap(){
    emit removeRect();
    QPrinter printer;
    QPrintDialog dlg( &printer );
    if( dlg.exec() == QDialog::Accepted )
    {
        QPainter painter( &printer );
        QPixmap pixmap=this->pixmap()->copy();
        pixmap=pixmap.scaled(printer.pageRect().width(), printer.pageRect().height(), Qt::KeepAspectRatio);
        painter.drawPixmap(0,0,pixmap);
        painter.end();

    }
    dlg.close();
}

void ImageLabelWidget::savePixmap(){
    emit removeRect();
    emit pixelValueChange("*.BMP *.png *.jpg *.jpeg *.xbm *.xpm are supported");
    QString filename = QFileDialog::getSaveFileName(this,"save file dialog",QDir::currentPath(),
                                                    "Images (*.BMP *.png *.jpg *.jpeg *.xbm *.xpm)");
    this->pixmap()->save(filename);
}

void ImageLabelWidget::zoomClicked(bool checked){
    if(checked) {
        this->zoomActive=true;
    }else{
        this->zoomActive=false;
        emit removeRect();
    }

}
void ImageLabelWidget::rotateClicked(bool checked){
    if(checked){
        this->zoomActive=false;
        emit removeRect();
    }else{
        this->zoomActive=true;
    }
}
void ImageLabelWidget::plotClicked(bool checked){
    this->plot=checked;
    emit reloadImage();
}


void ImageLabelWidget::mouseReleaseEvent(QMouseEvent *event){
    if(this->mousePressed){
        if(!this->image->isLocalObserver()){
            if(zoomActive){
                this->mousePressed=false;
                if(mouseMoved){
                    mouseMoved=false;
                    emit calculateZoom(first,last);
                }
            }else{
                if(mouseMoved){
                    QSize size=this->pixmap()->size();
                    last=event->pos();
                    int xlenghth=last.x()-first.x();
                    int ylength=last.y()-first.y();
                    size.setWidth(size.width()-2*coordDistanceToImage);
                    size.setHeight(size.height()-2*coordDistanceToImage);
                    if(xlenghth!=0){
                        inclination=inclination-(((ylength*1./size.height()))*120.);
                    }
                    if(ylength!=0){
                        if(phi*inclination/180>0){
                            phi=phi+(((xlenghth*1.0/size.width()))*120.);
                        }else{
                            phi=phi-(((xlenghth*1.0/size.width()))*120.);
                        }
                    }
                    mouseMoved=false;
                    emit phiAndinclinactionValueChanged(phi,inclination);
                }
            }
        }
    }
}
