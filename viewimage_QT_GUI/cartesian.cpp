#include "cartesian.h"
Cartesian *Cartesian::instance=0;

Cartesian::Cartesian(QObject *parent) :
    QObject(parent)
{

}

Cartesian *Cartesian::getInstance(QObject *parent){
    if(!instance){
        instance=new Cartesian(parent);
    }
    return instance;
}
void Cartesian::setX(QVector<double> values){
    this->x.clear();
    x=values;

}
void Cartesian::setY(QVector<double> values){
    this->y.clear();
   y=values;
}
void Cartesian::setZ(QVector<double> values){
    this->z.clear();
    z=values;
}
