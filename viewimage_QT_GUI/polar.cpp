#include "polar.h"
Polar *Polar::instance=0;
Polar::Polar(QObject *parent) :
    QObject(parent)
{
}
Polar *Polar::getInstance(QObject *parent){
    if(!instance){
        instance=new Polar(parent);
    }
return instance;
}
void Polar::setPhii(QVector<double> vector){
    this->phii.clear();
    this->phii=vector;
}
void Polar::setRi(QVector<double> vector){
    this->ri.clear();
    this->ri=vector;
}
void Polar::setThetai(QVector<double> vector){
    this->thetai.clear();
    this->thetai=vector;
}
void Polar::setPhi(QVector<double> vector){
    this->phi.clear();
    this->phi=vector;
}
void Polar::setR(QVector<double> vector){
    this->r.clear();
    this->r=vector;
}
void Polar::setTheta(QVector<double> vector){
    this->theta.clear();
    this->theta=vector;
}
