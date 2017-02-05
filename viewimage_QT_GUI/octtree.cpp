#include "octtree.h"
Octtree * Octtree::instance=0;

Octtree *Octtree::getInstance(){
    if(!instance){
        instance=new Octtree;
    }
    return instance;
}
void Octtree::setOcttree(QVector<int> tree){
    this->tree.clear();
    this->tree=tree;
}

QVector<int>  Octtree::getOcttree(){
        return this->tree;
}

Octtree::Octtree(QObject *parent) :
    QObject(parent)
{
}
