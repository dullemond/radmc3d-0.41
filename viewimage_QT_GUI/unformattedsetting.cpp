#include "unformattedsetting.h"

UnformattedSetting *UnformattedSetting::instance=0;

UnformattedSetting::UnformattedSetting(QObject *parent) :
    QObject(parent)
{
    if(QSysInfo::ByteOrder==QSysInfo::BigEndian)
        this->setByteOrder(QDataStream::BigEndian);
    else
        this->setByteOrder(QDataStream::LittleEndian);
    this->setRecordLengthSize(4);
}

UnformattedSetting  *UnformattedSetting::getInstance(QObject *parent){
    if(!instance){
        instance=new UnformattedSetting(parent);
    }
    return instance;
}
