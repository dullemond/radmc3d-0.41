#ifndef UNFORMATTEDSETTING_H
#define UNFORMATTEDSETTING_H

#if QT_VERSION < 0x050000
      #include <QtGui>
#else
       #include <QtWidgets>
#endif

class UnformattedSetting : public QObject
{
    Q_OBJECT
public:
    explicit UnformattedSetting(QObject *parent = 0);
    void setByteOrder(QDataStream::ByteOrder order){this->byteOrder=order;}
    QDataStream::ByteOrder getByteOrder(){return this->byteOrder;}
    void setRecordLengthSize(qint32 size){this->recordLengthSize=size;}
    qint32 getRecordLengthSize(){return this->recordLengthSize;}
    static UnformattedSetting* getInstance(QObject *parent = 0);
    static UnformattedSetting* instance;
private:
    QDataStream::ByteOrder byteOrder;
    qint32 recordLengthSize;
signals:
    
public slots:
    
};

#endif // UNFORMATTEDSETTING_H
