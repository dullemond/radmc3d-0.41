#ifndef CARTESIAN_H
#define CARTESIAN_H

#if QT_VERSION < 0x050000
      #include <QtGui>
#else
       #include <QtWidgets>
#endif

class Cartesian : public QObject
{
    Q_OBJECT
public:
    explicit Cartesian(QObject *parent = 0);
    void setX(QVector<double> values);
    void setY(QVector<double> values);
    void setZ(QVector<double> values);
    static Cartesian *getInstance(QObject *parent=0);
    QVector<double> getXi(){return this->xi;}
    QVector<double> getYi(){return this->yi;}
    QVector<double> getZi(){return this->zi;}
    void setXi( QVector<double> vec ){this->xi=vec ;}
    void setYi( QVector<double> vec ){this->yi=vec ;}
    void setZi( QVector<double> vec ){this->zi=vec ;}
private:
    QVector<double> xi,yi,zi;
    static Cartesian *instance;
    QVector<double> x,y,z;
    
signals:
    
public slots:
    
};

#endif // CARTESIAN_H
