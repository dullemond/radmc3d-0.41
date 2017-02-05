#ifndef POLAR_H
#define POLAR_H

#if QT_VERSION < 0x050000
      #include <QtGui>
#else
       #include <QtWidgets>
#endif

class Polar : public QObject
{
    Q_OBJECT
public:
    explicit Polar(QObject *parent = 0);
    static Polar *getInstance(QObject *parent=0);
    void setRi(QVector<double> vector);
    QVector<double>  getRi(){return this->ri;}
    void setThetai(QVector<double> vector);
    QVector<double>  getThetai(){return this->thetai;}
    void setPhii(QVector<double> vector);
    QVector<double>  getPhii(){return this->phii;}
    void setR(QVector<double> vector);
    QVector<double>  getR(){return this->r;}
    void setTheta(QVector<double> vector);
    QVector<double>  getTheta(){return this->theta;}
    void setPhi(QVector<double> vector);
    QVector<double>  getPhi(){return this->phi;}
    void setNumberOfR(int number){this->numberOfR=number;}
    void setNumberOfTheta(int number){this->numberOfTheta=number;}
    void setNumberOfPhi(int number){this->numberOfPhi=number;}
    int getNumberOfR(){return this->numberOfR;}
    int getNumberOfPhi(){return this->numberOfPhi;}
    int getNumberOfTheta(){return this->numberOfTheta;}
private:
    static Polar* instance;
    int numberOfR;
    int numberOfPhi;
    int numberOfTheta;
    QVector<double> thetai;
    QVector<double> phii;
    QVector<double> ri;
    QVector<double> theta;
    QVector<double> phi;
    QVector<double> r;
signals:
    
public slots:
    
};

#endif // POLAR_H
