#ifndef OCTTREE_H
#define OCTTREE_H
#if QT_VERSION < 0x050000
      #include <QtGui>
#else
       #include <QtWidgets>
#endif

class Octtree : public QObject
{
    Q_OBJECT
public:
    explicit Octtree(QObject *parent = 0);
    static Octtree* getInstance();
    void setOcttree(QVector<int>);
    QVector<int> getOcttree();
private:
    QVector<int> tree;
    static Octtree *instance;
    
signals:
    
public slots:
    
};

#endif // OCTTREE_H
