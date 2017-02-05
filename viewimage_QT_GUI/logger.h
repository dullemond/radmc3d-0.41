#ifndef LOGGER_H
#define LOGGER_H
#include <QtCore>
#include <QObject>
#include <QStringListModel>
#include <QDesktopServices>
class Logger:public QObject
{
 Q_OBJECT
public:
    explicit  Logger(QObject *parent=0);
    static Logger *getInstance(QObject *parent=0);
    void setLogFilename(QString filename);
    QString getLogFilename();
    void setProtocolState(bool state);
    bool isProtocolOn();
    void redirectToLogFile(bool redirectToLogFile);
    bool createLogFile();
    void writeToLogFile(QString logText);
    QStringListModel *getMessageModel();
    void stop(int exitCode=0);
private:
    QStringListModel *messageModel;
    static Logger *instance;
    bool protocolOn;
    bool redirectToFile;
    QString logFilename;


};

#endif // LOGGER_H
