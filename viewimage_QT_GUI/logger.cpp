#include "logger.h"

Logger* Logger::instance=0;

/**
*   @file logger.cpp
*   @brief  Primary function of the class Logger is protocol the RADMC3D-GUI. Every protocol-action should be called over this class
*   @author Farzin Sereshti
*   @version 1.0
*/


/**
 *  @brief This is constructor.
 *
 *  @details This constructor initializes the messageModel. The messageModel is a QStringListModel, that contains all program's messages, that will be generated in run-time.
 *  This QStringModel will be needed by listview to assure the Model-View's function.
 *  @attention Do not call this constructor directly. Instead of that, you should call getInstance
 *  @see getInstance
 */
Logger::Logger(QObject *parent):
    QObject(parent)
{
    this->protocolOn=false;
    this->redirectToFile=false;
      this->messageModel=new QStringListModel();
}

/**
 *
 */
void Logger::stop(int exitCode){
    this->redirectToLogFile(false);
    this->writeToLogFile(trUtf8("Program will end, see LOG-file exitCode"));
   QDesktopServices::openUrl(this->logFilename);
    QCoreApplication::exit(exitCode);
    exit(exitCode);
}

/**
 *  @brief This method gives a uniq instance of Logger class.
 *
 *
 */

Logger* Logger::getInstance(QObject *parent)
{
    if(!instance) {
        instance=new Logger(parent);
        instance->writeToLogFile(trUtf8("Logger is initialized now"));
    }
    return instance;
}
/**
 * @brief This method returns messageModel, that will be needed by gui's listview.
 *
 * @return messageModel.
 */
QStringListModel *Logger::getMessageModel(){
    return messageModel;
}

/**
 * @brief This method create a new Logfile.
 *
 * @details This method  create a new logfile, that is defiend by logFilename. All the informations,that are contained in logFile, will be removed.
 *
 * @return boolean value. If logFilename is set by setLogfilename and there is no IO-Error, this method returns a false, otherwise true.
 *
 * @see setLogFilename
 */
bool Logger::createLogFile(){
    redirectToLogFile(true);
    if( !this->logFilename.trimmed().isEmpty()){
        QFile logFile(this->logFilename);
        logFile.remove();
        if(logFile.open(QIODevice::Text|QIODevice::WriteOnly)){

        }else{
            this->writeToLogFile(trUtf8("The log can not open"));
        }
        logFile.close();
        return true;
    }else{
        this->writeToLogFile(trUtf8("The log's filename is empty"));
    }
    return false;
}
/**
 * @brief This method protocol a message on the standard output or in a logFile. All messages,that generates by RAMD3DGUI, should be protocol over this function.
 *
 * @param logText, contains the information, that should be logged.
 *
 * @details This method writes the given message logText in logFile or on STDOUT. This depends on redirectToFile's value. If redirectToFile's  value is true, the information
 * will be redirected to the Logfile, that's defiened by logFilename.
 *
 * @see setLogFilename
 * @see redirectToLogFile
 */
void  Logger::writeToLogFile(QString logText){
    QDateTime time=QDateTime::currentDateTime();
    QString message=time.toString(trUtf8("dd/MM hh:mm:ss  ->  "))+logText.at(0).toUpper()+logText.mid(1);
    if(!message.endsWith("."))message=message.trimmed().append(".");
    int row=this->messageModel->rowCount();
    this->messageModel->insertRows(row,1);
    this->messageModel->setData(this->messageModel->index(row),message);
    if(redirectToFile){
        if( !this->logFilename.trimmed().isEmpty()){
            QFile logFile(this->logFilename);
            if(logFile.open(QIODevice::Text|QIODevice::WriteOnly|QIODevice::Append)){
                QTextStream in(&logFile);
                in<<message<<endl;
                logFile.close();
            }
        }else{

        }
    }else{
        QTextStream in(stdout);
        in<<message<<endl;
    }



}
/**
 * @brief This method redirect protocol to STDOUT/LogFile depends on what is given by parameter redirectToLogFile.
 *
 * @param redirectToLogFile is a boolean-parameter. If it is true, the protocol will be redirected to logfile, otherwise to STDOUT
 *
 */
void Logger::redirectToLogFile(bool redirectToLogFile){
    if(redirectToLogFile){
       this->writeToLogFile(trUtf8("Logs are redirected to %1 now").arg(this->logFilename.trimmed()));
    }
    this->redirectToFile=redirectToLogFile;
}
/**
 * @brief This method set the logFilename.
 *
 * @param filename.
 *
 */
void Logger::setLogFilename(QString filename){
    this->logFilename=filename;
}
