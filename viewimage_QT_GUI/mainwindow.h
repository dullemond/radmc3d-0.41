#ifndef MAINWINDOW_H
#define MAINWINDOW_H



#if QT_VERSION < 0x050000
      #include <QtGui>
#else
       #include <QtWidgets>
#endif
#if QT_VERSION >= 0x050000
#include <QtConcurrent/QtConcurrent>
#endif
#include "setup.h"
#include "imagecontroller.h"
#include "physicalConstants.h"
#include "image.h"
#include "logger.h"
#include "makeimage.h"
#include "spectrum.h"
#include "makespectrum.h"
#include "spectrumlabelwidget.h"
#include "imagelabelwidget.h"
#include "unformattedsetting.h"
#include "spectrumsetting.h"


class MainWindow : public QMainWindow
{
    Q_OBJECT
public:
    MainWindow(QWidget *parent = 0);
    static int const EXIT_CODE_REBOOT;
    int exitCode;
private:
    struct Point{
        double x;
        double y;
        double x1;
        double x2;
        double y1;
        double y2;
    };

    UnformattedSetting *unformattedSetting;
    QSize lastsize;
    QAction* escAction;
    QAction *savedefaultsAction;
    QMenu *historyMenu;
    QAction *nativeMenuBarAction;
    Spectrum *spectrum;
    QMenu *helpMenu;
    QAction* loadPointsAction;
    QVector <spectrumSetting*> spectrumSettingVector;
    QImage scaleImage;
    QVector <QString> history;
    void startRadmc3d();
    void startRadmc3dAsChild();
    void radmc3dWriteGridFile(bool pipe);
    PhysicalConstants *physicalConstants;
    Qt::TransformationMode transformMode;
    bool radmc3dProcessIsLocalExisting;
    void setRadmc3dProcessLocalExisting(bool existing){this->radmc3dProcessIsLocalExisting=existing;}
    bool IsRadmc3dProcessLocalExisting(){return this->radmc3dProcessIsLocalExisting;}
    QImage lastImage;
    Logger *logger;
    const static int coordDistanceToImage=75;
    SetUp *setup;
    QVector<Point> loadedPoints;
    QMenu *imageMenu;
    QMenu *spectrumMenu;
    QLabel *label;
    QMenu *configMenu;
    QAction *localAction;
    QCheckBox *formatted;
    QTextEdit *textEdit;
    QAction *manualAction;
    QAction *guimanualAction;
    QWidgetAction *readImageAction;
    QWidgetAction *readSpectrumAction;
    QAction *userTransferAction;
    QAction *lineOptionAction;
    QAction *colorAction;
    QAction *smallAction;
    QAction *xLinAction;
    QActionGroup* sizeGroup;
    QActionGroup* freqGroup;
    QActionGroup* fluxGroup;
    QAction *yLinAction;
    QAction *changeCurrentDirectoryAction;
    QPushButton *closeLogPushButton;
    QAction *aboutAction;
    QAction *whatsthisAction;
    QGroupBox *textGroup;
    QString radmc3dStarterPATH;
    QVBoxLayout* texteditbox;
    QStringList PATHDirectories;
    void renameImageFile();
    void renameSpectrumFile();
    QString radmc3dPATH;
    QStackedWidget *leftWidget;
    QMessageBox *msgBox;
    QTabWidget *rightWidget;
    QBoxLayout::Direction currentDirection;
    QGroupBox *spectrumGroup;
    QGroupBox *imageGroup;
    QHBoxLayout *mainLayout;
    Image *image;
    QString imageControllerWidgetTabString;
    ImageControllerWidget *imageControllerWidget;
    ImageLabelWidget *imageLabel;
    SpectrumLabelWidget *spectrumLabel;
    QStringList freqUnitsList;
    QMenu *sizemenu;
    QMenu *freqmenu;
    QMenu *fluxmenu;
    QStringList fluxUnitsList;
    QStringList sizesList;
    QVector<QVector <QPoint > > contourPoints;
    QVector <QVector<QPoint> > Contour(QVector< QVector<double > > d,int ilb,int iub,int jlb,int jub,
                                       QVector<double >  x,QVector<double >  y,int nc,QVector<double >  z);
    void changeImageColorLookUpTable(bool invert=false);
    QVector<QPoint> PointInterpolation(const QVector<QPoint> &points);
    QProcess *radmc3dProcess;
    QVector<double> zoomBox;
    void calculateContour();
    void translate();
    void save();
    void saveAsDefaults();
    void setUpWhatsThis();
    void load();
    void loadAsdefaults();
    void setImageMenu();
    void checkLaststate();
    void setShortcuts();
    void setUp();
    void setIcons();
    void sendRadmc3dExitCommand();
    void updateImageMenu();
    void setMainWindowIcon(const QIcon& icon);
    QProgressDialog*progressDialog;
    void setSlots();
    QTimer *timer;
    bool findRadmc3dStarterPath();
    bool readRadmc3dPath();
    bool nativeMenuBar;
    bool zoomActive;
    bool initialized;
    bool firstTime;
    bool showLog;
    bool exit;
    bool loaded;
    bool setZoomBox;
    bool startnexttime;
    bool waitForFinishing;
    bool spectrumIsRead;
    bool imageIsRead;
    bool isErrorbar;
    MakeImage *makeImage;
    MakeSpectrum *makeSpectrum;
    QFuture<void>  future;
    QFutureWatcher<void> futureWacther;
    void setMainLayoutDirection(QBoxLayout::Direction dir);
    void pixelEdit();
    void imageEdit(bool changeColorTable,bool  invertPixel,bool calculateContour,bool pixelEdit);
    QThread  *thread;
    QMutex mutex;
    QVector <QAction *> historyActions;
    void initialize();
    void setupHistory(bool);
    int historyMax;
    QString readLastWorkingDirectory();
    void changeCurrentDirectory(QString path);
    void unzoomSpectrum();
protected:
    void closeEvent(QCloseEvent *event);
protected:
    void resizeEvent(QResizeEvent* event);
private slots:
    void drawRect(QPoint,QPoint);
    void removeZoom();
    void saveAsDefaultsSlot();
    void nativeMenuBarChanged(bool);
    void reloadImage();
    void reloadSpectrum();
    void checkFormatted(QString);
    void readImage();
    bool readSpectrum();
    void freqUnitChanged(QAction*);
    void fluxUnitChanged(QAction *action);
    void killRadmc3d(QCloseEvent *event);
    void killRadmc3d();
    void showProgressBar();
    void sendCommand(const QString objectName);
    void hideProgressBar();
    void hideProgressBarAndReload();
    void hideProgressBarAndShowRadmc3dLog();
    void hideProgressBarAndEditAndReload();
    void showAbout();
    void showManual();
    void showGUIManual();
    void changeIcon(QString);
    void showRadmc3dLog();
    void sizeUnitChanged();
    void closeMsg(int);
    void setWindowSize(bool);
    void hideLog();
    void showPixelValue(QString value);
    void calculateZoom(QPoint,QPoint);
    void radmc3dProcessFinished(int,QProcess::ExitStatus);
    void removeRect();
    void userTransferChanged(bool);
    void localObserverChanged(bool checked,bool reload=true);
    void colorChanged(bool,bool reload=true);
    void lineOptionChanged(bool,bool reload=true);
    void readMyAction();
    void setActionVisible(const QString objectName,bool visibility);
    void changeCurrentDirectory();
    void changeCurrentDirectoryByAction(QAction*);
    bool saveDialog();
    void closeMessage(int);
    void setCancelable(bool);
    void readAndReloadImage();
    void readAndReloadSpectrum();
    void drawCircle(int,int);
    void xLinearActive(bool);
    void yLinearActive(bool);
    void loadPointsSlot(bool);
signals:
    void setZoomOutVisible(bool);
    void readTransferFile();
    void activateLocalObserver(bool);
    void removeTransferTab();
};

#endif // MAINWINDOW_H
