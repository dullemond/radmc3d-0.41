#include <mainwindow.h>

#include "setup.h"



/**
*   @file main.cpp
*   @brief This is the main source file, that the startpoint of the program
*
*   @author Farzin Sereshti
*   @version 1.0
*/

/**
* @brief this the main, startpoint of the program.
*
* @param argc the number of arguments.
* @param argv the arguments adequate to number argc of the program, that can be given to the program.
*
*/

int main(int argc, char *argv[]){

    int currentExitCode = 0;

    MainWindow *mainWindow;
    QApplication application(argc,argv);
    do{


        mainWindow=new MainWindow();
        mainWindow->resize(800,600);


        mainWindow->show();

        currentExitCode=application.exec();

    }while( mainWindow->exitCode == MainWindow::EXIT_CODE_REBOOT );

    return currentExitCode;



}
