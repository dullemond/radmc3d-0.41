/*
* uname.c
*
*  Created on: 19.04.2013
*      Author: farzin sereshti
*/
#include <sys/utsname.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <utime.h>
#include <time.h>
#include <sys/types.h>
#include <unistd.h>


/*
* \brief this method returns the system name.
* \param name is a pointer of char array, that will contains the system-name.
*
*/
void get_system_name(char** name)
{
struct utsname utname;
uname(&utname);
*name= utname.sysname;
}
/*
* \brief this method set the last modification time for the given file-name
* \param file-name is the name of file.
* \param seconds is the number of seconds that should be added to current
* last modification time. It's value can be also negative.
*
*/
void change_last_modification_time(char *filename,int seconds)
{
struct utimbuf time;
struct stat buf;
stat(filename,&buf);
#if defined(__linux__)
	    time.modtime=buf.st_mtim.tv_sec+seconds;
        time.actime=buf.st_atim.tv_sec;
#else
	time.modtime=buf.st_mtimespec.tv_sec+seconds;
	time.actime=buf.st_atimespec.tv_sec;
#endif
utime(filename,&time);
}

/*
* \brief this method create a named pipe with the given file-name.
* \param filename is the name of pipe named, that should be created.
*
*/
int create_named_pipe(char *filename)
{
return mkfifo(filename,S_IRWXU);
}
/*
* \brief this checks the path, if the path is a directory
* \param foldername is a path.
* If the given path is a directory, this method checks the permission's modus.
* If the given directory is writeable, this method will return 1, otherwise -1
* If the given is not a directory, 0 will return.
*/
int is_directory(char *foldername){
 struct stat m;
 stat(foldername,&m);
 if(S_ISDIR(m.st_mode)){
     if(access(foldername,W_OK)!=-1){
        return 1;
     }
        return -1;
 }
 return 0;
}

