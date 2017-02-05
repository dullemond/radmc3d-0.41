#!/usr/bin/perl
#=======================================================
#           INSTALL ROUTINE FOR "RADMC-3D"
#=======================================================
use Config;

$pwd  = `pwd` ; 
chop($pwd) ;
$home = $ENV{"HOME"};
$bin  = $home . "/bin" ;
if(!(-e $bin)) {
    print "You must have a bin/ directory in your home directory\n" ;
    print "-----Shall I make $bin for you?\n" ;
    $input = <STDIN> ;
    print $input ;
    if($input=~/^[yY]/) {
	print "Creating $bin for you...\n" ;
	system("mkdir $bin") ;
    }
}
$path = $ENV{"PATH"};
if(!($path =~ /$bin/)) {
    print "The $bin directory exists, but it not in the PATH environment variable\n" ;
    print "You must put the \n" ;
    print "$bin directory \n" ;
    print "in the path yourself (in the .tcshrc file if you use the tcsh shell...)\n" ;
    print "If you do it now, don't forget to type 'rehash'.\n" ;
}
print "  Creating a link 'viewimage' in '$bin/'\n" ;
$radmc3dgui   = $pwd . "-build/RADMC3D_QT_GUI" ;
$radmc3dguilnk = $bin . "/viewimage" ;
if(!(-e $radmc3dguilnk)) {
    print "------ Warning: file $radmc3dguilnk did not exist previously. You might want to type 'rehash'\n" ;
}
 
open(FILE,">$radmc3dguilnk") || die "Could not open file\n" ;
print FILE "#!/usr/bin/perl\n" ;
print FILE "system(\"$radmc3dgui \@ARGV\");" ;
close (FILE) ;


`chmod u+rwx $radmc3dguilnk` ;

