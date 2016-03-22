#!/usr/bin/env bash

#Script installing tmux with dependencies in one go
#Check the target directory for tmux and the two dependencies libevent and ncurses
#And their version
# Shared by Mitra Barzine

#Target directory: $HOME/local (as /homes/yourlogin/local)
#tmux version: 1.8
#libevent version 2.0.21 stable
#ncurses version 5.9

# exit on error
Set -e
# create the target directories
mkdir $HOME/tmux_tmp
#move to the tempory directory
cd $HOME/tmux_tmp
 
# download source files for tmux, libevent, and ncurses ## change with the correct version if not the ones bellow
wget -O tmux-1.8.tar.gz http://sourceforge.net/projects/tmux/files/tmux/tmux-1.8/tmux-1.8.tar.gz/download
wget https://github.com/downloads/libevent/libevent/libevent-2.0.21-stable.tar.gz
wget ftp://ftp.gnu.org/gnu/ncurses/ncurses-5.9.tar.gz

#the dependencies will be installed first

#libevent 
##### don't forget to change the version to the correct one for the two following instructions
tar xvzf libevent-2.0.21-stable.tar.gz
cd libevent-2.0.21-stable

##### should change the /local by what you want in the following instruction but you will have to change them in the tmux part as well
##### and it would be better that you change accordingly the path for ncurses
##### DO NOT CHANGE the target directory if you're not at easy with the bash, NFS and how to move along the directory tree
#once the correct value given for the configure part, the installing part is as usual
./configure --prefix=$HOME/local --disable-shared ##that prevents quite a lot of the problem due to our low privileges
make
make install
cd ..
 
 
#ncurses
##### don't forget to change the version to the correct one for the two following instructions
tar xvzf ncurses-5.9.tar.gz
cd ncurses-5.9

##### should change the /local by what you want but you will have to change them in the tmux part as well
##### if you change it for libevent you should change it for the ncurses as well
##### once again, if not sure about what you're doing, don't touch that
./configure --prefix=$HOME/local 
make
make install
cd ..
 

#tmux
##### don't forget to change the version to the correct one for the two following instructions
tar xvzf tmux-1.8.tar.gz
cd tmux-1.8

##### if you didn't touch anything before, you're done :)
##### if not, be careful to the path you're changing. If you have just change $HOME/local by /some/path just replace every $HOME/local bellow by /some/path
##### if that not the case as well, well you know what you're doing, aren't you?

./configure CFLAGS="-I$HOME/local/include -I$HOME/local/include/ncurses" LDFLAGS="-L$HOME/local/lib -L$HOME/local/include/ncurses -L$HOME/local/include"
CPPFLAGS="-I$HOME/local/include -I$HOME/local/include/ncurses" LDFLAGS="-static -L$HOME/local/include -L$HOME/local/include/ncurses -L$HOME/local/lib" make

##### you can normally change this directory easily, just past wherever you want tmux in your path 
##### but then keep in mind that the last message is wrong and should been changed
cp tmux $HOME/local/bin 
cd ..
 
# the temporary directory is not needed anymore and then erased
rm -rf $HOME/tmux_tmp

# last message (has to been changed to the correct path of the directory into which you've pasted tmux):
echo "$HOME/local/bin/tmux is now installed and available. You should optionally add $HOME/local/bin to your PATH. (edit your PATH in $HOME/.bashrc or $HOME/.bash_profile)"
