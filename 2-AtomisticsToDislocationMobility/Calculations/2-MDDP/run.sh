#!/usr/bin/expect

if [fork]!=0 exit
disconnect
set timeout -1
set dir [lindex $argv 0]
spawn $dir/MDDP08-2008
expect " ok?, yes(y), No(n)"
send -- "y\r"

expect eof