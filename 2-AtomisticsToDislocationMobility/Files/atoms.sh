#!/usr/bin/expect

set crystal_structure [lindex $argv 0]
set material [lindex $argv 1]
set x [lindex $argv 2]
set y [lindex $argv 3]
set z [lindex $argv 4]
set dislocation_type [lindex $argv 5]
set structure_type [lindex $argv 6]

spawn ./atom-dislocation
expect "Which crystal structure do you want to consider?"
send -- "$crystal_structure\r"
expect "++   what kind of material do you want to consider      ++"
send -- "$material\r"
expect "number of unit cells along X, Y and Z"
send -- "$x $y $z\r"
expect "What kind of dislocation do you want?"
send -- "$dislocation_type\r"
expect "!** What kind of structure do you want? **!"
send -- "$structure_type\r"

expect eof