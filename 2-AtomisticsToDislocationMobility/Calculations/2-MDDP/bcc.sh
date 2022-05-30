#!/usr/bin/expect

set proceed [lindex $argv 0]
set data_type [lindex $argv 1]
set frs_pairs [lindex $argv 2]
set source_length_min [lindex $argv 3]
set source_length_max [lindex $argv 4]
set min_sep [lindex $argv 5]
set max_sep [lindex $argv 6]
set glide_plane [lindex $argv 7]
set random_seed [lindex $argv 8]
set year [lindex $argv 9]
set max_segment_length [lindex $argv 10]
set strain_rate [lindex $argv 11]
set loading_direction [lindex $argv 12]
set loading_condition [lindex $argv 13]

spawn ./BCCdata
expect " Yes=1 (continue), No=2 (Stop to go back to make it)"
send -- "$proceed\r"
expect " 12 = Lassila s  Case"
send -- "$data_type\r"
expect " Enter the number of Frank-Read source pairs"
send -- "$frs_pairs\r"
expect " Enter the minimum source length in units of b"
send -- "$source_length_min\r"
expect " Enter the maximum source length in units of b"
send -- "$source_length_max\r"
expect " a paired Frank-Read sources, in units of b"
send -- "$min_sep\r"
expect " a paired Frank-Read sources, in units of b"
send -- "$max_sep\r"
expect " 20) All glide planes are of type {110}"
send -- "$glide_plane\r"
expect " as a large integer say from 100001-200001:"
send -- "$random_seed\r"
expect " Are you done? 1=YES, 2=NO"
send -- "1\r"
expect " Enter Year"
send -- "$year\r"
expect " Max segment length ? (=100b)"
send -- "$max_segment_length\r"
expect " What is the strain rate?"
send -- "$strain_rate\r"
expect " what is the loading direction?"
send -- "$loading_direction\r"
expect " Loading condition? 1=creep, 0=constrant rate"
send -- "$loading_condition\r"

expect eof