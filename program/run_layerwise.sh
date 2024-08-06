#!/bin/csh
@ i=0
while ($i<25)
touch out/test_fd_$i.XYZ
touch logs/test_fd_$i.log
nohup ./ATAL_FD sources/test.XYZ out/test_fd_$i.XYZ 1 $i> logs/test_fd_$i.log &
@ i = $i + 1
end
