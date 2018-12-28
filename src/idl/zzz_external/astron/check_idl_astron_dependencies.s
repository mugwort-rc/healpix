#!/bin/csh -f
#
# =================================
# For HEALPix internal usage only
# =================================
#
#
# check that native astron routines are not called by those copied here
#
# created 2007-May 23
#

set dir = /home/soft/rsi/external_contributions/astron_2008Mar07/pro/

# full list of routine
set fulllist = `ls $dir/*/*pro | awk -F/ '{print $NF}' | awk -F. '{print $1}' | sort`

# # name of routines already copied
# set loclist = `ls *.pro  | awk -F. '{print $1}' | sort`

set missing = 0
foreach file ($fulllist)
	if (-e $file.pro) then
	   # nothing
	   # echo 'already here '$file
	else
 		set n = `grep -i $file *pro | grep -v ':;' | wc -c`
		if ($n > 0) then
		   @ missing ++
		   echo '--------'$file $missing'--------'
		   grep -i $file *pro | grep -v ':;' | grep -i $file
		   echo	
		endif	
	endif	
end

# 2007-05-23:
# ./astro/month_cnv.pro ./fits/fits_test_checksum.pro ./fits_table/ftaddcol.pro
# ./misc/blkshift.pro ./misc/xdispstr.pro ./misc/n_bytes.pro ./misc/wherenan.pro
#

exit

