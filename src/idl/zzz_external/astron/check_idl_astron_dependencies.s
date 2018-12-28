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
# 2009-10-08: deal with 2009Aug20 version
#             smarter search of name occurence
#

#set dir = /home/soft/rsi/external_contributions/astron_2008Mar07/pro/
#set dir = /Applications/itt/external_contributions/astron_2009Aug20/pro/
set dir = /Applications/itt/external_contributions/astron_2009Nov25/pro/

# full list of routine name (no path, no trailing .pro)
set fulllist = `ls $dir/*/*pro | awk -F/ '{print $NF}' | awk -F. '{print $1}' | sort`

# # name of routines already copied
# set loclist = `ls *.pro  | awk -F. '{print $1}' | sort`

set missing = 0
foreach name ($fulllist)
	if (-e $name.pro) then
	   # local copy exist:
	   # nothing
	   # echo 'already here '$name
	else
 		set n = `grep -i $name *pro | grep -v ':;' | wc -c`
		if ($n > 0) then
		   @ missing ++
		   echo '--------'$name $missing'--------'
		   grep -i "[ =(]${name}[ ,;(]" *pro | grep -v ':;' | grep -i $name
		   echo	
		endif	
	endif	
end

# 2007-05-23:
# ./astro/month_cnv.pro ./fits/fits_test_checksum.pro ./fits_table/ftaddcol.pro
# ./misc/blkshift.pro ./misc/xdispstr.pro ./misc/n_bytes.pro ./misc/wherenan.pro
#

exit

