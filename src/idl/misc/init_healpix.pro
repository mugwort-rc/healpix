; -----------------------------------------------------------------------------
;
;  Copyright (C) 1997-2008  Krzysztof M. Gorski, Eric Hivon, Anthony J. Banday
;
;
;
;
;
;  This file is part of HEALPix.
;
;  HEALPix is free software; you can redistribute it and/or modify
;  it under the terms of the GNU General Public License as published by
;  the Free Software Foundation; either version 2 of the License, or
;  (at your option) any later version.
;
;  HEALPix is distributed in the hope that it will be useful,
;  but WITHOUT ANY WARRANTY; without even the implied warranty of
;  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;  GNU General Public License for more details.
;
;  You should have received a copy of the GNU General Public License
;  along with HEALPix; if not, write to the Free Software
;  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
;
;  For more information about HEALPix see http://healpix.jpl.nasa.gov
;
; -----------------------------------------------------------------------------
pro init_healpix, verbose=verbose
;+
; defines the (structure) system variable Healpix
;
; 2006-05-22: print caution message in version 5.5a which has touchy expand_path
; 2006-oct : v 2.10, enabled nside > 8192
; 2008-oct: v2.11, use getenv
;-

; system variable name
healpix_sysvar = '!HEALPIX'

; Healpix version
version = '2.11'

; release data
date = '2008-11-13'

; Healpix directory
sv = 'HEALPIX'
;dsv = '$'+sv

; if (!version.release eq '5.5a') then message,'testing existence of '+sv+' environnement variable',/info
; directory = expand_path(dsv)
; directory = strtrim(directory,2)
; if (directory eq dsv or directory eq '') then begin
;     print,' system variable '+sv+' not found '
;     directory = ''
; endif
directory = getenv(sv)
if (strtrim(directory,2) eq '') then begin
    print,' system variable '+sv+' not found '
    directory = ''
endif

; list of possible Nside's
nside = 2L^lindgen(30)  ; 1, 2, 4, 8, ..., 8192, ..., 2^29 = 0.54e9

; flag for missing values
bad_value = -1.6375e30

comment = ['This system variable contains some information on Healpix :', $
           healpix_sysvar+'.VERSION   = current version number,', $
           healpix_sysvar+'.DATE      = date of release,',$
           healpix_sysvar+'.DIRECTORY = directory containing Healpix package,',$
           healpix_sysvar+'.NSIDE     = list of all valid values of Nside parameter,',$
           healpix_sysvar+'.BAD_VALUE = value of flag given to missing pixels in FITS files,',$
           healpix_sysvar+'.COMMENT   = this description.']

; create structure
stc = {version:version, date:date, directory:directory, nside:nside, bad_value:bad_value, comment:comment}

; fill variable out
defsysv, healpix_sysvar, exists = exists
if (exists) then begin
    !Healpix = stc
endif else begin
    defsysv, healpix_sysvar, stc
endelse

if (keyword_set(verbose)) then begin

    print,'Initializing '+healpix_sysvar+' system variable'
    print
    print,comment,form='(a)'
    print
;     help,/st,healpix_sysvar
;     print
endif



return
end

