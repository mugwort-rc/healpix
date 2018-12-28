; -----------------------------------------------------------------------------
;
;  Copyright (C) 1997-2005  Krzysztof M. Gorski, Eric Hivon, Anthony J. Banday
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
;-

; system variable name
healpix_sysvar = '!HEALPIX'

; Healpix version
version = 2.00

; release data
date = '2005-08-29'

; Healpix directory
sv = 'HEALPIX'
dsv = '$'+sv
directory = expand_path(dsv)
if (strtrim(directory,2) eq dsv) then begin
    print,' system variable '+sv+' not found '
    directory = ''
endif

; list of possible Nside's
nside = 2L^lindgen(14)  ; 1, 2, 4, 8, ..., 8192

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

