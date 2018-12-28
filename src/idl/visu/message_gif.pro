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
pro message_gif, code = code, error=error
;+
; message_gif : prints an error message about the non-availability of GIF
;
;
;-
error = 0

version = !version.release
vers_no_gif = 5.3 ; to be checked (smaller than 5.4 for sure)

scode = code+'> '
if (version ge vers_no_gif) then begin
    print,'=========================================='
    print,scode+'     IDL version : '+string(version,form='(f3.1)')
    print,scode+'ERROR : GIF is not supported anymore by IDL'
    print,scode+'ERROR : Use PNG instead'
    print,'=========================================='
    error = 1
endif


return
end
