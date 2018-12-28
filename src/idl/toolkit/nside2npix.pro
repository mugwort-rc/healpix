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
function nside2npix, nside, error=error
;+
; npix = nside2npix(nside, error=error)
;
; returns npix = 12*nside*nside
; number of pixels on a Healpix map of resolution nside
;
; if nside is not a power of 2 <= 8192,
; -1 is returned and the error flag is set to 1
;
; MODIFICATION HISTORY:
;
;     v1.0, EH, Caltech, 2000-02-11
;     v1.1, EH, Caltech, 2002-08-16 : uses !Healpix structure
;-

defsysv, '!healpix', exists = exists
if (exists ne 1) then init_healpix

error = 1
; is nside a power of 2 ?
;  junk = where(nside eq [1L,2,4,8,16,32,64,128,256,512,1024,2048,4096,8192], count)
junk = where(nside eq !healpix.nside, count)
if count ne 1 then return,-1

npix = 12L* long(nside)^2

error = 0
return, npix
end

