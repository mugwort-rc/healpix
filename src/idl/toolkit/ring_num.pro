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
function ring_num, nside, z
;+
; gives the ring number corresponding to z for the resolution nside
;
;-
twothird = 2.d0 /3.d0

;     ----- equatorial regime ---------
iring = NINT( nside*(2.d0-1.500d0*z))

;     ----- north cap ------
if (z gt twothird) then begin
    iring = NINT( nside* SQRT(3.d0*(1.d0-z)))
    if (iring eq 0) then iring = 1
endif

;     ----- south cap -----
if (z lt -twothird) then begin
    iring = NINT( nside* SQRT(3.d0*(1.d0+z)))
    if (iring eq 0) then iring = 1
    iring = 4*nside - iring
endif

return, iring
end
