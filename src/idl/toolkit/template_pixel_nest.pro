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
pro template_pixel_nest, nside, ipix, itplt, reflexion
;+
; NAME:
;         template_pixel_nest
;
; PURPOSE:
;         Returns template-pixel index corresponding to each (NESTED ordered) pixel provided
;         The reflexion (E-W, or N-W or both) to apply to the pixel to match the
;         template is optionally returned
;
; CATEGORY:
;
; CALLING SEQUENCE:
;         TEMPLATE_PIXEL_NEST, Nside, Ipix, Template [,Reflexion]
;
; INPUTS:
;         Nside : scalar,           integer, Healpix resolution parameter
;         Ipix  : vector or scalar, integer, list of nest ordered pixel indices
;
; OUTPUTS:
;         Template : vector or scalar, integer, template-pixel index
;         corresponding to Ipix (see Notes below). Same dimension as Ipix.
;
; OPTIONAL OUTPUTS:
;         Reflexion : vector or scalar, integer in {0,3} (see Notes below). 
;         Same dimension as Ipix.
;
; COMMON BLOCKS:
;         none
;
; SIDE EFFECTS:
;         none
;
; RESTRICTIONS:
;         ??
;
; PROCEDURE:
;
; EXAMPLE:
;         calls NEST2RING and TEMPLATE_PIXEL_RING
;
; RELATED ROUTINES:
;         NEST2RING, RING2NEST
;         SAME_SHAPE_PIXELS_NEST, SAME_SHAPE_PIXELS_RING
;         TEMPLATE_PIXEL_NEST,    TEMPLATE_PIXEL_RING
;
; MODIFICATION HISTORY:
;         2004-10, EH, v1.0
;
;
;-

routine = 'template_pixel_nest'
nparams = n_params()
if nparams lt 3 or nparams gt 4 then begin
    print,' syntax: '+routine+', Nside, Ipix, Template [,Reflexion]'
    message,'Abort'
    return
endif
do_sym = (nparams eq 4)

if (N_ELEMENTS(nside) GT 1) then message,'Nside should be a scalar in '+routine
npix = nside2npix(nside, error = error)
if (error ne 0) then message,'Invalid Nside '+string(nside)

lnside = long(nside)

min_pix = MIN(ipix)
max_pix = MAX(ipix)
IF (min_pix LT 0) THEN BEGIN
    PRINT,'pixel index : ',min_pix,FORMAT='(A,I10)'
    PRINT,'is out of range : ',0,npix-1,FORMAT='(A,I2,I8)'
    message,'Abort'
ENDIF
IF (max_pix GT npix-1) THEN BEGIN
    PRINT,'pixel index : ',max_pix,FORMAT='(A,I10)'
    PRINT,'is out of range : ',0,npix-1,FORMAT='(A,I2,I8)'
    message,'Abort'
ENDIF

nest2ring, nside, ipix, ip_ring
template_pixel_ring, nside, ip_ring, itplt, reflexion

return
end
