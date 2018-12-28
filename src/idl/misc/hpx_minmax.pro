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
function hpx_minmax, array,NAN=nan, DIMEN=dimen, miss=miss
;+
; NAME:
;      HPX_MINMAX
; PURPOSE:
;      Return a 2 element array giving the minimum and maximum of an array
; EXPLANATION:
;      Using MINMAX() is faster than doing a separate MAX and MIN.
;
; CALLING SEQUENCE:
;      value = minmax( array )
; INPUTS:
;      array - an IDL numeric scalar, vector or array.
;
; OUTPUTS:
;      value = a two element vector (if DIMEN is not supplied)
;            value[0] = minimum value of array
;            value[1] = maximum value of array
;
;            If the DIMEN keyword is supplied then value will be a 2 x N element
;            array where N is the number of elements in the specified
;            dimension
;              
; OPTIONAL INPUT KEYWORD:
;      /NAN   - Set this keyword to cause the routine to check for occurrences
;            of the IEEE floating-point value NaN in the input data.  Elements 
;            with the value NaN are treated as missing data.
;
;      DIMEN - (V5.5 or later) integer (either 1 or 2) specifying which 
;            dimension of a 2-d array to  take the minimum and maximum.   Note
;            that DIMEN is only valid for a 2-d array, larger dimensions are 
;            not supported.
; EXAMPLE:
;     (1)  Print the minimum and maximum of an image array, im
; 
;            IDL> print, minmax( im )
;
;     (2) Given a 2-dimension array of (echelle) wavelengths w, print the
;         minimum and maximum of each order (requires V5.5 or later)
;
;         print,minmax(w,dimen=1)
;
; PROCEDURE:
;      The MIN function is used with the MAX keyword
;
; REVISION HISTORY:
;      Written W. Landsman                January, 1990
;      Converted to IDL V5.0   W. Landsman   September 1997
;      Added NaN keyword.      M. Buie       June 1998
;      Added DIMENSION keyword    W. Landsman  January 2002
;-


dtype = size(array,/type)
t_flt = (dtype eq 4 or dtype eq 6) ; data is float (real or complex)
t_dbl = (dtype eq 5 or dtype eq 9) ; data is dbl (real or complex)

nbad = 0
my_nan = keyword_set(nan)

; Healpix-flagged values are replaced by NaN
; (only applies to float or double data)
if (t_flt or t_dbl) then begin
    defsysv, '!healpix', EXISTS = exists
    if (not exists) then init_healpix
    fbad = !healpix.bad_value

    ibad   = where(array eq fbad, nbad)
    if (nbad gt 0) then begin
        if (t_flt) then fnan = !values.f_nan
        if (t_dbl) then fnan = !values.d_nan
        array[ibad]= fnan
        my_nan     = 1
    endif
endif

On_error,2
if N_elements(DIMEN) GT 0 then begin
    amin = min(array, MAX = amax, NAN = my_nan, DIMEN = dimen) 
    if (nbad gt 0) then array[ibad] = fbad ; returns NaN to Healpix bad value
    return, transpose( [[amin], [amax] ])
endif else  begin 
    amin = min(array, MAX = amax, NAN = my_nan)
    if (nbad gt 0) then array[ibad] = fbad ; returns NaN to Healpix bad value
    return, [ amin, amax ]
endelse

end
