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
pro vec2ang, vector, theta, phi, astro = astro
;*****************************************************************************
;+
; NAME:
;       VEC2ANG
;
; PURPOSE:
;       converts  cartesian coordinates vector Vector 
;       into angular coordinates Theta, Phi
;
; CALLING SEQUENCE:
;       Vec2ang, Vector, Theta, Phi, [/Astro]
;
; INPUTS:
;       Vector : array of cartesian position of dimension (np, 3)
;        with north pole = (x=0,y=0,z=1) and origin at the center of
;        the sphere
;       and vector = (x(0), x(1), ... x(np-1), y(0), y(1), ... y(np-1),
;         z(0), z(1), ...z(np-1))
;
; KEYWORD PARAMETERS:
;       Astro : see below
;
; OUTPUTS:
;       Theta : vector of angular coordinate of size np
;       Phi   : vector of angular coordinate of size np
;
;
;        * if ASTRO is NOT set (default) : geometric coordinates
;       Theta : colatitude in RADIAN measured Southward from North pole
;       Phi   : longitude in RADIAN, increasing Eastward
;
;        * if ASTRO is set : astronomic coordinates
;       Theta : latitude in DEGREE measured Northward from Equator
;       Phi   : longitude in DEGREE, increasing Eastward
;
; MODIFICATION HISTORY:
;       March 6, 1999    Eric Hivon, Caltech, Version 1.0
;       March 22, 2002     use Message
;-
;*****************************************************************************

if N_params() lt 2 then begin
    message,' syntax = vec2ang, vec, theta, phi, astro=',/noprefix
endif

np1 = N_ELEMENTS(vector)
np = np1 / 3L
if (np1 NE 3 * np) then begin
    message,'inconsistent Vector in vec2ang',/noprefix
endif

twopi = 2. * !Pi
radeg = !RaDeg
if (DATATYPE(vector,2) EQ 5) then begin ; double precision
    twopi = 2.d0 *!DPi 
    radeg = 180.d0 / !DPi
endif

;---------------

Vector = REFORM(vector, np, 3, /OVER) ; condition the input vector

theta_rad = ACOS( vector(*,2)  /  SQRT(  TOTAL(vector^2, 2) ) )
phi_rad = ATAN( vector(*,1), vector(*,0) )  ; in [-Pi,Pi]

phi_rad = phi_rad + twopi * (phi_rad LT 0.)

IF KEYWORD_SET(astro) THEN BEGIN
    theta = 90. - theta_rad * RaDeg
    phi   = phi_rad * RaDeg
ENDIF ELSE BEGIN
    theta = theta_rad
    phi = phi_rad
ENDELSE

theta_rad = 0
phi_rad = 0


return
end

