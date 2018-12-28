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
PRO pix2ang_ring, nside, ipix, theta, phi
;****************************************************************************************
;+
; PIX2ANG_RING, nside, ipix, theta, phi
; 
;       renders Theta and Phi coordinates of the nominal pixel center
;       given the RING scheme pixel number Ipix and map resolution parameter Nside
;
; INPUT
;    Nside     : determines the resolution (Npix = 12* Nside^2)
;	SCALAR
;    Ipix  : pixel number in the RING scheme of Healpix pixelisation in [0,Npix-1]
;	can be an ARRAY 
;       pixels are numbered along parallels (ascending phi), 
;       and parallels are numbered from north pole to south pole (ascending theta)
;
; OUTPUT
;    Theta : angle (along meridian = co-latitude), in [0,Pi], theta=0 : north pole,
;	is an ARRAY of same size as Ipix
;    Phi   : angle (along parallel = azimut), in [0,2*Pi]
;	is an ARRAY of same size as Ipix
;
; SUBROUTINE
;    nside2npix
;
; HISTORY
;    June-October 1997,  Eric Hivon & Kris Gorski, TAC
;    Aug  1997 : treats correctly the case nside = 1
;    Feb 1999,           Eric Hivon,               Caltech
;         renamed pix2ang_ring
;    Sept 2000,          EH
;           free memory by discarding unused variables
;    June 2003,  EH, replaced STOPs by MESSAGEs
;
;-
;****************************************************************************************
routine = 'PIX2ANG_RING'
if N_params() lt 3 then begin
    message,' syntax: '+routine+', nside, ipix, theta, phi'
endif

if (N_ELEMENTS(nside) GT 1) then message,'Nside should be a scalar in '+routine
npix = nside2npix(nside, error = error)
if (error ne 0) then message,'Invalid Nside '+string(nside)

nside = LONG(nside)
nl2 = 2*nside
nl3 = 3*nside
nl4 = 4*nside
ncap = nl2*(nside-1L)
nsup = nl2*(5L*nside+1L)
fact1 = 1.5d0*nside
fact2 = (3.d0*nside)*nside
np = N_ELEMENTS(ipix)
theta = DBLARR(np)
phi   = DBLARR(np)

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

pix_np = WHERE(ipix LT ncap,   n_np)   ; north polar cap
IF (n_np GT 0) THEN BEGIN ; north polar cap ; ---------------------------------

   ip = long(ipix(pix_np)) + 1
   iring = LONG( SQRT( ip/2.d0 - SQRT(ip/2) ) ) + 1L ; counted from NORTH pole, starting at 1
   iphi  = ip - 2L*iring*(iring-1L)

   theta(pix_np) = ACOS( 1.d0 - iring^2 / fact2 )
   phi(pix_np)   = (iphi - 0.5d0) * !DPI/(2.d0*iring)
   ip = 0 & iring =0 & iphi = 0    ; free memory
   pix_np = 0                      ; free memory

ENDIF ; ------------------------------------------------------------------------

pix_eq = WHERE(ipix GE ncap AND ipix LT nsup,  n_eq) ; equatorial strip
IF (n_eq GT 0) THEN BEGIN ; equatorial strip ; ---------------------------------

   ip    = long(ipix(pix_eq)) - ncap
   iring = LONG( ip / nl4) + nside                        ; counted from NORTH pole
   iphi  = ( ip MOD nl4 )  + 1

   fodd  = 0.5d0 * (1 + ((iring+nside) MOD 2)) ; 1 if iring is odd, 1/2 otherwise

   theta(pix_eq) = ACOS( (nl2 - iring) / fact1 )
   phi(pix_eq)   = (iphi - fodd) * !DPI/(2.d0*nside)
   ip = 0 & iring =0 & iphi = 0    ; free memory
   pix_eq = 0                      ; free memory

ENDIF ; ------------------------------------------------------------------------

pix_sp = WHERE(ipix GE nsup,   n_sp)   ; south polar cap
IF (n_np + n_sp + n_eq) NE np THEN message,'error in '+routine
IF (n_sp GT 0) THEN BEGIN ; south polar cap ; ---------------------------------

   ip =  npix - long(ipix(pix_sp))
   iring = LONG( SQRT( ip/2.d0 - SQRT(ip/2) ) ) + 1      ; counted from SOUTH pole, starting at 1
   iphi  = 4*iring + 1 - (ip - 2L*iring*(iring-1L))
   
   theta(pix_sp) = ACOS( - 1.d0 + iring^2 / fact2 )
   phi(pix_sp)   = (iphi - 0.5d0) * !DPI/(2.d0*iring)
   ip = 0 & iring =0 & iphi = 0    ; free memory
   pix_sp = 0                      ; free memory

ENDIF ; ------------------------------------------------------------------------


RETURN
END

;=======================================================================
; The permission to use and copy this software and its documentation, 
; without fee or royalty is limited to non-commercial purposes related to 
; Microwave Anisotropy Probe (MAP) and
; PLANCK Surveyor projects and provided that you agree to comply with
; the following copyright notice and statements,
; and that the same appear on ALL copies of the software and documentation.
;
; An appropriate acknowledgement has to be included in any
; publications based on work where the package has been used
; and a reference to the homepage http://www.tac.dk/~healpix
; should be included
;
; Copyright 1997 by Eric Hivon and Kris Gorski.
;  All rights reserved.
;=======================================================================
