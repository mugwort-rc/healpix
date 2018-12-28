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
pro remove_dipole, map, weight, bad_data=bad_data, gal_cut = gal_cut, coord_in = coord_in, coord_out=coord_out, dipole = dipole, $
monopole=monopole, noremove=noremove, nside=nside_usr, onlymonopole=onlymonopole, $
ordering=ordering, pixel=pixel, $
units = units, help=help
;+
; NAME:
;   remove_dipole
;
; PURPOSE:
;   remove best fit monopole and dipole (simultaneous fit)
;
; CATEGORY:
;   Healpix map processing
;
; CALLING SEQUENCE:
;  remove_dipole, map, [weight, bad_data=, gal_cut=, coord_in=, coord_out=, dipole=, monopole=, $
;             noremove=, nside=, onlymonopole=, ordering=, pixel=, units=, help=]
; 
; INPUTS:
;   map : array on which monopole and dipole are to be removed
;      (also used for output)
;      assumed to be a full sky data set, unless pixel is set and has the same
;      size as map
;
; OPTIONAL INPUTS:
;   weight : array of same size as map,
;     describe weighting to apply to each pixel for the fit
;
; KEYWORD PARAMETERS:
;   bad_data : scalar float, value given on input to bad pixels
;         default = !healpix.bad_value = -1.63750000e+30
;   gal_cut : the pixels with galactic latitude |b|<gal_cut are not considered in the
;        fit. In Degrees, in [0, 90]
;   coord_in : coordinate system (either 'Q' or 'C' equatorial, 'G' galactic or 'E' ecliptic) of the input map
;         default = 'G' (galactic)
;   coord_out : coordinate system in which to output dipole vector in Dipole
;         default = same as coord_in
;   noremove : if set, the dipole and monopole are computed but not removed
;   nside : scalar integer, healpix resolution parameter
;   onlymonopole : fit and remove only the monopole
;   ordering : string, describe pixelisation (either 'RING' or 'NESTED')
;   pixel : vector, gives the actual list of pixel whose temperature is given in map
;         useful in case of very limited sky coverage
;   units : units of the input map
;   help:   displays this information header
;
; OUTPUTS:
;   map : contains the map minus monopole and dipole 
;    unless noremove is set
;
; OPTIONAL OUTPUTS:
;   dipole : 3-vector : amplitude and direction of dipole
;   monopole : amplitude of monopole, same units as map
;
; SIDE EFFECTS:
;
;
;
; PROCEDURE:
;   least square fit
;   directly hacked from A. J. Banday's Dipole
;   'optimized' to reduce memory requirement
;
; EXAMPLE:
;
;
;
; MODIFICATION HISTORY:
;      2000-02-16, EH, Caltech, 1.0
;      2002-08-16, EH, Caltech, 1.1 cosmetic editions
;      2006-06-23, EH, IAP: total() -> total(,/double) for improved accuracy
;      2008-08-21, EH, IAP: accept Nside>8192; added /HELP keyword; slight speed
;      increase
;      2009-10-30, EH, IAP: replaced obsolete SVD with SVDC (+ SVSOL)
;                           monopole and dipole now output in double precision
;-

defsysv, '!healpix', exists = exists
if (exists ne 1) then init_healpix

syntax = 'remove_dipole, map, [weight, bad_data=, gal_cut=, coord_in=, coord_out=, dipole=, monopole=, noremove=, nside=, onlymonopole=, ordering=, pixel=, units=, help=]'

code = 'remove_dipole'
if keyword_set(help) then begin
    doc_library,code
    return
endif

if n_params() lt 1 or n_params() gt 2 then begin
    print,syntax
    if n_params() gt 2 then message,'Abort' else return
endif
loadsky

obs_npix = n_elements(map)
obs_npix2 = n_elements(weight)
obs_npix3 = n_elements(pixel)
if undefined(bad_data) then bad_data = !healpix.bad_value ;-1.63750000e+30


dow8 = 0
if (obs_npix2 gt 1) then begin
    if (obs_npix2 ne obs_npix) then begin
        message,syntax,/info
        message,'inconsistent map and weight'
    endif else begin
        dow8 = 1
    endelse
endif

partial = 0
if (obs_npix3 gt 1) then begin
    if (obs_npix3 ne obs_npix) then begin
        message,syntax,/info
        message,'inconsistent map and pixel'
    endif else begin
        partial = 1
    endelse
endif

galcut = 0
scut = ' '
if defined(gal_cut) then begin
    if (gal_cut) gt 0. and gal_cut lt 90 then begin
        galcut = 1
        zcut = sin(gal_cut*!dtor)
        sgal_cut = '|b| > '+string(gal_cut,form='(f7.3)')+' Deg'
    endif
endif
if undefined(coord_in) then coord_in ='G'
if coord_in  EQ 'C' then coord_in =  'Q'
fullcoord = decode_coord(coord_in,error=error_coord1)
if undefined(coord_out) then coord_out = coord_in
if coord_out  EQ 'C' then coord_out =  'Q'
fullcoord_out = decode_coord(coord_out,error=error_coord2)

if (error_coord1+error_coord2) ne 0 then begin
    message,syntax,/info
    message,'invalid choice of coordinate'
endif

if undefined(units) then units = ' '

dodipole = 1
n = dblarr(4)
M = dblarr(4,4)
if keyword_set(onlymonopole) then begin
    dipole_out = 0
    dodipole = 0
    n = [0.d0]
    M = [0.d0]
endif
if (dodipole or galcut) then begin
    if undefined(nside_usr) or undefined(ordering) then begin
        print,code+': Nside and Ordering have to be known when '+$
          'removing dipole '
        return
    endif
endif
if defined(ordering)  then pix_type = STRUPCASE(STRMID(ordering,0,1))
if defined(nside_usr) then pix_param = nside_usr

npiece = 24
if (not dodipole and not galcut) then npiece = 12
stride = obs_npix/npiece > 1


;-------------- mean square equation ------------------------------
;-------------- construct matrix and vector ------------------------
; do a loop to minimize arrays size
for is = 0,npiece do begin
    imin = is*stride 
    imax = (imin+stride-1) < (obs_npix-1)
    if (imin gt imax) then goto, enough
    map_tmp = map[imin:imax]
    ; replace NaN by bad_data to allow calculations
    nan = where(finite(map_tmp,/nan), nnan)
    if (nnan gt 0) then map_tmp[nan] = bad_data
    nan = 0

    ; pixel weight
    flag = (abs(map_tmp/bad_data-1.) gt .01) ; byte array, either 1 for good pixels or 0 for bad ones
    if dow8 then flag = flag*weight[imin:imax]

    if (dodipole or galcut) then begin
        ; compute pixel position from pixel number
        if partial then begin
            id_pix = pixel[imin:imax]
        endif else begin
            id_pix = (pix_param gt 8192) ? lindgen64(imax-imin+1)+imin : lindgen(imax-imin+1)+imin
        endelse
        case pix_type of
            'R' : PIX2VEC_RING, pix_param, id_pix, vec ; Healpix ring
            'N' : PIX2VEC_NEST, pix_param, id_pix, vec ; Healpix nest
            'Q' : vec = PIX2UV(id_pix, pix_param) ; QuadCube (COBE cgis software)
            else : print,code+': error on pix_type'
        endcase
        if (galcut) then begin
            if coord_in ne 'G' then vecg = SKYCONV(vec, inco = coord_in, outco =  'G') else vecg=vec
            flagg = abs(vecg[*,2]) gt zcut ; 1 for pixels outside galactic cut, 0 otherwise
            flag = flag*flagg
            vecg = 0
        endif
    endif

    if dodipole then begin
        temp = map_tmp*flag
        n[0] += total(temp, /double)
        n[1] += total(temp*vec[*,0], /double)
        n[2] += total(temp*vec[*,1], /double)
        n[3] += total(temp*vec[*,2], /double)
    endif else begin
        n[0] += total(map_tmp*flag, /double)
    endelse
    M[0,0] += TOTAL(flag, /double)
    if dodipole then begin
        for i = 1,3 do begin
            M[i,0] += TOTAL(vec[*,i-1]*flag, /double)
            for j = 1,i do begin ; only fill one half of symmetric matrix
                M[i,j] += TOTAL(vec[*,i-1]*vec[*,j-1]*flag, /double)
            endfor
        endfor
    endif
endfor
enough:

print,'Best fit: '+scut
if not dodipole then begin
    ;monopole = float( (n/M)[0] )
    monopole = (n/M)[0]
    print,'Monopole ['+units+']: ',monopole
endif else begin
    ; fill other half of matrix by symmetry
    for i=1,3 do M[0:i-1,i] = M[i,0:i-1]

    ; check matrix diagonal
    err_on_mat = abs((M[1,1]+M[2,2]+M[3,3])/M[0,0]-1.d0)
    if (err_on_mat gt 1.d-10) then begin
        print, err_on_mat
        message,'Error in matrix construction'
    endif

   ; use SVD to invert matrix M and solve, A is inverse of Single Value
   ; Decomposition of M
    SVDC, M, w, U, V, /double ; w contains (positive) eigenvalues
    w_threshold = max(abs(w)) * 1.0d-06
    if (is_gdl()) then begin ; GDL: no SVSOL, no DIAG_MATRIX
        nw = n_elements(w)
        Wp = dblarr(nw,nw)
        large = where(abs(w) GT w_threshold, count) ; find large eigenvalues
        if (count gt 0) then Wp[large, large] = 1.d0/w[large] ; only invert those
        C = transpose(U) # ( Wp # (V # N))
    endif else begin
        w = w * ( abs(w) gt w_threshold) ; set low eigenvalues to zero
        C = SVSOL( U, w, V, N, /double)
    endelse

    ;;;C = float(C)
    Monopole = C[0]
    Dipole = C[1:3]

    print,'Monopole ['+units+']: ',monopole
    dip_ampli = sqrt(total(dipole^2, /double))
    print,'Dipole : amplitude: ',dip_ampli
    vec2ang,dipole,lat,long,/astro
    print,fullcoord+' coordinates [Deg]: ',long, lat,form='(a,f7.2,f7.2)'
    dipole_out = dipole
    if (coord_in ne coord_out) then begin
        dipole_out = SKYCONV(dipole, inco = coord_in, outco =  coord_out)
        vec2ang,dipole_out,lat,long,/astro
        print,fullcoord_out+' coordinates [Deg]: ',long, lat,form='(a,f7.2,f7.2)'
    endif
endelse

; correct map for monopole and dipole

if not keyword_set(noremove) then begin

; do a loop to minimize arrays size
    for is=0,npiece do begin
        imin = is*stride 
        imax = (imin+stride-1) < (obs_npix-1)
        if (imin gt imax) then goto, done
        map_tmp = map[imin:imax]

        ; replace NaN by bad_data to allow calculations
        nan = where(finite(map_tmp,/nan), nnan)
        if (nnan gt 0) then map_tmp[nan] = bad_data
        nan = 0

        good = imin + where( abs(map_tmp/bad_data-1.) gt .01, ngood )

        if (ngood gt 0) then begin
            if (dodipole) then begin
           ; compute pixel position from pixel number
                if partial then id_pix = pixel[good] $
                           else id_pix = good
                case pix_type of
                    'R' : PIX2VEC_RING, pix_param, id_pix, vec ; Healpix ring
                    'N' : PIX2VEC_NEST, pix_param, id_pix, vec ; Healpix nest
                    'Q' : vec = PIX2UV(id_pix, pix_param) ; QuadCube (COBE cgis software)
                    else : print,code+': error on pix_type'
                endcase
                map[good] = map[good] - (monopole + dipole ## vec )
            endif else begin
                map[good] = map[good] - monopole
            endelse
        endif

    endfor
done:
endif

dipole = dipole_out

return
end

