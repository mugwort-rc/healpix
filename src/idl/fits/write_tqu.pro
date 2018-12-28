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
pro write_tqu, fitsfile, TQU, Coordsys=coordsys, Nested=nested, Ring=ring, Ordering=ordering, Extension=extension_id, Hdr=hdr, Xhdr=xhdr, Units=units, Help=help, Error=error
;+
; NAME:
;    write_tqu
;
; PURPOSE:
;    writes a temperature+polarization Healpix map into a
;    binary table FITS file, 
;    2 formats are supported:
;    1) (T,Q,U) in extention 0, 
;       and optionally the error (dT,dQ,dU) in ext. 1
;       and correlation (dQU, dTU, dTQ) in extension 2
;    2) (T, Q, U, N_OBS) in extension 0
;       and optionally the counts (N_OBS, N_QQ, N_QU, N_UU) in ext 1
;
; CATEGORY:
;
;
; CALLING SEQUENCE:
;   WRITE_TQU, fitsfile, TQU, Coordsys=, Nested=, Ring=, Ordering=, $
;                      extension=, Hdr=, Xhdr=, Units=, Help=
; 
; INPUTS:
;   fitsfile : output filename
;
;   TQU : array of Healpix maps of size (npix,n_maps,n_ext) where npix is the total
;   number of Healpix pixels on the sky, n_maps=3 or 4, and n_ext <=3.
;     n_maps maps are written in each extension of the FITS file : 
;   either
;    1) (T,Q,U) in extention 0, 
;       and optionally the error (dT,dQ,dU) in ext. 1
;       and correlation (dQU, dTU, dTQ) in extension 2
;    2) (T, Q, U, N_OBS) in extension 0
;       and optionally the counts (N_OBS, N_QQ, N_QU, N_UU) in ext 1
;
;     it is also possible to write n_maps maps directly in a given
;     extension (provided the preceding extension, if any, is already
;     filled in)
;     by setting EXTENSION to the extension number in which to write
;     (0 based) and if n_ext + Extension <= 3
;
; OPTIONAL INPUTS:
;      
; KEYWORD
;    Ring-     if set, add 'ORDERING= RING' to the extension fits header
;    Nested-   if set, add 'ORDERING= NESTED' to the extension fits header
;    Ordering- if set to the string 'nested' or 'ring', set the
;              keyword 'ORDERING' to the respective value
;
;    one of them has to be present unless the ordering information is already
;    present in the fits header.
;
; OPTIONAL KEYWORD
;   Coordsys = if set to either 'C', 'E' or 'G' specifies that the
;     Healpix coordinate system is respectively Celestial=equatorial, Ecliptic or
;     Galactic
;   
;   Extension = extension unit a which to put the data (0 based)
;
;   Help= if set, an extensive help (this IDL header) is printed
;
;   Hdr= string containing the primary FITS header
;
;   Units = physical units of the maps
;
;   Xhdr= string containing the extension FITS header (it will be
;    repeated in each extension, except for TTYPE* and EXTNAME which
;    are generated by the routine and depend on the extension)
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;    Error = takes value 1 on output if error occurs
;
; COMMON BLOCKS:
;
; SIDE EFFECTS:
;
; RESTRICTIONS:
;
; PROCEDURE:
;   calls WRITE_FITS_SB
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
;    version 1.0, Eric Hivon, Nov-Dec 2002
;
;-

if (keyword_set(help)) then begin
    doc_library,'write_tqu'
    return
endif

syntax = ['SYNTAX : WRITE_TQU, fitsfile, TQU, Coordsys=, Nested=, Ring=, Ordering=, $',$
          '            Extension=, Hdr=, Xhdr=, Units=, Help=, Error=']

if n_params() lt 2 then begin
    print,syntax,form='(a)'
    if n_params() eq 0 then return
    message,'Not enough argument, no file written'
endif

if (datatype(fitsfile) ne 'STR' or datatype(TQU) eq 'STR') then begin
    print,syntax,form='(a)'
    message,'Error in argument type, no file written'
endif

sz = size(TQU)
ndim = sz[0]
npix = sz[1]
nmap = 1
n_ext = 1
if (ndim lt 2) then begin
;    message,' Expect a Npix * Nmaps * Next  TQU array'
    print,'WARNING: write_tqu: Expected Npix * Nmaps * Next array, getting Npix vector'
endif else begin
    nmap = sz[2]
    if (ndim ge 3) then n_ext = sz[3]
endelse
if (npix2nside(npix) lt 0) then message,' TQU 1st dimension (Npix) is not a valid Npix pixel number'
if (nmap ne 1 && nmap ne 3 &&  nmap ne 4) then message,' Expect a Npix map or Npix * [3 or 4]  TQU array'
if (nmap eq 1 && n_ext gt 1) then message, ' Expect Npix map'

i_ext0 = 0
if keyword_set(extension_id) then i_ext0 = extension_id
if (n_ext + i_ext0) gt 3 then begin
    print, i_ext0, n_ext
    message,' can not write more than 3 extensions'
endif

wmap_format = (nmap eq 4)
;-----------------
; create structure for primary unit
punit = 0
if defined(hdr) then punit = create_struct('HDR',hdr)

;------------------
; create names for extensions

if undefined(xhdr) then xhdr = [' ']


if (wmap_format) then begin
; name columns
    name =  [['TEMPERATURE',  'Q_POLARISATION',  'U_POLARISATION', 'N_OBS'], $
             ['N_OBS',        'N_QQ',            'N_QU',           'N_UU']]
; extension names
    xtname = ['Stokes Maps',  'Weight Arrays']
endif else begin
; name columns
;     name =  [['TEMPERATURE',  'Q_POLARIZATION',  'U_POLARIZATION'], $
    name =  [['TEMPERATURE',  'Q_POLARISATION',  'U_POLARISATION'], $
             ['dT',           'dQ',              'dU'], $
             ['dQU',          'dTU',             'dTQ']]
; extension names
    xtname = ['SIGNAL','ERROR','CORRELATION']
endelse

for i_ext=i_ext0, i_ext0+n_ext-1 do begin
; create structure for 1st/2nd/3rd extension (number 0,1,2)
    sxaddpar,xhdr,'EXTNAME',xtname[i_ext]
    if (nmap eq 1) then sxaddpar,xhdr,'POLAR','F' else sxaddpar,xhdr,'POLAR','T'

; add UNITS information
    if defined(units) then begin
        sunits = units
        if (i_ext eq 2) then sunits=strtrim(units,2)+'**2'
        add_units_fits,xhdr,units=sunits,error=error,colnum=[1,2,3]
        if error ne 0 then return
    endif

    ia = i_ext - i_ext0
    if (nmap eq 1) then begin
        xtns = create_struct('HDR', xhdr, $
                             name[0,i_ext], TQU[*,0,ia])
    endif else begin
        xtns = create_struct('HDR', xhdr, $
                             name[0,i_ext], TQU[*,0,ia],$
                             name[1,i_ext], TQU[*,1,ia],$
                             name[2,i_ext], TQU[*,2,ia])
    endelse
    if (wmap_format) then begin
        xtns = create_struct(xtns, name[3,i_ext], TQU[*,3,ia])
    endif

; write it
    ii_ext = i_ext
    write_fits_sb, fitsfile, punit, xtns, $
      Coordsys=coordsys, Nested=nested, Ring=ring, Ordering=ordering, Extension=ii_ext
endfor


return
end

