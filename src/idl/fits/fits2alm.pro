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
PRO FITS2ALM, index, alm_array, fitsfile, signal, HDR = hdr, XHDR = xhdr

;+
; NAME:
;       FITS2ALM
;
; PURPOSE:
;       read from a FITS file a binary table extension containing spherical
;       harmonic (scalar or tensor) coefficients together with their index 
;       and return them. Reads header information if required.
;
; CALLING SEQUENCE:
;       FITS2ALM, index, alm_array, fitsfile, [signal, HDR = , XHDR = ]
; 
; INPUTS:
;       fitsfile = String containing the name of the file to be written      
;
; OPTIONAL INPUTS:
;       signal   = String defining the signal coefficients to read
;                  Valid options: 'T' 'E' 'B' 'ALL' 
;
; OUTPUTS:
;       index    = Integer array containing the index for the corresponding 
;                  array of alm coefficients (and erralm if required). The
;                  index {i} is related to {l,m} by the relation
;                  i = l^2 + l + m + 1
;      alm_array = real/double array of alm coefficients read from the file.
;                  This has dimension (nl,nalm,nsig) -- corresponding to
;                  nl   = number of {l,m} indices
;                  nalm = 2 for real and imaginary parts of alm coefficients;
;                         4 for above plus corresponding error values
;                  nsig = number of signals extracted (1 for any of T E B
;                         or 3 if ALL extracted). Each signal is stored
;                         in a separate extension.
;
; OPTIONAL OUTPUT KEYWORDS:
;       HDR      = String array containing the header for the FITS file.
;       XHDR     = String array containing the extension header(s). If 
;                  ALL signals are required, then the three extension 
;                  headers are returned appended into one string array.
;
; EXAMPLE:
;       Read the B tensor harmonic coefficients into the real array, alm,
;       from the FITS file, coeffs.fits, with extension header information
;       contained in the string array, ext_text
;
;       IDL> fits2alm,alm,'coeffs.fits','B',XHDR=ext_txt
;
;
; PROCEDURES CALLED:
;       HEADFITS,  FITS_INFO, FITS_READ, TBINFO, TBGET
;
; MODIFICATION HISTORY:
;       May 1999: written by A.J. Banday (MPA)     
;       Dec 2004, EH: edited to avoid faulty /use_colnum keyword in MRDFITS  
;       Feb 2005, EH: replaced MRDFITS by faster FITS_READ+TBINFO+TBGET
;       May 2005, EH, replaces FINDFILE by FILE_TEST
;       Aug 2005, EH: make output alm_array of same precision as FITS file data
;  Jan 2008, EH: calls tbfree to remove heap pointer created by TBINFO
;
; requires the THE IDL ASTRONOMY USER'S LIBRARY 
; that can be found at http://idlastro.gsfc.nasa.gov/homepage.html
;
;-

if N_params() LT 3 or N_params() gt 6 then begin
      print,'Syntax : FITS2ALM, index, alm_array, fitsfile, [signal, HDR = , XHDR = ] '
      goto, Exit
endif

if (not file_test(fitsfile)) then message,'file '+fitsfile+' not found'

; run astrolib routine to set up non-standard system variables
defsysv, '!DEBUG', EXISTS = i  ; check if astrolib variables have been set-up
if (i ne 1) then astrolib      ; if not, run astrolib to do so

hdr  = HEADFITS(fitsfile)

; -------- extension -----------

if(undefined(signal))then signal = 'T'

signal = STRUPCASE(signal)

CASE signal OF
            'T'  : BEGIN
                     extension = 1 
                     nsig = 1
                   END
            'E'  : BEGIN
                     extension = 2
                     nsig = 1
                   END
            'B'  : BEGIN
                     extension = 3 
                     nsig = 1
                   END
            'ALL': BEGIN
                     extension = 1 
                     nsig = 3 ; extension value initialises read-loop
                   END

             else: BEGIN
                     print,' Incorrect signal selected'
                     goto, Exit
                   END
ENDCASE

; count extensions
fits_info,fitsfile, /silent, n_ext=n_ext

; simply read the extensions from the FITS file
savehdr = ''
nrows_old = -1
for i = 0,nsig-1 do begin
    exten = extension+i
    if (exten gt n_ext) then begin
        message,' Required extension does not exist in file'
    endif

    ; read data
    fits_read, fitsfile, tmpout, xhdr, /no_pdu, exten_no = exten
    nrows = sxpar(xhdr,'NAXIS2')
    savehdr = [savehdr,xhdr]
  
    if (i eq 0) then begin
        ; first extension, create arrays
        tbinfo, xhdr, tab_xhdr
        ncols = n_elements(tab_xhdr.(0))

        index = tbget(tab_xhdr, tmpout, 1)
        alm_array = make_array(nrows, ncols-1, nsig, type=(tab_xhdr.idltype)[1])
    endif else begin
        ; other extensions, check consistency
        if (nrows ne nrows_old) then begin
            print,exten-1,nrows_old,exten,nrows,   $
              form='("#",i1,":", i8,",  #",i1,":",i8)'
            message,' ERROR: Extensions have different sizes'
        endif
        indtmp = tbget(tab_xhdr, tmpout, 1)
        if (not array_equal(index, indtmp)) then begin
            message,' ERROR: Alm''s defined for different (l,m)'
        endif
    endelse
    nrows_old = nrows

    for jc=2,ncols do begin
        alm_array[*,jc-2,i] = tbget(tab_xhdr, tmpout, jc)
    endfor

endfor

alm_array = reform(alm_array)
xhdr = savehdr
tbfree, tab_xhdr

; Exit routine
Exit:
return
end



