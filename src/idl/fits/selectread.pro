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
pro selectread, file, array, polvec, header=exthdr, columns=columns, extension=extension_id, poltype=poltype, tonan=tonan, offset=offset, factor=factor, flip = flip, no_pdu= no_pdu
;+
; NAME:
;    selectread
;
; PURPOSE:
;    routine to read a FITS file containing temperature or temperature+
;    polarisation data. It has been optimized for speed and low memory usage
;
; CATEGORY:
;    healpix I/O
;
; CALLING SEQUENCE:
;    selectread, File, Array, [Polvec, HEADER=, COLUMNS=, EXTENSION=, POLTYPE=,
;    TONAN=, OFFSET=, FACTOR=, FLIP=, NO_PDU=]
;
; INPUTS:
;     File : string containing FITS file name
;
; OPTIONAL INPUTS:
;
;      
; KEYWORD PARAMETERS:
;     COLUMNS= list of numbers specifying columns to be read (1 based).
;       By default all columns are read, unless POLTYPE is set, in which case
;       only the relevant columns are read and processed
;
;     EXTENSION= extension to be read (0 based)
;
;     POLTYPE= processing to apply to polarised data before output
;        0: do nothing: output selected column(s) as read (after application of
;               OFFSET and FACTOR to each columns)
;        1: read the Q and U columns of the FITS file and outputs in Array the
;               norm = sqrt(U^2+Q^2) * factor
;        2: read the Q and U columns of the FITS file and output in Array the
;              angle psi = 1/2 atan(U/Q)
;        3: read the T,Q and U columns and outputs 
;              in Array (T + offset) * factor
;              polvec will contain norm and psi as described above
;
;     TONAN= set flagged pixel to NaN. Is set whn FACTOR and/or OFFSET are set
;
;     OFFSET= multiplicative factor applied to data, default = 1.
;
;     FACTOR= additive factor applied to data 
;      (ie, data_out = (data_fits + offset) * factor
;
;     FLIP= flip orientation of polarisation vector 
;      (ie, change the sign of psi, only applies when POLTYPE=2 or 3)
;
;     NO_PDU= if set, the primary data unit header is not included in the
;     returned header
;
; OUTPUTS:
;     Array: array of size (np, ncols) 
;        where np is the number of data points, ncols is the number of columns
;        read, and the array type is determined by the data read.
;
; OPTIONAL OUTPUTS:
;     Polvec: if POLTYPE=3, will contain the polarisation norm and direction in
;       an array of size (np, 2)
;
;     HEADER: contains the FITS header (combination of the primary unit header
;     and the one of the extension being read, unless NO_PDU is set)
;
; COMMON BLOCKS:
;
;
; SIDE EFFECTS:
;
;
; RESTRICTIONS:
;
;
; PROCEDURE:
;
;
; EXAMPLE:
;
;
; MODIFICATION HISTORY:
;         v1.0, EH, 2005-02-07
;                   2005-08-29: added no pdu
;
;-

if (n_params() lt 2) then begin
    print,'selectread, File, Array, [Polvec, HEADER=, COLUMNS=, EXTENSION=, POLTYPE=,'
    print,'                                  TONAN=, OFFSET=, FACTOR=, FLIP=, NO_PDU=]'
    return
endif

defsysv, '!healpix', exists = exists
if (exists ne 1) then init_healpix

do_rescale = (keyword_set(tonan) or keyword_set(offset) or keyword_set(factor))
xtn   = (keyword_set(extension_id)) ? (extension_id+1) : 1
polar = (keyword_set(poltype))      ? poltype          : 0
flipconv = (keyword_set(flip))      ? 1                : -1

; open file
fits_open, file, fcb

; find number of words in extension
if (xtn gt fcb.nextend) then begin
    message,'not enough extensions in '+file
endif

n_wpr  = (fcb.axis)[0,xtn] ; words per row in extension
n_rows = long64((fcb.axis)[1,xtn]) ; number of rows
n_words = n_wpr * n_rows

; read header for extension + optionally PDU
fits_read, fcb, void, exthdr, exten_no = xtn, /header_only, no_pdu = no_pdu
tfields  = round(float(sxpar(exthdr,'TFIELDS')))
bad_data =       float(sxpar(exthdr,'BAD_DATA', count=nbd))
if (nbd eq 0) then bad_data = !healpix.bad_value

; read image if file is in deprecated format
if (xtn eq 0 and fcb.nextend eq 0) then begin
    if polar gt 0 then message,'no polarisation information found in '+file
    fits_read, fcb, array, exten_no = xtn
    if (do_rescale) then begin
        bad_pixels = where(array le (bad_data*0.9) or finite(array,/nan), nbad)
        if (nbad gt 0)    then array[bad_pixels] = !values.f_nan
        if (factor ne 1.) then array = temporary(array) * factor 
        if (offset ne 0.) then array = temporary(array) + (factor*offset)
    endif
    goto, all_done
endif

; check out file for polarisation
if (tfields lt 3 and polar gt 0) then begin
    message,'no polarisation information found in '+file
endif

tbinfo, exthdr, tab_xhdr
nentry = max(tab_xhdr.numval)
npix = nentry * n_rows

; columns to be output in array
cols = indgen(tfields) + 1
if defined(columns) then cols = columns
if (max(cols) gt tfields or min(cols) le 0) then begin
    print,'columns = ',cols
    message,'invalid choice of columns for '+file
endif

nmaps = n_elements(cols)
if (polar eq 1 or polar eq 2) then begin
    cols = [2,3] ; read Q and U
    nmaps = 1
endif
if (polar eq 3) then begin
    cols = [1,2,3]
    nmaps = 1
endif
if (polar gt 3 or polar lt 0) then begin
    print,'Poltype = ',polar
    message,'invalid poltype for '+file
endif
types = (tab_xhdr.idltype)[cols-1]
type = max(types) ; find out type of output array to fit all data
junk = where(types eq 7, n_string)
if (n_string ge 1) then begin
    message,'Table in input file contains strings'
endif


; create array receiving final data
array = make_array(type=type, npix, nmaps, /nozero)
if (polar eq 3) then polvec = make_array(type=type, npix, 2, /nozero)

; read data piece by piece and process each piece individually
stride = 5.e6 ; 5 MB per piece
stride = FLOOR(stride / n_wpr) * n_wpr
w_start = long64(0)
pstart = long64(0)
while (w_start LE (n_words-1) ) do begin
    ; read one piece
    w_end = (w_start + stride - 1L) < (n_words-1)
    fits_read, fcb, data, exten_no = xtn, first=w_start, last=w_end
    nr = (w_end - w_start + 1) / n_wpr ; number of rows read
    np = nr * nentry ; number of pixels read
    data = reform(data, n_wpr, nr, /overwrite) ; required by tbget
    ; select useful columns and/or process data read
    if (polar eq 0) then begin ; standard case
        for i=0,nmaps-1 do begin
            x = (tbget(tab_xhdr, data, cols[i]))[*]
            if (do_rescale) then begin
                bad_pixels = where(x le (bad_data*0.9) or finite(x,/nan), nbad)
                if (nbad gt 0)    then x[bad_pixels] = !values.f_nan
                if (factor ne 1.) then x = x * factor 
                if (offset ne 0.) then x = x + (factor*offset)
            endif
            array[pstart,i] = x
        endfor
    endif else begin
        if (w_start eq 0) then begin
            tmp  = make_array(type=type,np,3,/nozero) 
            norm = make_array(type=type,np,  /nozero) & psi = norm
        endif
        for i=0,2 do begin
            tmp[0,i] = (tbget(tab_xhdr, data, i+1))[*]
            if (do_rescale) then begin
                bad_pixels = where(tmp[*,i] le (bad_data*0.9) or finite(tmp[*,i],/nan), nbad)
                if (nbad gt 0)    then tmp[bad_pixels,i] = !values.f_nan
                if (factor ne 1.) then tmp[*,i] = tmp[*,i] * factor 
                if (offset ne 0. and i eq 0) then tmp[*,i] = tmp[*,i] + (factor*offset) ; only offset temperature
            endif
        endfor
        if (polar eq 1 or polar eq 3) then begin
            norm[0:np-1] = sqrt(tmp[0:np-1,1]^2+tmp[0:np-1,2]^2)  ; (Q,U) --> P
        endif
        if (polar eq 2 or polar eq 3) then begin
            psi[0:np-1] = 0.5 * atan(tmp[0:np-1,2]*flipconv, tmp[0:np-1,1]) ; (Q,U) --> psi
        endif
        if (polar eq 1) then array[pstart] = norm[0:np-1]
        if (polar eq 2) then array[pstart] = psi[0:np-1]
        if (polar eq 3) then begin
            array[pstart,0] = tmp[0:np-1,0]
            polvec[pstart,0] = norm[0:np-1]
            polvec[pstart,1] = psi[0:np-1]
        endif
    endelse
    ; get ready for next piece
    w_start = w_end + 1
    pstart = pstart + np
endwhile
data = 0
tmp = 0 & norm = 0 & psi = 0

; close file
all_done:
fits_close, fcb



return
end
