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
;+
; NAME:
;       FITS2CL
;
; PURPOSE:
;       read from a FITS file an ascii or binary table extension containing 
;       power spectrum (C_l) or spherical harmonics coefficients (a_lm) and return the
;       corresponding power spectra (C_l). Reads header information if required.
;       It can also read a FITS file containing a (beam) window function (B_l or
;       W_l), and will return the raw field (B_l or W_l).
;
; CALLING SEQUENCE:
;       FITS2CL, cl_array, fitsfile, [HELP=, SILENT =, SHOW=, RSHOW= $
;                                     HDR =, MULTIPOLES=, XHDR =, INTERACTIVE= ]
; 
; INPUTS:
;       fitsfile = String containing the name of the file to be read     
;
; OPTIONAL INPUTS:
;       HELP = if set, this header is printed out, no file is read
;
;       SILENT = if set, no comments are written while the file is read
;
;       SHOW = if set, the power spectra read from the file and 
;              multiplied by l(l+1)/2Pi   are plotted
;
;       RSHOW = if set, the RAW power spectra read from the file are plotted
;
;       INTERACTIVE = if set, the plots generated by SHOW and RSHOW options 
;           are produced using iPlot routine, allowing 
;           for interactive cropping, zooming and annotation of the plots. This
;           requires IDL 6.4 or newer to work properly.
;
; OUTPUTS:
;       cl_array = real/double array of Cl coefficients read from the file.
;                  This has dimension (lmax+1) or (lmax+1,4) or
;                  (lmax+1,6) given in the sequence T E B TxE TxB BxE
;                  The convention for the power
;                  spectrum is 
;                  cl = Sum_m {a_lm^2}/(2l+1)
;                  ie. NOT normalised by the Harrison-Zeldovich spectrum.
;
; OPTIONAL OUTPUT KEYWORDS:
;       HDR      = String array containing the header for the FITS file.
;       XHDR     = String array containing the extension header for the 
;                  extension.
;       MULTIPOLES = array containing the multipoles l for which the
;       power spectra are provided
;           * either read from the file (1st column in the Planck format),
;           * or generated by the routine (assuming that all
;               multipoles from 0 to lmax included are provided).
;
; EXAMPLE:
;       Read the power spectrum coefficients into the real array, pwr,
;       from the FITS file, spectrum.fits, with extension header information
;       contained in the string array, ext_text
;
;       IDL> fits2cl,pwr,'spectrum.fits',XHDR=ext_txt
;
;
; PROCEDURES CALLED:
;       HEADFITS, MRDFITS
;
; MODIFICATION HISTORY:
;       May 1999: written by A.J. Banday (MPA)         
;       Oct 2001, EH : can now deal with 6 columns,
;                      added silent keyword
;       Dec 2002, EH : can deal with Planck format (1st column is
;       multipole),
;                      added multipoles output variable
;       Mar 2003, EH : added SHOW keyword
;       Aug 2004, EH : added HELP keyword, updated header
;       May 2005, EH, replaced FINDFILE by FILE_TEST, can read alm file
;       Feb 2007, EH, accept files generated by LevelS beam2alm
;       May 2007, EH, read (beam) window file
;       Oct 2007, EH, more effective /silent
;       Jan 2008, EH, addition of /interactive
;
; requires the THE IDL ASTRONOMY USER'S LIBRARY 
; that can be found at http://idlastro.gsfc.nasa.gov/homepage.html
;
;-
;***************************************************************************
pro alm2cl_sub, a1_r, a1_i, a2_r, a2_i, cl=cl
;
; ALM2CL_SUB
;
; computes (cross-spectra) from alm
;
; alm2cl_sub, a1_r, a1_i, a2_r, a2_i, cl=cl
;
; the alm are ordered with m varying the fastest (first index)
;
;

lmax1 = n_elements(a1_r[0,*]) -1
lmax2 = n_elements(a2_r[0,*]) -1

mmax1 = n_elements(a1_r[*,0]) -1
mmax2 = n_elements(a2_r[*,0]) -1

lmax = min([lmax1,lmax2])
mmax = min([mmax1,mmax2])
cl = dblarr(lmax+1)

for ll=0L,lmax do begin
    mm = min([ll,mmax])
    xx = a1_r[0,ll]*a2_r[0,ll]
    if (mm gt 0) then begin
        xx = xx + 2.0d0*TOTAL(a1_r[1:mm,ll]*a2_r[1:mm,ll] + a1_i[1:mm,ll]*a2_i[1:mm,ll], /double)
    endif
    cl[ll] = xx/(2.d0*ll+1.d0)
endfor

return
end
;***************************************************************************

PRO FITS2CL, cl_array, fitsfile, HDR = hdr, HELP = help, MULTIPOLES=multipoles, SILENT=silent, $
             SHOW=show, RSHOW=rshow, XHDR = xhdr, INTERACTIVE=interactive

code = 'FITS2CL'
syntax = ['Syntax : '+code+', cl_array, fitsfile, [/HELP, /SILENT, /SHOW, /RSHOW, ',$
          '                                       HDR=, MULTIPOLES=, XHDR=, /INTERACTIVE]' ]

if keyword_set(help) then begin
      doc_library,code
      return
endif

if N_params() NE 2  then begin
      print,syntax,form='(a)'
      print,'   for more details: '+code+',/help'
      print,'   file NOT read '
      goto, Exit
endif

if datatype(cl_array) eq 'STR' or datatype(fitsfile) ne 'STR' then begin
      print,syntax,form='(a)'
      print,'   the array comes first, then the file name '
      print,'   file NOT read '
      goto, Exit
endif

if (not file_test(fitsfile)) then message,'file '+fitsfile+' not found'

; run astrolib routine to set up non-standard system variables
defsysv, '!DEBUG', EXISTS = i  ; check if astrolib variables have been set-up
if (i ne 1) then astrolib      ; if not, run astrolib to do so

hdr  = HEADFITS(fitsfile)
xhdr = HEADFITS(fitsfile,EXTEN=1)
fits_info, fitsfile, /silent, n_ext=n_ext

ttype1 = sxpar(xhdr,'TTYPE1',count=nttype1)

;read_cl = (strtrim(strupcase(ttype1),2) ne 'INDEX')
read_cl = (strmid(strupcase(ttype1),0,5) ne 'INDEX') ; to accept files generated by beam2alm

sytit = 'C' ; generic y-title

if (read_cl) then begin
; ******** FITS file contains C(l) ********
; -------- Planck format ? (ie, first column is multipole)
    pdmtype = sxpar(xhdr,'PDMTYPE',count=nextra)
    if (nextra gt 0) then junk = where(strupcase(pdmtype) eq 'POWERSPEC', nextra)
    nextra = nextra < 1         ; either 0 or 1
    
    if (nextra eq 0) then begin
        if (strupcase(ttype1) eq 'L') then nextra = 1
    endif
    
    if ( (not keyword_set(silent) and (nextra eq 1))) then print,'Info : C(l) file '+fitsfile+' has Planck format'
; -------- extension -----------

; simply read the extension from the FITS file
    tmpout = MRDFITS(fitsfile,1,/use_colnum,silent=silent)

; get the dimensions of the input array
    info  = size(tmpout)
    nrows = info(1)             ; # of entries for l-range: nrows = lmax+1
    ncols = n_elements(tag_names(tmpout)) - nextra

    if ( (ncols ne 1) and (ncols ne 4) and (ncols ne 6)) then begin
;         print,' Input file does not conform to expected structure'
;         print,code+' expects either 1, 4 or 6 columns,      found ',ncols
;         goto, Exit
        print,'WARNING: Input file does not conform to structure expected for C(l) file'
        print,'WARNING: '+code+' expects either 1, 4 or 6 columns,      found ',ncols
;        goto, Exit
    endif

    extname = sxpar(xhdr,'EXTNAME', count=nextname)
;    extname = (nextname eq 1) ? strupcase(strmid(extname,0,5)) : 'POWER'
    extname = strupcase(strmid(extname,0,5))
    case extname of
        'POWER': sytit = 'C'
        'BEAM ': sytit = 'B'
        'WINDO': sytit = 'W'
        else: sytit = 'C'
    endcase

; output array
    if (ncols eq 1) then begin
        cl_array = tmpout.(0+nextra)
    endif else begin
        cl_array = make_array(nrows,ncols,type=datatype(tmpout.(nextra),2))
        for i=0,ncols-1 do begin
            cl_array[*,i] = tmpout.(i+nextra)
        endfor
    endelse

; generate l multipole 
    if nextra eq 1 then begin
        multipoles = tmpout.(0)
    endif else begin
        multipoles = findgen(nrows)
    endelse

    names = sxpar(xhdr,'TTYPE*')
    units = sxpar(xhdr,'TUNIT*',count=nunits)
    if (nunits eq 0) then begin
        units = replicate('',ncols)
    endif else begin
        if nunits eq 1 then units = replicate(units,ncols)
        units = '['+strtrim(units,2)+']'
    endelse
endif else begin
; ******** FITS file contains a_lm ********
    if (not keyword_set(silent)) then print,'Info : computes C(l) from alm'
    n_ext = n_ext < 3
    nfcl = long(n_ext * (n_ext + 1))/2
    nextra = 0 ; no l column
    ncols = nfcl

    type = ['T','E','B']
    for i=0,n_ext-1 do begin
        fits2alm, index, alm_array, fitsfile, type[i]

        index2lm, index, l, m
        lmax = max(l)
        mmax = max(m)
        nl = n_elements(index)
        index = 0
        if (not keyword_set(silent)) then print,'type, lmax, mmax: ',type[i],lmax,mmax

        if (i eq 0) then begin
            cl_array = dblarr(lmax+1, nfcl)
            lmax0 = lmax
        endif

        alm_real = dblarr(mmax+1, lmax+1)
        alm_imag = alm_real

        alm_real[ m, l] = alm_array[0:nl-1, 0]
        alm_imag[ m, l] = alm_array[0:nl-1, 1]

                                ; auto-spectra for TT, EE and BB
        alm2cl_sub, alm_real, alm_imag, alm_real, alm_imag, cl=cl
        cl_array[0,i] = cl[0:lmax]

        if (i eq 0) then begin
                                ; store T for latter cross-spectra calculation
            at_r = alm_real
            at_i = alm_imag
        endif

        if (i eq 1) then begin
                                ; computes TE and store E for latter calculation
            alm2cl_sub, alm_real, alm_imag, at_r, at_i, cl=cl
            cl_array[0,3] = cl[0:lmax]
            ae_r = alm_real
            ae_i = alm_imag
        endif
        
        if (i eq 2) then begin
                                ; computes TB
            alm2cl_sub, alm_real, alm_imag, at_r, at_i, cl=cl
            cl_array[0,4] = cl[0:lmax]
                                ; computes EB
            alm2cl_sub, alm_real, alm_imag, ae_r, ae_i, cl=cl
            cl_array[0,5] = cl[0:lmax]
        endif
    endfor

    multipoles = findgen(lmax0+1)

    names = ['TT','EE','BB','TxE','TxB','ExB']
    units = sxpar(xhdr,'TUNIT2',count=nunits)
    if (nunits eq 0) then begin
        units = replicate('',6)
    endif else begin
        units = '['+replicate(units,6)+'^2]'
    endelse
endelse


if (keyword_set(show) or keyword_set(rshow)) then begin
    ; plot size and lay out
    screen_size = get_screen_size()
    xs = 800<screen_size[0]
    ys = ((ncols+1)/2)*400 < screen_size[1]
    grid=[2<ncols,(ncols+1)/2]

    ; data to plot and labels
    l = multipoles
    if (keyword_set(show)) then begin ; show l(l+1)C(l)/2 Pi
        fl2 = l*(l+1.)/(2.*!pi)
        ytitle = '!12  l (l+1) !6'+sytit+'!12!dl!n /!6 2!7p!n!6 '
    endif else begin ;                  show C(l)
        fl2 = 1
        ytitle = ' !6'+sytit+'!12!dl!n !6 '        
    endelse
    xtitle = '!6 Multipole !12 l !6'

    ; actual plot
    if (keyword_set(interactive)) then begin
        if (float(!version.release) lt 6.4) then begin
            message,/info,code+',/interactive  requires IDL 6.4 or more to work properly'
            return
        endif
        view_grid = grid
        icol=0
        iplot,l,fl2*cl_array[*,icol], dimension=[xs,ys], title=fitsfile, $ 
              xtitle=xtitle,ytitle=ytitle+units[icol+nextra],view_title='!6'+names[icol+nextra], $
              view_grid=view_grid,view_zoom=1.0, anisotropic_scale_2d=1.0
        for icol=1,ncols-1 do begin
            iplot,l,fl2*cl_array[*,icol], dimension=[xs,ys], title=fitsfile, view_number = (icol+1), $ 
                  xtitle=xtitle,ytitle=ytitle+units[icol+nextra],view_title='!6'+names[icol+nextra], $
                  view_grid=view_grid,view_zoom=1.0, anisotropic_scale_2d=1.0
        endfor
    endif else begin
        window, /free, title=fitsfile, xs=xs, ys=ys
        !p.multi=[0, grid]
        for icol=0,ncols-1 do begin
            plot,l,fl2*cl_array[*,icol], $
                 xtitle=xtitle,ytitle=ytitle+units[icol+nextra],title='!6'+names[icol+nextra], $
                 charsize=1.3
        endfor
        !p.multi=0
    endelse
endif

; Exit routine
Exit:
return
end



