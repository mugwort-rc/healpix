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
pro oplot_sphere, u, v,  line_type=line_type, _extra = oplot_kw

if undefined(line_type) then line_type = 0

bad = where(abs(u-shift(u,1)) gt .1, nbad)

iw = index_word(tag_names(oplot_kw),'PSYM',err=errword)
if errword eq 0 then begin
    if (oplot_kw.psym gt 0) then nbad = 0
endif

if (nbad eq 0) then begin
    if (line_type lt 0) then begin
        oplot, u, v, _extra = oplot_kw, col=1 ; white background
        oplot, u, v, _extra = oplot_kw, col=0, lines=abs(line_type)
    endif else begin
        oplot, u, v, _extra = oplot_kw, lines=line_type
    endelse
endif else begin
    bad = [0,bad,n_elements(u)-1]
    for j=0,nbad do begin
        if (bad[j+1] gt bad[j]) then begin
            u1 = u[bad[j]:bad[j+1]-1]
            v1 = v[bad[j]:bad[j+1]-1]
            if (line_type lt 0) then begin
                oplot, u1, v1, _extra = oplot_kw, col=1 ; white background
                oplot, u1, v1, _extra = oplot_kw, col=0, lines=abs(line_type)
            endif else begin
                oplot, u1, v1, _extra = oplot_kw, lines=line_type
            endelse        
        endif
    endfor
endelse




return
end

