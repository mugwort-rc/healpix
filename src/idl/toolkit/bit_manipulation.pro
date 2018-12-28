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

function swapLSBMSB, i
; Returns i with even and odd bit positions interchanged.

oddbits = 89478485L
evenbits=178956970L

li = long(i)
swapLSBMSB = (li AND evenbits)/2 + (li AND oddbits)*2

return, swapLSBMSB
end

; -------

function invLSBMSB, i
; Returns NOT(i)

invLSBMSB = NOT (long(i))

return, invLSBMSB
end

; -------

function invswapLSBMSB, i
; Returns NOT(i) with even and odd bit positions interchanged.

invswapLSBMSB = NOT swapLSBMSB(i)

return, invswapLSBMSB
end

; -------

function invLSB, i
; Returns i with odd (1,3,5,...) bits inverted.

oddbits = 89478485L

li = long(i)
invLSB = (li XOR oddbits)

return, invLSB
end

; -------

function invMSB, i
; Returns i with even (0,2,4,...) bits inverted.

evenbits=178956970L

li = long(i)
invMSB = (li XOR evenbits)

return, invMSB
end

;--------------------

pro bit_manipulation

x =  swapLSBMSB(1)
x =  invLSBMSB(1)
x =  invswapLSBMSB(1)
x =  invLSB(1)
x =  invMSB(1)

return
end
