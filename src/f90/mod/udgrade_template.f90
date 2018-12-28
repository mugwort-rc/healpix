!-----------------------------------------------------------------------------
!
!  Copyright (C) 1997-2008 Krzysztof M. Gorski, Eric Hivon, 
!                          Benjamin D. Wandelt, Anthony J. Banday, 
!                          Matthias Bartelmann, Hans K. Eriksen, 
!                          Frode K. Hansen, Martin Reinecke
!
!
!  This file is part of HEALPix.
!
!  HEALPix is free software; you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation; either version 2 of the License, or
!  (at your option) any later version.
!
!  HEALPix is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with HEALPix; if not, write to the Free Software
!  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
!
!  For more information about HEALPix see http://healpix.jpl.nasa.gov
!
!-----------------------------------------------------------------------------
! template for routine SP/DP overloading for module udgrade_nr

! K M A P   : map kind                 either SP or DP
!
! K L O A D : suffixe of routine name, to be replaced by either s (SP) or d (DP)
  !=======================================================================
  subroutine sub_udgrade_nest_KLOAD(map_in, nside_in, map_out, nside_out, fmissval, pessimistic)
    !=======================================================================
    !
    !     SUB_UDGRADE_NEST(map_in, nside_in, map_out, nside_out)
    !
    !     map_in    : input map (REAL)
    !     nside_in  : resolution parameter of the input map (npix =12*nside**2)
    !                 must be a power of 2 (INTEGER)
    !     map_out   : output map (REAL)
    !     nside_out : resolution parameter of the output map (npix =12*nside**2)
    !                 must be a power of 2 (INTEGER)
    !
    !     upgrades or degrades a map in the NESTED numbering
    !
    !     if nside_out > nside_in : upgrades
    !     value(each new sub-pixel) = value(old big pixel)
    !
    !     if nside_in < nside_out : degrades
    !     value(new big pixel) = sum values(old small pixels) / number_of_small_pixels
    !
    !     version 1.0 : Eric Hivon, TAC, September 1997
    !     version 1.1 : Dec 1997, correction of a bug in the map indices
    !     version 1.2 : Dec 2001, addition of pessimistic and non-pessimistic behaviors
    !=======================================================================
    USE pix_tools, only : nside2npix
    INTEGER(I4B), INTENT(IN) :: nside_in, nside_out
    REAL(KMAP),     INTENT(IN),  DIMENSION(0:), target :: map_in
    REAL(KMAP),     INTENT(OUT), DIMENSION(0:) :: map_out
    REAL(KMAP),     INTENT(IN), OPTIONAL :: fmissval
    LOGICAL(LGT), INTENT(IN), OPTIONAL :: pessimistic

    INTEGER(I4B) :: npix_in, npix_out, npratio
    INTEGER(I4B) :: iu, ip, id
    LOGICAL(LGT), ALLOCATABLE, DIMENSION(:) :: good
    integer(I4B) :: nobs
    LOGICAL(LGT) :: do_pessimistic = .false.
    REAL(KMAP) :: bad_value
    REAL(KMAP), DIMENSION(:), POINTER :: p_in
    !-----------------------------------------------------------------------

    npix_out = nside2npix(nside_out)
    npix_in  = nside2npix(nside_in)

    bad_value = -1.6375e30_KMAP
    if (present(fmissval)) then
       bad_value = fmissval
    endif

    do ip=0, npix_out-1
       map_out(ip) = bad_value
    enddo

    if (nside_out < nside_in) then ! degrade
       if (present(pessimistic)) then
          do_pessimistic = pessimistic
       endif

       npratio = npix_in  / npix_out

!$OMP parallel default(none) &
!$OMP   shared(map_in, map_out, npix_out, npratio, do_pessimistic, bad_value) &
!$OMP   private(good, id, p_in, nobs)
       allocate(good(0:npratio-1))
!$OMP do schedule(dynamic,64)
       do id=0,npix_out-1
          p_in => map_in(id*npratio:(id+1)*npratio-1)
          good = (p_in /= bad_value)
          nobs = count(good)

          if (do_pessimistic) then
             if (nobs == npratio) then
                map_out(id) = sum(p_in) / nobs
             endif
          else
             if (nobs > 0) then
                map_out(id) = sum(p_in, mask=good) / nobs
             endif
          endif
       enddo
!$OMP end do
       deallocate(good)
!$OMP end parallel

    else ! upgrade
       npratio = npix_out  / npix_in
!$OMP parallel default(none) &
!$OMP   shared(map_in, map_out, npix_out, npratio) private(iu, ip)
!$OMP do schedule(dynamic,64)
       do iu = 0, npix_out - 1
          ip = iu/npratio
          map_out(iu) = map_in(ip)
       enddo
!$OMP end do
!$OMP end parallel
    endif

    return
  end subroutine sub_udgrade_nest_KLOAD
!=======================================================================
  !
  !     UDGRADE_RING(map_in, nside_in, map_out, nside_out)
  !
  !     map_in    : input map (REAL)
  !                 This map will be converted to nested format on output.
  !     nside_in  : resolution parameter of the input map (npix =12*nside**2)
  !                 must be a power of 2 (INTEGER)
  !     map_out   : output map (REAL)
  !     nside_out : resolution parameter of the output map (npix =12*nside**2)
  !                 must be a power of 2 (INTEGER)
  !
  !     upgrades or degrades a map in the RING numbering
  !
  !     if nside_out > nside_in : upgrades
  !     value(each new sub-pixel) = value(old big pixel)
  !
  !     if nside_in < nside_out : degrades
  !     value(new big pixel) = sum values(old small pixels) / number_of_small_pixels
  !
  !
  !     calls the ring2nest and nest2ring routines available in pix_tools module
  !
  !
  !
  !     NB : the speed is not the main concern of this routine,
  !     which not supposed to be used intensively for big maps
  !
  !     version 1.0 : Eric Hivon, TAC, September 1997
  !     version 1.1 : Dec 1997, correction of a bug in the map indices
  !     version 2.0 : Dec 2004-March 2005: SP/DP and 1D/ND overload
  !=======================================================================
  subroutine udgrade_ring_1d_KLOAD(map_in, nside_in, map_out, nside_out, fmissval, pessimistic)
    !=======================================================================
    !    1 dim implementation
    !=======================================================================
    USE pix_tools, only : nside2npix, ring2nest, convert_ring2nest, convert_nest2ring
    INTEGER(I4B), INTENT(IN) :: nside_in, nside_out
    REAL(KMAP), INTENT(INOUT), dimension(0:), target :: map_in
    REAL(KMAP), INTENT(OUT), dimension(0:),   target :: map_out
    REAL(KMAP), INTENT(IN), OPTIONAL :: fmissval
    LOGICAL(LGT) , INTENT(IN), OPTIONAL :: pessimistic

    INTEGER(I4B) :: npix_in, npix_out
    !-----------------------------------------------------------------------

    !     checks that the 2 nside are valid
    npix_out = nside2npix(nside_out)
    npix_in  = nside2npix(nside_in)
    if (npix_out < 0) then
       print*,"wrong nside_out in udgrade_ring : ", nside_out
       call fatal_error
    endif
    if (npix_in  < 0) then
       print*,"wrong nside_in  in udgrade_ring : ", nside_in
       call fatal_error
    endif

    call convert_ring2nest(nside_in, map_in)
    call sub_udgrade_nest(map_in, nside_in, map_out, nside_out, fmissval, pessimistic)
    call convert_nest2ring(nside_out, map_out)

    return
  end subroutine udgrade_ring_1d_KLOAD

  !=======================================================================
  subroutine udgrade_ring_nd_KLOAD(map_in, nside_in, map_out, nside_out, fmissval, pessimistic)
    !=======================================================================
    !    N dim implementation
    !=======================================================================
    USE pix_tools, only : nside2npix, ring2nest, convert_ring2nest, convert_nest2ring
    INTEGER(I4B), INTENT(IN) :: nside_in, nside_out
    REAL(KMAP), INTENT(INOUT), dimension(0:,1:), target :: map_in
    REAL(KMAP), INTENT(OUT),   dimension(0:,1:), target :: map_out
    REAL(KMAP), INTENT(IN), OPTIONAL :: fmissval
    LOGICAL(LGT) , INTENT(IN), OPTIONAL :: pessimistic

    INTEGER(I4B) :: npix_in, npix_out
    INTEGER(I4B) :: nd_in, nd_out, id
    real(KMAP), dimension(:), pointer :: p_in, p_out
    !-----------------------------------------------------------------------

    !    checks that the 2nd dimensions match
    nd_in  = size(map_in, 2)
    nd_out = size(map_out,2)
    if (nd_in /= nd_out) then
       print*,"UDGRADE_NEST: unconsistent dimension of input and output maps",nd_in,nd_out
       call fatal_error
    endif
       
    !     checks that the 2 nside are valid
    npix_out = nside2npix(nside_out)
    npix_in  = nside2npix(nside_in)
    if (npix_out < 0) then
       print*,"wrong nside_out in udgrade_ring : ", nside_out
       call fatal_error
    endif
    if (npix_in  < 0) then
       print*,"wrong nside_in  in udgrade_ring : ", nside_in
       call fatal_error
    endif

    call convert_ring2nest(nside_in, map_in)
    do id = 1, nd_in
       p_in  => map_in (0:npix_in -1, id) ! avoid actual replication in memory
       p_out => map_out(0:npix_out-1, id)
       call sub_udgrade_nest(p_in, nside_in, p_out, nside_out, fmissval, pessimistic)
    enddo
    call convert_nest2ring(nside_out, map_out)

    return
  end subroutine udgrade_ring_nd_KLOAD


  !=======================================================================
  !
  !     UDGRADE_NEST(map_in, nside_in, map_out, nside_out)
  !
  !     map_in    : input map (REAL)
  !     nside_in  : resolution parameter of the input map (npix =12*nside**2)
  !                 must be a power of 2 (INTEGER)
  !     map_out   : output map (REAL)
  !     nside_out : resolution parameter of the output map (npix =12*nside**2)
  !                 must be a power of 2 (INTEGER)
  !
  !     upgrades or degrades a map in the NESTED numbering
  !
  !     if nside_out > nside_in : upgrades
  !     value(each new sub-pixel) = value(old big pixel)
  !
  !     if nside_in < nside_out : degrades
  !     value(new big pixel) = sum values(old small pixels) / number_of_small_pixels
  !
  !     version 1.0 : Eric Hivon, TAC, September 1997
  !     version 1.1 : Dec 1997, correction of a bug in the map indices
  !     version 2.0 : Dec 2004-March 2005: SP/DP and 1D/ND overload
  !=======================================================================
  subroutine udgrade_nest_1d_KLOAD(map_in, nside_in, map_out, nside_out, fmissval, pessimistic)
    !=======================================================================
    !  1 dim. implementation
    !=======================================================================
    USE pix_tools, only : nside2npix
    INTEGER(I4B), INTENT(IN) :: nside_in, nside_out
    REAL(KMAP), INTENT(IN), dimension(0:)  :: map_in
    REAL(KMAP), INTENT(OUT), dimension(0:) :: map_out
    REAL(KMAP), INTENT(IN), OPTIONAL :: fmissval
    LOGICAL(LGT) , INTENT(IN), OPTIONAL :: pessimistic

    INTEGER(I4B) :: npix_in, npix_out
    !-----------------------------------------------------------------------

    !     checks that the 2 nside are valid
    npix_out = nside2npix(nside_out)
    npix_in  = nside2npix(nside_in)
    if (npix_out < 0) then
       print*,"wrong nside_out in udgrade_nest : ", nside_out
       call fatal_error
    endif
    if (npix_in  < 0) then
       print*,"wrong nside_in  in udgrade_nest : ", nside_in
       call fatal_error
    endif

    call sub_udgrade_nest(map_in, nside_in, map_out, nside_out, fmissval, pessimistic)

    return
  end subroutine udgrade_nest_1d_KLOAD
  !=======================================================================
  subroutine udgrade_nest_nd_KLOAD(map_in, nside_in, map_out, nside_out, fmissval, pessimistic)
    !=======================================================================
    !  N dim. implementation
    !=======================================================================
    USE pix_tools, only : nside2npix
    INTEGER(I4B), INTENT(IN) :: nside_in, nside_out
    REAL(KMAP), INTENT(IN),  dimension(0:,1:), target :: map_in
    REAL(KMAP), INTENT(OUT), dimension(0:,1:), target :: map_out
    REAL(KMAP), INTENT(IN), OPTIONAL :: fmissval
    LOGICAL(LGT) , INTENT(IN), OPTIONAL :: pessimistic

    INTEGER(I4B) :: npix_in, npix_out
    INTEGER(I4B) :: nd_in, nd_out, id
    real(KMAP), dimension(:), pointer :: p_in, p_out
    !-----------------------------------------------------------------------

    !    checks that the 2nd dimensions match
    nd_in  = size(map_in, 2)
    nd_out = size(map_out,2)
    if (nd_in /= nd_out) then
       print*,"UDGRADE_NEST: unconsistent dimension of input and output maps",nd_in,nd_out
       call fatal_error
    endif
       
    !     checks that the 2 nside are valid
    npix_out = nside2npix(nside_out)
    npix_in  = nside2npix(nside_in)
    if (npix_out < 0) then
       print*,"UDGRADE_NEST: wrong nside_out: ", nside_out
       call fatal_error
    endif
    if (npix_in  < 0) then
       print*,"UDGRADE_NEST: wrong nside_in: ", nside_in
       call fatal_error
    endif

    do id = 1, nd_in
       p_in  => map_in (0:npix_in -1, id) ! avoid actual replication in memory
       p_out => map_out(0:npix_out-1, id)
       call sub_udgrade_nest(p_in, nside_in, p_out, nside_out, fmissval, pessimistic)
    enddo


    return
  end subroutine udgrade_nest_nd_KLOAD
