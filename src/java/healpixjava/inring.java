  public void  in_ring (nside, iz, phi0, dphi, listir, nir, nest);
    //=======================================================================;
    //     returns the list of pixels in RING or NESTED scheme (listir);
    //     and their number (nir);
    //     with latitude in [phi0-dphi, phi0+dphi] on the ring ir;
    //     (in {1,4*nside-1})
    //     the pixel id-numbers are in {0,12*nside^2-1}
    //     the indexing is RING, unless NEST is set to 1;
    //=======================================================================;
    int  nside, iz;
    int  nir;
    double  phi0, dphi;
    int  listir;
    int  nest;

//     boolean conservative = true;
    boolean conservative = false;
    boolean take_all, to_top, do_ring;

    int  ip_low, ip_hi, i, in, inext, diff;
    int  npix, nr, nir1, nir2, ir, ipix1, ipix2, kshift, ncap;
    double  phi_low, phi_hi, shift;
    //=======================================================================;

    take_all = false;
    to_top   = false;
    do_ring  = true;
    if (present(nest)){
       do_ring = (nest == 0);
    }
    npix = 12 * nside * nside;
    ncap  = 2*nside*(nside-1) // number of pixels in the north polar cap;
    listir = -1;
    nir = 0;

    phi_low = (phi0 - dphi)%( twopi);
    phi_hi  = (phi0 + dphi)%( twopi);
    if (Math.abs(dphi-PI) < 1.0e-6_dp) take_all = true;

    //     ------------ identifies ring number --------------;
    if (iz >= nside .and. iz <= 3*nside){ // equatorial region
       ir = iz - nside + 1  // in {1, 2*nside + 1}
       ipix1 = ncap + 4*nside*(ir-1) //  lowest pixel number in the ring;
       ipix2 = ipix1 + 4*nside - 1   // highest pixel number in the ring;
       kshift = (ir)%(2);
       nr = nside*4;
    } else {
       if (iz < nside){       //    north pole
          ir = iz;
          ipix1 = 2*ir*(ir-1)        //  lowest pixel number in the ring;
          ipix2 = ipix1 + 4*ir - 1   // highest pixel number in the ring;
       } else {                          //    south pole
          ir = 4*nside - iz;
          ipix1 = npix - 2*ir*(ir+1) //  lowest pixel number in the ring;
          ipix2 = ipix1 + 4*ir - 1   // highest pixel number in the ring;
       }
       nr = ir*4;
       kshift = 1;
    }

    //     ----------- constructs the pixel list --------------;
    if (take_all){
       nir    = ipix2 - ipix1 + 1;
       if (do_ring){
          listir(0:nir-1) = (/ (i, i=ipix1,ipix2) /);
       } else {
           ring2nest(nside, ipix1, in);
          listir(0) = in;
          do i=1,nir-1;
              next_in_line_nest(nside, in, inext);
             in = inext;
             listir(i) = in;
          }//enddo
       }
       return;
    }

    shift = kshift * 0.5_dp;
    if (conservative){
       // conservative : include every intersected pixels,;
       // even if pixel CENTER is not in the range [phi_low, phi_hi];
       ip_low = nint (nr * phi_low / TWOPI - shift);
       ip_hi  = nint (nr * phi_hi  / TWOPI - shift);
       ip_low = modulo (ip_low, nr) // in {0,nr-1}
       ip_hi  = modulo (ip_hi , nr) // in {0,nr-1}
    } else {
       // strict : include only pixels whose CENTER is in [phi_low, phi_hi];
       ip_low = ceiling (nr * phi_low / TWOPI - shift);
       ip_hi  = floor   (nr * phi_hi  / TWOPI - shift);
//        if ((ip_low - ip_hi == 1) .and. (dphi*nr < PI)){ // EH, 2004-06-01
       diff = modulo(ip_low - ip_hi, nr) // in {-nr+1, nr-1} or {0,nr-1} ???
       if (diff < 0) diff = diff + nr    // in {0,nr-1}
       if ((diff == 1) .and. (dphi*nr < PI)){
          // the interval is so small (and away from pixel center);
          // that no pixel is included in it;
          nir = 0;
          return;
       }
//        ip_low = min(ip_low, nr-1) //  EH, 2004-05-28;
//        ip_hi  = max(ip_hi , 0   );
       if (ip_low >= nr) ip_low = ip_low - nr;
       if (ip_hi  <  0 ) ip_hi  = ip_hi  + nr;
    }
    //;
    if (ip_low > ip_hi) to_top = true;
    ip_low = ip_low + ipix1;
    ip_hi  = ip_hi  + ipix1;

    if (to_top){
       nir1 = ipix2 - ip_low + 1;
       nir2 = ip_hi - ipix1  + 1;
       nir  = nir1 + nir2;
       if (do_ring){
          listir(0:nir1-1)   = (/ (i, i=ip_low, ipix2) /);
          listir(nir1:nir-1) = (/ (i, i=ipix1, ip_hi) /);
       } else {
           ring2nest(nside, ip_low, in);
          listir(0) = in;
          do i=1,nir-1;
              next_in_line_nest(nside, in, inext);
             in = inext;
             listir(i) = in;
          }//enddo
       }
    } else {
       nir = ip_hi - ip_low + 1;
       if (do_ring){
          listir(0:nir-1) = (/ (i, i=ip_low, ip_hi) /);
       } else {
           ring2nest(nside, ip_low, in);
          listir(0) = in;
          do i=1,nir-1;
              next_in_line_nest(nside, in, inext);
             in = inext;
             listir(i) = in;
          }//enddo
       }
    }

    return;
  }// in_ring

