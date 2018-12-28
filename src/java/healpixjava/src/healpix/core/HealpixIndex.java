
/** Generic healpix routines but tied to a given NSIDE in the constructor 
   Java version of some healpix routines from DSRI
   in java everthing must be in a class - no functions floating about.
   Original algorithms Eric Hivon and Krzysztof M. Gorski.
   This code written by William O'Mullane
 */
package healpix.core;

import java.awt.*;
import edu.jhu.htm.core.Vector3d;



public class HealpixIndex {
    public static final double piover2 = 1.57079632679489661923;
    public static final double pi = 3.1415926535897932384626434;
    public static final double twopi = 6.283185307179586476925287;
    public static final int ns_max = 8192;
	public  static final double z0 =  0.66666666666666666666666; // 2/3
    private static final int x2pix[] = new int [128];
    private static final int y2pix[] = new int [128];
    private static final int pix2x[] = new int [1024];
    private static final int pix2y[] = new int [1024];
    private static int ix;
    private static int iy;

	public int nside = 1024;
	public int nl2 , nl3,nl4,npface,npix,ncap;
	public double fact1,fact2;
   
	protected void init() {
		if (pix2x[1023] <= 0) mkpix2xy();
		if (y2pix[127]==0) mkxy2pix();
		nl2 = 2*nside;
		nl3 = 3*nside;
		nl4 = 4*nside;
		npface = (int)Math.pow(nside,2);
	    ncap = 2*nside*(nside-1) ;// points in each polar cap, =0 for nside =1
		npix = (int)(12 * Math.pow(nside,2));
	   fact1 = 1.50*nside;
	   fact2 = 3.0*npix;
	}

	public HealpixIndex() {
		init();
	}
	public HealpixIndex(int nside ) throws Exception{
		if (nside > ns_max || nside < 1){ 
	   		throw new Exception("nsides must be between 1 and "+ns_max);
		}	
		this.nside = nside;
		init();
	} 
    /**
       Initialize pix2x and pix2y 
       constructs the array giving x and y in the face from pixel number
       for the nested (quad-cube like) ordering of pixels
       the bits corresponding to x and y are interleaved in the pixel number
       one breaks up the pixel number by even and odd bits 
     */
    protected static void mkpix2xy() {
	   int kpix, jpix, ix, iy, ip, id;

	   for (kpix=0;kpix<=1023;kpix++) {          // pixel number
	      jpix = kpix;
	      ix = 0;
	      iy = 0;
	      ip = 1               ;// bit position (in x and y)
	      while (jpix!=0) { // go through all the bits;
	    	id = jpix%2  ;// bit value (in kpix), goes in ix
	    	jpix = jpix/2;
	    	ix = id*ip+ix;

	    	id = jpix%2  ;// bit value (in kpix), goes in iy
	    	jpix = jpix/2;
	    	iy = id*ip+iy;

	    	ip = 2*ip         ;// next bit (in x and y)
	      };
	      pix2x[kpix] = ix     ;// in 0,31
	      pix2y[kpix] = iy     ;// in 0,31
	   };
	
    }
    
    /**
       Initialize x2pix and y2pix 
     */
    protected static void mkxy2pix() {
	int j,k,id,ip;
	for (int i = 0 ; i < 128 ; i++) {
	   j=i;
	   k =0;
	   ip=1;
	   while (j != 0 ) {
		id= j%2;	
		j = j/2;
		k = ip*id+k;
		ip = ip*4;
	   }
	   x2pix[i]=k;
	   y2pix[i]= 2*k;
	}
    }
    
    /**
       renders the pixel number ipix (NESTED scheme) for a pixel which contains
       a point on a sphere at coordinates theta and phi, given the map
       resolution parametr nside
       the computation is made to the highest resolution available (nside=8192)
       and then degraded to that required (by integer division)
       this doesn't cost more, and it makes sure
       that the treatement of round-off will be consistent
       for every resolution
     */
    public int ang2pix_nest(double theta, double phi) throws Exception {
	int ipix;
	double  z, za, tt, tp, tmp;
	int jp, jm, ifp, ifm, face_num,ix, iy; 
	int ix_low, ix_hi, iy_low, iy_hi, ipf, ntt;

      if (phi >= twopi) phi = phi - twopi;
      if (phi < 0.)    phi = phi + twopi;
	if (theta > pi || theta < 0) {
	   throw new Exception("theta must be between 0 and "+pi);
	}	
	if (phi > twopi || phi < 0) {
	   throw new Exception("phi must be between 0 and "+twopi);
	}	
	// Note excpetion thrown means method does not get further.

      z  = Math.cos(theta);
      za = Math.abs(z);
      tt = phi / piover2 ;// in [0,4]

  //System.out.println("Za:"+za +" z0:"+z0+" tt:"+tt+" z:"+z+" theta:"+theta+" phi:"+phi);
      if (za <= z0) { //Equatorial region
 // System.out.println("Equatorial !");
// (the index of edge lines increase when the longitude=phi goes up)
         jp =(int)Math.rint(ns_max*(0.50+tt-(z*0.750)));//ascending edge line index
         jm =(int)Math.rint(ns_max*(0.50+tt+(z*0.750)));//descending edge line index

//        finds the face
         ifp = jp / ns_max;  // in {0,4}
         ifm = jm / ns_max;
         if (ifp == ifm) {          // faces 4 to 7
            face_num = (ifp%4) + 4;
         } else { 
		 	if (ifp < ifm) {     // (half-)faces 0 to 3
            	   face_num = (ifp%4);
         	} else {                            // (half-)faces 8 to 11
            	   face_num = (ifm%4) + 8;
         	};
		};

         ix = (jm%ns_max);
         iy = ns_max - (jp%ns_max) - 1;
      } else { // polar region, za > 2/3

         ntt = (int)(tt);
         if (ntt >= 4) ntt = 3;
         tp = tt - ntt;
         tmp = Math.sqrt( 3.0*(1.0 - za) ) ; // in ]0,1]

//   (the index of edge lines increase when distance from the closest pole goes up)
         jp =(int)Math.rint(ns_max * tp * tmp);// line going toward the pole as phi increases
         jm =(int)Math.rint(ns_max * (1.0 - tp) * tmp ); // that one goes away of the closest pole
         jp = Math.min(ns_max-1, jp); // for points too close to the boundary
         jm = Math.min(ns_max-1, jm);

//        finds the face and pixel's (x,y)
         if (z >= 0) { //North Pole
//System.out.println("Polar z>=0 ntt:"+ntt+" tt:"+tt);
            face_num = ntt ; // in {0,3}
            ix = ns_max - jm - 1;
            iy = ns_max - jp - 1;
         } else {
//System.out.println("Polar z<0 ntt:"+ntt+" tt:"+tt);
            face_num = ntt + 8 ;// in {8,11}
            ix =  jp;
            iy =  jm;
         };
      };

      ix_low = ix%128;
      ix_hi  =     ix/128;
      iy_low = iy%128;
      iy_hi  =     iy/128;
      ipf =  (x2pix[ix_hi]+y2pix[iy_hi]) * (128 * 128)
           + (x2pix[ix_low]+y2pix[iy_low]);   // in {0, nside**2 - 1}

	  ipf = ipf / (int)Math.rint((Math.pow((ns_max/nside),2)));
 // System.out.println("ix_low:"+ix_low +" ix_hi:"+ix_hi+" iy_low:"+iy_low+" iy_hi:"+iy_hi+" ipf:"+ipf+ " face:"+face_num);

      ipix =(int) Math.rint(ipf + face_num* npface) ;   // in {0, 12*nside**2 - 1}

      return ipix;

    }
    
    /**
       Convert from pix number to angle 
       renders theta and phi coordinates of the nominal pixel center
       for the pixel number ipix (NESTED scheme)  
       given the map resolution parameter nside
     */
    public double[] pix2ang_nest(int ipix) throws Exception {

		int ipf, ip_low, ip_trunc, ip_med, ip_hi;
	 	int jrt, jr, nr, jpt, jp, kshift ;
		double z, fn, theta,phi;

		// cooordinate of the lowest corner of each face
		// add extra zero in front so array like in fortran
		int jrll[] = {0, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4};
		int jpll[] = {0, 1, 3, 5, 7, 0, 2, 4, 6, 1, 3, 5, 7};
//-----------------------------------------------------------------------
		if (ipix <0 || ipix>npix-1) throw new Exception("ipix out of range");


		fn = 1.0*nside;
		double fact1 = 1.0/(3.0*fn*fn);
		double fact2 = 2.0/(3.0*fn);

		//finds the face, and the number in the face

		int face_num = ipix/npface  ;// face number in {0,11}
		ipf = ipix%npface  ;// pixel number in the face {0,npface-1}

		//finds the x,y on the face (starting from the lowest corner)
		//from the pixel number
		ip_low = ipf%1024  ;     // content of the last 10 bits
		ip_trunc = ipf/1024 ;       // truncation of the last 10 bits
		ip_med = ip_trunc%1024;  // content of the next 10 bits
		ip_hi  =   ip_trunc/1024;   // content of the high weight 10 bits

		ix = 1024*pix2x[ip_hi] + 32*pix2x[ip_med] + pix2x[ip_low];
		iy = 1024*pix2y[ip_hi] + 32*pix2y[ip_med] + pix2y[ip_low];

		// transforms this in (horizontal, vertical) coordinates
		jrt = ix + iy  ;// 'vertical' in {0,2*(nside-1)}
		jpt = ix - iy  ;// 'horizontal' in {-nside+1,nside-1}

		// computes the z coordinate on the sphere
		jr =  jrll[face_num+1]*nside - jrt - 1	;// ring number in {1,4*nside-1}

		nr = nside	;// equatorial region (the most frequent)
		z  = (nl2-jr)*fact2;
		kshift = jr - nside% 2;
		if (jr < nside) {	  // north pole region
			nr = jr;
			z = 1.0 - nr*nr*fact1;
			kshift = 0;
		} else { 
			if (jr > nl3) { // south pole region
				nr = nl4 - jr;
				z = - 1.0 + nr*nr*fact1;
				kshift = 0;
			}
		};
		theta = Math.acos(z);

		// computes the phi coordinate on the sphere, in [0,2Pi]
		jp = (jpll[face_num+1]*nr + jpt + 1 + kshift)/2  ;// 'phi' number in the ring in {1,4*nr}
		if (jp > nl4) jp = jp - nl4;
		if (jp < 1)	jp = jp + nl4;

		phi = (jp - (kshift+1)*0.50) * (piover2 / nr);

		double[] ret = {theta,phi};
		return ret;

	
    }

    /**
       Convert from pix number to angle 
       renders theta and phi coordinates of the nominal pixel center
       for the pixel number ipix (RING scheme)  
       given the map resolution parameter nside
     */
    public double[] pix2ang_ring(int ipix) throws Exception {

	   double theta, phi;
	   int iring, iphi, ip, ipix1;
	   double fodd, hip, fihip;
	//-----------------------------------------------------------------------
	   if (ipix <0 || ipix>npix-1)  throw new Exception ("ipix out of range");

	   ipix1 = ipix + 1 ;// in {1, npix}

	   if (ipix1 <= ncap) { // North Polar cap -------------
         
	      hip   = ipix1/2.0;
	      fihip = (int)( hip );
	      iring = (int)( Math.sqrt( hip - Math.sqrt(fihip) ) ) + 1 ;// counted from North pole
	      iphi  = ipix1 - 2*iring*(iring - 1);

	      theta = Math.acos( 1.0 - Math.pow(iring,2) / fact2 );
	      phi   = ((double)(iphi) - 0.50) * pi/(2.0*iring);

	   } else {if (ipix1 <= nl2*(5*nside+1)) { // Equatorial region ------
	      ip    = ipix1 - ncap - 1;
	      iring = (int)( ip / nl4 ) + nside ;// counted from North pole
	      iphi  = (int) ip%nl4 + 1;

	      fodd  = 0.50 * (1 + ((iring+nside)%2)) ; // 1 if iring+nside is odd, 1/2 otherwise
	      theta = Math.acos( (nl2 - iring) / fact1 );
	      phi   = ((double)(iphi) - fodd) * pi /(double)nl2;

	   } else { // South Polar cap -----------------------------------
	      ip    = npix - ipix1 + 1;
	      hip   = ip/2.0;
	      fihip = (int)( hip );
	      iring = (int)( Math.sqrt( hip - Math.sqrt(fihip) ) ) + 1 ;// counted from South pole
	      iphi  = 4*iring + 1 - (ip - 2*iring*(iring-1));

	      theta = Math.acos( -1.0 + Math.pow(iring,2) / fact2 );
	      phi   = ((double)(iphi) - 0.50) * pi/(2.0*iring);

	   }};

		double[] ret = {theta,phi};
		return ret;
    }
    
    /**
       renders the pixel number ipix (RING scheme) for a pixel which contains
       a point on a sphere at coordinates theta and phi, given the map
       resolution parametr nside
       the computation is made to the highest resolution available (nside=8192)
       and then degraded to that required (by integer division)
       this doesn't cost more, and it makes sure
       that the treatement of round-off will be consistent
       for every resolution
     */
    public int ang2pix_ring(double theta, double phi) throws Exception {

	   int ipix;

	   int jp, jm, ipix1;
	   double z, za, tt, tp, tmp;
	   int ir, ip, kshift;

	//-----------------------------------------------------------------------
	   if (nside<1 || nside>ns_max)  throw new Exception ("nside out of range");
	   if (theta<0.0 || theta>pi)   throw new Exception ("theta out of range");

	   z = Math.cos( theta );
	   za = Math.abs(z);
	   if (phi >= twopi)  phi = phi - twopi;
	   if (phi < 0.)     phi = phi + twopi;
	   tt = phi / piover2  ;// in [0,4)


	   if ( za <= z0 ) { 

	      jp = (int)(nside*(0.50 + tt - z*0.750)) ;// index of  ascending edge line 
	      jm = (int)(nside*(0.50 + tt + z*0.750)) ;// index of descending edge line

	      ir = nside + 1 + jp - jm ;// in {1,2n+1} (ring number counted from z=2/3)
	      kshift = 0;
	      if (ir%2 == 0) kshift = 1 ;// kshift=1 if ir even, 0 otherwise

	      ip = (int)( ( jp+jm - nside + kshift + 1 ) / 2 ) + 1 ;// in {1,4n}
	      if (ip > nl4) ip = ip - nl4;

	      ipix1 = ncap + nl4*(ir-1) + ip ;

	   } else {

	      tp = tt - (int)(tt)      ;//MOD(tt,1.0)
	      tmp = Math.sqrt( 3.0*(1.0 - za) );

	      jp = (int)( nside * tp          * tmp ) ;// increasing edge line index
	      jm = (int)( nside * (1.0 - tp) * tmp ) ;// decreasing edge line index

	      ir = jp + jm + 1        ;// ring number counted from the closest pole
	      ip = (int)( tt * ir ) + 1 ;// in {1,4*ir}
	      if (ip > 4*ir) ip = ip - 4*ir;

	      ipix1 = 2*ir*(ir-1) + ip;
	      if (z <= 0.0) {
	         ipix1 = npix - 2*ir*(ir+1) + ip;
	      };

	   };

	   ipix = ipix1 - 1 ;// in {0, npix-1}

	   return ipix;

    }
    
    /**
       performs conversion from NESTED to RING pixel number
     */
    public int nest2ring(int ipnest) throws Exception {

	   int  ipring;

	   int face_num, n_before;
	   int ipf, ip_low, ip_trunc, ip_med, ip_hi;
	   int ix, iy, jrt, jr, nr, jpt, jp, kshift;

	// coordinate of the lowest corner of each face
	// 0 added in front because the java array is zero offset
	   int jrll[] = { 0 ,2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4};
	   int jpll[] = { 0, 1, 3, 5, 7, 0, 2, 4, 6, 1, 3, 5, 7};
	//-----------------------------------------------------------------------
	   if (ipnest<0 || ipnest>npix-1)  throw new Exception ("ipnest out of range");

		// initiates the array for the pixel number -> (x,y) mapping
	   if (pix2x[1023] <= 0) mkpix2xy();

		// finds the face, and the number in the face
	   face_num = ipnest/npface  ;// face number in {0,11}
	   ipf = ipnest%npface  ;// pixel number in the face {0,npface-1}

		// finds the x,y on the face (starting from the lowest corner)
		// from the pixel number
	   ip_low = ipf%1024       ;// content of the last 10 bits
	   ip_trunc =   ipf/1024        ;// truncation of the last 10 bits
	   ip_med = ip_trunc%1024  ;// content of the next 10 bits
	   ip_hi  =     ip_trunc/1024   ;// content of the high weight 10 bits

	   ix = 1024*pix2x[ip_hi] + 32*pix2x[ip_med] + pix2x[ip_low];
	   iy = 1024*pix2y[ip_hi] + 32*pix2y[ip_med] + pix2y[ip_low];

		// transforms this in (horizontal, vertical) coordinates
	   jrt = ix + iy  ;// 'vertical' in {0,2*(nside-1)}
	   jpt = ix - iy  ;// 'horizontal' in {-nside+1,nside-1}

		// computes the z coordinate on the sphere
	   jr =  jrll[face_num+1]*nside - jrt - 1   ;// ring number in {1,4*nside-1}

	   nr = nside                  ;// equatorial region (the most frequent)
	   n_before = ncap + nl4 * (jr - nside);
	   kshift = (jr - nside)% 2;
	   if (jr < nside) {     // north pole region
	      nr = jr;
	      n_before = 2 * nr * (nr - 1);
	      kshift = 0;
	   } else { if (jr > nl3) { // south pole region
	      nr = nl4 - jr;
	      n_before = npix - 2 * (nr + 1) * nr;
	      kshift = 0;
	   }};

		// computes the phi coordinate on the sphere, in [0,2Pi]
	   jp = (jpll[face_num+1]*nr + jpt + 1 + kshift)/2  ;// 'phi' number in the ring in {1,4*nr}

	   if (jp > nl4) jp = jp - nl4;
	   if (jp < 1)   jp = jp + nl4;

	   ipring = n_before + jp - 1 ;// in {0, npix-1}
	   return ipring;
    }
    
    /**
       performs conversion from  RING to NESTED pixel number
     */
    public int ring2nest(int ipring) throws Exception {
	   int  ipnest;

	   double fihip, hip;
	   int  ip, iphi, ipt, ipring1;
	   int kshift, face_num=0, nr;
	   int irn, ire, irm, irs, irt, ifm , ifp;
	   int ix, iy, ix_low, ix_hi, iy_low, iy_hi, ipf;

	// coordinate of the lowest corner of each face
	// 0 added in front because the java array is zero offset
	   int jrll[] = { 0, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4};
	   int jpll[] = { 0, 1, 3, 5, 7, 0, 2, 4, 6, 1, 3, 5, 7};
	//-----------------------------------------------------------------------
	   if (ipring <0 || ipring>npix-1) throw new Exception ("ipring out of range");
	   if (x2pix[127] <= 0) mkxy2pix();

	   ipring1 = ipring + 1;

		// finds the ring number, the position of the ring and the face number
	   if (ipring1 <= ncap) {  // North Pole
	      hip   = ipring1/2.0;
	      fihip = Math.rint ( hip );
	      irn   = (int)Math.floor(Math.sqrt(hip - Math.sqrt(fihip) ) ) + 1 ;// counted from North pole
	      iphi  = ipring1 - 2*irn*(irn - 1);
	      kshift = 0;
	      nr = irn                  ;// 1/4 of the number of points on the current ring
	      face_num = (iphi-1) / irn ;// in {0,3}

	   } else {if (ipring1 <= nl2*(5*nside+1)) { // Equatorial

	      ip    = ipring1 - ncap - 1;
	      irn   = (int)Math.floor( ip / nl4 ) + nside               ;// counted from North pole
	      iphi  = ip%nl4 + 1;

	      kshift  = (irn+nside)%2; // 1 if irn+nside is odd, 0 otherwise
	      nr = nside;
	      ire =  irn - nside + 1 ;// in {1, 2*nside +1}
	      irm =  nl2 + 2 - ire;
	      ifm = (iphi - ire/2 + nside -1) / nside ;// face boundary
	      ifp = (iphi - irm/2 + nside -1) / nside;
	      if (ifp == ifm) {          // faces 4 to 7
	         face_num = ifp%4 + 4;
	      } else { 
			if (ifp + 1 == ifm) { // (half-)faces 0 to 3
	         face_num = ifp;
	      	} else { 
				if (ifp - 1 == ifm) { // (half-)faces 8 to 11
	         		face_num = ifp + 7;
				}
	      	};
		 }

	   } else { // South

	      ip    = npix - ipring1 + 1;
	      hip   = ip/2.0;
	      fihip = Math.rint ( hip );
	      irs   = (int)Math.floor( Math.sqrt( hip - Math.sqrt(fihip) ) ) + 1  ;// counted from South pole
	      iphi  = 4*irs + 1 - (ip - 2*irs*(irs-1));

	      kshift = 0;
	      nr = irs;
	      irn   = nl4 - irs;
	      face_num = (iphi-1) / irs + 8 ;// in {8,11}

	   }};
      
		
		// finds the (x,y) on the face
	   irt =   irn  - jrll[face_num+1]*nside + 1       ;// in {-nside+1,0}
	   ipt = 2*iphi - jpll[face_num+1]*nr - kshift - 1 ;// in {-nside+1,nside-1}
	   if (ipt >= nl2) ipt = ipt - 8*nside;// for the face #4

	   ix =  (ipt - irt ) / 2;
	   iy = -(ipt + irt ) / 2;

	   //  System.out.println("face:"+face_num+" irt:"+irt+" ipt:"+ipt+" ix:"+ix+" iy:"+iy);
	   ix_low = ix%128;
	   ix_hi  = ix/128;
	   iy_low = iy%128;
	   iy_hi  = iy/128;

	   ipf =  (x2pix[ix_hi ]+y2pix[iy_hi ]) * (128 * 128)
	        + (x2pix[ix_low]+y2pix[iy_low])        ;// in {0, Math.pow(nside,2) - 1}
         

	   ipnest = (int)(ipf + face_num* npface);// in {0, 12*Math.pow(nside,2) - 1}

	   return ipnest;
    }
    
    /**
       Convert from pix number to x,y inside a given face.  0,0 is the lower right corner of the face.
     */
    public Point pix2xy_nest(int ipix) throws Exception {
		if (ipix <0 || ipix>=npface)
			throw new Exception ("ipix out of range");
		return pix2xy_nestface(ipix);

    }
    
    /**
       Convert from a x,y in a given face to a pix number.
     */
    public int xy2pix_nest(int ix, int iy, int face) throws Exception {
		if (nside<1 || nside>ns_max)  throw new Exception ("nside out of range");
       if (ix<0 || ix>(nside-1))      throw new Exception ("ix out of range");
       if (iy<0 || iy>(nside-1))      throw new Exception ("iy out of range");

		int ipf,ipix;
		ipf = xy2pix_nest(ix,iy);
       ipix = (int)(ipf + face* npface) ;// in {0, 12*(nside^2)-1}
       return ipix;
    }
    
    /**
       Convert from a point in a given face to a pix number. Convenience method just unpacks the point to x and y and calls the other xy2pix_nest method.
     */
    public int xy2pix_nest(Point p, int face) throws Exception {
		return xy2pix_nest(p.x,p.y,face);
    }
    
    /**
       Convert from a x,y in a given face to a pix number in a face withour offset.
     */
    public int xy2pix_nest(int ix, int iy) throws Exception {
		int ix_low, ix_hi, iy_low, iy_hi, ipf;
       if (x2pix[127] <= 0) mkxy2pix();

       ix_low = ix%128;
       ix_hi  =     ix/128;
       iy_low = iy%128;
       iy_hi  =     iy/128;

       ipf =  (x2pix[ix_hi ]+y2pix[iy_hi ]) * (128 * 128)
            + (x2pix[ix_low]+y2pix[iy_low]);
		return ipf;
    }
    
    /**
       Convert from pix number to x,y inside a given face.  0,0 is the lower right corner of the face.
     */
    public Point pix2xy_nestface(int ipix) throws Exception {
		if (pix2x[1023] <= 0) mkpix2xy();
		int ip_low, ip_trunc, ip_med, ip_hi,ix,iy;
		ip_low = ipix%1024       ;// content of the last 10 bits
		ip_trunc =   ipix/1024        ;// truncation of the last 10 bits
		ip_med = ip_trunc%1024  ;// content of the next 10 bits
		ip_hi  =     ip_trunc/1024   ;// content of the high weight 10 bits
		ix = 1024*pix2x[ip_hi] + 32*pix2x[ip_med] + pix2x[ip_low];
		iy = 1024*pix2y[ip_hi] + 32*pix2y[ip_med] + pix2y[ip_low];
//System.out.println("ip_low:"+ip_low+" iptrun:"+ip_trunc+" ip_med"+ip_med+" ip_hi"+ip_hi+" ix:"+ix+" iy:"+iy);
		return new Point(ix,iy);
    }

	public double[] integration_limits_in_costh(int i_th){
	
	   double a, ab, b, r_n_side;

//  integration limits in cos(theta) for a given ring i_th
// i > 0 !!!

 	r_n_side = 1.0 * nside;
	if( i_th <= nside ) { 
  		ab = 1.0 - (Math.pow(i_th ,2.0)/3.0) / (double)npface;
  		b =  1.0 - (Math.pow((i_th-1),2.0)/ 3.0) / (double)npface;
		if( i_th == nside ) {
  			a =   2.0*(nside - 1.0)/3.0/r_n_side;
		} else {
  			a =  1.0 - Math.pow((i_th+1),2)/3.0/ (double)npface;
		};

	} else {
		if( i_th < nl3 ){
  			ab =  2.0*(2*nside - i_th    )/3.0/r_n_side;
  			b =   2.0*(2*nside - i_th + 1)/3.0/r_n_side;
  			a =   2.0*(2*nside - i_th - 1)/3.0/r_n_side;
 		} else {
			if( i_th == nl3 ){
  				b = 2.0*(-nside + 1)/3.0/r_n_side;
			} else {
				b = -1.0 + Math.pow((4*nside - i_th + 1),2)/3.0/(double)npface;
			}

			a =  -1.0 + Math.pow((nl4 - i_th - 1),2)/3.0/(double)npface;
  			ab = -1.0 + Math.pow((nl4 - i_th),2)/3.0/(double)npface;
		}

 	}
	//    END integration limits in cos(theta)
	double[] ret ={b,ab,a}; 
	return ret; 
}
	/** calculate the points of crosing for a given theata on the boundaries of the pixel - returns the left and right phi crosings */
	public double[] pixel_boundaries(double i_th,double  i_phi, int i_zone, double cos_theta){
        double sq3th,factor,jd,ju,ku,kd,phi_l, phi_r;
        double r_n_side= 1.0*nside;

	
	//	HALF a pixel away from both poles
	if( Math.abs(cos_theta) >= 1.0-1.0/3.0/(double)npface ){
		phi_l = i_zone  * piover2;
		phi_r = (i_zone+1) * piover2;
		double[] ret = {phi_l,phi_r};
		return ret;
	}		
//-------
//	NORTH POLAR CAP	
	if( 1.50*cos_theta >= 1.0){
        sq3th = Math.sqrt(3.0* (1.0 - cos_theta));
        factor = 1.0/r_n_side/sq3th;
        jd = (float)(i_phi);
        ju = jd - 1;
        ku = (float)(i_th - i_phi);
        kd = ku + 1;
//System.out.println(" cos_theta:"+cos_theta+" sq3th:"+sq3th+" factor:"+factor+" jd:"+jd+" ju:"+ju+" ku:"+ku+" kd:"+kd+ " izone:"+i_zone);
        phi_l = piover2 * (Math.max((ju*factor),
									(1.0 -  (kd*factor)))  
							+ i_zone);
        phi_r = piover2 * (Math.min((1.0 - (ku*factor)),
									(jd*factor))
							+ i_zone);

	} else {
		if( -1.0 < 1.50*cos_theta ){
//-------
//-------
//	EQUATORIAL ZONE
	    double cth34   = 0.50 * ( 1.0 - 1.50*cos_theta );
		double cth34_1 = cth34 + 1.0;
 		int modfactor = (int)(nside + (i_th% 2));
	
		jd =         i_phi - ( modfactor - i_th )/2.0;
		ju = jd - 1;
		ku =                 ( modfactor + i_th )/2.0 - i_phi;
		kd = ku + 1;

		phi_l = piover2 * (  Math.max( (cth34_1  - (kd/r_n_side)), 
	                                    (-cth34 + (ju/r_n_side)))+ i_zone)   ;

		phi_r = piover2 * (  Math.min( (cth34_1 - (ku/r_n_side)),  
	                                    (-cth34 + (jd/r_n_side))) +i_zone )  ;
//-------
//-------
//	SOUTH POLAR CAP

		} else {
		sq3th = Math.sqrt(3.0*(1.0+cos_theta));
		factor = 1.0/r_n_side/sq3th;
		int ns2 = 2*nside;

		jd =            i_th - ns2 + i_phi ;
		ju = jd -1;
		ku =                   ns2 - i_phi ;
		kd = ku + 1;

		phi_l = piover2 * (Math.max((1.0 - (ns2-ju)*factor), 
	                                    ((ns2-kd)*factor )) +i_zone);

		phi_r = piover2 * (Math.min((1.0 - (ns2-jd)*factor),  
	                                     ((ns2-ku)*factor)) +i_zone); 
	    }// of SOUTH POLAR CAP
	}
//    and that's it
//System.out.println(" nside:"+nside+" i_th:"+i_th+" i_phi:"+i_phi+" izone:"+i_zone+" cos_theta:"+cos_theta+" phi_l:"+phi_l+" phi_r:"+phi_r);

		double[] ret = {phi_l,phi_r};
		return ret;
		
	}
    
	/** return ring number for given pix in ring scheme */ 
    public int ring(int ipix) throws Exception {
	   int iring=0;
	   int ipix1 = ipix + 1 ;// in {1, npix}
		int ip;
		double hip,fihip=0;
	   if (ipix1 <= ncap) { // North Polar cap -------------
	      hip   = ipix1/2.0;
	      fihip = (int)( hip );
	      iring = (int)( Math.sqrt( hip - Math.sqrt(fihip) ) ) + 1 ;// counted from North pole
	   } else {
			if (ipix1 <= nl2*(5*nside+1)) { // Equatorial region ------
	      		ip    = ipix1 - ncap - 1;
	      		iring = (int)( ip / nl4 ) + nside ;// counted from North pole
	   		} else { // South Polar cap -----------------------------------
	      		ip    = npix - ipix1 + 1;
	     		hip   = ip/2.0;
	      		fihip = (int)( hip );
	      		iring = (int)( Math.sqrt( hip - Math.sqrt(fihip) ) ) +1 ;// counted from South pole
				iring = nl4 - iring;
			}
		};
		return iring;
	}

	public static  Vector3d vector(double theta, double phi) {
        double x,y,z;
        x = 1 * Math.sin(theta) * Math.cos(phi);
        y = 1 * Math.sin(theta) * Math.sin(phi);
        z = 1 * Math.cos(theta);
        return new Vector3d(x,y,z);
	}

	public static double[] ang(Vector3d vec) {
		double theta = Math.acos(vec.z());	
		double phi = Math.asin(vec.y()/Math.sin(theta));	
		double[] ret = { theta,phi};
		return ret;
		
	}

	/** unit vector to pix number */
    public int vec2pix_nest(Vector3d vec) throws Exception   {
		double[] angs = ang(vec);
		return ang2pix_nest(angs[0],angs[1]);
	}

    public int vec2pix_ring(Vector3d vec) throws Exception {
		double[] angs = ang(vec);
		return ang2pix_ring(angs[0],angs[1]);
	}

	public Vector3d pix2vec_nest(int pix) throws Exception {
		double[] angs = pix2ang_nest(pix);
		return vector(angs[0],angs[1]);
	}

	public Vector3d pix2vec_ring(int pix) throws Exception {
		double[] angs = pix2ang_ring(pix);
		return vector(angs[0],angs[1]);
	}


	public Vector3d[] corners_nest(int pix, int step ) throws Exception{
        int pixr = nest2ring(pix);
        return corners_ring(pixr,step);
    }

	/** return set of points along the boundary of the given pixel step 1
		gives 4 points on the corners */

	    public Vector3d[] corners_ring(int pix, int step)  throws Exception{
		int nPoints = step*2 +2;
		Vector3d[] points = new Vector3d[nPoints];
		double[] p0 = pix2ang_ring(pix);
		double cos_theta = Math.cos(p0[0]);
		double theta = p0[0];
		double phi=p0[1];

		int i_zone = (int)(phi/piover2);
		int ringno = ring(pix);
		int i_phi_count = Math.min(ringno, Math.min(nside, (nl4) -ringno));
		int i_phi = 0;
		double phifac = piover2/i_phi_count ;
		if (ringno >=nside && ringno <= nl3) {
		    // adjust by 0.5 for odd numbered rings in equatorial since
		    // they start out of phase by half phifac.
		    i_phi=(int)(phi/phifac + ((ringno%2)/2.0))+1;
		} else {
		    i_phi=(int)(phi/phifac )+1;
		}
		// adjust for zone offset
		i_phi=i_phi - (i_zone * i_phi_count);
		int spoint = (int)(nPoints/2);
		// get north south middle - mddle should match theata !
		double[] nms = integration_limits_in_costh(ringno);
		double ntheta = Math.acos(nms[0]);
		double stheta = Math.acos(nms[2]);
		double[] philr=pixel_boundaries(ringno,i_phi,i_zone,nms[0]);
		if (i_phi > (i_phi_count/2)) {
		    points[0] = vector(ntheta,philr[1]);
		}else {
		    points[0] = vector(ntheta,philr[0]);
		}
		philr=pixel_boundaries(ringno,i_phi,i_zone,nms[2]);
		if (i_phi > (i_phi_count/2)) {
		    points[spoint] = vector(stheta,philr[1]);
		}else {
		    points[spoint] = vector(stheta,philr[0]);
        }
	if (step ==1) {
	    double mtheta = Math.acos(nms[1]);
            philr=pixel_boundaries(ringno,i_phi,i_zone,nms[1]);
            points[1] = vector(mtheta,philr[0]);
            points[3] = vector(mtheta,philr[1]);
	} else {
     	   double cosThetaLen = nms[2] - nms[0];
       	   double cosThetaStep = (cosThetaLen/(step+1)); // skip North and south
       	   for (int p=1; p <= step; p++) {
            /* Integrate points along the sides */
            cos_theta = nms[0] + (cosThetaStep*p);
            theta=Math.acos(cos_theta);
            philr=pixel_boundaries(ringno,i_phi,i_zone,cos_theta);
            points[p] = vector(theta,philr[0]);
            points[nPoints-p] = vector(theta,philr[1]);
          }
	}

        return points;
    }

}
