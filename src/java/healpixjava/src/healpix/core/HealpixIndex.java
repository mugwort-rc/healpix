/*
 * HEALPix Java code supported by the Gaia project.
 * Copyright (C) 2006-2011 Gaia Data Processing and Analysis Consortium
 * 
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */
package healpix.core;

import healpix.core.base.BitManipulation;
import healpix.core.base.HealpixException;
import healpix.tools.Constants;
import healpix.tools.SpatialVector;

import java.awt.Point;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Vector;

import javax.vecmath.Vector3d;

/**
 * Generic healpix routines but tied to a given NSIDE in the constructor Java
 * version of some healpix routines from DSRI in java everthing must be in a
 * class - no functions floating about. Original algorithms Eric Hivon and
 * Krzysztof M. Gorski. This code written by William O'Mullane extended by
 * Emmanuel Joliet with some methods added from pix_tools F90 code port to Java.
 * 
 * @author William O'Mullane, extended by Emmanuel Joliet
 * @version $Id: HealpixIndex.java,v 1.1.2.2 2009/08/03 16:25:20 healpix Exp $
 */

public class HealpixIndex implements Serializable {
	/**
	 * Default serial version
	 */
	private static final long serialVersionUID = 1L;
    public static final String REVISION =
        "$Id: HealpixIndex.java,v 1.1.2.2 2009/08/03 16:25:20 healpix Exp $";
	/** The Constant ns_max. */
	public static final int ns_max = 8192;
	
	public static long[] nsidelist = { 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048,
			4096, 8192, 16384, 32768, 65563, 131072, 262144, 524288,
			1048576 };


	/** The Constant z0. */
	public static final double z0 = Constants.twothird; // 2/3

	/** The Constant x2pix. */
	private static final int x2pix[] = new int[128];

	/** The Constant y2pix. */
	private static final int y2pix[] = new int[128];

	/** The Constant pix2x. */
	private static final int pix2x[] = new int[1024];

	/** The Constant pix2y. */
	private static final int pix2y[] = new int[1024];

	/** The ix. */
	private static int ix;

	/** The iy. */
	private static int iy;

	/** The nside. */
	public int nside = 1024;

	/** The ncap. */
	public int nl2, nl3, nl4, npface, npix, ncap;

	/** The fact2. */
	public double fact1, fact2;

	/** The bm. */
	transient private BitManipulation bm = new BitManipulation();

	/**
	 * Inits the.
	 */
	protected void init() {
		if ( pix2x[1023] <= 0 )
			mkpix2xy();
		if ( y2pix[127] == 0 )
			mkxy2pix();
		nl2 = 2 * nside;
		nl3 = 3 * nside;
		nl4 = 4 * nside;
		npface = (int) Math.pow(nside, 2);
		ncap = 2 * nside * ( nside - 1 );// points in each polar cap, =0 for
		// nside =1
		npix = (int) ( 12 * Math.pow(nside, 2) );
		fact1 = 1.50 * nside;
		fact2 = 3.0 * npface;
	}

	/**
	 * Default constructor nside = 1024.
	 */
	public HealpixIndex() {
		init();
	}

	/**
	 * Construct healpix routines tied to a given nside
	 * 
	 * @param nside
	 *            resolution number
	 * @throws Exception
	 */
	public HealpixIndex(int nside) throws Exception {
		if ( nside > ns_max || nside < 1 ) {
			throw new Exception("nsides must be between 1 and " + ns_max);
		}
		this.nside = nside;
		init();
	}

	/**
	 * Initialize pix2x and pix2y constructs the array giving x and y in the
	 * face from pixel number for the nested (quad-cube like) ordering of pixels
	 * the bits corresponding to x and y are interleaved in the pixel number one
	 * breaks up the pixel number by even and odd bits
	 */
	protected static void mkpix2xy() {
		int kpix, jpix, ix, iy, ip, id;

		for ( kpix = 0; kpix <= 1023; kpix++ ) { // pixel number
			jpix = kpix;
			ix = 0;
			iy = 0;
			ip = 1;// bit position (in x and y)
			while ( jpix != 0 ) { // go through all the bits;
				id = jpix % 2;// bit value (in kpix), goes in ix
				jpix = jpix / 2;
				ix = id * ip + ix;

				id = jpix % 2;// bit value (in kpix), goes in iy
				jpix = jpix / 2;
				iy = id * ip + iy;

				ip = 2 * ip;// next bit (in x and y)
			};
			pix2x[kpix] = ix;// in 0,31
			pix2y[kpix] = iy;// in 0,31
		};

	}

	/**
	 * Initialize x2pix and y2pix
	 */
	protected static void mkxy2pix() {
		int j, k, id, ip;
		for ( int i = 0; i < 128; i++ ) {
			j = i;
			k = 0;
			ip = 1;
			while ( j != 0 ) {
				id = j % 2;
				j = j / 2;
				k = ip * id + k;
				ip = ip * 4;
			}
			x2pix[i] = k;
			y2pix[i] = 2 * k;
		}
	}

	/**
	 * renders the pixel number ipix (NESTED scheme) for a pixel which contains
	 * a point on a sphere at coordinates theta and phi, given the map
	 * resolution parametr nside the computation is made to the highest
	 * resolution available (nside=8192) and then degraded to that required (by
	 * integer division) this doesn't cost more, and it makes sure that the
	 * treatement of round-off will be consistent for every resolution
	 * 
	 * @param theta
	 *            angle (along meridian), in [0,Pi], theta=0 : north pole
	 * @param phi
	 *            angle (along parallel), in [0,2*Pi]
	 * @return pixel index number
	 * @throws Exception
	 */
	public int ang2pix_nest(double theta, double phi) throws Exception {
		int ipix;
		double z, za, tt, tp, tmp;
		int jp, jm, ifp, ifm, face_num, ix, iy;
		int ix_low, ix_hi, iy_low, iy_hi, ipf, ntt;

		if ( phi >= Constants.twopi )
			phi = phi - Constants.twopi;
		if ( phi < 0. )
			phi = phi + Constants.twopi;
		if ( theta > Constants.PI || theta < 0 ) {
			throw new Exception("theta must be between 0 and " + Constants.PI);
		}
		if ( phi > Constants.twopi || phi < 0 ) {
			throw new Exception("phi must be between 0 and " + Constants.twopi);
		}
		// Note excpetion thrown means method does not get further.

		z = Math.cos(theta);
		za = Math.abs(z);
		tt = phi / Constants.piover2;// in [0,4]

		// System.out.println("Za:"+za +" z0:"+z0+" tt:"+tt+" z:"+z+"
		// theta:"+theta+" phi:"+phi);
		if ( za <= z0 ) { // Equatorial region
			// System.out.println("Equatorial !");
			// (the index of edge lines increase when the longitude=phi goes up)
			jp = (int) Math.rint(ns_max * ( 0.50 + tt - ( z * 0.750 ) ));// ascending
			// edge
			// line
			// index
			jm = (int) Math.rint(ns_max * ( 0.50 + tt + ( z * 0.750 ) ));// descending
			// edge
			// line
			// index

			// finds the face
			ifp = jp / ns_max; // in {0,4}
			ifm = jm / ns_max;
			if ( ifp == ifm ) { // faces 4 to 7
				face_num = ( ifp % 4 ) + 4;
			} else {
				if ( ifp < ifm ) { // (half-)faces 0 to 3
					face_num = ( ifp % 4 );
				} else { // (half-)faces 8 to 11
					face_num = ( ifm % 4 ) + 8;
				};
			};

			ix = ( jm % ns_max );
			iy = ns_max - ( jp % ns_max ) - 1;
		} else { // polar region, za > 2/3

			ntt = (int) ( tt );
			if ( ntt >= 4 )
				ntt = 3;
			tp = tt - ntt;
			tmp = Math.sqrt(3.0 * ( 1.0 - za )); // in ]0,1]

			// (the index of edge lines increase when distance from the closest
			// pole goes up)
			jp = (int) Math.rint(ns_max * tp * tmp);// line going toward the
			// pole as phi increases
			jm = (int) Math.rint(ns_max * ( 1.0 - tp ) * tmp); // that one goes
			// away of the
			// closest pole
			jp = Math.min(ns_max - 1, jp); // for points too close to the
			// boundary
			jm = Math.min(ns_max - 1, jm);

			// finds the face and pixel's (x,y)
			if ( z >= 0 ) { // North Pole
				// System.out.println("Polar z>=0 ntt:"+ntt+" tt:"+tt);
				face_num = ntt; // in {0,3}
				ix = ns_max - jm - 1;
				iy = ns_max - jp - 1;
			} else {
				// System.out.println("Polar z<0 ntt:"+ntt+" tt:"+tt);
				face_num = ntt + 8;// in {8,11}
				ix = jp;
				iy = jm;
			};
		};

		ix_low = ix % 128;
		ix_hi = ix / 128;
		iy_low = iy % 128;
		iy_hi = iy / 128;
		ipf = ( x2pix[ix_hi] + y2pix[iy_hi] ) * ( 128 * 128 )
				+ ( x2pix[ix_low] + y2pix[iy_low] ); // in {0, nside**2 - 1}

		ipf = ipf / (int) Math.rint(( Math.pow(( ns_max / nside ), 2) ));
		// System.out.println("ix_low:"+ix_low +" ix_hi:"+ix_hi+"
		// iy_low:"+iy_low+" iy_hi:"+iy_hi+" ipf:"+ipf+ " face:"+face_num);

		ipix = (int) Math.rint(ipf + face_num * npface); // in {0,
		// 12*nside**2 - 1}

		return ipix;

	}

	/**
	 * Convert from pix number to angle renders theta and phi coordinates of the
	 * nominal pixel center for the pixel number ipix (NESTED scheme) given the
	 * map resolution parameter nside
	 * 
	 * @param ipix
	 *            pixel index number
	 * @return double array of [theta, phi] angles
	 * @throws Exception
	 */
	public double[] pix2ang_nest(int ipix) throws Exception {

		int ipf, ip_low, ip_trunc, ip_med, ip_hi;
		int jrt, jr, nr, jpt, jp, kshift;
		double z, fn, theta, phi;

		// cooordinate of the lowest corner of each face
		// add extra zero in front so array like in fortran
		int jrll[] = { 0, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4 };
		int jpll[] = { 0, 1, 3, 5, 7, 0, 2, 4, 6, 1, 3, 5, 7 };
		// -----------------------------------------------------------------------
		if ( ipix < 0 || ipix > npix - 1 )
			throw new Exception("ipix out of range");

		fn = 1.0 * nside;
		double fact1 = 1.0 / ( 3.0 * fn * fn );
		double fact2 = 2.0 / ( 3.0 * fn );

		// finds the face, and the number in the face

		int face_num = ipix / npface;// face number in {0,11}
		ipf = ipix % npface;// pixel number in the face {0,npface-1}

		// finds the x,y on the face (starting from the lowest corner)
		// from the pixel number
		ip_low = ipf % 1024; // content of the last 10 bits
		ip_trunc = ipf / 1024; // truncation of the last 10 bits
		ip_med = ip_trunc % 1024; // content of the next 10 bits
		ip_hi = ip_trunc / 1024; // content of the high weight 10 bits

		ix = 1024 * pix2x[ip_hi] + 32 * pix2x[ip_med] + pix2x[ip_low];
		iy = 1024 * pix2y[ip_hi] + 32 * pix2y[ip_med] + pix2y[ip_low];

		// transforms this in (horizontal, vertical) coordinates
		jrt = ix + iy;// 'vertical' in {0,2*(nside-1)}
		jpt = ix - iy;// 'horizontal' in {-nside+1,nside-1}

		// computes the z coordinate on the sphere
		jr = jrll[face_num + 1] * nside - jrt - 1;// ring number in
		// {1,4*nside-1}

		nr = nside;// equatorial region (the most frequent)
		z = ( nl2 - jr ) * fact2;
		kshift = jr - nside % 2;
		if ( jr < nside ) { // north pole region
			nr = jr;
			z = 1.0 - nr * nr * fact1;
			kshift = 0;
		} else {
			if ( jr > nl3 ) { // south pole region
				nr = nl4 - jr;
				z = -1.0 + nr * nr * fact1;
				kshift = 0;
			}
		};
		theta = Math.acos(z);

		// computes the phi coordinate on the sphere, in [0,2Pi]
		jp = ( jpll[face_num + 1] * nr + jpt + 1 + kshift ) / 2;// 'phi' number
		// in the ring
		// in {1,4*nr}
		if ( jp > nl4 )
			jp = jp - nl4;
		if ( jp < 1 )
			jp = jp + nl4;

		phi = ( jp - ( kshift + 1 ) * 0.50 ) * ( Constants.piover2 / nr );

		// if (phi < 0)
		// phi += 2.0 * Math.PI; // phi in [0, 2pi]

		double[] ret = { theta, phi };
		return ret;

	}

	/**
	 * Convert from pix number to angle renders theta and phi coordinates of the
	 * nominal pixel center for the pixel number ipix (RING scheme) given the
	 * map resolution parameter nside
	 * 
	 * @param ipix
	 *            pixel index number
	 * @return double array of [theta, phi] angles
	 * @throws Exception
	 */
	public double[] pix2ang_ring(int ipix) throws Exception {

		double theta, phi;
		int iring, iphi, ip, ipix1;
		double fodd, hip, fihip;
		// -----------------------------------------------------------------------
		if ( ipix < 0 || ipix > npix - 1 )
			throw new Exception("ipix out of range");

		ipix1 = ipix + 1;// in {1, npix}

		if ( ipix1 <= ncap ) { // North Polar cap -------------

			hip = ipix1 / 2.0;
			fihip = (int) ( hip );
			iring = (int) ( Math.sqrt(hip - Math.sqrt(fihip)) ) + 1;// counted
			// from
			// North
			// pole
			iphi = ipix1 - 2 * iring * ( iring - 1 );

			theta = Math.acos(1.0 - Math.pow(iring, 2) / fact2);
			phi = ( (double) ( iphi ) - 0.50 ) * Constants.PI / ( 2.0 * iring );

		} else {
			if ( ipix1 <= nl2 * ( 5 * nside + 1 ) ) { // Equatorial region
				// ------
				ip = ipix1 - ncap - 1;
				iring = (int) ( ip / nl4 ) + nside;// counted from North pole
				iphi = (int) ip % nl4 + 1;

				fodd = 0.50 * ( 1 + ( ( iring + nside ) % 2 ) ); // 1 if
				// iring+nside
				// is odd, 1/2
				// otherwise
				theta = Math.acos(( nl2 - iring ) / fact1);
				phi = ( (double) ( iphi ) - fodd ) * Constants.PI
						/ (double) nl2;

			} else { // South Polar cap -----------------------------------
				ip = npix - ipix1 + 1;
				hip = ip / 2.0;
				fihip = (int) ( hip );
				iring = (int) ( Math.sqrt(hip - Math.sqrt(fihip)) ) + 1;// counted
				// from
				// South
				// pole
				iphi = 4 * iring + 1 - ( ip - 2 * iring * ( iring - 1 ) );

				theta = Math.acos(-1.0 + Math.pow(iring, 2) / fact2);
				phi = ( (double) ( iphi ) - 0.50 ) * Constants.PI
						/ ( 2.0 * iring );

			}
		};

		double[] ret = { theta, phi };
		return ret;
	}

	/**
	 * renders the pixel number ipix (RING scheme) for a pixel which contains a
	 * point on a sphere at coordinates theta and phi, given the map resolution
	 * parametr nside the computation is made to the highest resolution
	 * available (nside=8192) and then degraded to that required (by integer
	 * division) this doesn't cost more, and it makes sure that the treatement
	 * of round-off will be consistent for every resolution
	 * 
	 * @param theta
	 *            angle (along meridian), in [0,Pi], theta=0 : north pole
	 * @param phi
	 *            angle (along parallel), in [0,2*Pi]
	 * @return pixel index number
	 * @throws Exception
	 */
	public int ang2pix_ring(double theta, double phi) throws Exception {

		int ipix;

		int jp, jm, ipix1;
		double z, za, tt, tp, tmp;
		int ir, ip, kshift;

		// -----------------------------------------------------------------------
		if ( nside < 1 || nside > ns_max )
			throw new Exception("nside out of range");
		if ( theta < 0.0 || theta > Constants.PI )
			throw new Exception("theta out of range");

		z = Math.cos(theta);
		za = Math.abs(z);
		if ( phi >= Constants.twopi )
			phi = phi - Constants.twopi;
		if ( phi < 0. )
			phi = phi + Constants.twopi;
		tt = phi / Constants.piover2;// in [0,4)

		if ( za <= z0 ) {

			jp = (int) ( nside * ( 0.50 + tt - z * 0.750 ) );// index of
			// ascending edge
			// line
			jm = (int) ( nside * ( 0.50 + tt + z * 0.750 ) );// index of
			// descending edge
			// line

			ir = nside + 1 + jp - jm;// in {1,2n+1} (ring number counted from
			// z=2/3)
			kshift = 0;
			if ( ir % 2 == 0 )
				kshift = 1;// kshift=1 if ir even, 0 otherwise

			ip = (int) ( ( jp + jm - nside + kshift + 1 ) / 2 ) + 1;// in {1,4n}
			if ( ip > nl4 )
				ip = ip - nl4;

			ipix1 = ncap + nl4 * ( ir - 1 ) + ip;

		} else {

			tp = tt - (int) ( tt );// MOD(tt,1.0)
			tmp = Math.sqrt(3.0 * ( 1.0 - za ));

			jp = (int) ( nside * tp * tmp );// increasing edge line index
			jm = (int) ( nside * ( 1.0 - tp ) * tmp );// decreasing edge line
			// index

			ir = jp + jm + 1;// ring number counted from the closest pole
			ip = (int) ( tt * ir ) + 1;// in {1,4*ir}
			if ( ip > 4 * ir )
				ip = ip - 4 * ir;

			ipix1 = 2 * ir * ( ir - 1 ) + ip;
			if ( z <= 0.0 ) {
				ipix1 = npix - 2 * ir * ( ir + 1 ) + ip;
			};

		};

		ipix = ipix1 - 1;// in {0, npix-1}

		return ipix;

	}

	/**
	 * performs conversion from NESTED to RING pixel number
	 * 
	 * @param ipnest
	 *            pixel NEST index number
	 * @return RING pixel index number
	 * @throws Exception
	 */
	public int nest2ring(int ipnest) throws Exception {

		int ipring;

		int face_num, n_before;
		int ipf, ip_low, ip_trunc, ip_med, ip_hi;
		int ix, iy, jrt, jr, nr, jpt, jp, kshift;

		// coordinate of the lowest corner of each face
		// 0 added in front because the java array is zero offset
		int jrll[] = { 0, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4 };
		int jpll[] = { 0, 1, 3, 5, 7, 0, 2, 4, 6, 1, 3, 5, 7 };
		// -----------------------------------------------------------------------
		if ( ipnest < 0 || ipnest > npix - 1 )
			throw new Exception("ipnest out of range");

		// initiates the array for the pixel number -> (x,y) mapping
		if ( pix2x[1023] <= 0 )
			mkpix2xy();

		// finds the face, and the number in the face
		face_num = ipnest / npface;// face number in {0,11}
		ipf = ipnest % npface;// pixel number in the face {0,npface-1}

		// finds the x,y on the face (starting from the lowest corner)
		// from the pixel number
		ip_low = ipf % 1024;// content of the last 10 bits
		ip_trunc = ipf / 1024;// truncation of the last 10 bits
		ip_med = ip_trunc % 1024;// content of the next 10 bits
		ip_hi = ip_trunc / 1024;// content of the high weight 10 bits

		ix = 1024 * pix2x[ip_hi] + 32 * pix2x[ip_med] + pix2x[ip_low];
		iy = 1024 * pix2y[ip_hi] + 32 * pix2y[ip_med] + pix2y[ip_low];

		// transforms this in (horizontal, vertical) coordinates
		jrt = ix + iy;// 'vertical' in {0,2*(nside-1)}
		jpt = ix - iy;// 'horizontal' in {-nside+1,nside-1}

		// computes the z coordinate on the sphere
		jr = jrll[face_num + 1] * nside - jrt - 1;// ring number in
		// {1,4*nside-1}

		nr = nside;// equatorial region (the most frequent)
		n_before = ncap + nl4 * ( jr - nside );
		kshift = ( jr - nside ) % 2;
		if ( jr < nside ) { // north pole region
			nr = jr;
			n_before = 2 * nr * ( nr - 1 );
			kshift = 0;
		} else {
			if ( jr > nl3 ) { // south pole region
				nr = nl4 - jr;
				n_before = npix - 2 * ( nr + 1 ) * nr;
				kshift = 0;
			}
		};

		// computes the phi coordinate on the sphere, in [0,2Pi]
		jp = ( jpll[face_num + 1] * nr + jpt + 1 + kshift ) / 2;// 'phi' number
		// in the ring
		// in {1,4*nr}

		if ( jp > nl4 )
			jp = jp - nl4;
		if ( jp < 1 )
			jp = jp + nl4;

		ipring = n_before + jp - 1;// in {0, npix-1}
		return ipring;
	}

	/**
	 * performs conversion from RING to NESTED pixel number
	 * 
	 * @param ipring
	 *            pixel RING index number
	 * @return NEST pixel index number
	 * @throws Exception
	 */
	public int ring2nest(int ipring) throws Exception {
		int ipnest;

		double fihip, hip;
		int ip, iphi, ipt, ipring1;
		int kshift, face_num = 0, nr;
		int irn, ire, irm, irs, irt, ifm, ifp;
		int ix, iy, ix_low, ix_hi, iy_low, iy_hi, ipf;

		// coordinate of the lowest corner of each face
		// 0 added in front because the java array is zero offset
		int jrll[] = { 0, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4 };
		int jpll[] = { 0, 1, 3, 5, 7, 0, 2, 4, 6, 1, 3, 5, 7 };
		// -----------------------------------------------------------------------
		if ( ipring < 0 || ipring > npix - 1 )
			throw new Exception("ipring out of range");
		if ( x2pix[127] <= 0 )
			mkxy2pix();

		ipring1 = ipring + 1;

		// finds the ring number, the position of the ring and the face number
		if ( ipring1 <= ncap ) { // North Pole
			hip = ipring1 / 2.0;
			fihip = Math.rint(hip);
			irn = (int) Math.floor(Math.sqrt(hip - Math.sqrt(fihip))) + 1;// counted
			// from
			// North
			// pole
			iphi = ipring1 - 2 * irn * ( irn - 1 );
			kshift = 0;
			nr = irn;// 1/4 of the number of points on the current ring
			face_num = ( iphi - 1 ) / irn;// in {0,3}

		} else {
			if ( ipring1 <= nl2 * ( 5 * nside + 1 ) ) { // Equatorial

				ip = ipring1 - ncap - 1;
				irn = (int) Math.floor(ip / nl4) + nside;// counted from
				// North pole
				iphi = ip % nl4 + 1;

				kshift = ( irn + nside ) % 2; // 1 if irn+nside is odd, 0
				// otherwise
				nr = nside;
				ire = irn - nside + 1;// in {1, 2*nside +1}
				irm = nl2 + 2 - ire;
				ifm = ( iphi - ire / 2 + nside - 1 ) / nside;// face boundary
				ifp = ( iphi - irm / 2 + nside - 1 ) / nside;
				if ( ifp == ifm ) { // faces 4 to 7
					face_num = ifp % 4 + 4;
				} else {
					if ( ifp + 1 == ifm ) { // (half-)faces 0 to 3
						face_num = ifp;
					} else {
						if ( ifp - 1 == ifm ) { // (half-)faces 8 to 11
							face_num = ifp + 7;
						}
					};
				}

			} else { // South

				ip = npix - ipring1 + 1;
				hip = ip / 2.0;
				fihip = Math.rint(hip);
				irs = (int) Math.floor(Math.sqrt(hip - Math.sqrt(fihip))) + 1;// counted
				// from
				// South
				// pole
				iphi = 4 * irs + 1 - ( ip - 2 * irs * ( irs - 1 ) );

				kshift = 0;
				nr = irs;
				irn = nl4 - irs;
				face_num = ( iphi - 1 ) / irs + 8;// in {8,11}

			}
		};

		// finds the (x,y) on the face
		irt = irn - jrll[face_num + 1] * nside + 1;// in {-nside+1,0}
		ipt = 2 * iphi - jpll[face_num + 1] * nr - kshift - 1;// in
		// {-nside+1,nside-1}
		if ( ipt >= nl2 )
			ipt = ipt - 8 * nside;// for the face #4

		ix = ( ipt - irt ) / 2;
		iy = -( ipt + irt ) / 2;

		// System.out.println("face:"+face_num+" irt:"+irt+" ipt:"+ipt+"
		// ix:"+ix+" iy:"+iy);
		ix_low = ix % 128;
		ix_hi = ix / 128;
		iy_low = iy % 128;
		iy_hi = iy / 128;

		ipf = ( x2pix[ix_hi] + y2pix[iy_hi] ) * ( 128 * 128 )
				+ ( x2pix[ix_low] + y2pix[iy_low] );// in {0, Math.pow(nside,2)
		// - 1}

		ipnest = (int) (ipf + face_num * npface);// in {0,
		// 12*Math.pow(nside,2) - 1}

		return ipnest;
	}

	/**
	 * Convert from pix number to x,y inside a given face. 0,0 is the lower
	 * right corner of the face.
	 * 
	 * @param ipix
	 *            pixel index NEST number
	 * @return {@link Point} coordinate x,y
	 * @throws Exception
	 */
	public Point pix2xy_nest(int ipix) throws Exception {
		if ( ipix < 0 || ipix >= npface )
			throw new Exception("ipix out of range");
		return pix2xy_nestface(ipix);

	}

	/**
	 * Convert from a x,y in a given face to a pix number.
	 * 
	 * @param ix
	 *            x coordinate
	 * @param iy
	 *            y coordinate
	 * @param face
	 *            face nested number
	 * @return pixel index number in NEST scheme
	 * @throws Exception
	 */
	public int xy2pix_nest(int ix, int iy, int face) throws Exception {
		if ( nside < 1 || nside > ns_max )
			throw new Exception("nside out of range");
		if ( ix < 0 || ix > ( nside - 1 ) )
			throw new Exception("ix out of range");
		if ( iy < 0 || iy > ( nside - 1 ) )
			throw new Exception("iy out of range");

		int ipf, ipix;
		ipf = xy2pix_nest(ix, iy);
		ipix = (int) ( ipf + face * npface );// in {0, 12*(nside^2)-1}
		return ipix;
	}

	/**
	 * Convert from a point in a given face to a pix number. Convenience method
	 * just unpacks the point to x and y and calls the other xy2pix_nest method.
	 * 
	 * @param p
	 *            {@link Point} coordinate x,y
	 * @param face
	 *            face nested number
	 * @return pixel index nested number
	 * @throws Exception
	 */
	public int xy2pix_nest(Point p, int face) throws Exception {
		return xy2pix_nest(p.x, p.y, face);
	}

	/**
	 * Convert from a x,y in a given face to a pix number in a face without
	 * offset.
	 * 
	 * @param ix
	 *            x coordinate
	 * @param iy
	 *            y coordinate
	 * @return pixel index number
	 * @throws Exception
	 */
	public int xy2pix_nest(int ix, int iy) throws Exception {
		int ix_low, ix_hi, iy_low, iy_hi, ipf;
		if ( x2pix[127] <= 0 )
			mkxy2pix();

		ix_low = ix % 128;
		ix_hi = ix / 128;
		iy_low = iy % 128;
		iy_hi = iy / 128;

		ipf = ( x2pix[ix_hi] + y2pix[iy_hi] ) * ( 128 * 128 )
				+ ( x2pix[ix_low] + y2pix[iy_low] );
		return ipf;
	}

	/**
	 * Convert from pix number to x,y inside a given face. 0,0 is the lower
	 * right corner of the face.
	 * 
	 * @param ipix
	 *            pixel index number in that face
	 * @return {@link Point} coordinate x,y
	 * @throws Exception
	 */
	public Point pix2xy_nestface(int ipix) throws Exception {
		if ( pix2x[1023] <= 0 )
			mkpix2xy();
		int ip_low, ip_trunc, ip_med, ip_hi, ix, iy;
		ip_low = ipix % 1024;// content of the last 10 bits
		ip_trunc = ipix / 1024;// truncation of the last 10 bits
		ip_med = ip_trunc % 1024;// content of the next 10 bits
		ip_hi = ip_trunc / 1024;// content of the high weight 10 bits
		ix = 1024 * pix2x[ip_hi] + 32 * pix2x[ip_med] + pix2x[ip_low];
		iy = 1024 * pix2y[ip_hi] + 32 * pix2y[ip_med] + pix2y[ip_low];
		// System.out.println("ip_low:"+ip_low+" iptrun:"+ip_trunc+"
		// ip_med"+ip_med+" ip_hi"+ip_hi+" ix:"+ix+" iy:"+iy);
		return new Point(ix, iy);
	}

	/**
	 * integration limits in cos(theta) for a given ring i_th, i_th > 0
	 * 
	 * @param i_th
	 *            ith ring
	 * @return limits
	 */
	public double[] integration_limits_in_costh(int i_th) {

		double a, ab, b, r_n_side;

		// integration limits in cos(theta) for a given ring i_th
		// i > 0 !!!

		r_n_side = 1.0 * nside;
		if ( i_th <= nside ) {
			ab = 1.0 - ( Math.pow(i_th, 2.0) / 3.0 ) / (double) npface;
			b = 1.0 - ( Math.pow(( i_th - 1 ), 2.0) / 3.0 ) / (double) npface;
			if ( i_th == nside ) {
				a = 2.0 * ( nside - 1.0 ) / 3.0 / r_n_side;
			} else {
				a = 1.0 - Math.pow(( i_th + 1 ), 2) / 3.0 / (double) npface;
			};

		} else {
			if ( i_th < nl3 ) {
				ab = 2.0 * ( 2 * nside - i_th ) / 3.0 / r_n_side;
				b = 2.0 * ( 2 * nside - i_th + 1 ) / 3.0 / r_n_side;
				a = 2.0 * ( 2 * nside - i_th - 1 ) / 3.0 / r_n_side;
			} else {
				if ( i_th == nl3 ) {
					b = 2.0 * ( -nside + 1 ) / 3.0 / r_n_side;
				} else {
					b = -1.0 + Math.pow(( 4 * nside - i_th + 1 ), 2) / 3.0
							/ (double) npface;
				}

				a = -1.0 + Math.pow(( nl4 - i_th - 1 ), 2) / 3.0
						/ (double) npface;
				ab = -1.0 + Math.pow(( nl4 - i_th ), 2) / 3.0 / (double) npface;
			}

		}
		// END integration limits in cos(theta)
		double[] ret = { b, ab, a };
		return ret;
	}

	/**
	 * calculate the points of crosing for a given theata on the boundaries of
	 * the pixel - returns the left and right phi crosings
	 * 
	 * @param i_th
	 *            ith pixel
	 * @param i_phi
	 *            phi angle
	 * @param i_zone
	 *            ith zone (0,...,3), a quarter of sphere
	 * @param cos_theta
	 *            theta cosinus
	 * @return the left and right phi crossings
	 */
	public double[] pixel_boundaries(double i_th, double i_phi, int i_zone,
			double cos_theta) {
		double sq3th, factor, jd, ju, ku, kd, phi_l, phi_r;
		double r_n_side = 1.0 * nside;

		// HALF a pixel away from both poles
		if ( Math.abs(cos_theta) >= 1.0 - 1.0 / 3.0 / (double) npface ) {
			phi_l = i_zone * Constants.piover2;
			phi_r = ( i_zone + 1 ) * Constants.piover2;
			double[] ret = { phi_l, phi_r };
			return ret;
		}
		// -------
		// NORTH POLAR CAP
		if ( 1.50 * cos_theta >= 1.0 ) {
			sq3th = Math.sqrt(3.0 * ( 1.0 - cos_theta ));
			factor = 1.0 / r_n_side / sq3th;
			jd = (double) ( i_phi );
			ju = jd - 1;
			ku = (double) ( i_th - i_phi );
			kd = ku + 1;
			// System.out.println(" cos_theta:"+cos_theta+" sq3th:"+sq3th+"
			// factor:"+factor+" jd:"+jd+" ju:"+ju+" ku:"+ku+" kd:"+kd+ "
			// izone:"+i_zone);
			phi_l = Constants.piover2
					* ( Math.max(( ju * factor ), ( 1.0 - ( kd * factor ) )) + i_zone );
			phi_r = Constants.piover2
					* ( Math.min(( 1.0 - ( ku * factor ) ), ( jd * factor )) + i_zone );

		} else {
			if ( -1.0 < 1.50 * cos_theta ) {
				// -------
				// -------
				// EQUATORIAL ZONE
				double cth34 = 0.50 * ( 1.0 - 1.50 * cos_theta );
				double cth34_1 = cth34 + 1.0;
				int modfactor = (int) ( nside + ( i_th % 2 ) );

				jd = i_phi - ( modfactor - i_th ) / 2.0;
				ju = jd - 1;
				ku = ( modfactor + i_th ) / 2.0 - i_phi;
				kd = ku + 1;

				phi_l = Constants.piover2
						* ( Math.max(( cth34_1 - ( kd / r_n_side ) ),
								( -cth34 + ( ju / r_n_side ) )) + i_zone );

				phi_r = Constants.piover2
						* ( Math.min(( cth34_1 - ( ku / r_n_side ) ),
								( -cth34 + ( jd / r_n_side ) )) + i_zone );
				// -------
				// -------
				// SOUTH POLAR CAP

			} else {
				sq3th = Math.sqrt(3.0 * ( 1.0 + cos_theta ));
				factor = 1.0 / r_n_side / sq3th;
				int ns2 = 2 * nside;

				jd = i_th - ns2 + i_phi;
				ju = jd - 1;
				ku = ns2 - i_phi;
				kd = ku + 1;

				phi_l = Constants.piover2
						* ( Math.max(( 1.0 - ( ns2 - ju ) * factor ),
								( ( ns2 - kd ) * factor )) + i_zone );

				phi_r = Constants.piover2
						* ( Math.min(( 1.0 - ( ns2 - jd ) * factor ),
								( ( ns2 - ku ) * factor )) + i_zone );
			}// of SOUTH POLAR CAP
		}
		// and that's it
		// System.out.println(" nside:"+nside+" i_th:"+i_th+" i_phi:"+i_phi+"
		// izone:"+i_zone+" cos_theta:"+cos_theta+" phi_l:"+phi_l+"
		// phi_r:"+phi_r);

		double[] ret = { phi_l, phi_r };
		return ret;

	}

	/**
	 * return ring number for given pix in ring scheme
	 * 
	 * @param ipix
	 *            pixel index number in ring scheme
	 * @return ring number
	 * @throws Exception
	 */
	public int ring(int ipix) throws Exception {
		int iring = 0;
		int ipix1 = ipix + 1;// in {1, npix}
		int ip;
		double hip, fihip = 0;
		if ( ipix1 <= ncap ) { // North Polar cap -------------
			hip = ipix1 / 2.0;
			fihip = (int) ( hip );
			iring = (int) ( Math.sqrt(hip - Math.sqrt(fihip)) ) + 1;// counted
			// from
			// North
			// pole
		} else {
			if ( ipix1 <= nl2 * ( 5 * nside + 1 ) ) { // Equatorial region
				// ------
				ip = ipix1 - ncap - 1;
				iring = (int) ( ip / nl4 ) + nside;// counted from North pole
			} else { // South Polar cap -----------------------------------
				ip = npix - ipix1 + 1;
				hip = ip / 2.0;
				fihip = (int) ( hip );
				iring = (int) ( Math.sqrt(hip - Math.sqrt(fihip)) ) + 1;// counted
				// from
				// South
				// pole
				iring = nl4 - iring;
			}
		};
		return iring;
	}

	/**
	 * Construct a {@link SpatialVector} from the angle (theta,phi)
	 * 
	 * @param theta
	 *            angle (along meridian), in [0,Pi], theta=0 : north pole
	 * @param phi
	 *            angle (along parallel), in [0,2*Pi]
	 * @return vector {@link SpatialVector}
	 */
	public SpatialVector vector(double theta, double phi) {
		double x, y, z;
		x = 1 * Math.sin(theta) * Math.cos(phi);
		y = 1 * Math.sin(theta) * Math.sin(phi);
		z = 1 * Math.cos(theta);
		return new SpatialVector(x, y, z);
	}

	/**
	 * Converts the unit vector to pix number in NEST scheme
	 * 
	 * @param vec
	 *            {@link SpatialVector}
	 * @return pixel index number in nest scheme
	 * @throws Exception
	 */
	public int vec2pix_nest(SpatialVector vec) throws Exception {
		double[] angs = Vect2Ang(vec);//ang(vec);
		return ang2pix_nest(angs[0], angs[1]);
	}

	/**
	 * Converts the unit vector to pix number in RING scheme
	 * 
	 * @param vec
	 *            {@link SpatialVector}
	 * @return pixel index number in ring scheme
	 * @throws Exception
	 */
	public int vec2pix_ring(SpatialVector vec) throws Exception {
		double[] angs = Vect2Ang(vec);
		return ang2pix_ring(angs[0], angs[1]);
	}

	/**
	 * Converts pix number in NEST scheme to the unit vector
	 * 
	 * @param pix
	 *            pixel index number in nest scheme
	 * @return {@link SpatialVector}
	 * @throws Exception
	 */
	public SpatialVector pix2vec_nest(int pix) throws Exception {
		double[] angs = pix2ang_nest(pix);
		return vector(angs[0], angs[1]);
	}

	/**
	 * Converts pix number in RING scheme to the unit vector
	 * 
	 * @param pix
	 *            pixel index number in ring scheme
	 * @return {@link SpatialVector}
	 * @throws Exception
	 */
	public SpatialVector pix2vec_ring(int pix) throws Exception {
		double[] angs = pix2ang_ring(pix);
		return vector(angs[0], angs[1]);
	}

	/**
	 * Returns set of points along the boundary of the given pixel in NEST
	 * scheme. Step 1 gives 4 points on the corners.
	 * 
	 * @param pix
	 *            pixel index number in nest scheme
	 * @param step
	 * @return {@link SpatialVector} for each points
	 * @throws Exception
	 */
	public SpatialVector[] corners_nest(int pix, int step) throws Exception {
		int pixr = nest2ring(pix);
		return corners_ring(pixr, step);
	}

	/**
	 * Returns set of points along the boundary of the given pixel in RING
	 * scheme. Step 1 gives 4 points on the corners.
	 * 
	 * @param pix
	 *            pixel index number in ring scheme
	 * @param step
	 * @return {@link SpatialVector} for each points
	 * @throws Exception
	 */
	public SpatialVector[] corners_ring(int pix, int step) throws Exception {
		int nPoints = step * 2 + 2;
		SpatialVector[] points = new SpatialVector[nPoints];
		double[] p0 = pix2ang_ring(pix);
		double cos_theta = Math.cos(p0[0]);
		double theta = p0[0];
		double phi = p0[1];

		int i_zone = (int) ( phi / Constants.piover2 );
		int ringno = ring(pix);
		int i_phi_count = Math.min(ringno, Math.min(nside, ( nl4 ) - ringno));
		int i_phi = 0;
		double phifac = Constants.piover2 / i_phi_count;
		if ( ringno >= nside && ringno <= nl3 ) {
			// adjust by 0.5 for odd numbered rings in equatorial since
			// they start out of phase by half phifac.
			i_phi = (int) ( phi / phifac + ( ( ringno % 2 ) / 2.0 ) ) + 1;
		} else {
			i_phi = (int) ( phi / phifac ) + 1;
		}
		// adjust for zone offset
		i_phi = i_phi - ( i_zone * i_phi_count );
		int spoint = (int) ( nPoints / 2 );
		// get north south middle - middle should match theta !
		double[] nms = integration_limits_in_costh(ringno);
		double ntheta = Math.acos(nms[0]);
		double stheta = Math.acos(nms[2]);
		double[] philr = pixel_boundaries(ringno, i_phi, i_zone, nms[0]);
		if ( i_phi > ( i_phi_count / 2 ) ) {
			points[0] = vector(ntheta, philr[1]);
		} else {
			points[0] = vector(ntheta, philr[0]);
		}
		philr = pixel_boundaries(ringno, i_phi, i_zone, nms[2]);
		if ( i_phi > ( i_phi_count / 2 ) ) {
			points[spoint] = vector(stheta, philr[1]);
		} else {
			points[spoint] = vector(stheta, philr[0]);
		}
		if ( step == 1 ) {
			double mtheta = Math.acos(nms[1]);
			philr = pixel_boundaries(ringno, i_phi, i_zone, nms[1]);
			points[1] = vector(mtheta, philr[0]);
			points[3] = vector(mtheta, philr[1]);
		} else {
			double cosThetaLen = nms[2] - nms[0];
			double cosThetaStep = ( cosThetaLen / ( step + 1 ) ); // skip
			// North
			// and south
			for ( int p = 1; p <= step; p++ ) {
				/* Integrate points along the sides */
				cos_theta = nms[0] + ( cosThetaStep * p );
				theta = Math.acos(cos_theta);
				philr = pixel_boundaries(ringno, i_phi, i_zone, cos_theta);
				points[p] = vector(theta, philr[0]);
				points[nPoints - p] = vector(theta, philr[1]);
			}
		}
		return points;
	}

	/**
	 * calculates angular resolution of the pixel map in arc seconds.
	 * 
	 * @param nside
	 * @return double resolution in arcsec
	 */
	static public double getPixRes(long nside) {
		double res = 0.;
		double degrad = Math.toDegrees(1.0);
		double skyArea = 4. * Constants.PI * degrad * degrad; // 4PI steredian
		// in deg^2
		double arcSecArea = skyArea * 3600. * 3600.; // 4PI steredian in
		// (arcSec^2)
		long npixels = 12 * nside * nside;
		res = arcSecArea / npixels; // area per pixel
		res = Math.sqrt(res); // angular size of the pixel arcsec
		return res;
	}

	/**
	 * calculate required nside given pixel size in arcsec
	 * 
	 * @param pixsize
	 *            in arcsec
	 * @return long nside parameter
	 */
	static public long calculateNSide(double pixsize) {
		long res = 0;
		double pixelArea = pixsize * pixsize;
		double degrad = Math.toDegrees(1.);
		double skyArea = 4. * Constants.PI * degrad * degrad * 3600. * 3600.;
		long npixels = (long) ( skyArea / pixelArea );
		long nsidesq = npixels / 12;
		long nside_req = (long) Math.sqrt(nsidesq);
		long mindiff = ns_max;
		int indmin = 0;
		for ( int i = 0; i < nsidelist.length; i++ ) {
			if ( Math.abs(nside_req - nsidelist[i]) <= mindiff ) {
				mindiff = Math.abs(nside_req - nsidelist[i]);
				res = nsidelist[i];
				indmin = i;
			}
			if ( ( nside_req > res ) && ( nside_req < ns_max ) )
				res = nsidelist[indmin + 1];
			if ( nside_req > ns_max ) {
				System.out.println("nside cannot be bigger than " + ns_max);
				return ns_max;
			}

		}
		return res;
	}

	/**
	 * calculates vector corresponding to angles theta (co-latitude measured
	 * from North pole, in [0,pi] radians) phi (longitude measured eastward in
	 * [0,2pi] radians) North pole is (x,y,z) = (0, 0, 1)
	 * 
	 * @param theta
	 *            angle (along meridian), in [0,Pi], theta=0 : north pole
	 * @param phi
	 *            angle (along parallel), in [0,2*Pi]
	 * @return SpatialVector
	 * @throws IllegalArgumentException
	 */
	public Vector3d Ang2Vec(double theta, double phi) {
		double PI = Math.PI;
		String SID = "Ang2Vec:";
		SpatialVector v;
		if ( ( theta < 0.0 ) || ( theta > PI ) ) {
			throw new IllegalArgumentException(SID
					+ " theta out of range [0.,PI]");
		}
		double stheta = Math.sin(theta);
		double x = stheta * Math.cos(phi);
		double y = stheta * Math.sin(phi);
		double z = Math.cos(theta);
		v = new SpatialVector(x, y, z);
		return (Vector3d) v;
	}

	/**
	 * converts a SpatialVector in a tuple of angles tup[0] = theta co-latitude
	 * measured from North pole, in [0,PI] radians, tup[1] = phi longitude
	 * measured eastward, in [0,2PI] radians
	 * 
	 * @param v
	 *            SpatialVector
	 * @return double[] out_tup out_tup[0] = theta out_tup[1] = phi
	 */
	public double[] Vect2Ang(SpatialVector v) {
		double[] out_tup = new double[2];
		double norm = v.length();
		double z = v.z() / norm;
		double theta = Math.acos(z);
		double phi = 0.;
		if ( ( v.x() != 0. ) || ( v.y() != 0 ) ) {
			phi = Math.atan2(v.y(), v.x()); // phi in [-pi,pi]
		}
		if ( phi < 0 )
			phi += 2.0 * Math.PI; // phi in [0, 2pi]
		out_tup[0] = theta;
		out_tup[1] = phi;
		return out_tup;
	}

	/**
	 * returns nside such that npix = 12*nside^2 nside should by power of 2 and
	 * smaller than ns_max if not return -1
	 * 
	 * @param npix
	 *            long the number of pixels in the map
	 * @return nside long the map resolution parameter
	 */
	public long Npix2Nside(long npix) {
		long nside = 0;
		long npixmax = 12 * (long) ns_max * (long) ns_max;
		System.out.println("ns_max=" + ns_max + "  npixmax=" + npixmax);
		String SID = "Npix2Nside:";
		nside = (long) Math.rint(Math.sqrt(npix / 12));
		if ( npix < 12 ) {
			throw new IllegalArgumentException(SID
					+ " npix is too small should be > 12");
		}
		if ( npix > npixmax ) {
			throw new IllegalArgumentException(SID
					+ " npix is too large > 12 * ns_max^2");
		}
		double fnpix = 12.0 * nside * nside;
		if ( Math.abs(fnpix - npix) > 1.0e-2 ) {
			throw new IllegalArgumentException(SID
					+ "  npix is not 12*nside*nside");
		}
		double flog = Math.log((double) nside) / Math.log(2.0);
		double ilog = Math.rint(flog);
		if ( Math.abs(flog - ilog) > 1.0e-6 ) {
			throw new IllegalArgumentException(SID
					+ "  nside is not power of 2");
		}
		return nside;
	}

	/**
	 * calculates npix such that npix = 12*nside^2 nside should be a power of 2,
	 * and smaller than ns_max otherwise return -1
	 * 
	 * @param nside
	 *            long the map resolution
	 * @return npix long the number of pixels in the map
	 */
	public long Nside2Npix(long nside) {

		long[] nsidelist = { 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048,
				4096, 8192, 16384, 32768, 65563, 131072, 262144, 524288,
				1048576, 2097152, 4194304 };

		long res = 0;
		String SID = "Nside2Npix:";
		if ( Arrays.binarySearch(nsidelist, nside) < 0 ) {
			throw new IllegalArgumentException(SID
					+ " nside should be >0, power of 2, <" + ns_max);
		}
		res = 12 * nside * nside;
		return res;
	}

	/**
	 * calculates the surface of spherical triangle defined by vertices v1,v2,v3
	 * Algorithm: finds triangle sides and uses l'Huilier formula to compute
	 * "spherical excess" = surface area of triangle on a sphere of radius one
	 * see, eg Bronshtein, Semendyayev Eq 2.86 half perimeter hp =
	 * 0.5*(side1+side2+side3) l'Huilier formula x0 = tan( hp/2.) x1 = tan((hp -
	 * side1)/2.) x2 = tan((hp - side2)/2.) x3 = tan((hp - side3)/2.)
	 * 
	 * @param v1
	 *            SpatialVector
	 * @param v2
	 *            SpatialVector
	 * @param v3
	 *            SpatialVector vertices of the triangle
	 * @return double the triangle surface in steradians of the spherical
	 *         triangle with vertices vec1, vec2, vec3
	 * @throws Exception
	 */
	public double SurfaceTriangle(SpatialVector v1, SpatialVector v2,
			SpatialVector v3) throws Exception {
		double res = 0.;
		double side1 = AngDist(v2, v3) / 4.0;
		double side2 = AngDist(v3, v1) / 4.0;
		double side3 = AngDist(v1, v2) / 4.0;
		double x0 = Math.tan(side1 + side2 + side3);
		double x1 = Math.tan(side2 + side3 - side1);
		double x2 = Math.tan(side1 + side3 - side2);
		double x3 = Math.tan(side1 + side2 - side3);
		res = 4.0 * Math.atan(Math.sqrt(x0 * x1 * x2 * x3));

		return res;
	}

	/**
	 * calculates angular distance (in radians) between 2 Vectors v1 and v2 In
	 * general dist = acos(v1.v2) except if the vectors are almost aligned
	 * 
	 * @param v1
	 *            SpatialVector
	 * @param v2
	 *            SpatialVector
	 * @return double dist
	 * @throws Exception
	 */
	public double AngDist(SpatialVector v1, SpatialVector v2) throws Exception {
		double dist = 0.;
		double aligned = 0.999;
		/* Normalize both vectors */
		SpatialVector r1 = v1;
		SpatialVector r2 = v2;
		r1.normalize();
		r2.normalize();
		double sprod = r1.dot(r2);
		/* This takes care about the bug in vecmath method from java3d project */
		if ( sprod > aligned ) { // almost aligned
			r1.sub(r2);
			double diff = r1.length();
			dist = 2.0 * Math.asin(diff / 2.0);

		} else if ( sprod < -aligned ) {
			r1.add(r2);
			double diff = r1.length();
			dist = Math.PI - 2.0 * Math.asin(diff / 2.0);
		} else {
			// javax.vecmath.Vector3d r3d1 = new Vector3d(r1.x(), r1.y(),
			// r1.z());
			// javax.vecmath.Vector3d r3d2 = new Vector3d(r2.x(), r2.y(),
			// r2.z());
			dist = Math.acos(sprod);// r3d1.angle(r3d2);
		}
		return dist;
	}

	//
	// public static double[] ang(Vector3d vec) {
	// double theta = Math.acos(vec.getZ());
	// double phi = Math.asin(vec.getY() / Math.sin(theta));
	// double[] ret = { theta, phi };
	// return ret;
	//
	// }

	/**
	 * calculates a vector production of two vectors.
	 * 
	 * @param v1
	 *            Vectror containing 3 elements of Number type
	 * @param v2
	 *            Vector containing 3 elements of Number type
	 * @return Vector of 3 Objects of Double type
	 * @throws Exception
	 */
	public Vector<Object> VectProd(Vector v1, Vector v2) throws Exception {
		Vector<Object> res = new Vector<Object>();
		double[] v1_element = new double[3];
		double[] v2_element = new double[3];
		for ( int i = 0; i < 3; i++ ) {
			if ( v1.get(i).getClass().isInstance(Number.class) ) {
				v1_element[i] = ( (Number) v1.get(i) ).doubleValue();
			} else {
				throw new Exception();
			}
			if ( v2.get(i).getClass().isInstance(Number.class) ) {
				v2_element[i] = ( (Number) v2.get(i) ).doubleValue();
			} else {
				throw new Exception();
			}

		}

		Double value = new Double(v1_element[1] * v2_element[2] - v1_element[2]
				* v2_element[1]);
		res.add((Object) value);
		value = new Double(v1_element[1] * v2_element[2] - v1_element[2]
				* v2_element[1]);
		res.add((Object) value);
		value = new Double(v1_element[1] * v2_element[2] - v1_element[2]
				* v2_element[1]);
		res.add((Object) value);
		return res;
	}

	/**
	 * calculates a dot product (inner product) of two 3D vectors the result is
	 * double
	 * 
	 * @param v1
	 *            3d Vector of Number Objects (Double, long .. )
	 * @param v2
	 *            3d Vector
	 * @return double
	 * @throws Exception
	 */
	public double dotProduct(SpatialVector v1, SpatialVector v2)
			throws Exception {
		double prod = v1.x() * v2.x() + v1.y() * v2.y() + v1.z() * v2.z();

		return prod;
	}

	/**
	 * calculate cross product of two vectors
	 * 
	 * @param v1
	 *            SpatialVector
	 * @param v2
	 *            SpatialVector
	 * @return SpatialVector result of the product
	 */
	public SpatialVector crossProduct(SpatialVector v1, SpatialVector v2) {
		SpatialVector res = new SpatialVector(0., 0., 0.);
		// double x = v1.y() * v2.z() - v1.z() * v2.y();
		// double y = v1.z() * v2.x() - v1.x() * v2.z();
		// double z = v1.x() * v2.y() - v1.y() * v2.x();
		res.cross(v1, v2);
		return res;
	}

	/**
	 * generates in the RING or NESTED scheme all pixels that lies within an
	 * angular distance Radius of the center.
	 * 
	 * @param nside
	 *            long map resolution
	 * @param vector
	 *            Vector3d pointing to the disc center
	 * @param radius
	 *            double angular radius of the disk (in RADIAN )
	 * @param nest
	 *            int 0 (default) if output is in RING scheme, if set to 1
	 *            output is in NESTED
	 * @param inclusive
	 *            int 0 (default) only pixsels whose center lie in the triangle
	 *            are listed, if set to 1, all pixels overlapping the triangle
	 *            are listed
	 * @return ArrayList of pixel numbers calls: RingNum(nside, ir)
	 *         InRing(nside, iz, phi0, dphi,nest)
	 */
	public ArrayList<Long> query_disc(long nside, Vector3d vector,
			double radius, int nest, int inclusive) {
		ArrayList<Long> res = new ArrayList<Long>();
		long irmin, irmax, iz;
		double x0, y0, z0, radius_eff;
		double a, b, c, cosang;
		double dth1, dth2;
		double phi0, cosphi0, cosdphi, dphi;
		double rlat0, rlat1, rlat2, zmin, zmax, z;
		boolean do_inclusive = false;
		boolean do_nest = false;
		String SID = "QUERY_DISC";

		if ( radius < 0.0 || radius > Constants.PI ) {
			throw new IllegalArgumentException(SID
					+ ": angular radius is in RADIAN and should be in [0,pi]");
		}
		if ( inclusive == 1 )
			do_inclusive = true;
		if ( nest == 1 )
			do_nest = true;

		dth1 = 1.0 / ( 3.0 * nside * nside );
		dth2 = 2.0 / ( 3.0 * nside );

		radius_eff = radius;
		if ( do_inclusive )
			radius_eff += Constants.PI / ( 4.0 * nside ); // increas radius by
		// half pixel
		cosang = Math.cos(radius_eff);
		/* disc center */
		vector.normalize();
		x0 = vector.x; // norm_vect0;
		y0 = vector.y; // norm_vect0;
		z0 = vector.z; // norm_vect0;

		phi0 = 0.0;
		dphi = 0.0;
		if ( x0 != 0. || y0 != 0. )
			phi0 = bm.MODULO(Math.atan2(y0, x0) + Constants.twopi,
					Constants.twopi); // in [0,
		// 2pi]
		cosphi0 = Math.cos(phi0);
		a = x0 * x0 + y0 * y0;
		/* coordinate z of highest and lowest points in the disc */
		rlat0 = Math.asin(z0); // latitude in RAD of the center
		rlat1 = rlat0 + radius_eff;
		rlat2 = rlat0 - radius_eff;
		//
		if ( rlat1 >= Constants.piover2 ) {
			zmax = 1.0;
		} else {
			zmax = Math.sin(rlat1);
		}
		irmin = RingNum(nside, zmax);
		irmin = Math.max(1, irmin - 1); // start from a higher point to be safe
		if ( rlat2 <= -Constants.piover2 ) {
			zmin = -1.0;
		} else {
			zmin = Math.sin(rlat2);
		}
		irmax = RingNum(nside, zmin);
		irmax = Math.min(4 * nside - 1, irmax + 1); // go down to a lower point

		/* loop on ring number */
		for ( iz = irmin; iz <= irmax; iz++ ) {
			if ( iz <= nside - 1 ) { // north polar cap
				z = 1.0 - iz * iz * dth1;
			} else if ( iz <= 3 * nside ) { // tropical band + equator
				z = ( 2.0 * nside - iz ) * dth2;
			} else {
				z = -1.0 + ( 4.0 * nside - iz ) * ( 4.0 * nside - iz ) * dth1;
			}
			/* find phi range in the disc for each z */
			b = cosang - z * z0;
			c = 1.0 - z * z;
			cosdphi = b / Math.sqrt(a * c);
			long done = 0;

			if ( Math.abs(x0) <= 1.0e-12 && Math.abs(y0) <= 1.0e-12 ) {
				cosdphi = -1.0;
				dphi = Constants.PI;
				done = 1;
			}
			if ( done == 0 ) {
				if ( Math.abs(cosdphi) <= 1.0 ) {
					dphi = Math.acos(cosdphi); // in [0,PI]
				} else {
					if ( cosphi0 >= cosdphi ) {
						dphi = Constants.PI; // all the pixels at this
						// elevation are in
						// the disc
					} else {
						done = 2; // out of the disc
					}
				}

			}
			if ( done < 2 ) { // pixels in disc
				/* find pixels in the disc */

				ArrayList<Long> listir = InRing(nside, iz, phi0, dphi, do_nest);//InRing(nside, iz, phi0, dphi, do_nest)
				res.addAll(listir);
			}
		}
		return res;
	}

	/**
	 * returns the ring number in {1, 4*nside - 1} calculated from z coordinate
	 * 
	 * @param nside
	 *            long resolution
	 * @param z
	 *            double z coordinate
	 * @return long ring number
	 */
	public long RingNum(long nside, double z) {
		long iring = 0;
		/* equatorial region */

		iring = (long) Math.round(nside * ( 2.0 - 1.5 * z ));
		/* north cap */
		if ( z > Constants.twothird ) {
			iring = (long) Math.round(nside * Math.sqrt(3.0 * ( 1.0 - z )));
			if ( iring == 0 )
				iring = 1;
		}
		/* south cap */
		if ( z < -Constants.twothird ) {
			iring = (long) Math.round(nside * Math.sqrt(3.0 * ( 1.0 + z )));
			if ( iring == 0 )
				iring = 1;
			iring = 4 * nside - iring;
		}
		return iring;
	}
	/**
	 * returns the list of pixels in RING or NEST scheme with latitude in [phi0 -
	 * dpi, phi0 + dphi] on the ring iz in [1, 4*nside -1 ] The pixel id numbers
	 * are in [0, 12*nside^2 - 1] the indexing is in RING, unless nest is set to
	 * 1
	 * 
	 * @param nside
	 *            long the map resolution
	 * @param iz
	 *            long ring number
	 * @param phi0
	 *            double
	 * @param dphi
	 *            double
	 * @param nest
	 *            boolean format flag
	 * @return ArrayList of pixels
	 * @throws IllegalArgumentException *
	 */
	public ArrayList<Long> InRing(long nside, long iz, double phi0,
			double dphi, boolean nest) {
		return InRing(nside,iz,phi0,dphi,nest,false);
	}
	/**
	 * returns the list of pixels in RING or NEST scheme with latitude in [phi0 -
	 * dpi, phi0 + dphi] on the ring iz in [1, 4*nside -1 ] The pixel id numbers
	 * are in [0, 12*nside^2 - 1] the indexing is in RING, unless nest is set to
	 * 1
	 * 
	 * @param nside
	 *            long the map resolution
	 * @param iz
	 *            long ring number
	 * @param phi0
	 *            double
	 * @param dphi
	 *            double
	 * @param nest
	 *            boolean format flag
	 * @param conservative
	 *            if true, include every intersected pixels, even if pixel
	 *            CENTER is not in the range [phi_low, phi_hi]. If not, strict :
	 *            include only pixels whose CENTER is in [phi_low, phi_hi]
	 * @return ArrayList of pixels
	 * @throws IllegalArgumentException *
	 */
	public ArrayList<Long> InRing(long nside, long iz, double phi0,
			double dphi, boolean nest, boolean conservative) {
		boolean take_all = false;
		boolean to_top = false;
		boolean do_ring = true;
		// String SID = "InRing:";
		double epsilon = 1.0e-13; // the constant to eliminate
		// java calculation jitter
		if ( nest )
			do_ring = false;
		double shift = 0.;
		long ir = 0;
		long kshift, nr, ipix1, ipix2, nir1, nir2, ncap, npix;
		long ip_low = 0, ip_hi = 0, in, inext, nir;
		ArrayList<Long> res = new ArrayList<Long>();
		npix = 12 * nside * nside; // total number of pixels
		ncap = 2 * nside * ( nside - 1 ); // number of pixels in the north
		// polar
		// cap
		double phi_low = bm.MODULO(phi0 - dphi, Constants.twopi) - epsilon; // phi
		// min,
		// excluding
		// 2pi period
		double phi_hi = bm.MODULO(phi0 + dphi, Constants.twopi) + epsilon;

		if ( Math.abs(dphi - Constants.PI) < 1.0e-6 )
			take_all = true;
		/* identifies ring number */
		if ( ( iz >= nside ) && ( iz <= 3 * nside ) ) { // equatorial region
			ir = iz - nside + 1; // in [1, 2*nside + 1]
			ipix1 = ncap + 4 * nside * ( ir - 1 ); // lowest pixel number in
			// the
			// ring
			ipix2 = ipix1 + 4 * nside - 1; // highest pixel number in the ring
			kshift = (long) bm.MODULO(ir, 2.);

			nr = nside * 4;
		} else {
			if ( iz < nside ) { // north pole
				ir = iz;
				ipix1 = 2 * ir * ( ir - 1 ); // lowest pixel number
				ipix2 = ipix1 + 4 * ir - 1; // highest pixel number
			} else { // south pole
				ir = 4 * nside - iz;

				ipix1 = npix - 2 * ir * ( ir + 1 ); // lowest pixel number
				ipix2 = ipix1 + 4 * ir - 1; // highest pixel number
			}
			nr = ir * 4;
			kshift = 1;
		}
		// Construct the pixel list
		if ( take_all ) {
			nir = ipix2 - ipix1 + 1;
			if ( do_ring ) {
				long ind = 0;
				for ( long i = ipix1; i <= ipix2; i++ ) {
					if(i>-1) {
						res.add((int) ind, new Long(i));					
						ind++;
					}
				}
			} else {
				try {
					in = ring2nest((int) ipix1);
					if(in>-1)
						res.add(0, new Long(in));
					for ( int i = 1; i < nir; i++ ) {
						inext = next_in_line_nest(nside, in);
						in = inext;
						if(in>-1)
							res.add(i, new Long(in));
					}
				} catch ( Exception e ) {
					e.printStackTrace();
				}
			}
			return res;
		}
		shift = kshift / 2.0;

		// conservative : include every intersected pixel, even if the
		// pixel center is out of the [phi_low, phi_hi] region
		if ( conservative ) {
			ip_low = (long) Math.round(( nr * phi_low ) / Constants.twopi
					- shift);
			ip_hi = (long) Math
					.round(( nr * phi_hi ) / Constants.twopi - shift);
			ip_low = (long) bm.MODULO(ip_low, nr); // in [0, nr - 1]
			ip_hi = (long) bm.MODULO(ip_hi, nr); // in [0, nr - 1]
		} else { // strict: includes only pixels whose center is in
			// [phi_low,phi_hi]

			ip_low = (long) Math.round(( nr * phi_low ) / Constants.twopi
					- shift);
			ip_hi = (long) Math
					.round(( nr * phi_hi ) / Constants.twopi - shift);
			if ( ip_low == ip_hi + 1 )
				ip_low = ip_hi;

			if ( ( ip_low - ip_hi == 1 ) && ( dphi * nr < Constants.PI ) ) {
				// the interval is too small ( and away from pixel center)
				// so no pixels is included in the list
				System.out
						.println("the longerval is too small and avay from center");
				return res; // return empty list
			}
			ip_low = Math.min(ip_low, nr - 1);
			ip_hi = Math.max(ip_hi, 0);
		}
		//
		if ( ip_low > ip_hi )
			to_top = true;
		ip_low += ipix1;
		ip_hi += ipix1;
		if ( to_top ) {
			nir1 = ipix2 - ip_low + 1;
			nir2 = ip_hi - ipix1 + 1;
			nir = nir1 + nir2;
			if ( do_ring ) {
				int ind = 0;
				for ( long i = ip_low; i <= ipix2; i++ ) {
					if(i>-1) {
						res.add(ind, new Long(i));
						ind++;
					}
				}
				// ind = nir1;
				for ( long i = ipix1; i <= ip_hi; i++ ) {
					if(i>-1) {
						res.add(ind, new Long(i));
						ind++;
					}					
				}
			} else {
				try {
					in = ring2nest((int) ip_low);
					res.add(0, new Long(in));
					for ( long i = 1; i <= nir - 1; i++ ) {
						inext = next_in_line_nest(nside, in);
						in = inext;
						if(i>-1)
							res.add((int) i, new Long(in));
					}
				} catch ( Exception e ) {
					e.printStackTrace();
				}
			}
		} else {
			nir = ip_hi - ip_low + 1;
			if ( do_ring ) {
				int ind = 0;
				for ( long i = ip_low; i <= ip_hi; i++ ) {
					if(i>-1) {
						res.add(ind, new Long(i));					
						ind++;
					}
				}
			} else {
				try {
					in = ring2nest((int) ip_low);
					res.add(0, new Long(in));
					for ( int i = 1; i <= nir - 1; i++ ) {
						inext = next_in_line_nest(nside, in);
						in = inext;
						if(i>-1)
							res.add(i, new Long(in));
					}
				} catch ( Exception e ) {
					e.printStackTrace();
				}
			}
		}
		return res;
	}

	/**
	 * returns the list of pixels in RING or NEST scheme with latitude in [phi0 -
	 * dpi, phi0 + dphi] on the ring iz in [1, 4*nside -1 ] The pixel id numbers
	 * are in [0, 12*nside^2 - 1] the indexing is in RING, unless nest is set to
	 * 1
	 * NOTE: this is the f90 code 'in_ring' method ported to java with 'conservative' flag to false
	 * 
	 * @param nside
	 *            long the map resolution
	 * @param iz
	 *            long ring number
	 * @param phi0
	 *            double
	 * @param dphi
	 *            double
	 * @param nest
	 *            boolean format flag
	 * @return ArrayList of pixels
	 * @throws IllegalArgumentException
	 */
	public ArrayList<Long> InRingCxx(long nside, long iz, double phi0,
			double dphi, boolean nest) {
		long nr, ir, ipix1;

		double shift = 0.5;

		if ( iz < nside ) // north pole
		{
			ir = iz;
			nr = ir * 4;
			ipix1 = 2 * ir * ( ir - 1 ); // lowest pixel number in the ring
		} else if ( iz > ( 3 * nside ) ) // south pole
		{
			ir = 4 * nside - iz;
			nr = ir * 4;
			ipix1 = npix - 2 * ir * ( ir + 1 ); // lowest pixel number in the
			// ring
		} else // equatorial region
		{
			ir = iz - nside + 1; // within {1, 2*nside + 1}
			nr = nside * 4;
			if ( ( ir & 1 ) == 0 )
				shift = 0.;
			ipix1 = ncap + ( ir - 1 ) * nr; // lowest pixel number in the ring
		}

		long ipix2 = ipix1 + nr - 1; // highest pixel number in the ring
		ArrayList<Long> listir = new ArrayList<Long>();
		// ----------- constructs the pixel list --------------
		if ( dphi > ( Constants.PI - 1e-7 ) )
			for ( Long i = ipix1; i <= ipix2; ++i )
				listir.add(i);
		else {
			int ip_lo = (int) ( Math.floor(nr * ( 1 / Constants.twopi )
					* ( phi0 - dphi ) - shift) + 1 );
			int ip_hi = (int) ( Math.floor(nr * 1 / Constants.twopi
					* ( phi0 + dphi ) - shift) );
			long pixnum = (int) ( ip_lo + ipix1 );
			if ( pixnum < ipix1 )
				pixnum += nr;
			for ( int i = ip_lo; i <= ip_hi; ++i, ++pixnum ) {
				if ( pixnum > ipix2 )
					pixnum -= nr;
				listir.add(pixnum);
			}
		}
		 ArrayList<Long> listirnest = new ArrayList<Long>();
		listir.trimToSize();
		if ( nest ) {
			int i = 0;
			while ( i < listir.size() ) {
				long ipring = listir.get(i);
				try {
					long ipnest = ring2nest((int) ipring);
					listirnest.add(ipnest);
					i++;
				} catch ( Exception ex ) {
					ex.printStackTrace();
					break;// Very bad!
				}
			}
			return listirnest;
		}
		return listir;

	}

	/**
	 * calculates the pixel that lies on the East side (and the same latitude)
	 * as the given NESTED pixel number - ipix
	 * 
	 * @param nside
	 *            long resolution
	 * @param ipix
	 *            long pixel number
	 * @return long next pixel in line
	 * @throws Exception
	 * @throws IllegalArgumentException
	 */
	public long next_in_line_nest(long nside, long ipix) throws Exception {
		long npix, ipf, ipo, ix, ixp, iy, iym, ixo, iyo, face_num, other_face;
		@SuppressWarnings("unused")
		long ia, ib, ibp, ibm, ib2, nsidesq;
		int icase;
		long local_magic1, local_magic2;
		long[] ixiy = new long[2];
		long inext = 0; // next in line pixel in Nest scheme
		String SID = "next_in_line:";
		if ( ( nside < 1 ) || ( nside > ns_max ) ) {
			throw new IllegalArgumentException(SID
					+ " nside should be power of 2 >0 and < " + ns_max);
		}
		nsidesq = nside * nside;
		npix = 12 * nsidesq; // total number of pixels
		if ( ( ipix < 0 ) || ( ipix > npix - 1 ) ) {
			throw new IllegalArgumentException(SID
					+ " ipix out of range defined by nside");
		}
		// initiates array for (x,y) -> pixel number -> (x,y) mapping
		if ( x2pix[127] <= 0 ) // xmax-1
			mkxy2pix();
		local_magic1 = ( nsidesq - 1 ) / 3;
		local_magic2 = 2 * local_magic1;
		face_num = ipix / nsidesq;
		ipf = (long) bm.MODULO(ipix, nsidesq); // Pixel number in face
		ixiy[0] = pix2xy_nest((int) ipf).x;
		ixiy[1] = pix2xy_nest((int) ipf).y;
		ix = ixiy[0];
		iy = ixiy[1];
		ixp = ix + 1;
		iym = iy - 1;
		boolean sel = false;
		icase = -1; // iside the nest flag
		// Exclude corners
		if ( ipf == local_magic2 ) { // west coirner
			inext = ipix - 1;
			return inext;
		}
		if ( ( ipf == nsidesq - 1 ) && !sel ) { // North corner
			icase = 6;
			sel = true;
		}
		if ( ( ipf == 0 ) && !sel ) { // Siuth corner
			icase = 7;
			sel = true;
		}
		if ( ( ipf == local_magic1 ) && !sel ) { // East corner
			icase = 8;
			sel = true;
		}
		// Detect edges
		if ( ( ( ipf & local_magic1 ) == local_magic1 ) && !sel ) { // North-East
			icase = 1;
			sel = true;
		}
		if ( ( ( ipf & local_magic2 ) == 0 ) && !sel ) { // South-East
			icase = 4;
			sel = true;
		}
		if ( !sel ) { // iside a face
			inext = xy2pix_nest((int) ixp, (int) iym, (int) face_num);
			return inext;
		}
		//
		ia = face_num / 4; // in [0,2]
		ib = (long) bm.MODULO(face_num, 4); // in [0,3]
		ibp = (long) bm.MODULO(ib + 1, 4);
		ibm = (long) bm.MODULO(ib + 4 - 1, 4);
		ib2 = (long) bm.MODULO(ib + 2, 4);

		if ( ia == 0 ) { // North pole region
			switch ( icase ) {
				case 1:
					other_face = 0 + ibp;
					ipo = (long) bm.MODULO(bm.swapLSBMSB(ipf), nsidesq);
					inext = other_face * nsidesq + ipo;
					break;
				case 4:
					other_face = 4 + ibp;
					ipo = (long) bm.MODULO(bm.invMSB(ipf), nsidesq); // SE-NW
					// flip

					ixiy[0] = pix2xy_nest((int) ipo).x;
					ixiy[1] = pix2xy_nest((int) ipo).y;
					ixo = ixiy[0];
					iyo = ixiy[1];

					inext = xy2pix_nest((int) ixo + 1, (int) iyo,
							(int) other_face);

					break;
				case 6: // North corner
					other_face = 0 + ibp;
					inext = other_face * nsidesq + nsidesq - 1;
					break;
				case 7:
					other_face = 4 + ibp;
					inext = other_face * nsidesq + local_magic2 + 1;
					break;
				case 8:
					other_face = 0 + ibp;
					inext = other_face * nsidesq + local_magic2;
					break;
			}

		} else if ( ia == 1 ) { // Equatorial region
			switch ( icase ) {
				case 1: // NorthEast edge
					other_face = 0 + ib;
					// System.out.println("ipf="+ipf+" nsidesq="+nsidesq+"
					// invLSB="+bm.invLSB(ipf));
					ipo = (long) bm.MODULO((double) bm.invLSB(ipf),
							(double) nsidesq); // NE-SW flip
					// System.out.println(" ipo="+ipo);

					ixiy[0] = pix2xy_nest((int) ipo).x;
					ixiy[1] = pix2xy_nest((int) ipo).y;

					ixo = ixiy[0];
					iyo = ixiy[1];
					inext = xy2pix_nest((int) ixo, (int) iyo - 1,
							(int) other_face);
					break;
				case 4: // SouthEast edge
					other_face = 8 + ib;
					ipo = (long) bm.MODULO(bm.invMSB(ipf), nsidesq);
					ixiy[0] = pix2xy_nest((int) ipo).x;
					ixiy[1] = pix2xy_nest((int) ipo).y;
					inext = xy2pix_nest((int) ixiy[0] + 1, (int) ixiy[1],
							(int) other_face);
					break;
				case 6: // Northy corner
					other_face = 0 + ib;
					inext = other_face * nsidesq + local_magic2 - 2;
					break;
				case 7: // South corner
					other_face = 8 + ib;
					inext = other_face * nsidesq + local_magic2 + 1;
					break;
				case 8: // East corner
					other_face = 4 + ibp;
					inext = other_face * nsidesq + local_magic2;
					break;

			}
		} else { // South pole region
			switch ( icase ) {
				case 1: // NorthEast edge
					other_face = 4 + ibp;
					ipo = (long) bm.MODULO(bm.invLSB(ipf), nsidesq); // NE-SW
					// flip
					ixiy[0] = pix2xy_nest((int) ipo).x;
					ixiy[1] = pix2xy_nest((int) ipo).y;
					inext = xy2pix_nest((int) ixiy[0], (int) ixiy[1] - 1,
							(int) other_face);
					break;
				case 4: // SouthEast edge
					other_face = 8 + ibp;
					ipo = (long) bm.MODULO(bm.swapLSBMSB(ipf), nsidesq); // E-W
					// flip
					inext = other_face * nsidesq + ipo; // (8)
					break;
				case 6: // North corner
					other_face = 4 + ibp;
					inext = other_face * nsidesq + local_magic2 - 2;
					break;
				case 7: // South corner
					other_face = 8 + ibp;
					inext = other_face * nsidesq;
					break;
				case 8: // East corner
					other_face = 8 + ibp;
					inext = other_face * nsidesq + local_magic2;
					break;
			}
		}
		return inext;
	}

	/**
	 * finds pixels that lie within a CONVEX polygon defined by its vertex on
	 * sphere
	 * 
	 * @param nside
	 *            the map resolution
	 * @param vlist
	 *            ArrayList of vectors defining the polygon vertices
	 * @param nest
	 *            if set to 1 use NESTED scheme
	 * @param inclusive
	 *            if set 1 returns all pixels crossed by polygon boundaries
	 * @return ArrayList of pixels algorithm: the polygon is divided into
	 *         triangles vertex 0 belongs to all triangles
	 * @throws IllegalArgumentException
	 */
	public ArrayList<Long> query_polygon(long nside, ArrayList<Object> vlist,
			long nest, long inclusive) throws Exception {
		ArrayList<Long> res = new ArrayList<Long>();
		int nv = vlist.size();
		healpix.tools.SpatialVector vp0, vp1, vp2;
		healpix.tools.SpatialVector vo;
		ArrayList<Long> vvlist = new ArrayList<Long>();
		double hand;
		double[] ss = new double[nv];
		long npix;
		int ix = 0;

		int n_remain, np, nm, nlow;
		String SID = "QUERY_POLYGON";

		// Start polygon
		for ( int k = 0; k < nv; k++ )
			ss[k] = 0.;
		/* -------------------------------------- */
		n_remain = nv;
		if ( n_remain < 3 ) {
			throw new IllegalArgumentException(SID
					+ " Number of vertices should be >= 3");
		}
		/*---------------------------------------------------------------- */
		/* Check that the poligon is convex or has only one concave vertex */
		/*---------------------------------------------------------------- */
		int i0 = 0;
		int i2 = 0;
		if ( n_remain > 3 ) { // a triangle is always convex
			for ( int i1 = 1; i1 <= n_remain - 1; i1++ ) { // in [0,n_remain-1]
				i0 = (int) bm.MODULO(i1 - 1, n_remain);
				i2 = (int) bm.MODULO(i1 + 1, n_remain);
				vp0 = (SpatialVector) vlist.get(i0); // select vertices by 3
				// neighbour
				vp1 = (SpatialVector) vlist.get(i1);
				vp2 = (SpatialVector) vlist.get(i2);
				// computes handedness (v0 x v2) . v1 for each vertex v1
				vo = vp0.cross(vp2);
				hand = dotProduct(vo, vp1);
				if ( hand >= 0. ) {
					ss[i1] = 1.0;
				} else {
					ss[i1] = -1.0;
				}

			}
			np = 0; // number of vert. with positive handedness
			for ( int i = 0; i < nv; i++ ) {
				if ( ss[i] > 0. )
					np++;
			}
			nm = n_remain - np;

			nlow = Math.min(np, nm);

			if ( nlow != 0 ) {
				if ( nlow == 1 ) { // only one concave vertex
					if ( np == 1 ) { // ix index of the vertex in the list
						for ( int k = 0; k < nv - 1; k++ ) {
							if ( Math.abs(ss[k] - 1.0) <= 1.e-12 ) {
								ix = k;
								break;
							}
						}
					} else {
						for ( int k = 0; k < nv - 1; k++ ) {
							if ( Math.abs(ss[k] + 1.0) <= 1.e-12 ) {
								ix = k;
								break;
							}
						}
					}

					// rotate pixel list to put that vertex in #0
					int n_rot = vlist.size() - ix;
					int ilast = vlist.size() - 1;
					for ( int k = 0; k < n_rot; k++ ) {
						SpatialVector temp = (SpatialVector) vlist.get(ilast);
						vlist.remove(ilast);
						vlist.add(0, (Object) temp);
					}
				}
				if ( nlow > 1 ) { // more than 1concave vertex
					System.out
							.println(" The polygon has more than one concave vertex");
					System.out.println(" The result is unpredictable");
				}
			}
		}
		/* fill the poligon, one triangle at a time */
		npix = (long) Nside2Npix(nside);
		while ( n_remain >= 3 ) {
			vp0 = (SpatialVector) vlist.get(0);
			vp1 = (SpatialVector) vlist.get(n_remain - 2);
			vp2 = (SpatialVector) vlist.get(n_remain - 1);

			/* find pixels within the triangle */
			ArrayList<Long> templist = new ArrayList<Long>();
			templist = query_triangle(nside, vp0, vp1, vp2, nest, inclusive);

			vvlist.addAll(templist);
			n_remain--;
		}
		/* make final pixel list */
		npix = vvlist.size();
		long[] pixels = new long[(int) npix];
		for ( int i = 0; i < npix; i++ ) {
			pixels[i] = vvlist.get(i).longValue();
		}
		Arrays.sort(pixels);
		int k = 0;
		res.add(k, new Long(pixels[0]));
		for ( int i = 1; i < pixels.length; i++ ) {
			if ( pixels[i] > pixels[i - 1] ) {
				k++;
				res.add(k, new Long(pixels[i]));
			}
		}

		return res;
	}
	
	/**
	 * Prints the vec.
	 * 
	 * @param vec the vec
	 */
	public void printVec(double[] vec) {
		System.out.print("[");
		for ( int i = 0; i < vec.length; i++ ) {
			System.out.print(vec[i] + " ");
		}
		System.out.println("]");
	}
	/**
	 * generates a list of pixels that lie inside a triangle defined by the
	 * three vertex vectors
	 * 
	 * @param nside
	 *            long map resolution parameter
	 * @param v1
	 *            Vector3d defines one vertex of the triangle
	 * @param v2
	 *            Vector3d another vertex
	 * @param v3
	 *            Vector3d yet another one
	 * @param nest
	 *            long 0 (default) RING numbering scheme, if set to 1 the NESTED
	 *            scheme will be used.
	 * @param inclusive
	 *            long 0 (default) only pixels whose centers are inside the
	 *            triangle will be listed, if set to 1 all pixels overlaping the
	 *            triangle will be listed
	 * @return ArrayList with pixel numbers
	 * @throws healpixException
	 *             if the triangle is degenerated
	 * @throws IllegalArgumentException
	 */
	public ArrayList<Long> query_triangle(long nside, SpatialVector v1,
			SpatialVector v2, SpatialVector v3, long nest, long inclusive)
			throws Exception {
		ArrayList<Long> res;
		res = new ArrayList<Long>();
		ArrayList<Long> listir;
		long npix, iz, irmin, irmax, n12, n123a, n123b, ndom = 0;
		boolean test1, test2, test3;
		double dth1, dth2, determ, sdet;
		double zmax, zmin, z1max, z1min, z2max, z2min, z3max, z3min;
		double z, tgth, st, offset, sin_off;
		double phi_pos, phi_neg;
		SpatialVector[] vv = new SpatialVector[3];
		SpatialVector[] vo = new SpatialVector[3];
		double[] sprod = new double[3];
		double[] sto = new double[3];
		double[] phi0i = new double[3];
		double[] tgthi = new double[3];
		double[] dc = new double[3];
		double[][] dom = new double[3][2];
		double[] dom12 = new double[4];
		double[] dom123a = new double[4];
		double[] dom123b = new double[4];
		double[] alldom = new double[6];
		double a_i, b_i, phi0, dphiring;
		long idom;
		boolean do_inclusive = false;
		boolean do_nest = false;
		String SID = "QUERY_TRIANGLE";
		long nsidesq = nside * nside;
		/*                                       */

		// System.out.println("in query_triangle");
		npix = Nside2Npix(nside);
		if ( npix < 0 ) {
			throw new IllegalArgumentException(SID
					+ " Nside should be power of 2 >0 and < " + ns_max);
		}
		if ( inclusive == 1 )
			do_inclusive = true;
		if ( nest == 1 )
			do_nest = true;
		vv[0] = new SpatialVector(v1);
//		vv[0].normalize();
		vv[1] = new SpatialVector(v2);
//		vv[1].normalize();
		vv[2] = new SpatialVector(v3);
//		vv[2].normalize();
//		printVec(vv[0].get());
//		printVec(vv[1].get());
//		printVec(vv[2].get());
		/*                                  */
		dth1 = 1.0 / ( 3.0 * nsidesq );
		dth2 = 2.0 / ( 3.0 * nside );
		/*
		 * determ = (v1 X v2) . v3 determines the left ( <0) or right (>0)
		 * handedness of the triangle
		 */
		SpatialVector vt = new SpatialVector(0., 0., 0.);
		vt = crossProduct(vv[0], vv[1]);
		determ = dotProduct(vt, vv[2]);

		if ( Math.abs(determ) < 1.0e-20 ) {
			throw new HealpixException(
					SID
							+ ": the triangle is degenerated - query cannot be performed");
		}
		if ( determ >= 0. ) { // The sign of determ
			sdet = 1.0;
		} else {
			sdet = -1.0;
		}

		sprod[0] = dotProduct(vv[1], vv[2]);
		sprod[1] = dotProduct(vv[2], vv[0]);
		sprod[2] = dotProduct(vv[0], vv[1]);
		/* vector ortogonal to the great circle containing the vertex doublet */

		vo[0] = crossProduct(vv[1], vv[2]);
		vo[1] = crossProduct(vv[2], vv[0]);
		vo[2] = crossProduct(vv[0], vv[1]);
		vo[0].normalize();
		vo[1].normalize();
		vo[2].normalize();
//		System.out.println("Orthogonal vectors:");
		
//		printVec(vo[0].get());
//		printVec(vo[1].get());
//		printVec(vo[2].get());
		/* test presence of poles in the triangle */
		zmax = -1.0;
		zmin = 1.0;
		test1 = ( vo[0].z() * sdet >= 0.0 ); // north pole in hemisphere
		// defined
		// by
		// 2-3
		test2 = ( vo[1].z() * sdet >= 0.0 ); // north pole in the hemisphere
		// defined
		// by 1-2
		test3 = ( vo[2].z() * sdet >= 0.0 ); // north pole in hemisphere
		// defined
		// by
		// 1-3
		if ( test1 && test2 && test3 )
			zmax = 1.0; // north pole in the triangle
		if ( ( !test1 ) && ( !test2 ) && ( !test3 ) )
			zmin = -1.0; // south pole in the triangle
		/* look for northenest and southernest points in the triangle */
		// ! look for northernest and southernest points in the triangle
		// ! node(1,2) = vector of norm=1, in the plane defined by (1,2) and
		// with z=0
		
		 boolean test1a = ((vv[2].z() - sprod[0] * vv[1].z()) >= 0.0); //
		 boolean test1b = ((vv[1].z() - sprod[0] * vv[2].z()) >= 0.0);
		 boolean test2a = ((vv[2].z() - sprod[1] * vv[0].z()) >= 0.0); //
		 boolean test2b = ((vv[0].z() - sprod[1] * vv[2].z()) >= 0.0);
		 boolean test3a = ((vv[1].z() - sprod[2] * vv[0].z()) >= 0.0); //
		 boolean test3b = ((vv[0].z() - sprod[2] * vv[1].z()) >= 0.0);
		/* sin of theta for orthogonal vector */
		for ( int i = 0; i < 3; i++ ) {
			sto[i] = Math.sqrt(( 1.0 - vo[i].z() ) * ( 1.0 + vo[i].z() ));
		}
		/*
		 * for each segment ( side of the triangle ) the extrema are either -
		 * -the 2 vertices 
		 * - one of the vertices and a point within the segment
		 */
		z1max = vv[1].z();
		z1min = vv[2].z();
		double zz;
//		segment 2-3
		if ( test1a == test1b ) {
			zz = sto[0];
			if ( ( vv[1].z() + vv[2].z() ) >= 0.0 ) {
				z1max = zz;
			} else {
				z1min = -zz;
			}
		}
		// segment 1-3
		z2max = vv[2].z();
		z2min = vv[0].z();
		if ( test2a == test2b ) {
			zz = sto[1];
			if ( ( vv[0].z() + vv[2].z() ) >= 0.0 ) {
				z2max = zz;
			} else {
				z2min = -zz;
			}
		}
		// segment 1-2
		z3max = vv[0].z();
		z3min = vv[1].z();
		if ( test3a == test3b ) {
			zz = sto[2];
			if ( ( vv[0].z() + vv[1].z() ) >= 0.0 ) {
				z3max = zz;
			} else {
				z3min = -zz;
			}
		}

		zmax = Math.max(Math.max(z1max, z2max), Math.max(z3max, zmax));
		zmin = Math.min(Math.min(z1min, z2min), Math.min(z3min, zmin));
		/*
		 * if we are inclusive, move upper point up, and lower point down, by a
		 * half pixel size
		 */
		offset = 0.0;
		sin_off = 0.0;
		if ( do_inclusive ) {
			offset = Constants.PI / ( 4.0 * nside ); // half pixel size
			sin_off = Math.sin(offset);
			zmax = Math.min(1.0, Math.cos(Math.acos(zmax) - offset));
			zmin = Math.max(-1.0, Math.cos(Math.acos(zmin) + offset));
		}

		irmin = RingNum(nside, zmax);
		irmax = RingNum(nside, zmin);

//		System.out.println("irmin = " + irmin + " irmax =" + irmax);

		/* loop on the rings */
		for ( int i = 0; i < 3; i++ ) {
			tgthi[i] = -1.0e30 * vo[i].z();
			phi0i[i] = 0.0;
		}
		for ( int j = 0; j < 3; j++ ) {
			if ( sto[j] > 1.0e-10 ) {
				tgthi[j] = -vo[j].z() / sto[j]; // - cotan(theta_orth)

				phi0i[j] = Math.atan2(vo[j].y(), vo[j].x()); // Should make
				// it
				// 0-2pi
				// ?
				/* Bring the phi0i to the [0,2pi] domain if need */

				 if ( phi0i[j] < 0.) {
					phi0i[j] = bm
							.MODULO(
									( Math.atan2(vo[j].y(), vo[j].x()) + Constants.twopi ),
									Constants.twopi); // [0-2pi]
				}
//				System.out.println("phi0i = " + phi0i[j] + " tgthi = "
//						+ tgthi[j]);
			}
		}
		//MOD(ATAN2(X,Y) + TWOPI, TWOPI) : ATAN2 in 0-2pi
		/*
		 * the triangle boundaries are geodesics: intersection of the sphere
		 * with plans going through (0,0,0) if we are inclusive, the boundaries
		 * are the intersection of the sphere with plains pushed outward by
		 * sin(offset)
		 */
		boolean found = false;
		for ( iz = irmin; iz <= irmax; iz++ ) {
			found = false;
			if ( iz <= nside - 1 ) { // North polar cap
				z = 1.0 - iz * iz * dth1;
			} else if ( iz <= 3 * nside ) { // tropical band + equator
				z = ( 2.0 * nside - iz ) * dth2;
			} else {
				z = -1.0 + ( 4.0 * nside - iz ) * ( 4.0 * nside - iz ) * dth1;
			}

			/* computes the 3 intervals described by the 3 great circles */
			st = Math.sqrt(( 1.0 - z ) * ( 1.0 + z ));
			tgth = z / st; // cotan(theta_ring)
			for ( int j = 0; j < 3; j++ ) {
				dc[j] = tgthi[j] * tgth - sdet * sin_off
						/ ( ( sto[j] + 1.0e-30 ) * st ) ;

			}
			for ( int k = 0; k < 3; k++ ) {
				if ( dc[k] * sdet <= -1.0 ) { // the whole iso-latitude ring
					// is on
					// right side of the great circle
					dom[k][0] = 0.0;
					dom[k][1] = Constants.twopi;
				} else if ( dc[k] * sdet >= 1.0 ) { // all on the wrong side
					dom[k][0] = -1.000001 * ( k + 1 );
					dom[k][1] = -1.0 * ( k + 1 );
				} else { // some is good some is bad
					phi_neg = phi0i[k] - ( Math.acos(dc[k]) * sdet );
					phi_pos = phi0i[k] + ( Math.acos(dc[k]) * sdet );
					//					
					 if ( phi_pos < 0. )
						phi_pos += Constants.twopi;
					if ( phi_neg < 0. )
						phi_neg += Constants.twopi;
					//
					dom[k][0] = bm.MODULO(phi_neg, Constants.twopi);
					dom[k][1] = bm.MODULO(phi_pos, Constants.twopi);
//					double domk0 = (phi0i[k] - ( Math.acos(dc[k]) * sdet )) % Constants.twopi;
//					double domk1 = (phi0i[k] + ( Math.acos(dc[k]) * sdet )) %  Constants.twopi;
				}
//				System.out.println("dom["+k+"][0] = " + dom[k][0] + " [1]= "
//						+ dom[k][1]);
				//

			}
			/* identify the intersections (0,1,2 or 3) of the 3 intervals */

			dom12 = intrs_intrv(dom[0], dom[1]);
			n12 = dom12.length / 2;
			if ( n12 != 0 ) {
				if ( n12 == 1 ) {
					dom123a = intrs_intrv(dom[2], dom12);
					n123a = dom123a.length / 2;

					if ( n123a == 0 )
						found = true;
					if ( !found ) {
						for ( int l = 0; l < dom123a.length; l++ ) {
							alldom[l] = dom123a[l];
						}

						ndom = n123a; // 1 or 2
					}
				}
				if ( !found ) {
					if ( n12 == 2 ) {
						double[] tmp = { dom12[0], dom12[1] };
						dom123a = intrs_intrv(dom[2], tmp);
						double[] tmp1 = { dom12[2], dom12[3] };
						dom123b = intrs_intrv(dom[2], tmp1);
						n123a = dom123a.length / 2;
						n123b = dom123b.length / 2;
						ndom = n123a + n123b; // 0, 1, 2 or 3

						if ( ndom == 0 )
							found = true;
						if ( !found ) {
							if ( n123a != 0 ) {
								for ( int l = 0; l < 2 * n123a; l++ ) {
									alldom[l] = dom123a[l];
								}
							}
							if ( n123b != 0 ) {
								for ( int l = 0; l < 2 * n123b; l++ ) {
									alldom[(int) ( l + 2 * n123a )] = dom123b[l];
								}
							}
							if ( ndom > 3 ) {
								throw new HealpixException(SID
										+ ": too many intervals found");
							}
						}
					}
				}
				if ( !found ) {
					for ( idom = 0; idom < ndom; idom++ ) {
						a_i = alldom[(int) ( 2 * idom )];
						b_i = alldom[(int) ( 2 * idom + 1 )];
						phi0 = ( a_i + b_i ) * 0.5;
						dphiring = (b_i - a_i) * 0.5;
						if ( dphiring < 0.0 ) {
							phi0 += Constants.PI;
							dphiring += Constants.PI;
						}
						/* finds pixels in the triangle on that ring */
						listir = InRing(nside, iz, phi0, dphiring, do_nest);
//						ArrayList<Long> listir2 = InRing(nside, iz, phi0, dphiring, do_nest);						
						res.addAll(listir);

					}
				}
			}

		}
		return res;
	}

	/**
	 * computes the intersection di of 2 intervals d1 (= [a1,b1]) and d2 (=
	 * [a2,b2]) on the periodic domain (=[A,B] where A and B arbitrary) ni is
	 * the resulting number of intervals (0,1, or 2) if a1 <b1 then d1 = {x |a1 <=
	 * x <= b1} if a1>b1 then d1 = {x | a1 <=x <= B U A <=x <=b1}
	 * 
	 * @param d1
	 *            double[] first interval
	 * @param d2
	 *            double[] second interval
	 * @return double[] one or two intervals intersections
	 */
	public double[] intrs_intrv(double[] d1, double[] d2) {
		double[] res;
		double epsilon = 1.0e-10;
		double[] dk;
		double[] di = { 0. };
		int ik = 0;
		boolean tr12, tr21, tr34, tr43, tr13, tr31, tr24, tr42, tr14, tr32;
		/*                                             */

		tr12 = ( d1[0] < d1[1] + epsilon );
		tr21 = !tr12; // d1[0] >= d1[1]
		tr34 = ( d2[0] < d2[1] + epsilon );
		tr43 = !tr34; // d2[0]>d2[1]
		tr13 = ( d1[0] < d2[0] + epsilon ); // d2[0] can be in interval
		tr31 = !tr13; // d1[0] in longerval
		tr24 = ( d1[1] < d2[1] + epsilon ); // d1[1] upper limit
		tr42 = !tr24; // d2[1] upper limit
		tr14 = ( d1[0] < d2[1] + epsilon ); // d1[0] in interval
		tr32 = ( d2[0] < d1[1] + epsilon ); // d2[0] in interval

		ik = 0;
		dk = new double[] { -1.0e9, -1.0e9, -1.0e9, -1.0e9 };
		/* d1[0] lower limit case 1 */
		if ( ( tr34 && tr31 && tr14 ) || ( tr43 && ( tr31 || tr14 ) ) ) {
			ik++; // ik = 1;
			dk[ik - 1] = d1[0]; // a1

		}
		/* d2[0] lower limit case 1 */
		if ( ( tr12 && tr13 && tr32 ) || ( tr21 && ( tr13 || tr32 ) ) ) {
			ik++; // ik = 1
			dk[ik - 1] = d2[0]; // a2

		}
		/* d1[1] upper limit case 2 */
		if ( ( tr34 && tr32 && tr24 ) || ( tr43 && ( tr32 || tr24 ) ) ) {
			ik++; // ik = 2
			dk[ik - 1] = d1[1]; // b1

		}
		/* d2[1] upper limit case 2 */
		if ( ( tr12 && tr14 && tr42 ) || ( tr21 && ( tr14 || tr42 ) ) ) {
			ik++; // ik = 2
			dk[ik - 1] = d2[1]; // b2

		}
		di = new double[1];
		di[0] = 0.;
		switch ( ik ) {

			case 2:
				di = new double[2];

				di[0] = dk[0] - epsilon;
				di[1] = dk[1] + epsilon;
				break;
			case 4:

				di = new double[4];
				di[0] = dk[0] - epsilon;
				di[1] = dk[3] + epsilon;
				di[2] = dk[1] - epsilon;
				di[3] = dk[2] + epsilon;
				break;
		}
		res = di;

		return res;
	}

	/**
	 * finds pixels having a colatitude (measured from North pole) : theta1 <
	 * colatitude < theta2 with o <= theta1 < theta2 <= Pi if theta2 < theta1
	 * then pixels with 0 <= colatitude < theta2 or theta1 < colatitude < Pi are
	 * returned
	 * 
	 * @param nside
	 *            long the map resolution parameter
	 * @param theta1
	 *            lower edge of the colatitude
	 * @param theta2
	 *            upper edge of the colatitude
	 * @param nest
	 *            long if = 1 result is in NESTED scheme
	 * @return ArrayList of pixel numbers (long)
	 * @throws IllegalArgumentException
	 */
	public ArrayList<Long> query_strip(long nside, double theta1,
			double theta2, long nest) throws Exception {
		ArrayList<Long> res = new ArrayList<Long>();
		ArrayList<Long> listir = new ArrayList<Long>();
		long npix, nstrip;
		long iz, irmin, irmax;
		int is;
		double phi0, dphi;
		double[] colrange = new double[4];
		boolean nest_flag = false;
		String SID = " QUERY_STRIP";
		/* ---------------------------------------- */
		npix = Nside2Npix(nside);
		if ( nest == 1 )
			nest_flag = true;
		if ( npix < 0 ) {
			throw new IllegalArgumentException(SID
					+ " Nside should be power of 2");
		}
		if ( ( theta1 < 0.0 || theta1 > Constants.PI )
				|| ( theta2 < 0.0 || theta2 > Constants.PI ) ) {
			throw new IllegalArgumentException(SID
					+ " Illegal value of theta1, theta2");
		}
		if ( theta1 <= theta2 ) {
			nstrip = 1;
			colrange[0] = theta1;
			colrange[1] = theta2;
		} else {
			nstrip = 2;
			colrange[0] = 0.0;
			colrange[1] = theta2;
			colrange[2] = theta1;
			colrange[3] = Constants.PI;
		}
		/* loops on strips */
		for ( is = 0; is < nstrip; is++ ) {
			irmin = RingNum(nside, Math.cos(colrange[2 * is]));
			irmax = RingNum(nside, Math.cos(colrange[2 * is + 1]));
			/* loop on ring number */
			for ( iz = irmin; iz <= irmax; iz++ ) {
				phi0 = 0.;
				dphi = Constants.PI;
				listir = InRing(nside, iz, phi0, dphi, nest_flag);//InRing(nside, iz, phi0, dphi, nest_flag);
				res.addAll(listir);
			}
		}
		return res;
	}

	/**
	 * Spatial vector2 vector3d.
	 * 
	 * @param vec the vec
	 * 
	 * @return the vector3d
	 */
	Vector3d SpatialVector2Vector3d(SpatialVector vec) {
		Vector3d res = new Vector3d();
		res.set(new double[] { vec.x(), vec.y(), vec.z() });
		return res;
	}

	/**
	 * returns 7 or 8 neighbours of any pixel in the nested scheme The neighbours
	 * are ordered in the following way: First pixel is the one to the south (
	 * the one west of the south direction is taken for pixels that don't have a
	 * southern neighbour). From then on the neighbors are ordered in the
	 * clockwise direction.
	 * 
	 * @param nside the map resolution
	 * @param ipix long pixel number
	 * @return ArrayList
	 * @throws Exception 
	 * @throws IllegalArgumentException
	 */
	public ArrayList<Long> neighbours_nest(int nside, long ipix) throws Exception  {
		ArrayList<Long> res = new ArrayList<Long>(8);
		int ipf, ipo, ix, ixm, ixp, iy, iym, iyp, ixo, iyo;
		int face_num, other_face;
		int ia, ib, ibp, ibm, ib2,  nsidesq;
        int icase;
		long local_magic1, local_magic2;
		long arb_const = 0;
		Point ixiy = new Point();
		Point ixoiyo = new Point();
		/* fill the pixel list with 0 */
		res.add(0, new Long(0));
		res.add(1, new Long(0));
		res.add(2, new Long(0));
		res.add(3, new Long(0));
		res.add(4, new Long(0));
		res.add(5, new Long(0));
		res.add(6, new Long(0));
		res.add(7, new Long(0));
		icase = 0;
		nsidesq = nside * nside;
		
		local_magic1 = (nsidesq - 1) / 3;
		local_magic2 = 2 * local_magic1;
		face_num = (int) ( ipix / nsidesq );
		ipf = (int) bm.MODULO(ipix, nsidesq); // Pixel number in face
		ixiy = pix2xy_nest((int) ipf);
		ix = ixiy.x;
		iy = ixiy.y;
		//
		ixm = ix - 1;
		ixp = ix + 1;
		iym = iy - 1;
		iyp = iy + 1;

		icase = 0; // inside the face

		/* exclude corners */
		if (ipf == local_magic2 && icase == 0)
			icase = 5; // West corner
		if (ipf == (nsidesq - 1) && icase == 0)
			icase = 6; // North corner
		if (ipf == 0 && icase == 0)
			icase = 7; // South corner
		if (ipf == local_magic1 && icase == 0)
			icase = 8; // East corner

		/* detect edges */
		if ((ipf & local_magic1) == local_magic1 && icase == 0)
			icase = 1; // NorthEast
		if ((ipf & local_magic1) == 0 && icase == 0)
			icase = 2; // SouthWest
		if ((ipf & local_magic2) == local_magic2 && icase == 0)
			icase = 3; // NorthWest
		if ((ipf & local_magic2) == 0 && icase == 0)
			icase = 4; // SouthEast

		/* inside a face */
		if (icase == 0) {
			res.set(0, new Long( xy2pix_nest(ixm, iym, face_num)));
			res.set(1, new Long( xy2pix_nest(ixm, iy, face_num)));
			res.set(2, new Long( xy2pix_nest(ixm, iyp, face_num)));
			res.set(3, new Long( xy2pix_nest(ix, iyp, face_num)));
			res.set(4, new Long( xy2pix_nest(ixp, iyp, face_num)));
			res.set(5, new Long( xy2pix_nest(ixp, iy, face_num)));
			res.set(6, new Long( xy2pix_nest(ixp, iym, face_num)));
			res.set(7, new Long( xy2pix_nest(ix, iym, face_num)));
			return res;
		}
		/*                 */
		ia = face_num / 4; // in [0,2]
		ib = (int) bm.MODULO(face_num, 4); // in [0,3]
		ibp = (int) bm.MODULO(ib + 1, 4);
		ibm = (int) bm.MODULO(ib + 4 - 1, 4);
		ib2 = (int) bm.MODULO(ib + 2, 4);

		if (ia == 0) { // North pole region
			switch (icase) {
			case 1: // north-east edge
				other_face = 0 + ibp;
				res.set(0, new Long( xy2pix_nest(ixm, iym, face_num)));
				res.set(1, new Long( xy2pix_nest(ixm, iy, face_num)));
				res.set(2, new Long( xy2pix_nest(ixm, iyp, face_num)));
				res.set(3, new Long( xy2pix_nest(ix, iyp, face_num)));
				res.set(7, new Long( xy2pix_nest(ix, iym, face_num)));
				ipo = (int) bm.MODULO(bm.swapLSBMSB( ipf), nsidesq);
				ixoiyo = pix2xy_nest((int) ipo);
				ixo = ixoiyo.x;
				iyo = ixoiyo.y;
				res.set(4, new Long( xy2pix_nest(ixo + 1, iyo,
						other_face)));
				res.set(5, new Long( (other_face * nsidesq + ipo)));
				res.set(6, new Long( xy2pix_nest(ixo - 1, iyo,
						other_face)));
				break;
			case 2: // SouthWest edge
				other_face = 4 + ib;
				ipo = (int) bm.MODULO(bm.invLSB( ipf), nsidesq); // SW-NE flip
				ixoiyo = pix2xy_nest((int) ipo);
				ixo = ixoiyo.x;
				iyo = ixoiyo.y;
				res.set(0, new Long( xy2pix_nest(ixo, iyo - 1,
						other_face)));
				res.set(1, new Long( (other_face * nsidesq + ipo)));
				res.set(2, new Long( xy2pix_nest(ixo, iyo + 1,
						other_face)));
				res.set(3, new Long( xy2pix_nest(ix, iyp, face_num)));
				res.set(4, new Long( xy2pix_nest(ixp, iyp, face_num)));
				res.set(5, new Long( xy2pix_nest(ixp, iy, face_num)));
				res.set(6, new Long( xy2pix_nest(ixp, iym, face_num)));
				res.set(7, new Long( xy2pix_nest(ix, iym, face_num)));
				break;
			case 3: // NorthWest edge
				other_face = 0 + ibm;
				ipo = (int) bm.MODULO(bm.swapLSBMSB( ipf), nsidesq); // E-W flip
				ixoiyo = pix2xy_nest((int) ipo);
				ixo = ixoiyo.x;
				iyo = ixoiyo.y;
				res.set(0, new Long( xy2pix_nest(ixm, iym, face_num)));
				res.set(1, new Long( xy2pix_nest( ixm, iy, face_num)));
				res.set(2, new Long( xy2pix_nest( ixo, iyo - 1,
						other_face)));
				res.set(3, new Long( (other_face * nsidesq + ipo)));
				res.set(4, new Long( xy2pix_nest( ixo, iyo + 1,
						other_face)));
				res.set(5, new Long( xy2pix_nest( ixp, iy, face_num)));
				res.set(6, new Long( xy2pix_nest(ixp, iym, face_num)));
				res.set(7, new Long( xy2pix_nest( ix, iym, face_num)));
				break;
			case 4: // SouthEast edge
				other_face = 4 + ibp;
				ipo = (int) bm.MODULO(bm.invMSB( ipf), nsidesq); // SE-NW flip
				ixoiyo = pix2xy_nest((int) ipo);
				ixo = ixoiyo.x;
				iyo = ixoiyo.y;
				res.set(0, new Long( xy2pix_nest(ixo - 1, iyo,
						other_face)));
				res.set(1, new Long( xy2pix_nest(ixm, iy, face_num)));
				res.set(2, new Long( xy2pix_nest(ixm, iyp, face_num)));
				res.set(3, new Long( xy2pix_nest(ix, iyp, face_num)));
				res.set(4, new Long( xy2pix_nest(ixp, iyp, face_num)));
				res.set(5, new Long( xy2pix_nest(ixp, iy, face_num)));
				res.set(6, new Long( xy2pix_nest(ixo + 1, iyo,
						other_face)));
				res.set(7, new Long( (other_face * nsidesq + ipo)));
				break;
			case 5: // West corner
				other_face = 4 + ib;
				arb_const = other_face * nsidesq + nsidesq - 1;
				res.set(0, new Long( (arb_const - 2)));
				res.set(1, new Long( arb_const));
				other_face = 0 + ibm;
				arb_const = other_face * nsidesq + local_magic1;
				res.set(2, new Long( arb_const));
				res.set(3, new Long( (arb_const + 2)));
				res.set(4, new Long( (ipix + 1)));
				res.set(5, new Long( (ipix - 1)));
				res.set(6, new Long( (ipix - 2)));
				res.remove(7);
				break;
			case 6: //  North corner
				other_face = 0 + ibm;
				res.set(0, new Long( (ipix - 3)));
				res.set(1, new Long( (ipix - 1)));
				arb_const = other_face * nsidesq + nsidesq - 1;
				res.set(2, new Long( (arb_const - 2)));
				res.set(3, new Long( arb_const));
				other_face = 0 + ib2;
				res.set(4, new Long( (other_face * nsidesq + nsidesq - 1)));
				other_face = 0 + ibp;
				arb_const = other_face * nsidesq + nsidesq - 1;
				res.set(5, new Long( arb_const));
				res.set(6, new Long( (arb_const - 1)));
				res.set(7, new Long( (ipix - 2)));
				break;
			case 7: // South corner
				other_face = 8 + ib;
				res.set(0, new Long( (other_face * nsidesq + nsidesq - 1)));
				other_face = 4 + ib;
				arb_const = other_face * nsidesq + local_magic1;
				res.set(1, new Long( arb_const));
				res.set(2, new Long( (arb_const + 2)));
				res.set(3, new Long( (ipix + 2)));
				res.set(4, new Long( (ipix + 3)));
				res.set(5, new Long( (ipix + 1)));
				other_face = 4 + ibp;
				arb_const = other_face * nsidesq + local_magic2;
				res.set(6, new Long( (arb_const + 1)));
				res.set(7, new Long( arb_const));
				break;
			case 8: // East corner
				other_face = 0 + ibp;
				res.set(1, new Long( (ipix - 1)));
				res.set(2, new Long( (ipix + 1)));
				res.set(3, new Long( (ipix + 2)));
				arb_const = other_face * nsidesq + local_magic2;
				res.set(4, new Long( (arb_const + 1)));
				res.set(5, new Long(( arb_const)));
				other_face = 4 + ibp;
				arb_const = other_face * nsidesq + nsidesq - 1;
				res.set(0, new Long( (arb_const - 1)));
				res.set(6, new Long( arb_const));
				res.remove(7);
				break;
			}
		} else if (ia == 1) { // Equatorial region
			switch (icase) {
			case 1: // north-east edge
				other_face = 0 + ib;
				res.set(0, new Long( xy2pix_nest( ixm, iym, face_num)));
				res.set(1, new Long( xy2pix_nest(  ixm, iy, face_num)));
				res.set(2, new Long( xy2pix_nest(  ixm, iyp, face_num)));
				res.set(3, new Long( xy2pix_nest(  ix, iyp, face_num)));
				res.set(7, new Long( xy2pix_nest(  ix, iym, face_num)));
				ipo = (int) bm.MODULO(bm.invLSB( ipf), nsidesq); // NE-SW flip
				ixoiyo = pix2xy_nest(  ipo);
				ixo = ixoiyo.x;
				iyo = ixoiyo.y;
				res.set(4, new Long( xy2pix_nest(  ixo, iyo + 1,
						other_face)));
				res.set(5, new Long( (other_face * nsidesq + ipo)));
				res.set(6, new Long( xy2pix_nest(  ixo, iyo - 1,
						other_face)));
				break;
			case 2: // SouthWest edge
				other_face = 8 + ibm;
				ipo = (int) bm.MODULO(bm.invLSB( ipf), nsidesq); // SW-NE flip
				ixoiyo = pix2xy_nest(  ipo);
				ixo = ixoiyo.x;
				iyo = ixoiyo.y;
				res.set(0, new Long( xy2pix_nest(  ixo, iyo - 1,
						other_face)));
				res.set(1, new Long((other_face * nsidesq + ipo)));
				res.set(2, new Long( xy2pix_nest( ixo, iyo + 1,
						other_face)));
				res.set(3, new Long( xy2pix_nest( ix, iyp, face_num)));
				res.set(4, new Long( xy2pix_nest( ixp, iyp, face_num)));
				res.set(5, new Long( xy2pix_nest( ixp, iy, face_num)));
				res.set(6, new Long( xy2pix_nest(  ixp, iym, face_num)));
				res.set(7, new Long( xy2pix_nest(  ix, iym, face_num)));
				break;
			case 3: // NortWest edge
				other_face = 0 + ibm;
				ipo = (int) bm.MODULO(bm.invMSB( ipf), nsidesq); // NW-SE flip
				ixoiyo = pix2xy_nest(  ipo);
				ixo = ixoiyo.x;
				iyo = ixoiyo.y;
				res.set(2, new Long( xy2pix_nest( ixo - 1, iyo,
						other_face)));
				res.set(3, new Long( (other_face * nsidesq + ipo)));
				res.set(4, new Long( xy2pix_nest( ixo + 1, iyo,
						other_face)));
				res.set(0, new Long( xy2pix_nest(  ixm, iym, face_num)));
				res.set(1, new Long( xy2pix_nest(  ixm, iy, face_num)));
				res.set(5, new Long( xy2pix_nest(  ixp, iy, face_num)));
				res.set(6, new Long( xy2pix_nest(  ixp, iym, face_num)));
				res.set(7, new Long(xy2pix_nest(  ix, iym, face_num)));
				break;
			case 4: // SouthEast edge
				other_face = 8 + ib;
				ipo = (int) bm.MODULO(bm.invMSB( ipf), nsidesq); // SE-NW flip
				ixoiyo = pix2xy_nest(  ipo);
				ixo = ixoiyo.x;
				iyo = ixoiyo.y;
				res.set(0, new Long( xy2pix_nest(  ixo - 1, iyo,
						other_face)));
				res.set(1, new Long( xy2pix_nest( ixm, iy, face_num)));
				res.set(2, new Long( xy2pix_nest(  ixm, iyp, face_num)));
				res.set(3, new Long( xy2pix_nest(  ix, iyp, face_num)));
				res.set(4, new Long( xy2pix_nest( ixp, iyp, face_num)));
				res.set(5, new Long( xy2pix_nest(  ixp, iy, face_num)));
				res.set(6, new Long( xy2pix_nest(  ixo + 1, iyo,
						other_face)));
				res.set(7, new Long( (other_face * nsidesq + ipo)));
				break;
			case 5: // West corner
				other_face = 8 + ibm;
				arb_const = other_face * nsidesq + nsidesq - 1;
				res.set(0, new Long( (arb_const - 2)));
				res.set(1, new Long( arb_const));
				other_face = 4 + ibm;
				res.set(2, new Long( (other_face * nsidesq + local_magic1)));
				other_face = 0 + ibm;
				arb_const = other_face * nsidesq;
				res.set(3, new Long( arb_const));
				res.set(4, new Long( (arb_const + 1)));
				res.set(5, new Long( (ipix + 1)));
				res.set(6, new Long( (ipix - 1)));
				res.set(7, new Long( (ipix - 2)));
				break;
			case 6: //  North corner
				other_face = 0 + ibm;
				res.set(0, new Long( (ipix - 3)));
				res.set(1, new Long( (ipix - 1)));
				arb_const = other_face * nsidesq + local_magic1;
				res.set(2, new Long( (arb_const - 1)));
				res.set(3, new Long( arb_const));
				other_face = 0 + ib;
				arb_const = other_face * nsidesq + local_magic2;
				res.set(4, new Long( arb_const));
				res.set(5, new Long( (arb_const - 2)));
				res.set(6, new Long( (ipix - 2)));
				res.remove(7);
				break;
			case 7: // South corner
				other_face = 8 + ibm;
				arb_const = other_face * nsidesq + local_magic1;
				res.set(0, new Long( arb_const));
				res.set(1, new Long( (arb_const + 2)));
				res.set(2, new Long( (ipix + 2)));
				res.set(3, new Long( (ipix + 3)));
				res.set(4, new Long( (ipix + 1)));
				other_face = 8 + ib;
				arb_const = other_face * nsidesq + local_magic2;
				res.set(5, new Long( (arb_const + 1)));
				res.set(6, new Long( arb_const));
				res.remove(7);
				break;
			case 8: // East corner
				other_face = 8 + ib;
				arb_const = other_face * nsidesq + nsidesq - 1;
				res.set(0, new Long( (arb_const - 1)));
				res.set(1, new Long( (ipix - 1)));
				res.set(2, new Long( (ipix + 1)));
				res.set(3, new Long( (ipix + 2)));
				res.set(7, new Long( arb_const));
				other_face = 0 + ib;
				arb_const = other_face * nsidesq;
				res.set(4, new Long( (arb_const + 2)));
				res.set(5, new Long( arb_const));
				other_face = 4 + ibp;
				res.set(6, new Long( (other_face * nsidesq + local_magic2)));
				break;
			}
		} else { // South pole region
			switch (icase) {
			case 1: // North-East edge
				other_face = 4 + ibp;
				res.set(0, new Long( xy2pix_nest( ixm, iym, face_num)));
				res.set(1, new Long( xy2pix_nest( ixm, iy, face_num)));
				res.set(2, new Long( xy2pix_nest(  ixm, iyp, face_num)));
				res.set(3, new Long( xy2pix_nest(  ix, iyp, face_num)));
				res.set(7, new Long( xy2pix_nest(  ix, iym, face_num)));
				ipo = (int) bm.MODULO(bm.invLSB( ipf), nsidesq); // NE-SW flip
				ixoiyo = pix2xy_nest( ipo);
				ixo = ixoiyo.x;
				iyo = ixoiyo.y;
				res.set(4, new Long( xy2pix_nest(  ixo, iyo + 1,
						other_face)));
				res.set(5, new Long( (other_face * nsidesq + ipo)));
				res.set(6, new Long( xy2pix_nest( ixo, iyo - 1,
						other_face)));
				break;
			case 2: // SouthWest edge
				other_face = 8 + ibm;
				ipo = (int) bm.MODULO(bm.swapLSBMSB( ipf), nsidesq); // W-E flip
				ixoiyo = pix2xy_nest(  ipo);
				ixo = ixoiyo.x;
				iyo = ixoiyo.y;
				res.set(0, new Long( xy2pix_nest(  ixo - 1, iyo,
						other_face)));
				res.set(1, new Long( (other_face * nsidesq + ipo)));
				res.set(2, new Long( xy2pix_nest(  ixo + 1, iyo,
						other_face)));
				res.set(3, new Long( xy2pix_nest(  ix, iyp, face_num)));
				res.set(4, new Long( xy2pix_nest(  ixp, iyp, face_num)));
				res.set(5, new Long( xy2pix_nest(  ixp, iy, face_num)));
				res.set(6, new Long(xy2pix_nest(  ixp, iym, face_num)));
				res.set(7, new Long( xy2pix_nest(  ix, iym, face_num)));
				break;
			case 3: // NorthWest edge
				other_face = 4 + ib;
				ipo = (int) bm.MODULO(bm.invMSB( ipf), nsidesq);
				ixoiyo = pix2xy_nest(  ipo);
				ixo = ixoiyo.x;
				iyo = ixoiyo.y;
				res.set(0, new Long( xy2pix_nest(  ixm, iym, face_num)));
				res.set(1, new Long( xy2pix_nest(  ixm, iy, face_num)));
				res.set(2, new Long( xy2pix_nest(  ixo - 1, iyo,
						other_face)));
				res.set(3, new Long( (other_face * nsidesq + ipo)));
				res.set(4, new Long( xy2pix_nest(  ixo + 1, iyo,
						other_face)));
				res.set(5, new Long( xy2pix_nest(  ixp, iy, face_num)));
				res.set(6, new Long( xy2pix_nest(  ixp, iym, face_num)));
				res.set(7, new Long( xy2pix_nest(  ix, iym, face_num)));
				break;
			case 4: // SouthEast edge
				other_face = 8 + ibp;
				ipo = (int) bm.MODULO(bm.swapLSBMSB( ipf), nsidesq); // SE-NW
				// flip
				ixoiyo = pix2xy_nest(  ipo);
				ixo = ixoiyo.x;
				iyo = ixoiyo.y;
				res.set(0, new Long( xy2pix_nest(  ixo, iyo - 1,
						other_face)));
				res.set(1, new Long( xy2pix_nest(  ixm, iy, face_num)));
				res.set(2, new Long( xy2pix_nest(  ixm, iyp, face_num)));
				res.set(3, new Long( xy2pix_nest(  ix, iyp, face_num)));
				res.set(4, new Long( xy2pix_nest(  ixp, iyp, face_num)));
				res.set(5, new Long( xy2pix_nest(  ixp, iy, face_num)));
				res.set(6, new Long( xy2pix_nest(  ixo, iyo + 1,
						other_face)));
				res.set(7, new Long( (other_face * nsidesq + ipo)));
				break;
			case 5: // West corner
				other_face = 8 + ibm;
				arb_const = other_face * nsidesq + local_magic1;
				res.set(0, new Long( (arb_const - 2)));
				res.set(1, new Long( arb_const));
				other_face = 4 + ib;
				res.set(2, new Long( (other_face * nsidesq)));
				res.set(3, new Long( (other_face * nsidesq + 1)));
				res.set(4, new Long( (ipix + 1)));
				res.set(5, new Long( (ipix - 1)));
				res.set(6, new Long( (ipix - 2)));
				res.remove(7);
				break;
			case 6: //  North corner
				other_face = 4 + ib;
				res.set(0, new Long( (ipix - 3)));
				res.set(1, new Long((ipix - 1)));
				arb_const = other_face * nsidesq + local_magic1;
				res.set(2, new Long( (arb_const - 1)));
				res.set(3, new Long( arb_const));
				other_face = 0 + ib;
				res.set(4, new Long( (other_face * nsidesq)));
				other_face = 4 + ibp;
				arb_const = other_face * nsidesq + local_magic2;
				res.set(5, new Long( arb_const));
				res.set(6, new Long( (arb_const - 2)));
				res.set(7, new Long( (ipix - 2)));
				break;
			case 7: // South corner
				other_face = 8 + ib2;
				res.set(0, new Long( (other_face * nsidesq)));
				other_face = 8 + ibm;
				arb_const = other_face * nsidesq;
				res.set(1, new Long( arb_const));
				res.set(2, new Long( (arb_const + 1)));
				res.set(3, new Long( (ipix + 2)));
				res.set(4, new Long( (ipix + 3)));
				res.set(5, new Long( (ipix + 1)));
				other_face = 8 + ibp;
				arb_const = other_face * nsidesq;
				res.set(6, new Long( (arb_const + 2)));
				res.set(7, new Long( arb_const));
				break;
			case 8: // East corner
				other_face = 8 + ibp;
				res.set(1, new Long( (ipix - 1)));
				res.set(2, new Long( (ipix + 1)));
				res.set(3, new Long( (ipix + 2)));
				arb_const = other_face * nsidesq + local_magic2;
				res.set(6, new Long( arb_const));
				res.set(0, new Long( (arb_const - 2)));
				other_face = 4 + ibp;
				arb_const = other_face * nsidesq;
				res.set(4, new Long( (arb_const + 2)));
				res.set(5, new Long( arb_const));
				res.remove(7);
				break;
			}
		}
		return res;
	}

	/*
	 * return the parent PIXEL of a given pixel at some higher NSIDE. 
	 * One must also provide the nsode of the given pixel as otherwise it
	 * can not be known.
	 * 
	 * This only makes sense for Nested Scheme.
	 * This is basically a simple bit shift in the difference
	 * of number of bits between the two NSIDEs. 
	 * 
	 * @param child  the pixel 
	 * @param childnside nside of the pixel
	 * @param requirednside nside to upgrade to
	 * 
	 * @return the new pixel number
 	 */
	static public long parentAt(long child, int childnside, int requirednside) throws Exception{
	    // nside is the number of bits .. 
		if (childnside < requirednside) {
			throw new Exception ("Parent ("+requirednside+
					") should have smaller NSIDE than Child("+childnside+")");
		}
		long ppix =0;
		
		// number of bits in aid depdens on the depth of the nside

		int bitdiff = bitdiff(requirednside, childnside); 
	    ppix = child >> bitdiff;
    	return ppix;	 		
	}

	/**
	 * return difference of number of bits in pixels of two nsides.
	 * @param nside1
	 * @param nside2
	 * @return  number of bits difference between the pixel ids.
	 */
	public static int bitdiff(int nside1, int nside2){
		int pbits = 2;
		int childnside=nside2;
		int parentnside=nside1;
		if (nside1>=nside2){
			childnside=nside1;
			parentnside=nside2;
		}
		int tnside = 2;
		while (tnside < parentnside) {
			pbits+=2;
			tnside=tnside<<1 ;// next power of 2
		}
		// child is deeper 
		int cbits = pbits;
		while (tnside < childnside) {
			cbits+=2;
			tnside=tnside<<1 ;// next power of 2
		}
		return (cbits- pbits);//  
		
	}
	/**
	 * for a given pixel list all children pixels for it. 
	 * This is simply a matter of shifting the pixel number left by
	 * the difference in NSIDE bits and then listing all numbers 
	 * which fill the empty bits. 
	 * 
	 * BEWARE - not chcking you ar enot trying to go too DEEP. 
	 * 
	 * @param nside  nside of pix
	 * @param pix  the pixel 
	 * @param requiredNside  the nside you want the children at
	 * @return
	 */
	public static long[] getChildrenAt(int nside, long pix, int requiredNside) throws Exception{
	 
		if (nside >= requiredNside){
			throw new Exception("The requirend NSIDE should be greater than the pix NSIDE");
		}
		int bitdiff=bitdiff(nside,requiredNside);
		int numpix = bitdiff<<1;// square num bits is count of pix
		long[] pixlist= new long[numpix];
		long ppix=pix<<bitdiff; // shift the current pix over
		// nopw just keep adding to it ..
		for (int i=0;i < numpix; i++){
			pixlist[i]=ppix+i;
		}
		return pixlist;
	}
}
