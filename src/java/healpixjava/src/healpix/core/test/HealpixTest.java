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
package healpix.core.test;

import healpix.core.AngularPosition;
import healpix.core.HealpixIndex;
import healpix.core.base.BitManipulation;
import healpix.tools.Constants;
import healpix.tools.SpatialVector;

import java.util.ArrayList;

import junit.framework.TestCase;

/**
 * Test suite for Healpix
 * 
 * @author ejoliet
 */

public class HealpixTest extends TestCase {
	public static final double RADIAN_ARCSECOND = 1/Constants.ARCSECOND_RADIAN;//2.06264806247096e+5;
	/**
	 * Test modulo.
	 */
	public void testMODULO() {
		BitManipulation bm = new BitManipulation();
		double A = 8.;
		double B = 5.;
		double res = bm.MODULO(A, B);
		assertEquals("modulo = " + res, 3., res, 1e-10);
		A = -8.;
		B = 5.;
		res = bm.MODULO(A, B);
		assertEquals("modulo = " + res, 2., res, 1e-10);
		A = 8.;
		B = -5.;
		res = bm.MODULO(A, B);
		assertEquals("modulo = " + res, -2., res, 1e-10);
		A = -8.;
		B = -5.;
		res = bm.MODULO(A, B);
		assertEquals("modulo = " + res, -3., res, 1e-10);
		System.out.println(" test MODULO is done");
	}

	/**
	 * Test ang dist.
	 * 
	 * @throws Exception the exception
	 */
	public void testAngDist() throws Exception {
		SpatialVector v1 = new SpatialVector(1.0, 0., 0.);
		SpatialVector v2 = new SpatialVector(0., 1., 0.);
		HealpixIndex pt = new HealpixIndex();
		double res1 = pt.AngDist(v1, v2);
		double res2 = Math.PI / 2;
		System.out.println("res1 = " + res1 + " res2=" + res2);
		assertEquals("angular Distance=" + res2, 1.0, res1 / res2, 1e-10);

		// SpatialVector v5 = new SpatialVectorImp(1.5, 1.6, 0.);
		// SpatialVector v6 = new SpatialVectorImp(-1.5, -1.75, 0.);
		// Vector3d v55 = new Vector3d(1.5, 1.6, 0.);
		// Vector3d v66 = new Vector3d(-1.5, -1.75, 0.);
		// double res5 = pt.AngDist(v5, v6);
		// double res6 = v55.angle(v66);
		// System.out.println("res5 = " + res5 + " res6=" + res6);
		// assertEquals("angular Distance=" + res6, 2.0, res5 / res6, 1e-10);
		//		
		// /* Check known problem with vecmath for small vector differences */
		//
		// SpatialVector v3 = new SpatialVector(1.5, 1.6, 0.);
		// SpatialVector v4 = new SpatialVector(1.5, 1.601, 0.);
		// double res3 = pt.AngDist(v3, v4);
		// Vector3d v33 = new Vector3d(1.5, 1.6, 0.);
		// Vector3d v44 = new Vector3d(1.5, 1.601, 0.);
		// double res4 = v33.angle(v44);
		// double expected = res4 - Math.PI / 2.;
		// System.out.println("res3 = " + res3 + " res4=" + res4);
		// assertEquals("angular Distance=" + res4, 1., res3
		// / expected, 1e-3);
		System.out.println(" test of AngDist is done");
	}

	/**
	 * Test surface triangle.
	 * 
	 * @throws Exception the exception
	 */
	public void testSurfaceTriangle() throws Exception {
		SpatialVector v1 = new SpatialVector(1.0, 0.0, 0.0);
		SpatialVector v2 = new SpatialVector(0.0, 1.0, 0.0);
		SpatialVector v3 = new SpatialVector(0.0, 0.0, 1.0);
		HealpixIndex pt = new HealpixIndex();
		double res = pt.SurfaceTriangle(v1, v2, v3);
		System.out.println("Triangle surface is=" + res / Math.PI
				+ " steredians");
		assertEquals("Triangle surface=" + res, 0.5, res / Math.PI, 1e-10);
		System.out.println(" test of SurfaceTriangle is done");
	}

	/**
	 * Test nside2 npix.
	 * 
	 * @throws Exception the exception
	 */
	public void testNside2Npix() throws Exception {
		int nside = 1;
		int npix = 0;
		HealpixIndex pt = new HealpixIndex(nside);
		npix = (int) pt.Nside2Npix(nside);
		assertEquals("Npix=" + npix, 12, npix, 1e-10);
		nside = 2;
		pt = new HealpixIndex(nside);
		npix = (int) pt.Nside2Npix(nside);
		assertEquals("Npix=" + npix, 48, npix, 1e-10);
	}

	/**
	 * Test npix2 nside.
	 * 
	 * @throws Exception the exception
	 */
	public void testNpix2Nside() throws Exception {
		int npix = 48;
		int nside = 2;
		HealpixIndex pt = new HealpixIndex(nside);
		int nside1 = (int) pt.Npix2Nside(npix);
		assertEquals("Nside=" + nside1, nside, nside1, 1e-10);

	}

	/**
	 * Test vec2 ang.
	 */
	public void testVec2Ang() {
		double PI = Math.PI;
		SpatialVector v = new SpatialVector(0.0, 1.0, 0.0);
		HealpixIndex pt = new HealpixIndex();
		double[] ang_tup = { 0., 0. };
		ang_tup = pt.Vect2Ang(v);
		System.out.println(" Theta=" + ang_tup[0] / PI + " Phi=" + ang_tup[1]
				/ PI);
		assertEquals("Theta=" + ang_tup[0], 0.5, ang_tup[0] / PI, 1e-10);
		assertEquals("Phi=" + ang_tup[1], 0.5, ang_tup[1] / PI, 1e-10);
		v = new SpatialVector(1.0, 0.0, 0.0);
		ang_tup = pt.Vect2Ang(v);
		assertEquals("phi=" + ang_tup[1], 0., ang_tup[1] / PI, 1e-10);
		System.out.println(" test Vect2Ang is done");
	}

	/**
	 * Test ang2 pix.
	 * 
	 * @throws Exception the exception
	 */
	public void testAng2Pix() throws Exception {
		System.out.println(" Test ang2pix_ring ___________________");
		double PI = Math.PI;
		long pix = -1;
		double theta = PI / 2. - 0.2;
		double phi = PI / 2.;
		long nside = 4;
		HealpixIndex pt = new HealpixIndex((int) nside);
		try {
			pix = pt.ang2pix_ring(theta, phi);
		} catch ( Exception e ) {
			e.printStackTrace();
		}
		SpatialVector v = (SpatialVector) pt.Ang2Vec(theta, phi);
		long pix1 = pt.vec2pix_ring(v);
		assertEquals("pix=" + pix, pix1, pix, 1e-10);
		assertEquals("pix=" + pix, 76, pix, 1e-10);

		double[] radec = pt.pix2ang_ring(76);
		assertEquals("theta=" + theta, theta, radec[0], 4e-2);
		assertEquals("phi=" + phi, radec[1], phi, 1e-2);
		double[] radec1 = pt.pix2ang_nest(92);
		System.out.println("theta=" + radec1[0] + " phi=" + radec1[1]);
		assertEquals("theta=" + theta, theta, radec1[0], 4e-2);
		assertEquals("phi=" + phi, radec1[1], phi, 1e-2);
		System.out.println(" test Ang2Pix is done");
	}

	/**
	 * Test ang2 vect.
	 * 
	 * @throws Exception the exception
	 */
	public void testAng2Vect() throws Exception {
		System.out.println(" Start test Ang2Vect----------------");
		double PI = Math.PI;
		double theta = PI / 2.;
		double phi = PI / 2;
		HealpixIndex pt = new HealpixIndex();
		SpatialVector v = (SpatialVector) pt.Ang2Vec(theta, phi);
		System.out.println("Vector x=" + v.x() + " y=" + v.y() + " z=" + v.z());
		assertEquals("x=" + v.x(), 0., v.x(), 1e-10);
		assertEquals("y=" + v.y(), 1., v.y(), 1e-10);
		assertEquals("z=" + v.z(), 0., v.z(), 1e-10);
		System.out.println(" test Ang2Vect is done");
	}

	/**
	 * Test ring num.
	 * 
	 * @throws Exception the exception
	 */
	public void testRingNum() throws Exception {
		double z = 0.25;
		int nside = 1;
		System.out.println("Start test RingNum !!!!!!!!!!!!!!!!!!!!");
		HealpixIndex pt = new HealpixIndex(nside);
		int nring = (int) pt.RingNum(nside, z);
		System.out.println("z=" + z + " ring number =" + nring);
		assertEquals("z=" + z, 2, nring, 1e-10);
		z = -0.25;
		nring = (int) pt.RingNum(nside, z);
		assertEquals("z=" + z, 2, nring, 1e-10);
		z = 0.8;
		nring = (int) pt.RingNum(nside, z);
		assertEquals("z=" + z, 1, nring, 1e-10);
		z = -0.8;
		nring = (int) pt.RingNum(nside, z);
		assertEquals("z=" + z, 3, nring, 1e-10);
		System.out.println(" test RingNum is done");
		nside = 4;
		pt = new HealpixIndex(nside);
		int pix = 3;
		SpatialVector v = pt.pix2vec_ring(pix);
		z = v.z();
		nring = (int) pt.RingNum(nside, z);
		assertEquals("z=" + z, 1, nring, 1e-10);
		pix = 11;
		v = pt.pix2vec_ring(pix);
		z = v.z();
		nring = (int) pt.RingNum(nside, z);
		assertEquals("z=" + z, 2, nring, 1e-10);
		pix = 23;
		v = pt.pix2vec_ring(pix);
		z = v.z();
		nring = (int) pt.RingNum(nside, z);
		assertEquals("z=" + z, 3, nring, 1e-10);
		pix = 39;
		v = pt.pix2vec_ring(pix);
		z = v.z();
		nring = (int) pt.RingNum(nside, z);
		assertEquals("z=" + z, 4, nring, 1e-10);
		pix = 55;
		v = pt.pix2vec_ring(pix);
		z = v.z();
		nring = (int) pt.RingNum(nside, z);
		assertEquals("z=" + z, 5, nring, 1e-10);
		pix = 71;
		v = pt.pix2vec_ring(pix);
		z = v.z();
		nring = (int) pt.RingNum(nside, z);
		assertEquals("z=" + z, 6, nring, 1e-10);
		pix = 87;
		v = pt.pix2vec_ring(pix);
		z = v.z();
		nring = (int) pt.RingNum(nside, z);
		assertEquals("z=" + z, 7, nring, 1e-10);
		pix = 103;
		v = pt.pix2vec_ring(pix);
		z = v.z();
		nring = (int) pt.RingNum(nside, z);
		assertEquals("z=" + z, 8, nring, 1e-10);
		pix = 119;
		v = pt.pix2vec_ring(pix);
		z = v.z();
		nring = (int) pt.RingNum(nside, z);
		assertEquals("z=" + z, 9, nring, 1e-10);
		pix = 135;
		v = pt.pix2vec_ring(pix);
		z = v.z();
		nring = (int) pt.RingNum(nside, z);
		assertEquals("z=" + z, 10, nring, 1e-10);
		pix = 151;
		v = pt.pix2vec_ring(pix);
		z = v.z();
		nring = (int) pt.RingNum(nside, z);
		assertEquals("z=" + z, 11, nring, 1e-10);
		pix = 167;
		v = pt.pix2vec_ring(pix);
		z = v.z();
		nring = (int) pt.RingNum(nside, z);
		assertEquals("z=" + z, 12, nring, 1e-10);
		pix = 169;
		v = pt.pix2vec_ring(pix);
		z = v.z();
		nring = (int) pt.RingNum(nside, z);
		assertEquals("z=" + z, 13, nring, 1e-10);
		pix = 180;
		v = pt.pix2vec_ring(pix);
		z = v.z();
		nring = (int) pt.RingNum(nside, z);
		assertEquals("z=" + z, 14, nring, 1e-10);
		System.out.println("End test RingNum");
	}

	/**
	 * Test nest2 ring.
	 * 
	 * @throws Exception the exception
	 */
	public void testNest2Ring() throws Exception {
		int ipnest = 3;
		int nside = 2;
		HealpixIndex pt = new HealpixIndex(nside);
		int ipring = (int) pt.nest2ring(ipnest);
		assertEquals("ipring=" + ipring, 0, ipring, 1e-10);
		ipnest = 0;
		nside = 2;
		ipring = (int) pt.nest2ring(ipnest);
		assertEquals("ipring=" + ipring, 13, ipring, 1e-10);
		ipnest = 18;
		nside = 2;
		ipring = (int) pt.nest2ring(ipnest);
		assertEquals("ipring=" + ipring, 27, ipring, 1e-10);
		ipnest = 23;
		nside = 2;
		ipring = (int) pt.nest2ring(ipnest);
		assertEquals("ipring=" + ipring, 14, ipring, 1e-10);
		ipnest = 5;
		nside = 4;
		pt = new HealpixIndex(nside);
		ipring = (int) pt.nest2ring(ipnest);
		assertEquals("ipring = " + ipring, 27, ipring, 1e-10);
		System.out.println(" test Nest2Ring is done");
	}

	/**
	 * Test ring2 nest.
	 * 
	 * @throws Exception the exception
	 */
	public void testRing2Nest() throws Exception {
		System.out.println(" start test Ring2Nest !!!!!!!!!!!!!!!!!!!!!!");
		int ipring = 0;
		int nside = 2;
		HealpixIndex pt = new HealpixIndex(nside);
		int ipnest = (int) pt.ring2nest(ipring);
		assertEquals("ipnest=" + ipnest, 3, ipnest, 1e-10);
		ipring = 13;
		nside = 2;
		ipnest = (int) pt.ring2nest(ipring);
		assertEquals("ipnest=" + ipnest, 0, ipnest, 1e-10);
		ipring = 27;
		nside = 2;
		ipnest = (int) pt.ring2nest(ipring);
		assertEquals("ipnest=" + ipnest, 18, ipnest, 1e-10);
		ipring = 14;
		nside = 2;
		ipnest = (int) pt.ring2nest(ipring);
		assertEquals("ipnest=" + ipnest, 23, ipnest, 1e-10);
		ipring = 27;
		nside = 4;
		pt = new HealpixIndex(4);
		ipnest = (int) pt.ring2nest(ipring);
		assertEquals("ipnest = " + ipnest, 5, ipnest, 1e-10);
		ipring = 83;
		ipnest = (int) pt.ring2nest(ipring);
		assertEquals("ipnest = " + ipnest, 123, ipnest, 1e-10);
		System.out.println(" test Ring2Nest is done");
		pt = new HealpixIndex(64);
		ipnest = (int) pt.ring2nest(0);
		System.out.println(" ipnest (64,0)=" + ipnest);
	}
	/**
	 * Test next_ in_ line_ nest.
	 * 
	 * @throws Exception the exception
	 */
	public void testNext_In_Line_Nest() throws Exception {
		int ipix = 0;
		int nside = 2;
		HealpixIndex pt = new HealpixIndex(nside);
		int ipnext = (int) pt.next_in_line_nest(nside, ipix);
		assertEquals("ipnext=" + ipnext, 23, ipnext, 1e-10);
		ipix = 1;
		nside = 2;
		ipnext = (int) pt.next_in_line_nest(nside, ipix);
		assertEquals("ipnext=" + ipnext, 6, ipnext, 1e-10);
		ipix = 4;
		nside = 2;
		ipnext = (int) pt.next_in_line_nest(nside, ipix);
		assertEquals("ipnext=" + ipnext, 27, ipnext, 1e-10);
		ipix = 27;
		nside = 2;
		ipnext = (int) pt.next_in_line_nest(nside, ipix);
		assertEquals("ipnext=" + ipnext, 8, ipnext, 1e-10);
		ipix = 12;
		nside = 2;
		ipnext = (int) pt.next_in_line_nest(nside, ipix);
		assertEquals("ipnext=" + ipnext, 19, ipnext, 1e-10);
		ipix = 118;
		nside = 4;
		pt = new HealpixIndex(nside);
		ipnext = (int) pt.next_in_line_nest(nside, ipix);
		assertEquals("ipnext = " + ipnext, 117, ipnext, 1e-10);
		System.out.println(" test next_in_line_nest is done");
	}
	
	/**
	 * Test in ring cxx.
	 * 
	 * @throws Exception the exception
	 */
	public void testInRingCxx() throws Exception {
		System.out.println(" Start test InRing !!!!!!!!!!!!!!!!!!!!!!!!!");
		int[] nestComp = { 19, 0, 23, 4, 27, 8, 31, 12 };
		double PI = Math.PI;
		long nside = 2;
		HealpixIndex pt = new HealpixIndex((int) nside);
		boolean nest = false;
		int iz = 3;
		double phi = PI;
		double dphi = PI;
		ArrayList ring = pt.InRingCxx(nside, iz, phi, dphi, nest);
		for ( int i = 0; i < ring.size(); i++ ) {
			assertEquals("ipnext = " + ( (Long) ring.get(i) ).longValue(),
					i + 12, ( (Long) ring.get(i) ).longValue(), 1e-10);
		}
		nest = true;
		ring = pt.InRingCxx(nside, iz, phi, dphi, nest);
		for ( int i = 0; i < ring.size(); i++ ) {
			assertEquals("ipnext = " + ( (Long) ring.get(i) ).longValue(),
					nestComp[i], ( (Long) ring.get(i) ).longValue(), 1e-10);
		}
		nest = false;
		nside = 4;
		pt = new HealpixIndex((int) nside);
		phi = 2.1598449493429825;
		iz = 8;
		dphi = 0.5890486225480867;
		// System.out.println(" iz="+iz+" phi="+phi+" dphi="+dphi);
		ring = pt.InRingCxx(nside, iz, phi, dphi, nest);
		// for (int i = 0; i<ring.size(); i++) {
		// System.out.println("ipnext = "+ ((Integer)ring.get(i)).intValue());
		// }
		nest = false;
		dphi = 0. * PI;
		iz = 8;
		phi = 2.1598449493429825;
		// System.out.println(" iz="+iz+" phi="+phi+" dphi="+dphi);
		ring = pt.InRingCxx(nside, iz, phi, dphi, nest);
		// for (int i = 0; i<ring.size(); i++) {
		// System.out.println("ipnext = "+ ((Integer)ring.get(i)).intValue());
		// }
		System.out.println(" test InRing is done");
	}
	
	/**
	 * Test in ring.
	 * 
	 * @throws Exception the exception
	 */
	public void testInRing() throws Exception {
		System.out.println(" Start test InRing !!!!!!!!!!!!!!!!!!!!!!!!!");
		int[] nestComp = { 19, 0, 23, 4, 27, 8, 31, 12 };
		double PI = Math.PI;
		long nside = 2;
		HealpixIndex pt = new HealpixIndex((int) nside);
		boolean nest = false;
		int iz = 3;
		double phi = PI;
		double dphi = PI;
		ArrayList ring = pt.InRing(nside, iz, phi, dphi, nest);
		for ( int i = 0; i < ring.size(); i++ ) {
			assertEquals("ipnext = " + ( (Long) ring.get(i) ).longValue(),
					i + 12, ( (Long) ring.get(i) ).longValue(), 1e-10);
		}
		nest = true;
		ring = pt.InRing(nside, iz, phi, dphi, nest);
		for ( int i = 0; i < ring.size(); i++ ) {
			assertEquals("ipnext = " + ( (Long) ring.get(i) ).longValue(),
					nestComp[i], ( (Long) ring.get(i) ).longValue(), 1e-10);
		}
		nest = false;
		nside = 4;
		pt = new HealpixIndex((int) nside);
		phi = 2.1598449493429825;
		iz = 8;
		dphi = 0.5890486225480867;
		// System.out.println(" iz="+iz+" phi="+phi+" dphi="+dphi);
		ring = pt.InRing(nside, iz, phi, dphi, nest);
		// for (int i = 0; i<ring.size(); i++) {
		// System.out.println("ipnext = "+ ((Integer)ring.get(i)).intValue());
		// }
		nest = false;
		dphi = 0. * PI;
		iz = 8;
		phi = 2.1598449493429825;
		// System.out.println(" iz="+iz+" phi="+phi+" dphi="+dphi);
		ring = pt.InRing(nside, iz, phi, dphi, nest);
		// for (int i = 0; i<ring.size(); i++) {
		// System.out.println("ipnext = "+ ((Integer)ring.get(i)).intValue());
		// }
		System.out.println(" test InRing is done");
	}

	/**
	 * Test intrs_ intrv.
	 * 
	 * @throws Exception the exception
	 */
	public void testIntrs_Intrv() throws Exception {
		System.out.println(" test intrs_intrv !!!!!!!!!!!!!!!!!!!!!!!!!!!");
		HealpixIndex pt = new HealpixIndex();
		double[] d1 = { 1.0, 9.0 };
		double[] d2 = { 3.0, 16.0 };
		double[] di;
		// System.out.println("Case "+d1[0]+" "+d1[1]+" | "+d2[0]+" "+d2[1]);
		di = pt.intrs_intrv(d1, d2);
		// System.out.println("Result "+di[0]+" - "+di[1]);
		int n12 = di.length / 2;
		assertEquals("n12 = " + n12, 1, n12, 1e-6);
		assertEquals("di[0] = " + di[0], 3.0, di[0], 1e-6);
		assertEquals("di[1] = " + di[1], 9.0, di[1], 1e-6);
		d1 = new double[] { 0.537, 4.356 };
		d2 = new double[] { 3.356, 0.8 };
		// System.out.println("Case "+d1[0]+" "+d1[1]+" | "+d2[0]+" "+d2[1]);
		di = pt.intrs_intrv(d1, d2);
		n12 = di.length / 2;
		assertEquals("n12 = " + n12, 2, n12, 1e-6);
		assertEquals("di[0] = " + di[0], 0.537, di[0], 1e-6);
		assertEquals("di[1] = " + di[1], 0.8, di[1], 1e-6);
		assertEquals("di[2] = " + di[2], 3.356, di[2], 1e-6);
		assertEquals("di[1] = " + di[3], 4.356, di[3], 1e-6);

		d1 = new double[] { 2.356194490092345, 2.356194490292345 };
		d2 = new double[] { 1.251567, 4.17 };
		// System.out.println("Case "+d1[0]+" "+d1[1]+" | "+d2[0]+" "+d2[1]);
		di = pt.intrs_intrv(d1, d2);
		n12 = di.length / 2;
		assertEquals("n12 = " + n12, 1, n12, 1e-6);
		assertEquals("di[0] = " + di[0], 2.35619449009, di[0], 1e-6);
		assertEquals("di[1] = " + di[1], 2.35619449029, di[1], 1e-6);

		System.out.println(" test intrs_intrv is done");
	}

	/**
	 * Test pix2 vect_ring.
	 * 
	 * @throws Exception the exception
	 */
	public void testPix2Vect_ring() throws Exception {
		System.out.println("Start test pix2vec_ring !!!!!!!!!!!!!!!!!!!");
		double TWOPI = 2.0 * Math.PI;
		int nside = 2;
		int ipix = 0;
		HealpixIndex pt = new HealpixIndex(nside);
		SpatialVector v1 = new SpatialVector(0., 0., 0.);
		v1 = pt.pix2vec_ring(ipix);
		assertEquals("v1.z = " + v1.z(), 1.0, v1.z(), 1e-1);

		ipix = 20;
		SpatialVector v2 = new SpatialVector(0., 0., 0.);
		v2 = pt.pix2vec_ring(ipix);
		assertEquals("v2.x = " + v2.x(), 1.0, v2.x(), 1e-1);
		assertEquals("v2.z = " + v2.z(), 0.0, v2.z(), 1e-1);
		ipix = 22;
		SpatialVector v3 = new SpatialVector();
		v3 = pt.pix2vec_ring(ipix);
		assertEquals("v3.y = " + v3.y(), 1.0, v3.y(), 1e-1);
		assertEquals("v3.z = " + v3.z(), 0.0, v3.z(), 1e-1);
		// System.out.println("Vector3 x="+v3.x+" y="+v3.y+" z="+v3.z);
		ipix = 95;
		nside = 4;
		pt = new HealpixIndex(nside);
		v1 = pt.pix2vec_ring(ipix);
		v1.normalize();
		double phi1 = Math.atan2(v1.y(), v1.x());
		double[] tetphi = new double[2];
		tetphi = pt.pix2ang_ring(ipix);
		assertEquals("phi = " + phi1, 0.0, Math.abs(phi1 - tetphi[1]), 1e-10);
		ipix = 26;
		v1 = pt.pix2vec_ring(ipix);
		v1.normalize();
		phi1 = Math.atan2(v1.y(), v1.x());
		if ( phi1 < 0. )
			phi1 += TWOPI;
		tetphi = new double[2];
		tetphi = pt.pix2ang_ring(ipix);
		assertEquals("phi = " + phi1, 0.0, Math.abs(phi1 - tetphi[1]), 1e-10);
		System.out.println("------------------------------------------");
		System.out.println(" test pix2vec_ring is done");
	}

	/**
	 * Test pix2 vect_nest.
	 * 
	 * @throws Exception the exception
	 */
	public void testPix2Vect_nest() throws Exception {
		double TWOPI = 2.0 * Math.PI;
		int nside = 2;
		int ipix = 3;
		System.out.println(" Start test pix2vec_nest !!!!!!!!!!!!!!");
		HealpixIndex pt = new HealpixIndex(nside);
		SpatialVector v1 = new SpatialVector(0., 0., 0.);
		v1 = pt.pix2vec_nest(ipix);
		assertEquals("v1.z = " + v1.z(), 1.0, v1.z(), 1e-1);
		ipix = 17;
		SpatialVector v2 = new SpatialVector(0., 0., 0.);
		v2 = pt.pix2vec_nest(ipix);
		assertEquals("v2.x = " + v2.x(), 1.0, v2.x(), 1e-1);
		assertEquals("v2.z = " + v2.z(), 0.0, v2.z(), 1e-1);

		ipix = 21;
		SpatialVector v3 = new SpatialVector();
		v3 = pt.pix2vec_nest(ipix);
		assertEquals("v3.y = " + v3.y(), 1.0, v3.y(), 1e-1);
		assertEquals("v3.z = " + v3.z(), 0.0, v3.z(), 1e-1);
		nside = 4;
		pt = new HealpixIndex(nside);
		ipix = 105;
		v1 = pt.pix2vec_nest(ipix);
		v1.normalize();
		double phi1 = Math.atan2(v1.y(), v1.x());
		if ( phi1 < 0. )
			phi1 += TWOPI;
		double[] tetphi = new double[2];
		tetphi = pt.pix2ang_nest(ipix);
		assertEquals("phi = " + phi1, 0.0, Math.abs(phi1 - tetphi[1]), 1e-10);
		nside = 4;
		ipix = 84;
		v1 = pt.pix2vec_nest(ipix);
		v1.normalize();
		phi1 = Math.atan2(v1.y(), v1.x());
		if ( phi1 < 0. )
			phi1 += TWOPI;
		// tetphi = new double[2];
		tetphi = pt.pix2ang_nest(ipix);
		assertEquals("phi = " + phi1, 0.0, Math.abs(phi1 - tetphi[1]), 1e-10);

		System.out.println(" test pix2vec_nest is done");
		System.out.println("-------------------------------------------");
	}

	/**
	 * Test vect2 pix_ring.
	 * 
	 * @throws Exception the exception
	 */
	public void testVect2Pix_ring() throws Exception {
		System.out.println("Start test Vect2Pix_ring !!!!!!!!!!!!!!!!!!!");
		int nside = 4;
		int ipix = 84;
		int respix = 0;
		HealpixIndex pt = new HealpixIndex(nside);
		SpatialVector v1 = pt.pix2vec_ring(ipix);
		respix = pt.vec2pix_ring(v1);
		assertEquals("respix = " + respix, ipix, respix, 1e-10);

		System.out.println("------------------------------------------");
		System.out.println(" test vec2pix_ring is done");
	}

	/**
	 * Test vect2 pix_nest.
	 * 
	 * @throws Exception the exception
	 */
	public void testVect2Pix_nest() throws Exception {
		System.out.println("Start test Vect2Pix_nest !!!!!!!!!!!!!!!!!!!");
		int nside = 4;
		int ipix = 88;
		int respix = 0;
		HealpixIndex pt = new HealpixIndex(nside);
		// SpatialVector v1 = new SpatialVectorImp(0., 0., 0.);
		SpatialVector v1 = pt.pix2vec_nest(ipix);
		respix = pt.vec2pix_nest(v1);
		assertEquals("respix = " + respix, ipix, respix, 1e-10);

		System.out.println("------------------------------------------");
		System.out.println(" test vec2pix_nest is done");
	}

	/**
	 * Test query_ strip.
	 * 
	 * @throws Exception the exception
	 */
	public void testQuery_Strip() throws Exception {
		System.out.println(" Start test query Strip !!!!!!!!!!!!!!!!");
		int[] pixel1 = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
				16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31,
				32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47,
				48, 49, 50, 51, 52, 53, 54, 55 };
		int nside = 4;
		HealpixIndex pt = new HealpixIndex(nside);
		int nest = 0;
		double theta1 = 0.0;
		double theta2 = Math.PI / 4.0 + 0.2;
		ArrayList pixlist;
		pixlist = pt.query_strip(nside, theta1, theta2, nest);
		int nlist = pixlist.size();
		for ( int i = 0; i < nlist; i++ ) {
			long ipix = ( (Long) pixlist.get(i) ).longValue();
			assertEquals("pixel = " + ipix, pixel1[i], ipix, 1e-10);
		}

		System.out.println(" test query_strip is done");

	}

    /**
     * Converts radians to arcsecs.
     * 
     * @param rad
     *            radians
     * @return arcseconds
     */
    public double radToAs(final double rad) {
        return rad * RADIAN_ARCSECOND;
    }
    /**
     * Converts radians to arcmin.
     * 
     * @param rad
     *            radians
     * @return arcmin
     */
    public double radToAm(final double rad) {
        return rad * RADIAN_ARCSECOND/60.;
    }
	  /**
  	 * Test query_ disc.
  	 * 
  	 * @throws Exception the exception
  	 */
  	public void testQuery_Disc() throws Exception {
		System.out.println(" Start test query_disc !!!!!!!!!!!!!!!!!!!!!");
		int nside = 4;
		HealpixIndex pt = new HealpixIndex(nside);
		int nest = 0;
		long ipix = 0;
		int[] pixel1 = { 45, 46, 60, 61, 62, 76,77, 78, 79, 91, 92, 93, 94, 95, 108,
				109, 110, 111, 124, 125, 126, 141, 142 };//ring
		int[] pixel2 = { 24, 19, 93, 18, 17, 92, 87, 16, 107, 89, 86, 85, 106,
				105, 83, 84, 159, 104, 81, 158, 157, 155, 156 };//nest
		int inclusive = 1;
		double radius = Math.PI / 8.0;
		SpatialVector v = pt.pix2vec_ring(93);
		ArrayList pixlist;
		pixlist = pt.query_disc(nside, v, radius, nest, inclusive);
		int nlist = pixlist.size();
		for ( int i = 0; i < nlist; i++ ) {
			ipix = ( (Long) pixlist.get(i) ).longValue();
//			System.out.println(" i=" + i + "pixel=" + ipix);
			assertEquals("pixel = " + ipix, pixel1[i], ipix, 1e-10);
		}
//		v = pt.pix2vec_ring(103);
//		pixlist = pt.query_disc(nside, v, radius, nest, inclusive);
//		nlist = pixlist.size();
//		for ( int i = 0; i < nlist; i++ ) {
//			ipix = ( (Long) pixlist.get(i) ).longValue();
			// assertEquals("pixel = " +ipix, pixel1[i], ipix, 1e-10);
//			System.out.println(" i=" + i + "pixel=" + ipix);
//		}
		 v = pt.pix2vec_nest( (int) pt.ring2nest(93)); 
		nest = 1; 
		inclusive = 1;
		pixlist = pt.query_disc(nside, v, radius, nest, inclusive);
		nlist = pixlist.size();
		for ( int i = 0; i < nlist; i++ ) {
			ipix = ( (Long) pixlist.get(i) ).longValue();
//			 System.out.println("ipix="+ ipix);
			assertEquals("pixel = " + ipix, pixel2[i], ipix, 1e-10);
		}
		System.out.println(" test query_disc is done");
		
		pt = new HealpixIndex(4096);
		v = pt.pix2vec_nest( 175 );
		//1 arcmin
		radius = 1*Constants.ARCSECOND_RADIAN;
		pixlist = pt.query_disc(4096, v, radius, nest, inclusive);
		nlist = pixlist.size();
		for ( int i = 0; i < nlist; i++ ) {
			ipix = ( (Long) pixlist.get(i) ).longValue();
			System.out.println(" i=" + i + " pixel=" + ipix);
		}
	}
	 
	/**
	 * Test query_ triangle360 nest.
	 * 
	 * @throws Exception the exception
	 */
	public void testQuery_Triangle360Nest() throws Exception {
		int nside = 2;
		HealpixIndex pt = new HealpixIndex(nside);
		int nest = 1;
		long ipix = 0;
		int inclusive = 0;
		
		int triang[] = { 12, 0, 16 }; // crossing 360
		int pixels[] = { 12,19,0,18,17,16 };
		System.out.println("Start test Query Triangle");
		SpatialVector v[] = new SpatialVector[3];

		for ( int vi = 0; vi < v.length; vi++ ) {

			v[vi] = pt.pix2vec_nest(triang[vi]);
//			System.out.println("ipnest = "+triang[vi]+" ipring = "+pt.nest2ring(triang[vi]));
//			printVec(v[vi].get());
		}

		ArrayList pixlist;
		pixlist = pt.query_triangle(nside, v[0], v[1], v[2], nest, inclusive);

		int nlist = pixlist.size();

		for ( int i = 0; i < nlist; i++ ) {
			ipix =  ((Long) pixlist.get(i) ).longValue();
//			ipix = pt.ring2nest((int) ( (Long) pixlist.get(i) ).longValue());
//			System.out.println(ipix);
			 assertEquals("pixel = " + ipix + " wrong.", (long) pixels[i], ipix);

		}
		
		nside=4;
		pt = new HealpixIndex(nside);
		int triang4[] = { 78, 66, 68 }; // crossing 360
		int pixels4[] = { 78,75,76,73,70,72,67,68,66, };
		v = new SpatialVector[3];

		for ( int vi = 0; vi < v.length; vi++ ) {

			v[vi] = pt.pix2vec_ring(pt.nest2ring(triang4[vi]));
//			System.out.println("ipnest = "+triang4[vi]+" ipring = "+pt.nest2ring(triang4[vi]));
//			printVec(v[vi].get());
		}

		pixlist = pt.query_triangle(nside, v[0], v[1], v[2], nest, inclusive);

		nlist = pixlist.size();

		for ( int i = 0; i < nlist; i++ ) {
			ipix = ( (Long) pixlist.get(i) ).longValue();
//			System.out.println(ipix);
			 assertEquals("pixel = " + ipix + " wrong.", (long) pixels4[i], ipix);

		}
		System.out.println(" test query_triangle is done");
		
	}

	/**
	 * Test query_ triangle360 ring.
	 * 
	 * @throws Exception the exception
	 */
	public void testQuery_Triangle360Ring() throws Exception {
		int nside = 2;
		HealpixIndex pt = new HealpixIndex(nside);
		int nest = 0;
		long ipix = 0;
		int inclusive = 0;
		int triang[] = { 19, 13, 28 }; // crossing 360
		int pixels[] = { 19, 12, 13, 27, 20, 28, };
		System.out.println("Start test Query Triangle");
		SpatialVector v[] = new SpatialVector[3];

		for ( int vi = 0; vi < v.length; vi++ ) {

			v[vi] = pt.pix2vec_ring(triang[vi]);
//			System.out.println("iptested = "+triang[vi]);
//			printVec(v[vi].get());
		}
		ArrayList pixlist;
		pixlist = pt.query_triangle(nside, v[0], v[1], v[2], nest, inclusive);

		int nlist = pixlist.size();

		for ( int i = 0; i < nlist; i++ ) {
			ipix = ( (Long) pixlist.get(i) ).longValue();
//			System.out.println(ipix);
			 assertEquals("pixel = " + ipix+ " wrong.", (long)pixels[i],
			 ipix);

		}
		nside = 4;
		inclusive = 1;
		pt = new HealpixIndex(nside);
		int triang4[] = { 71, 135, 105 }; // crossing 360
		int pixels4[] = { 71, 87, 72, 103, 88, 119, 104, 105, 135, };//120 is missing because center outside the triangle (inclusive=0)
		System.out.println("Start test Query Triangle");
		v = new SpatialVector[3];

		for ( int vi = 0; vi < v.length; vi++ ) {

			v[vi] = pt.pix2vec_ring(triang4[vi]);
//			System.out.println("iptested = "+triang4[vi]);
		}

		pixlist = pt.query_triangle(nside, v[0], v[1], v[2], nest, inclusive);

		nlist = pixlist.size();
		if ( inclusive == 1 )
			pixels4 = new int[] { 55, 71, 87, 72, 103, 88, 119, 104, 105, 135,
					120, 151 };

		for ( int i = 0; i < nlist; i++ ) {
			ipix = ( (Long) pixlist.get(i) ).longValue();
//			System.out.println(ipix);
//			assertEquals("pixel = " + ipix + " wrong.", (long) pixels4[i], ipix);

		}
		System.out.println(" test query_triangle is done");
	}
	
	/**
	 * Test query_ triangle not crossing360 ring.
	 * 
	 * @throws Exception the exception
	 */
	public void testQuery_TriangleNotCrossing360Ring() throws Exception {
		int nside = 4;
		int triang[] = { 44, 46, 77 }; //not crossing 360
		int pixels[] = { 44,45,46,60,61,77, };//120 is missing because center outside the triangle (inclusive=0)
		
//		nside = 2;
//		int triang[] = { 14, 16, 31 }; // not crossing 360
//		int pixels[] = { 14,15,16,22,23,31, };
		
		HealpixIndex pt = new HealpixIndex(nside);
		int nest = 0;
		long ipix = 0;
		int inclusive = 0;
		
		System.out.println("Start test Query Triangle");
		SpatialVector v[] = new SpatialVector[3];

		for ( int vi = 0; vi < v.length; vi++ ) {

			v[vi] = pt.pix2vec_ring(triang[vi]);

		}

		ArrayList pixlist;
		pixlist = pt.query_triangle(nside, v[0], v[1], v[2], nest, inclusive);

		int nlist = pixlist.size();

		for ( int i = 0; i < nlist; i++ ) {
			ipix = ( (Long) pixlist.get(i) ).longValue();
//			System.out.println(ipix);
			assertEquals("pixel = " + ipix + " wrong.", (long) pixels[i], ipix);

		}
		System.out.println(" test query_triangle is done");
	}
	
	/**
	 * Test query_ triangle not crossing360 nest.
	 * 
	 * @throws Exception the exception
	 */
	public void testQuery_TriangleNotCrossing360Nest() throws Exception {
		int nside = 4;
		int nest = 1;
		int inclusive = 0;
		
		int triang[] = { 95, 20, 85 }; // not crossing 360
		int pixels[] = { 95,24,19,20,93,18,17,87,16,85, };//120 is missing because center outside the triangle (inclusive=0)
//		nside = 2;
//		nest = 1;
//		inclusive = 0;
//		 int triang[] = {23,27,39} ; // not crossing 360
//		 int pixels[] = { 23,4,27,21,26,39 };
		HealpixIndex pt = new HealpixIndex(nside);
		long ipix = 0;
		System.out.println("Start test Query Triangle");
		SpatialVector v[] = new SpatialVector[3];

		for ( int vi = 0; vi < v.length; vi++ ) {

			v[vi] = pt.pix2vec_nest(triang[vi]);

		}

		ArrayList pixlist;
		pixlist = pt.query_triangle(nside, v[0], v[1], v[2], nest, inclusive);

		int nlist = pixlist.size();

		for ( int i = 0; i < nlist; i++ ) {
			ipix = ( (Long) pixlist.get(i) ).longValue();
//			System.out.println(ipix);
			assertEquals("pixel = " + ipix + " wrong.", (long) pixels[i], (int) ipix);

		}
		System.out.println(" test query_triangle is done");
	}
	
	/**
	 * Test query_ polygon.
	 * 
	 * @throws Exception the exception
	 */
	@SuppressWarnings("unchecked")
	public void testQuery_Polygon() throws Exception {
		long nside = 4;
		HealpixIndex pt = new HealpixIndex((int) nside);
		int nest = 0;
		long ipix = 0;
		int inclusive = 0;
		int[] result = { 51, 52, 53, 66, 67, 68, 69, 82, 83, 84, 85, 86, 98,
				99, 100, 101, 115, 116, 117 };
		int[] result1 = { 55, 70, 71, 87 };
		int[] result2 = { 137, 152, 153, 168 };
		int[] result3 = { 27, 43, 44, 58, 59, 60, 74, 75, 76, 77, 89, 90, 91,
				92, 93, 105, 106, 107, 108, 109, 110, 121, 122, 123, 124, 125,
				138, 139, 140, 141, 154, 156 };

		System.out.println("Start test query_polygon !!!!!!!!!!!!!!!!!!!!!!");
		ArrayList vlist = new ArrayList();
		SpatialVector v = pt.pix2vec_ring(53);
		vlist.add((Object) v);
		v = pt.pix2vec_ring(51);
		vlist.add((Object) v);
		v = pt.pix2vec_ring(82);
		vlist.add((Object) v);
		v = pt.pix2vec_ring(115);
		vlist.add((Object) v);
		v = pt.pix2vec_ring(117);
		vlist.add((Object) v);
		v = pt.pix2vec_ring(86);
		vlist.add((Object) v);

		ArrayList pixlist;
		pixlist = pt.query_polygon(nside, vlist, nest, inclusive);
		// System.out.println(" List size="+pixlist.size());
		int nlist = pixlist.size();
		// System.out.println(" Pixel list:");
		for ( int i = 0; i < nlist; i++ ) {
			ipix = ( (Long) pixlist.get(i) ).longValue();
			assertEquals("pixel = " + ipix, result[i], ipix, 1e-10);
			// System.out.println("i="+i+" pixel # "+ipix);
		}

		/* Yet another test */

		ArrayList vlist1 = new ArrayList();
		v = pt.pix2vec_ring(71);
		vlist1.add((Object) v);
		v = pt.pix2vec_ring(55);
		vlist1.add((Object) v);
		v = pt.pix2vec_ring(70);
		vlist1.add((Object) v);
		v = pt.pix2vec_ring(87);
		vlist1.add((Object) v);
		pixlist = pt.query_polygon(nside, vlist1, nest, inclusive);

		nlist = pixlist.size();

		for ( int i = 0; i < nlist; i++ ) {
			ipix = ( (Long) pixlist.get(i) ).longValue();
			// System.out.println("i="+i+" pixel # "+ipix);
			assertEquals("pixel = " + ipix, result1[i], ipix, 1e-10);
		}

		/* Yet another test */
		ArrayList vlist2 = new ArrayList();
		v = pt.pix2vec_ring(153);
		vlist2.add((Object) v);
		v = pt.pix2vec_ring(137);
		vlist2.add((Object) v);
		v = pt.pix2vec_ring(152);
		vlist2.add((Object) v);
		v = pt.pix2vec_ring(168);
		vlist2.add((Object) v);
		pixlist = pt.query_polygon(nside, vlist2, nest, inclusive);

		nlist = pixlist.size();

		for ( int i = 0; i < nlist; i++ ) {
			ipix = ( (Long) pixlist.get(i) ).longValue();
			assertEquals("pixel = " + ipix, result2[i], ipix, 1e-10);
			// System.out.println("i="+i+" pixel # "+ipix);
		}
		/* Yet another test */

		ArrayList vlist3 = new ArrayList();
		v = pt.pix2vec_ring(110);
		vlist3.add((Object) v);
		v = pt.pix2vec_ring(27);
		vlist3.add((Object) v);
		v = pt.pix2vec_ring(105);
		vlist3.add((Object) v);
		v = pt.pix2vec_ring(154);
		vlist3.add((Object) v);
		v = pt.pix2vec_ring(123);
		vlist3.add((Object) v);
		v = pt.pix2vec_ring(156);
		vlist3.add((Object) v);
		pixlist = pt.query_polygon(nside, vlist3, nest, inclusive);

		nlist = pixlist.size();

		for ( int i = 0; i < nlist; i++ ) {
			ipix = ( (Long) pixlist.get(i) ).longValue();
			assertEquals("pixel = " + ipix, result3[i], ipix, 1e-10);
			// System.out.println("i="+i+" pixel # "+ipix);
		}
		System.out.println(" test query_polygon is done");

	}

	/**
	 * Test max resolution.
	 * 
	 * @throws Exception the exception
	 */
	public void testMaxResolution() throws Exception {
		System.out.println(" Start test MaxRes !!!!!!!!!!!!!!!!!!!!!");

		// long nside = 1048576;
		long nside = 8192;
		HealpixIndex pt = new HealpixIndex((int) nside);
		double res = pt.getPixRes(nside);
		System.out.println("Minimum size of the pixel side is " + res
				+ " arcsec.");
		assertEquals("res = " + res, 25.76, res, 1e-1);
		System.out.println(" End of MaxRes test _______________________");
	}

	/*
	 * public void testQueryDiskRes() throws Exception { System.out.println("
	 * Start test DiskRes !!!!!!!!!!!!!!!!!!!!!"); HealpixIndex pt = new
	 * HealpixIndex(); int nest = 0; long ipix = 0; int inclusive = 0; double
	 * theta= Math.PI; double phi = Math.PI; double radius =
	 * Math.toRadians(0.2/3600.); // One arcse long nside = pt.getNSide(radius);
	 * System.out.println(" calculated nside="+nside); long cpix =
	 * pt.ang2pix_ring(theta,phi); SpatialVector vc = pt.pix2vec_ring((int)
	 * cpix); ArrayList pixlist; pixlist = pt.query_disc(nside, vc, radius,
	 * nest, inclusive); int nlist = pixlist.size(); for (int i = 0; i < nlist;
	 * i++) { ipix = ((Long) pixlist.get(i)).longValue(); SpatialVector v =
	 * pt.pix2vec_ring((int) ipix); double dist = pt.AngDist(v,vc);
	 * assertTrue(dist<=2.*radius); } cpix = pt.ang2pix_nest(theta,phi);
	 * SpatialVector vc1 = pt.pix2vec_nest( (int) cpix); long cpixnest=
	 * pt.vec2pix_nest(vc1); ArrayList pixlist1; nest = 1; radius *=4; pixlist1 =
	 * pt.query_disc(nside, vc1, radius, nest, inclusive); int nlist1 =
	 * pixlist1.size(); for (int i = 0; i < nlist1; i++) { ipix = ((Long)
	 * pixlist1.get(i)).longValue(); SpatialVector v = pt.pix2vec_nest((int)
	 * ipix); double dist = pt.AngDist(v,vc1); assertTrue(dist<=2.*radius); }
	 * System.out.println(" test query disc Hi Res is done
	 * -------------------"); }
	 */
	/**
	 * Test get nside.
	 */
	public void testGetNside() {
		System.out.println(" Start test GetNside !!!!!!!!!!!!!!!!!!!!!");

		double pixsize = 25;
		HealpixIndex pt = new HealpixIndex();
		long nside = HealpixIndex.calculateNSide(pixsize);
		System.out.println("Required nside is " + nside);
		assertEquals("nside = " + nside, 8192, nside, 1e-1);
		System.out.println(" End of GetNSide test _______________________");
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
	 * Test keith cover.
	 * 
	 * @throws Exception the exception
	 */
	public void testKeithCover() throws Exception {
		SpatialVector galUnitVec = new SpatialVector(-0.8359528269,
				-0.4237775489, -0.3487055694);

		double[] orig = galUnitVec.get();

		System.err.println(galUnitVec.ra() + ", " + galUnitVec.dec());
		System.err.println(galUnitVec);
		final HealpixIndex healpixIndex = new HealpixIndex(512);

		double[] pos1a = healpixIndex.Vect2Ang(galUnitVec);
		AngularPosition pos1 = new AngularPosition(pos1a[0], pos1a[1]);

		int pixelIdxSel = healpixIndex.vec2pix_nest(galUnitVec);
		int pix2 = healpixIndex.ang2pix_nest(pos1a[0], pos1a[1]);

		assertEquals("pixels dont match", pixelIdxSel, pix2);
		double[] vec = healpixIndex.pix2vec_nest(pixelIdxSel).get();

		double[] posa = healpixIndex.pix2ang_nest(pixelIdxSel);
		AngularPosition pos = new AngularPosition(posa[0], posa[1]);

		System.err.println("orig:" + pos1 + "  ret:" + pos);

		assertEquals("Angular pos theta incorrect ", pos1.theta(), pos.theta(),
				0.01);
		assertEquals("Angular pos phi incorrect ", pos1.phi(), pos.phi(), 0.01);

		SpatialVector svec = new SpatialVector(vec[0], vec[1], vec[2]);

		printVec(vec);
		for ( int i = 0; i < vec.length; i++ ) {
			assertEquals(" original vector not returned position:" + i,
					orig[i], vec[i], 0.01);
		}

		System.err.println(galUnitVec);

	}
	/**
	 * tests Neighbour's method for nest schema of the pixelization
	 * 
	 * @throws Exception
	 */
	public void testNeighbours_Nest() throws Exception {
		System.out.println(" Start test Neighbours_Nest !!!!!!!!!!!!!!!!!");		
		int nside = 2;
		HealpixIndex pt = new HealpixIndex(nside);
		long ipix = 25;
		long[] test17 = { 34, 16, 18, 19, 2, 0, 22, 35 };
		long[] test32 = { 40, 44, 45, 34, 35, 33, 38, 36 };
		long[] test3 = { 0, 2, 13, 15, 11, 7, 6, 1 };
		long[] test25 = { 42, 24, 26, 27, 10, 8, 30, 43 };
		//
		ArrayList npixList = new ArrayList();
		npixList = pt.neighbours_nest(nside, ipix);
		for ( int i = 0; i < npixList.size(); i++ ) {
			assertEquals("ip = " + ( (Long) npixList.get(i) ).longValue(),
					test25[i], ( (Long) npixList.get(i) ).longValue(), 1e-10);
		}
		ipix = 17;

		npixList = pt.neighbours_nest(nside, ipix);
		for ( int i = 0; i < npixList.size(); i++ ) {
			assertEquals("ip = " + ( (Long) npixList.get(i) ).longValue(),
					test17[i], ( (Long) npixList.get(i) ).longValue(), 1e-10);
		}
		ipix = 32;

		npixList = pt.neighbours_nest(nside, ipix);
		for ( int i = 0; i < npixList.size(); i++ ) {
			assertEquals("ip = " + ( (Long) npixList.get(i) ).longValue(),
					test32[i], ( (Long) npixList.get(i) ).longValue(), 1e-10);
		}
		ipix = 3;

		npixList = pt.neighbours_nest(nside, ipix);
		for ( int i = 0; i < npixList.size(); i++ ) {
			assertEquals("ip = " + ( (Long) npixList.get(i) ).longValue(),
					test3[i], ( (Long) npixList.get(i) ).longValue(), 1e-10);
		}
		pt = new HealpixIndex(4096);
		npixList = pt.neighbours_nest(4096, 175);
		for ( int i = 0; i < npixList.size(); i++ ) {
			System.err.println("ip = " + npixList.get(i));
		}
		System.out.println(" test NeighboursNest is done");
	}
	/** run through angles ang vectors converting back and forth * */
	/**
	 * public void testAllAngVec() throws Exception { HealpixIndex heal = new
	 * HealpixIndex(512); for (double theta=0; theta < Constants.PI; theta+=0.1) {
	 * for (double phi=0; phi < 2*Constants.PI; phi+=0.1) { AngularPosition pos =
	 * new AngularPosition(theta,phi); // SpatialVector vec2 =
	 * heal.vector(theta, phi); SpatialVector vec =
	 * (SpatialVector)heal.Ang2Vec(theta, phi); double[] ang =
	 * heal.Vect2Ang(vec); assertEquals("Ang theta is incorrect",theta,ang[0] );
	 * assertEquals("Ang phi is incorrect",phi,ang[1] ); } } }
	 */
	
	public void testGetParentAt(int childnside, int requirednside) throws Exception {
		HealpixIndex c = new HealpixIndex(childnside);
		HealpixIndex p = new HealpixIndex(requirednside);
		
		for (double theta =0.1;  theta < 3.16; theta+=0.2){
			for (double phi=0.1; phi < 6.28 ; phi+=0.2){
				long child = c.ang2pix_nest(theta, phi);
				long ppix = p.ang2pix_nest(theta, phi);	
				long tppix = HealpixIndex.parentAt(child, childnside, requirednside);
				
				assertEquals("Parent is not as expexted",tppix,ppix) ;
	   }}
	}
	
	public void testGetParentAt() throws Exception{
		testGetParentAt(4,2);
		testGetParentAt(8,2);
		testGetParentAt(1024,128);
		
	}
	
	
	public void testGetChildrenAt() throws Exception {
		long[] pix = HealpixIndex.getChildrenAt(2, 1, 4);
		assertEquals("Incorrect number of pixels ",4, pix.length);
		int fpix = 4 ; // see the primer
		for (int i=0; i < pix.length; i++){
			assertEquals("Incorrect pixel number", fpix+i, pix[i]);
		}
	
	}
}
