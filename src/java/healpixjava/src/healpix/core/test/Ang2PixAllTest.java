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
import healpix.core.Healpix;
import healpix.core.HealpixIndex;
import healpix.tools.Constants;

import java.text.DecimalFormat;

import junit.framework.TestCase;

/**
 * Test the healpix pixel and angle related methods.
 * 
 * @author ejoliet
 * @version $Id: Ang2PixAllTest.java,v 1.1 2008/04/25 14:44:51 healpix Exp $
 */
public class Ang2PixAllTest extends TestCase {

	/**
	 * Test in ring scheme with nside = 128
	 * 
	 * @throws Exception
	 */
	public void testAllRing128() throws Exception {
		testAll(128, true);
	}

	/**
	 * Test in ring scheme with nside = 64
	 * 
	 * @throws Exception
	 */
	public void testAllRing64() throws Exception {
		testAll(64, true);
	}

	/**
	 * Test in nest scheme with nside = 128
	 * 
	 * @throws Exception
	 */
	public void testAllNest128() throws Exception {
		testAll(128, false);
	}

	/**
	 * Test pix to angle and angle to pixel tied to a resolution number nside
	 * 
	 * @param nside
	 *            resolution number
	 * @param ring
	 *            if true, ring scheme is selected
	 * @throws Exception
	 */
	public void testAll(int nside, boolean ring) throws Exception {

		HealpixIndex hi = new HealpixIndex(nside);

		int length = 12 * nside * nside;
		DecimalFormat form = new DecimalFormat("#.###");

		for (int i = 0; i < length; i++) {
			AngularPosition pos = null;
			double[] posHi = null;

			if (ring) {
				pos = Healpix.pix2ang_ring(nside, i);
				posHi = hi.pix2ang_ring(i);
			} else {
				pos = Healpix.pix2ang_nest(nside, i);
				posHi = hi.pix2ang_nest(i);
			}
			int pix = 0;
			if (ring)
				pix = Healpix.ang2pix_ring(nside, pos.theta(), pos.phi());
			else
				pix = Healpix.ang2pix_nest(nside, pos.theta(), pos.phi());
			assertEquals("i incorrect for  theta "
					+ form.format(Math.cos(pos.theta())) + " phi/pi "
					+ form.format(pos.phi() / Constants.PI), i, pix);

			assertEquals("Healpix and HealpixIndex disagree on theta for " + i,
					pos.theta(), posHi[0]);
			assertEquals("Healpix and HealpixIndex disagree on phi", pos.phi(),
					posHi[1]);

		}

	}

	/** Uwe's problem index* */
	public void test8063() throws Exception {

		int nside = 64;
		HealpixIndex hi = new HealpixIndex(nside);
		int pix = 8062;
		AngularPosition pos = Healpix.pix2ang_ring(nside, pix);
		double[] posHi = hi.pix2ang_ring(pix);
		System.out.println(pix + " -> " + pos + " HealpixIndex -> " + posHi[0]
				+ " ," + posHi[1]);
		assertEquals("Healpix and HealpixIndex disagree on theta for " + pix,
				pos.theta(), posHi[0]);
		assertEquals("Healpix and HealpixIndex disagree on phi for " + pix, pos
				.phi(), posHi[1]);
		pix++;

		pos = Healpix.pix2ang_ring(nside, pix);
		posHi = hi.pix2ang_ring(pix);
		System.out.println(pix + " -> " + pos + " HealpixIndex -> " + posHi[0]
				+ " ," + posHi[1]);
		assertEquals("Healpix and HealpixIndex disagree on theta for " + pix,
				pos.theta(), posHi[0]);
		assertEquals("Healpix and HealpixIndex disagree on phi for " + pix, pos
				.phi(), posHi[1]);
		pix++;

		pos = Healpix.pix2ang_ring(nside, pix);
		posHi = hi.pix2ang_ring(pix);
		System.out.println(pix + " -> " + pos + " HealpixIndex -> " + posHi[0]
				+ " ," + posHi[1]);
		assertEquals("Healpix and HealpixIndex disagree on theta for " + pix,
				pos.theta(), posHi[0]);
		assertEquals("Healpix and HealpixIndex disagree on phi for " + pix, pos
				.phi(), posHi[1]);
		pix++;

	}

	/** Uwe's problem index* */
	public void test41087() throws Exception {

		int nside = 64;
		AngularPosition pos0 = Healpix.pix2ang_ring(nside, 41086);
		AngularPosition pos1 = Healpix.pix2ang_ring(nside, 41087);
		AngularPosition pos2 = Healpix.pix2ang_ring(nside, 41088);

		System.out.println("41086 -> " + pos0);
		System.out.println("41087 -> " + pos1);
		System.out.println("41088 -> " + pos2);

		assertEquals("theta jump too big  ", pos1.theta(), pos2.theta(), 0.1);

	}

}
