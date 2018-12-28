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
package healpix.tools.test;

import healpix.core.dm.HealpixMap;
import healpix.fits.Fits2HealpixMap;
import healpix.fits.Fits2HealpixMapImp;
import healpix.plot3d.gui.view.MapView3d;
import healpix.tools.HealpixMapCreator;
import junit.framework.TestCase;

public class HealpixMapCreatorTest extends TestCase {

	HealpixMap map;

	private void viewMap() {
		MapView3d mview = new MapView3d();
		mview.setMap(getMap());
		mview.setVisible(true);
	}

	public void testGenerateMapFromFits() {
		Fits2HealpixMap fmap = new Fits2HealpixMapImp();
		try {
			fmap.fits2map("data/test/test_2.ds");
			setMap(fmap.getMap());
			viewMap();
			assertEquals(" Checks the nSide : ", fmap.getMap().nside(), 2);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public void testGenerateMapCreator() {
		HealpixMapCreator cr = new HealpixMapCreator(2);
		assertEquals(" Checks the nSide : ", cr.getMap().nside(), 2);
	}
	/**
	 * Create a copy of a map from his data column(s) and string name(s).
	 * Checks max val, nside and value test from a particular pixel index
	 */
	public void testGenMapCreator() {
		HealpixMapCreator cr = new HealpixMapCreator(64);
		double val2test = 3.9d;
		int ipix2test = 5;
		cr.getMap().setValueCell(0,ipix2test, val2test);
		double[][] mapData = cr.getMap().getMapItemData();
		HealpixMapCreator cr2 = new HealpixMapCreator(cr.getMap().getName(),mapData);
		assertEquals(" Checks the nSide : ", cr2.getMap().nside(), 64);
		assertEquals(" Checks max : ", cr2.getMap().getMax(0), (double)val2test);
		assertEquals(" Checks ipix val : ", cr2.getMap().get(0,ipix2test), val2test);
		
	}
	public HealpixMap getMap() {
		return map;
	}

	protected void setMap(HealpixMap m) {
		map = m;
	}
}
