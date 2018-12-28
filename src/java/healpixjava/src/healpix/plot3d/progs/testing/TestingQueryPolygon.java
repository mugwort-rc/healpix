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
package healpix.plot3d.progs.testing;

import healpix.core.AngularPosition;
import healpix.core.HealpixIndex;
import healpix.core.dm.HealpixMap;
import healpix.core.dm.AbstractHealpixMap.Scheme;
import healpix.plot3d.gui.view.MapView3d;
import healpix.tools.HealpixMapCreator;
import healpix.tools.SpatialVector;

import java.util.ArrayList;

/**
 * Testing query polygon to detect prime meridian black zone defect in Planck
 * Usage: java -cp jhealpix.jar healpix.plot3d.progs.testing.TestingQueryPolygon
 * 
 * @author ejoliet
 * @version $Id: TestingQueryPolygon.java,v 1.1 2008/04/25 14:44:51 healpix Exp $
 */
public class TestingQueryPolygon {
	private static ArrayList<Object> vlist;

	public static void main(String[] args) throws Exception {
		try {
			MapView3d mview = new MapView3d(false);
			vlist = new ArrayList<Object>();
			HealpixMap map = getMap3PixelsRing();// getMapWithPixRingTriangle();
													// // getMap();
			mview.setMap(map);
			System.out.println("Map min/max: " + map.getMin(0) + "/"
					+ map.getMax(0));
			mview.setSize(800, 800);
			mview.setVisible(true);
		} catch ( Exception e ) {
			e.printStackTrace();
		}
	}

	public static HealpixMap getMap() throws Exception {
		int nside = 4;
		HealpixMapCreator cr = new HealpixMapCreator(nside, true);
		HealpixMap map = cr.getMap();
		map.setScheme(Scheme.NEST);
		HealpixIndex hi = new HealpixIndex(nside);
		SpatialVector vec1 = new SpatialVector();
		double offsetra = 0.0;
		double offsetdec = 0.0;
		// make left polygon
		vec1.set(( 345.0 + offsetra ), 10.0 + offsetdec);// north
		SpatialVector vec2 = new SpatialVector();
		System.out.println("1:" + hi.vec2pix_nest(vec1));
		vec2.set(( 350.0 + offsetra ), -20.0 + offsetdec);// south
		System.out.println("2:" + hi.vec2pix_nest(vec2));
		// hemisphere
		SpatialVector vec3 = new SpatialVector();
		vec3.set(( 20.0 + offsetra ), -20.0 + offsetdec);// South hemisphere
		System.out.println("3:" + hi.vec2pix_nest(vec3));
		SpatialVector vec4 = new SpatialVector();
		vec4.set(( 20.0 + offsetra ), 10.0 + offsetdec);
		System.out.println("4:" + hi.vec2pix_nest(vec4));

		vlist.add(vec1);
		vlist.add(vec2);
		vlist.add(vec3);
		vlist.add(vec4);

		ArrayList pixlist = new HealpixIndex(map.nside()).query_polygon(map
				.nside(), vlist, 1, 0);
		pixlist = new HealpixIndex(map.nside()).query_triangle(map.nside(),
				vec1, vec2, vec3, 1, 0);
		int nlist = pixlist.size();
		for ( int i = 0; i < nlist; i++ ) {
			long ip = (int) ( (Long) pixlist.get(i) ).longValue();
			map.setValueCell((int) ip, 0.5);
			// System.out.println(ip);
		}

		addVec(vec1, map, 1);
		addVec(vec2, map, 2);
		addVec(vec3, map, 3);
		addVec(vec4, map, 4);
		return map;
	}

	public static HealpixMap getMapWithPixRingTriangle() throws Exception {
		int nside = 4;
		HealpixMapCreator cr = new HealpixMapCreator(nside, true);

		HealpixMap map = cr.getMap();
		map.setScheme(Scheme.RING);

		ArrayList vlist1 = new ArrayList();
		HealpixIndex pt = new HealpixIndex(nside);
		int nest = 0;
		long ipix = 0;
		int inclusive = 1;
		int triang[] = { 71, 135, 105 }; // crossing 360
		int pixels[] = { 71, 72, 88, 103, 104, 135, };// 120 is missing
		// because center
		// outside the triangle
		// (inclusive=0)
		System.out.println("Start test Query Triangle");
		SpatialVector v[] = new SpatialVector[3];

		for ( int vi = 0; vi < v.length; vi++ ) {
			map.add((int) triang[vi], 10);
			v[vi] = pt.pix2vec_ring(triang[vi]);

		}

		ArrayList pixlist;
		pixlist = pt.query_triangle(nside, v[0], v[1], v[2], nest, inclusive);

		int nlist = pixlist.size();
		for ( int i = 0; i < nlist; i++ ) {
			long ip = (int) ( (Long) pixlist.get(i) ).longValue();
			map.add((int) pt.ring2nest((int) ip), 5);
			System.out.println(ip);
		}

		return map;

	}

	public static HealpixMap getMap3PixelsRing() throws Exception {
		int nside = 4;
		HealpixMapCreator cr = new HealpixMapCreator(nside, true);
		HealpixMap map = cr.getMap();
		map.setScheme(Scheme.NEST);
		// int[] ipringtest = {19,13,28};//crossing
		int[] ipringtest = { 71, 135, 105 };// crossing (nside = 4)
		ArrayList vlist1 = new ArrayList();
		HealpixIndex pt = new HealpixIndex(nside);
		System.out.println("Start test Query Triangle");
		SpatialVector v[] = new SpatialVector[3];
		for ( int i = 0; i < ipringtest.length; i++ ) {
			int iptest = (int) pt.ring2nest(ipringtest[i]);
			v[i] = pt.pix2vec_ring(ipringtest[i]);
			double[] ang = pt.pix2ang_nest(iptest);
			System.out.println(ang[0] + "," + ang[1]);
			map.add(iptest, 10);
		}

		ArrayList pixlist;
		pixlist = pt.query_triangle(nside, v[0], v[1], v[2], 0, 0);

		int nlist = pixlist.size();
		for ( int i = 0; i < nlist; i++ ) {
			long ip = (int) ( (Long) pixlist.get(i) ).longValue();
			map.add((int) pt.ring2nest((int) ip), 5);
			System.out.println(ip);
		}
		return map;
	}

	public static HealpixMap getMapWithPixRing() throws Exception {
		int nside = 4;
		HealpixMapCreator cr = new HealpixMapCreator(nside, true);

		HealpixMap map = cr.getMap();
		map.setScheme(Scheme.RING);

		ArrayList vlist1 = new ArrayList();
		SpatialVector v = null;
		HealpixIndex pt = new HealpixIndex(nside);

		double pv = 3;
		v = pt.pix2vec_ring(1);
		vlist1.add((Object) v);
		map.setValueCell((int) 1, 9);

		addVec(v, map, pv++);
		v = pt.pix2vec_ring(48);
		vlist1.add((Object) v);
		addVec(v, map, pv++);
		v = pt.pix2vec_ring(94);
		vlist1.add((Object) v);
		addVec(v, map, pv++);
		v = pt.pix2vec_ring(112);
		vlist1.add((Object) v);
		addVec(v, map, pv++);
		v = pt.pix2vec_ring(81);
		vlist1.add((Object) v);
		addVec(v, map, pv++);

		ArrayList pixlist = pt.query_polygon(nside, vlist1, 0, 0);
		int nlist = pixlist.size();
		for ( int i = 0; i < nlist; i++ ) {
			long ip = (int) ( (Long) pixlist.get(i) ).longValue();
			map.add((int) ip, 0.5);
			System.out.println(ip);
		}

		return map;

	}

	public static HealpixMap getMapWithPixNest() throws Exception {
		int nside = 4;
		HealpixMapCreator cr = new HealpixMapCreator(nside, true);
		HealpixMap map = cr.getMap();
		map.setScheme(Scheme.NEST);

		ArrayList vlist1 = new ArrayList();
		SpatialVector v = null;
		HealpixIndex pt = new HealpixIndex(nside);

		double pv = 3;
		v = pt.pix2vec_nest(21);
		vlist1.add((Object) v);
		map.setValueCell((int) 21, pv);
		addVec(v, map, pv++);
		v = pt.pix2vec_nest(16);
		vlist1.add((Object) v);
		map.setValueCell(16, pv);
		addVec(v, map, pv++);
		v = pt.pix2vec_nest(104);
		vlist1.add((Object) v);
		map.setValueCell(104, pv);
		addVec(v, map, pv++);
		v = pt.pix2vec_nest(109);
		map.setValueCell(109, pv);
		vlist1.add((Object) v);
		addVec(v, map, pv++);

		addVec(v, map, pv++);
		ArrayList pixlist = pt.query_polygon(nside, vlist1, 1, 0);
		int nlist = pixlist.size();
		for ( int i = 0; i < nlist; i++ ) {
			long ip = (int) ( (Long) pixlist.get(i) ).longValue();
			map.add((int) ip, 0.5);
			System.out.println(ip);
		}

		return map;

	}

	/**
	 * @param vec1
	 * @param vec2
	 * @return
	 */
	private static double[] diffVec(SpatialVector vec1, SpatialVector vec2) {
		double x, y, z;
		x = -vec1.getX() + vec2.getX();
		y = -vec1.getY() + vec2.getY();
		z = -vec1.getZ() + vec2.getZ();
		return new double[] { x, y, z };
	}

	public static void addVec(SpatialVector vec, HealpixMap map, double v)
			throws Exception {

		HealpixIndex hi = new HealpixIndex(4);

		// double angs[] = HealpixIndex.ang(vec);
		double angs[] = hi.Vect2Ang(vec);

		AngularPosition ang2 = new AngularPosition(angs[0], angs[1]);

		map.add(ang2, v);

	}
}