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
package healpix.plot3d.gui.healpix3d;

import healpix.core.dm.HealpixMap;
import healpix.tools.SpatialVector;

import javax.media.j3d.Geometry;
import javax.media.j3d.GeometryArray;
import javax.vecmath.Color3f;
import javax.vecmath.Point3d;

/**
 * Despite name represents a single Healpix face. Uses a coloured quadrilateral
 * to indicate a data value for each pixel. DataSphere deals with different map
 * inside a HealpixMap object - e.g. read from fits file-.
 * 
 * @version $Id: HealpixPixelDataSphere.java,v 1.1 2008/04/25 14:44:51 healpix Exp $
 */
public class HealpixPixelDataSphere extends HealSphere {

	private static final double BASE_BLUE_FACTOR = 0.19f;

	private static final double BASE_GREEN_FACTOR = 0.29f;

	private static final double BASE_RED_FACTOR = 0.36f;

	private static final double BASE_BLUE_OFFSET = 0.029f;

	private static final double BASE_GREEN_OFFSET = -0.18f;

	private static final double BASE_RED_OFFSET = -0.46f;

	protected int face; // = 0, default initialization.

	protected int imap;

	int q;

	protected HealpixMap ch;

	protected double min, max;

	// --------------------------------------------------------------------------
	/** Default constructor. */
	public HealpixPixelDataSphere() {
		super();

	}

	/**
	 * Used to get the data sphere from a ith map
	 * 
	 * @param ch
	 *            the map data
	 * @param imap
	 *            the column data from this map ch
	 * @param iQuads
	 *            ith quad array
	 */
	public HealpixPixelDataSphere(HealpixMap ch, int imap, int iQuads) {
		super(ch.nside());
		this.imap = imap;
		this.ch = ch;
		this.q = iQuads;
		// System.out.println("Map col name"+ch.getName()[imap]);
		this.min = ch.getMin(imap);
		this.max = ch.getMax(imap);
		this.setGeometry(createGeometry());
	}

	protected Geometry createGeometry() {
		int ppq = 4; // points per quad
		// QuadArrayExt[] quads2 = new QuadArrayExt[nQuads];
		QuadArrayExt quads = new QuadArrayExt(ppq * 12
				* (int) Math.pow(nside, 2), GeometryArray.COORDINATES
				| GeometryArray.COLOR_3);

		// Specific colour scaling factors and offsets.
		double scaleColor = (double) Math.abs(max - min);
		double bc = BASE_BLUE_FACTOR * scaleColor;
		double gc = BASE_GREEN_FACTOR * scaleColor;
		double rc = BASE_RED_FACTOR * scaleColor;
		double bof = BASE_BLUE_OFFSET * scaleColor;
		double gof = BASE_GREEN_OFFSET * scaleColor;
		double rof = BASE_RED_OFFSET * scaleColor;

		int offset;
		Color3f c;

		try {
			// quads2[q] = new QuadArrayExt(4,GeometryArray.COORDINATES
			// | GeometryArray.COLOR_3);
			int pixindex = q;
			SpatialVector[] points = index.corners_nest(pixindex, 1);
			double val = (double) ch.get(imap, pixindex);// ch.getPixAsFloat(pixindex);
			// System.out.println("********** val(ipix=" +
			// pixindex+"):"+val);
			if (Double.isNaN(val)) {
				c = new Color3f(128, 128, 128);
				// System.out.println("Val isNaN, color:"+c.toString());
			} else {
				// double blue = (double) Math.sin((val - min + bof) / bc);
				// double green = (double) Math.sin((val - min + gof) / gc);
				// double red = (double) Math.sin((val - min + rof) / rc);
				// IndexedColorMap cmap = new IndexedColorMap();
				float blue = (float) Math.sin((val - min + bof) / bc);
				float green = (float) Math.sin((val - min + gof) / gc);
				float red = (float) Math.sin((val - min + rof) / rc);
				if (red < 0)
					red = 0;
				if (blue < 0)
					blue = 0;
				if (green < 0)
					green = 0;
				c = new Color3f(red, green, blue);
				// System.out.println(q+".-val="+val+"/min="+min+":color::"+c.toString());
			}
			offset = q * ppq;
			// System.out.println("Points length:"+points.length);
			for (int v = 0; v < points.length; v++) {
				Point3d p3d = new Point3d(points[v].x(), points[v].y(),
						points[v].z());
				System.out.println(q + ".- point=" + v + "offset=" + offset);
				System.out.println("point setCoord(offset+v,...)="
						+ (offset + v));
				System.out.println("Points[v]=" + points[v].toString());
				quads.setCoordinate(offset + v, p3d);
				quads.setColor(offset + v, c);
			}
			quads.setAngle(ch.pix2ang(pixindex));
			quads.setIpix(pixindex);
			quads.setValue(ch.get(imap, pixindex));
			// }
			// System.out.println("********** End ipix=" + (faceoff + q - 1));
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(0);
		}

		return quads;
	}
}