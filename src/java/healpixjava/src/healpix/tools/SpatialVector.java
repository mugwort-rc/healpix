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
package healpix.tools;

import javax.vecmath.Vector3d;

/**
 * The SpatialVector is a standard 3D vector class with the addition that each
 * coordinate (x,y,z) is also kept in ra,dec since we expect the vector to live
 * on the surface of the unit sphere, i.e.
 * 
 * <pre>
 *  2   2   2
 *  x + y + z  = 1
 * </pre>
 * 
 * This is not enforced, so you can specify a vector that has not unit length.
 * If you request the ra/dec of such a vector, it will be automatically
 * normalized to length 1 and you get the ra/dec of that vector (the
 * intersection of the vector's direction with the unit sphere.
 * 
 * This code comes originally from the HTM library of Peter Kunst during his
 * time at JHU.
 */

public class SpatialVector extends Vector3d {

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;


	/** The ra_. */
	double ra_;

	/** The dec_. */
	double dec_;

	/** The ok ra dec_. */
	boolean okRaDec_;

	/**
	 * Default constructor constructs (1,0,0), ra=0, dec=0.
	 */
	public SpatialVector() {
		super(1, 0, 0);
		ra_ = 0;
		dec_ = 0;
		okRaDec_ = true;
	}

	/**
	 * Constructor from three coordinates
	 * 
	 * @param x
	 * @param y
	 * @param z
	 */
	public SpatialVector(double x, double y, double z) {
		super(x, y, z);
		ra_ = 0;
		dec_ = 0;
		okRaDec_ = false;
	}

	/**
	 * Construct from ra/dec in degrees
	 * 
	 * @param ra
	 *            RA in degrees
	 * @param dec
	 *            DEC in degrees
	 */
	public SpatialVector(double ra, double dec) {
		ra_ = ra;
		dec_ = dec;
		okRaDec_ = true;
		updateXYZ();
	}

	/**
	 * Copy constructor - be aware this only copies x,y,z
	 * 
	 * @param copy
	 *            the vector to copy
	 */
	public SpatialVector(SpatialVector copy) {
		super(copy.x(), copy.y(), copy.z());
		normalize();
		updateRaDec();
	}

	/**
	 * Sets the ra and dec angles in degrees
	 * 
	 * @param ra
	 *            right ascension angle in degrees
	 * @param dec
	 *            declination angle in degrees
	 * 
	 */
	public void set(double ra, double dec) {
		ra_ = ra;
		dec_ = dec;
		okRaDec_ = true;
		updateXYZ();
	}

	/**
	 * Get the coordinates in a 3 elements 1D array
	 * 
	 * @return coordinates [x,y,z]
	 */
	public double[] get() {
		double ret[] = new double[3];
		ret[0] = x;
		ret[1] = y;
		ret[2] = x;
		return ret;
	}
	
	
	/** return x (only as rvalue) */
	public double x() {
		return x;
	}
	

	/** return y (only as rvalue) */
	public double y() {
		return y;
	}


	/** return z (only as rvalue) */
	public double z() {
		return z;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see javax.vecmath.Tuple3d#toString()
	 */
	public String toString() {
		return "" + x() + " " + y() + " " + z();
	}

	/**
	 * vector cross product
	 * 
	 * @param v
	 *            the vector to cross
	 * @return the vector cross product
	 */
	public SpatialVector cross(SpatialVector v) {

		return new SpatialVector(y() * v.z() - v.y() * z(), z() * v.x() - v.z()
				* x(), x() * v.y() - v.x() * y());
	}

	/**
	 * Compare vectors if coordinates are equals
	 * 
	 * @param v
	 *            the vector to be compared with
	 * @return true if both coordinates of vectors are equal
	 */
	public boolean equal(SpatialVector v) {
		return ((x() == v.x() && y() == v.y() && z() == v.z()) ? true : false);
	}

	/**
	 * multiply with a number
	 * 
	 * @param n
	 *            the scale number to be multiply to the coordinates x,y,z
	 * @return the vector with coordinates multiplied by n
	 */
	public SpatialVector mul(double n) {
		return new SpatialVector((n * x()), (n * y()), (n * z()));
	}

	/**
	 * vector addition
	 * 
	 * @param v
	 *            the vector to be added
	 * @return vector result by addition
	 */
	public SpatialVector add(SpatialVector v) {
		return new SpatialVector(x() + v.x(), y() + v.y(), z() + v.z());
	}

	/**
	 * vector subtraction
	 * 
	 * @param v
	 *            the vector to be substracted
	 * @return vector result by substraction
	 */
	public SpatialVector sub(SpatialVector v) {
		return new SpatialVector(x() - v.x(), y() - v.y(), z() - v.z());
	}

	/**
	 * Get the dec angle in degrees
	 * 
	 * @return declination angle
	 */
	public double dec() {
		if (!okRaDec_) {
			normalize();
			updateRaDec();
		}
		return dec_;
	}

	/**
	 * Get the ra angle in degrees
	 * 
	 * @return right ascension
	 */
	public double ra() {
		if (!okRaDec_) {
			normalize();
			updateRaDec();
		}
		return ra_;
	}

	/**
	 * Update x_ y_ z_ from ra_ and dec_ variables
	 */
	protected void updateXYZ() {
		double cd = Math.cos(dec_ * Constants.cPr);
		x = Math.cos(ra_ * Constants.cPr) * cd;
		y = Math.sin(ra_ * Constants.cPr) * cd;
		z = Math.sin(dec_ * Constants.cPr);
		set(x, y, z);
	}

	/**
	 * Update ra_ and dec_ from x_ y_ z_ variables
	 */
	protected void updateRaDec() {
		dec_ = Math.asin(z()) / Constants.cPr; // easy.
		double cd = Math.cos(dec_ * Constants.cPr);
		if (cd > Constants.EPS || cd < -Constants.EPS)
			if (y() > Constants.EPS || y() < -Constants.EPS) {
				if (y() < 0.0)
					ra_ = 360 - Math.acos(x() / cd) / Constants.cPr;
				else
					ra_ = Math.acos(x() / cd) / Constants.cPr;
			} else {
				ra_ = (x() < 0.0 ? 180 : 0.0);
				// ra_ = (x_ < 0.0 ? -1.0 : 1.0) *
				// Math.asin(y_/cd)/Constants.cPr;
			}
		else
			ra_ = 0.0;
		okRaDec_ = true;
	}

};
