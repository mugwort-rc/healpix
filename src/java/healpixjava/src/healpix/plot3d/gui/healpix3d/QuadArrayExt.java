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

import healpix.core.AngularPosition;

import javax.media.j3d.QuadArray;

/**
 * Quadrilatere Array with tooltip fonctionality ready to be used by
 * {@link DataSphere}.
 * 
 * @author ejoliet
 * @version $Id: QuadArrayExt.java,v 1.1 2008/04/25 14:44:51 healpix Exp $
 */
public class QuadArrayExt extends QuadArray {
	String text;

	int ipix;

	AngularPosition angle;

	double value;

	int index;

	HealpixDataIndex data[];

	QuadArrayExt(int nPoints, int color) {
		super(nPoints, color);
		data = new HealpixDataIndex[nPoints];
		init();
	}

	private void init() {
		for (int u = 0; u < data.length; u++) {
			data[u] = new HealpixDataIndex();
		}
	}

	public void setText(String txt) {
		this.text = txt;
	}

	public String getText() {
		return text;
	}

	public void setIpix(int ind, int ipix) {
		data[ind].ipix = ipix;
	}

	public int getIpix(int ind) {
		return data[ind].ipix;
	}

	public void setIpix(int ipix) {
		this.ipix = ipix;
	}

	public int getIpix() {
		return ipix;
	}

	public void setAngle(int ind, AngularPosition ang) {
		data[ind].angle = ang;
	}

	public void setAngle(AngularPosition ang) {
		this.angle = ang;
	}

	public AngularPosition getAngle(int ind) {
		return data[ind].angle;
	}

	public AngularPosition getAngle() {
		return angle;
	}

	public void setValue(int ind, double d) {
		data[ind].value = d;
	}

	public void setValue(double d) {
		this.value = d;
	}

	public String getToolTipTxt(int ind) {
		return "<html>" + "Value=" + data[ind].value + "<br>"
				+ data[ind].angle.toString() + "<br>" + "Healpix pixel:"
				+ data[ind].ipix + "<html>";
	}

	public String getToolTipTxt() {
		return "<html>" + "Value=" + value + "<br>" + angle.toString() + "<br>"
				+ "Healpix pixel:" + ipix + "<html>";
	}
}
