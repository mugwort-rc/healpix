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
package healpix.plot3d.gui;

import healpix.plot3d.canvas3d.SineColorTransform;

import java.awt.Canvas;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Image;
import java.awt.image.BufferedImage;

/**
 * This class is basicaly a plotting canvas for a box of pixels extracted from
 * an HealpixMap. It implements ExtBoxDisplayer.
 * 
 * @author ejoliet
 * @version $Id: BoxDisplay.java,v 1.1 2008/04/25 14:44:51 healpix Exp $
 */
public class BoxDisplay extends Canvas implements ExtBoxDisplayer {
	private static final long serialVersionUID = 1L;

	protected double[][] pix;

	protected double scalex = 0;

	protected double scaley = 0;

	protected BufferedImage theImage;

	protected Image theScaledImage;

	protected boolean newbox;

	public void setPixBox(double[][] pix) {
		System.out.println("setting the box..");
		this.pix = pix;
		newbox = true;
		theImage = new BufferedImage(pix[0].length, pix.length,
				BufferedImage.TYPE_INT_RGB);
		// finding max and min pixel value for color scale
		double max = Double.MIN_VALUE;
		double min = Double.MAX_VALUE;
		for (int y = 0; y < pix.length; y++) {
			for (int x = 0; x < pix[0].length; x++) {
				if (pix[y][x] < min)
					min = pix[y][x];
				if (pix[y][x] > max)
					max = pix[y][x];
			}
		}
		System.out.println("Box min " + min + " and box max " + max);
		SineColorTransform ctrans = new SineColorTransform((double) min,
				(double) max);
		for (int y = 0; y < pix.length; y++) {
			for (int x = 0; x < pix[0].length; x++) {
				if (Double.isNaN(pix[y][x])) {
					// blank... we'll set them to white for the moment
					Color c = new Color(128, 128, 128);
					theImage.setRGB(x, y, c.getRGB());
				} else {
					Color c = ctrans.getColor(pix[y][x]);// new Color(red,
					// green, blue);
					theImage.setRGB(x, y, c.getRGB());
				}
			}
		}

		rescale();
		paint(this.getGraphics());
	}

	public BoxDisplay() {
		setBackground(Color.black);
	}

	public void rescale() {
		setBackground(Color.black);
		boolean newImage = false;
		double temp = scalex;
		this.scalex = this.getWidth() / (pix[0].length * 2);
		newImage = temp != scalex;
		temp = scaley;
		this.scaley = this.getHeight() / (pix.length * 2);
		newImage = newImage || (temp != scaley);
		if (theImage != null && (newImage || newbox)) {

			newbox = false;
			theScaledImage = theImage.getScaledInstance(this.getWidth(), this
					.getHeight(), 0);
			System.out.println("Scaled image " + theScaledImage + " "
					+ this.getWidth());

		}
	}

	public void paint(Graphics g) {
		this.setBackground(Color.black);
		if (pix != null) {
			System.out.println("Plotting image ..");
			rescale();
			g.drawImage(theScaledImage, 0, 0, this);
		}

	}

	// Properties
	public double[][] getPixBox() {
		return pix;
	}

	public Dimension getPreferedSize() {
		return new Dimension(100, 100);
	}
}
