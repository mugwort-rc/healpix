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
package healpix.fits;

import healpix.core.AngularPosition;
import healpix.core.dm.HealpixMap;
import healpix.core.dm.HealpixMapImp;
import healpix.core.dm.AbstractHealpixMap.Scheme;

import java.util.Iterator;

import net.ivoa.fits.Fits;
import net.ivoa.fits.Header;
import net.ivoa.fits.HeaderCard;
import net.ivoa.fits.data.AsciiTable;
import net.ivoa.fits.data.BinaryTable;
import net.ivoa.fits.data.Data;
import net.ivoa.fits.hdu.AsciiTableHDU;
import net.ivoa.fits.hdu.BasicHDU;
import net.ivoa.fits.hdu.BinaryTableHDU;
import net.ivoa.fits.hdu.TableHDU;
import net.ivoa.util.ColumnTable;

/**
 * Converts fits file into healpix map. Read the fits and create the
 * {@link HealpixMap} object.
 * 
 * @version $Id: Fits2HealpixMapImp.java,v 1.1 2008/04/25 14:44:51 healpix Exp $
 * @author ejoliet
 */
public class Fits2HealpixMapImp implements Fits2HealpixMap {
	private int nrow;

	private int ncol;

	private static HealpixMap map;

	private String schemeFits;

	public String[] names;

	private ColumnTable objData;

	private static double[][] data;

	/**
	 * Default constructor.
	 */
	public Fits2HealpixMapImp() {
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see healpix.fits.Fits2HealpixMap#readFitsBinaryTable(java.lang.String)
	 */
	@SuppressWarnings("deprecation")
	public double[][] readFitsBinaryTable(String filename) throws Exception {
		Fits fits = new Fits(filename);
		BasicHDU[] bhdu = fits.read();
		TableHDU thdu = null;
		Header head = null;
		int headerInd = 0;
		Data dat = null;
		for ( int i = 0; i < fits.getNumberOfHDUs(); i++ ) {
			dat = bhdu[i].getData();
			System.out.println("Data:" + dat.toString());
			if ( dat instanceof BinaryTable ) {
				headerInd = i;
				thdu = (BinaryTableHDU) fits.getHDU(headerInd);
				BinaryTable data = (BinaryTable) dat;
				ncol = data.getNCols();
				nrow = data.getNRows();
				objData = (ColumnTable) data.getData();
				System.out.println("Ncols, Nrows:" + ncol + ", " + nrow);
				// char[] types = data.getTypes();
				break;
			}
			if ( dat instanceof AsciiTable ) {
				headerInd = i;
				thdu = (AsciiTableHDU) fits.getHDU(headerInd);
				AsciiTable data = (AsciiTable) dat;
				ncol = data.getNCols();
				nrow = data.getNRows();
				objData = (ColumnTable) data.getData();
				System.out.println("Ncols, Nrows:" + ncol + ", " + nrow);
				break;
			}
		}
		System.out.println(objData.getColumn(0).toString());
		head = thdu.getHeader();
		int ncard = head.getNumberOfCards();
		Iterator hci = head.iterator();
		while ( hci.hasNext() ) {
			HeaderCard hc = (HeaderCard) hci.next();
			// System.out.println("Key card read from fits:" + head.getKey(n));
			// System.out.println("-> " + head.getValue(head.getKey(n)));
			if (hc.getKey().equals("NSIDE") ) {
				System.out.println("NSIDE read from fits:"
						+ hc.getValue());
			}
			if ( hc.getKey().equals("ORDERING") ) {
				schemeFits = hc.getValue();
				System.out.println("Scheme found:" + schemeFits);
			}
			if (hc.getKey().contains("UNIT")) {
				System.out.println("UNIT read from fits:"
						+ hc.getValue());
			}
		}
		double[][] array = new double[ncol][nrow];
		names = new String[ncol];
		for ( int i = 0; i < ncol; i++ ) {
			String type = thdu.getColumnFormat(i);
			int dim = getDimension(type);
			if ( getDataType(type) != 'A' && getDataType(type) != 'B' ) {
				Object obj = thdu.getColumn(i);
				if ( getDataType(type) == 'F' ) {
					if ( dim > 1 ) {
						// for (int u = 0; u < dim; u++) {
						// float[][] dati = (float[][]) obj;
						// float[] tdoub = dati[u];
						// }
					} else {
						double[] tdoub = (double[]) obj;
						for ( int j = 0; j < nrow; j++ ) {
							array[i][j] = tdoub[j];
						}
					}
				} else {
					if ( dim > 1 ) {
						// for (int u = 0; u < dim; u++) {
						// double[][] dati = (double[][]) obj;
						// double[] tdoub = dati[u];
						// }
					} else {
						double[] tdoub = (double[]) obj;
						for ( int j = 0; j < nrow; j++ ) {
							array[i][j] = tdoub[j];
						}
					}
				}
				names[i] = thdu.getColumnName(i);
				System.out.println("Name:" + names[i]);
			} else {
				System.err
						.println("Format Byte or String couldn't be read, not implemented!");
				return null;
				// float[][] tdoub = (float[][]) thdu.getColumn(i);
				// int nr = tdoub.length;
				// int nc = tdoub[0].length;
				// names[i] = thdu.getColumnName(i);
				// System.out.println("Name:" + names[i]);
				// for (int m = 0; m < nc; m++) {
				// for (int j = 0; j < nr; j++) {
				// array[m][j] = (double) tdoub[m][j];
				// }
				// }
			}
		}
		setColname(names);
		data = array;
		printTable();
		return array;
	}

	public int getDimension(String ftype) {

		if ( ftype.length() == 1 ) // No multiplicity present, value is a
			// scalar
			return 0;
		if ( ftype.indexOf("J") >= 0 )
			return Short.parseShort(ftype.substring(0, ftype.indexOf("J")));
		if ( ftype.indexOf("K") >= 0 )
			return Integer.parseInt(ftype.substring(0, ftype.indexOf("K")));
		if ( ftype.indexOf("E") >= 0 )
			return Integer.parseInt(ftype.substring(0, ftype.indexOf("E")));
		if ( ftype.indexOf("D") >= 0 )
			return Integer.parseInt(ftype.substring(0, ftype.indexOf("D")));
		if ( ftype.indexOf("I") >= 0 )
			return Integer.parseInt(ftype.substring(0, ftype.indexOf("I")));
		if ( ftype.indexOf("B") >= 0 )
			return Integer.parseInt(ftype.substring(0, ftype.indexOf("B")));
		if ( ftype.indexOf("A") >= 0 )
			return Integer.parseInt(ftype.substring(0, ftype.indexOf("A")));
		return 0;
	}

	public char getDataType(String ftype) {

		if ( ftype.indexOf("J") >= 0 )
			return 'I';
		if ( ftype.indexOf("K") >= 0 )
			return 'L';
		if ( ftype.indexOf("E") >= 0 )
			return 'F';
		if ( ftype.indexOf("D") >= 0 )
			return 'D';
		if ( ftype.indexOf("I") >= 0 )
			return 'S';
		if ( ftype.indexOf("B") >= 0 )
			return 'B';
		if ( ftype.indexOf("A") >= 0 )
			return 'A';
		return ' ';
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see healpix.fits.Fits2HealpixMap#fits2map(java.lang.String)
	 */
	public HealpixMap fits2map(String filename) throws Exception {
		double[][] data = readFitsBinaryTable(filename);
		HealpixMapImp map = new HealpixMapImp((long) nrow, getColname());
		map.setScheme(Scheme.NEST);
		for ( int i = 0; i < ncol; i++ ) {
			for ( int j = 0; j < nrow; j++ ) {
				double val = data[i][j];
				if ( schemeFits.equals("RING") ) {
					map.setValueCell(i, (int) map.ring2nest(j), val);
				} else {
					map.setValueCell(i, j, val);
				}
			}
		}
		setMap(map);
		return getMap();
	}

	/**
	 * Print out the data table.
	 * 
	 * @throws Exception
	 */
	public void printTable() throws Exception {
		for ( int n = 0; n < ncol; n++ ) {
			int nr = 0;
			// the 10 first data...
			while ( nr < 10 ) {
				if ( getColname()[n].equals("RA") ) {
					System.out.println("*** Row/Index = " + nr);
					AngularPosition ang = map.pix2ang(nr);
					System.out.println(getColname()[n] + " = " + data[n][nr]);// RA
					System.out.println(getColname()[n + 1] + " = "
							+ data[n + 1][nr]);
					System.out.println("HealpixMap Theta/Phi ? = "
							+ ang.theta() + "/" + ang.phi());
					System.out.println("Pix number for DEC/RA above:"
							+ map.ang2pix(data[n + 1][nr], data[n][nr]));
				}
				nr++;
			}
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see healpix.fits.Fits2HealpixMap#getMap()
	 */
	public HealpixMap getMap() {
		return map;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see healpix.fits.Fits2HealpixMap#setMap(healpix.core.dm.HealpixMap)
	 */
	public void setMap(HealpixMap m) {
		map = m;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see healpix.fits.Fits2HealpixMap#setColname(java.lang.String[])
	 */
	public void setColname(String[] nam) {
		names = nam;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see healpix.fits.Fits2HealpixMap#getColname()
	 */
	public String[] getColname() {
		return names;
	}
}