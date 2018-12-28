package healpix.fits;

import healpix.core.dm.HealpixMap;

/**
 * Interface to read fits file into healpix map
 * 
 * @author ejoliet
 * @version $Id: Fits2HealpixMap.java,v 1.1.2.2 2009/08/03 16:25:20 healpix Exp $
 */
public interface Fits2HealpixMap {

	/**
	 * Reads a binary table from fits file
	 * 
	 * @param filename
	 *            the fits file name
	 * @return double[][] Array of doubles containing the values
	 * @throws Exception
	 */
	@SuppressWarnings("deprecation")
	public abstract double[][] readFitsBinaryTable(String filename)
			throws Exception;

	/**
	 * Read fits binary table from FITS file and set the {@link HealpixMap} with
	 * the data from table read
	 * 
	 * @param filename
	 *            FITS file name.
	 * @return {@link HealpixMap} map.
	 * @throws Exception
	 */
	public abstract HealpixMap fits2map(String filename) throws Exception;

	/**
	 * Getting the map.
	 * 
	 * @return {@link HealpixMap} map
	 */
	public abstract HealpixMap getMap();

	/**
	 * Set the {@link HealpixMap}
	 * 
	 * @param m
	 *            {@link HealpixMap}
	 */
	public abstract void setMap(HealpixMap m);

	/**
	 * Set the maps name
	 * 
	 * @param nam
	 *            names
	 */
	public abstract void setColname(String[] nam);

	/**
	 * Get maps name
	 * 
	 * @return names
	 */
	public abstract String[] getColname();

}