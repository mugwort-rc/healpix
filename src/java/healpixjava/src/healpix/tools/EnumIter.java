package healpix.tools;


import java.util.Iterator;
import java.util.Enumeration;
import java.lang.reflect.Method;

/**
	Masqueard an Enumeration as an Iteration
 */
public class EnumIter implements Iterator {
	protected Enumeration dataIter;
	public EnumIter(Enumeration e ) {
		this.dataIter=e;
    }
	/**
		Take next data object and Return
	 */
	public Object next() {
		return dataIter.nextElement();
    }

	/**
	 */
	public boolean hasNext() {
		return	dataIter.hasMoreElements();
    }
	
	/** optional - needed to satisfy interface Iterator */
	public void remove() {};

}
