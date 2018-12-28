//Source file: /usr4/users/womullan/Planck/healpix.corej/java/healpix.tools/Tools.java

package healpix.tools;

import nom.tam.fits.*;

public class Tools {
	
	public Tools() {}
	
	/**
	 */
	public static final void yieldToGC() {
		        // get the current thread and priority
        Thread myThread = Thread.currentThread ();
        int myPriority = myThread.getPriority ();

        // set prioroty lower so GC can run
        myThread.setPriority (Thread.MIN_PRIORITY);

        // call the system garbage collector
        System.gc ();

        // yield this thread so GC runs
        myThread.yield ();

        // after GC, set priority back to original state
        myThread.setPriority (myPriority);

	}
	
	/**
	   Check a string for T,t TRUE  or true and return true if any of these are found.
	 */
	public static boolean parseBool(String s) {
		return (s.toUpperCase().startsWith("T"));
	}

	public static String extractKey(String key, Header hd) throws Exception {
		String card =hd.findKey(key);
		if (card==null) throw new Exception();
		int pos1 = card.indexOf('=') +1;
		int pos2 = card.indexOf('/') -1;
		if (pos2 < 0) pos2= card.length() -1;

		return card.substring(pos1,pos2).trim();
	}
}
