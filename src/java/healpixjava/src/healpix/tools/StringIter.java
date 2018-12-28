//Source file: /usr4/users/womullan/Planck/healpix.corej/java/healpix.tools/StringIter.java

package healpix.tools;


import java.util.Iterator;
import java.util.Enumeration;
import java.lang.reflect.Method;

/**
   Wraps an iterator up and manufactures a string from each object returned. The String id compsed of the return values of the Get Methods provided in the COnstructor.
 */
public class StringIter implements Iterator, Enumeration {
	protected Iterator dataIter;
	protected Method getms[];
	protected String gets[];
	protected String delim = ",";
	
	/**
	   Need the list of gets and an iter Need classname to get the getMs
	 */
	public StringIter(Iterator dataIter, String[] gets) throws Exception {
		this.dataIter=dataIter;
		this.gets=gets;
	}
	
	/**
	 */
	public StringIter(Enumeration e, String[] gets) throws Exception {
		this.dataIter=new EnumIter(e);
		this.gets=gets;
	}
	
	/**
	   Need the list of gets and an iter- also take delim
	 */
	public StringIter(Iterator dataIter, String[] gets, String d) throws Exception {
		this.dataIter=dataIter;
		delim=d;
		this.gets=gets;
    }
	
	/**
	   change delimeter between fields
	 */
	public void setDelimiter(String d) {
		this.delim=d;
	}
	
	/**
	   Take next data object and create the String  and return it
	 */
	public Object next() {
		try {
			Object obj = dataIter.next();
			if (getms==null) init(obj.getClass());
			StringBuffer ret = new StringBuffer(20);
			for (int i = 0; i < getms.length;i++) {
				try {
					ret.append(getms[i].invoke(obj,null));
				} catch (Exception e) {
					e.printStackTrace();
					ret.append("ERROR");
				}
				if (i <  getms.length-1) ret.append(delim);
			}
            return  ret;
		} catch (Exception e) {
			e.printStackTrace();
		}
		return null;
    }
	
	/**
	 */
	public boolean hasNext() {
		return	dataIter.hasNext();
    }
	
	/**
	   optional - needed to satisfy interface Iterator
	 */
	public void remove() {}
	
	/**
	   for enumeration
	 */
	public Object nextElement() { 
		return next();
	}
	
	/**
	   for enumeration
	 */
	public boolean hasMoreElements() {
		return  dataIter.hasNext();
	}
	
	/**
	 */
	public void init(Class cls) throws Exception {
         getms = new Method[gets.length];
         for (int i = 0; i < gets.length; i++) {
                getms[i] = cls.getMethod(gets[i],null);
         }
    }
}
