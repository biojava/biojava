/*
 *					BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 *
 *	  http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *	  http://www.biojava.org/
 *
 */
/**
 * 
 */
package org.biojava.nbio.core.sequence.features;

/**
 * @author mckeee1
 *
 */
public class Qualifier {
	private String[] values;
	private String name;

	boolean needsQuotes;

	/**
	 * initialize qualifier with name (=key) and value 
	 */
	public Qualifier(String name, String value) {
		// TODO Auto-generated constructor stub
		this.name=name;
		this.values=new String[1];
		this.values[0]=value;
		needsQuotes = false;
	}
	/**
	 * initialize Qualifier with name (=key) and value array
	 * @param name
	 * @param values
	 */
	public Qualifier(String name, String[] values) {
		this.name=name;
		this.values=values;
	}
	/**
	 * 
	 * @param name
	 * @param value
	 * @param needsQuotes
	 */
	public Qualifier(String name, String value, boolean needsQuotes) {
		// TODO Auto-generated constructor stub
		this.name=name;
		this.values=new String[1];
		this.values[0]=value;
		this.needsQuotes = needsQuotes;
	}
	/**
	 * 
	 * @param name
	 * @param values
	 * @param needsQuotes
	 */
	public Qualifier(String name, String[] values, boolean needsQuotes) {
		this.name=name;
		this.values=values;
		this.needsQuotes = needsQuotes;
	}
	

	/**
	 * @return the name (=key)
	 */
	public String getName() {
		return name;
	}

	/**
	 * @return the first value of qualifier (which most often contains only one value)
	 */
	public String getFirstValue() {
		return values[0];
	}
	/**
	 * 
	 * @param i
	 * @return value i if exists, else return null
	 */
	public String getValue(int i) {
		if(i<values.length) return values[i];
		return null;
	}
	
	/**
	 * return values as array
	 * @return
	 */
	public String[] getValues() {
		return values;
	}

	/**
	 * @return the needsQuotes
	 */
	public boolean needsQuotes() {
		return needsQuotes;
	}
	/**
	 * @param set the name (=key) of the qualifier
	 **/
	public void setName(String name) {
		this.name = name;
	}
	/**
	 * @param needsQuotes the needsQuotes to set
	 */
	public void setNeedsQuotes(boolean needsQuotes) {
		this.needsQuotes = needsQuotes;
	}

	/**
	 * @param value the value to set
	 */
	public void setValue(String value) {
		this.values = new String[1];
		this.values[0]=value;
	}
	/**
	 * set value array
	 * @param values
	 */
	public void setValues(String[] values) {
		this.values=values;
	}
	
	/**
	 * add value to values
	 * @param value
	 */
	public void addValue(String newValue) {
		String[] tmp=values;
		this.values=new String[tmp.length+1];
		int i=0;
		while(i<tmp.length) {
			this.values[i]=tmp[i];
			i++;
		}
		this.values[i]=newValue;
	}
		
	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append("Qualifier[ name='" + name +"' value='" + values[0]);
		if(values.length>1) for(int i=0;i<values.length;i++) sb.append(", "+values[i]);
		sb.append("]");
		return sb.toString();	   	
	}
	
	/**
	 * return the name (=key) of the qualifier
	 * @return
	 */
	public String getKey() {
		return getName();
	}
	/**
	 * get size of string array of qualifier
	 * @return
	 */
	public int valueSize() {
		return values.length;
	}
	/**
	 * add array of values to qualifier
	 * @param newValues
	 */
	public void addValues(String[] newValues) {
		String[] tmp=values;
		this.values=new String[tmp.length+newValues.length];
		int i=0;
		while(i<tmp.length) {
			this.values[i]=tmp[i];
			i++;
		}
		int j=0;
		while(i<tmp.length+values.length) {
			this.values[i]=newValues[j];
			i++;
			j++;
		}
	}
	
	/**
	 * 
	 * @param value
	 * @return true when qualifier values contain value
	 */
	public boolean containsValue(String value) {
		for(int i=0;i<values.length;i++) if(values[i].equals(value)) return true;
		return false;
	}
}
