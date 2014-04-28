/**
 * 
 */
package org.biojava3.core.sequence.features;

/**
 * @author mckeee1
 *
 */
public class Qualifier {
	String value;
	String name;

	boolean needsQuotes;

	/**
	 * 
	 */
	public Qualifier(String name, String value) {
		// TODO Auto-generated constructor stub
		this.name=name;
		this.value=value;
		needsQuotes = false;
	}

	/**
	 * 
	 */
	public Qualifier(String name, String value, boolean needsQuotes) {
		// TODO Auto-generated constructor stub
		this.name=name;
		this.value=value;
		this.needsQuotes = needsQuotes;
	}

	/**
	 * @return the name
	 */
	public String getName() {
		return name;
	}

	/**
	 * @return the value
	 */
	public String getValue() {
		return value;
	}

	/**
	 * @return the needsQuotes
	 */
	public boolean needsQuotes() {
		return needsQuotes;
	}
	/**
	 * @param name the name to set
	 */
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
		this.value = value;
	}
}
