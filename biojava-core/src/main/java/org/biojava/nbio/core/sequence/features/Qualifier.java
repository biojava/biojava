/*
 *                    BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 *
 *      http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *      http://www.biojava.org/
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

	@Override
	public String toString() {
		return "Qualifier[ name='" + name +"' value='"+ value + "' ]";
	}
}
