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
package org.biojava.nbio.aaproperties.xml;

import javax.xml.bind.annotation.*;
import java.util.List;

@XmlRootElement(name = "compoundcomposition", namespace ="http://biojava.org")
@XmlAccessorType(XmlAccessType.NONE)
public class AminoAcidComposition {
	/**
	 * Amino acid symbol - single character
	 */
	@XmlAttribute(name = "symbol", required = true)
	private String symbol;
	/**
	 * Amino acid short name - three characters
	 */
	@XmlAttribute(name = "shortname")
	private String shortName;
	/**
	 * Amino acid full name
	 */
	@XmlAttribute(name = "name")
	private String name;
	/**
	 * Amino acid mass based on MOD ID
	 */
	@XmlAttribute(name = "id")
	String id;
	/**
	 * Stores the name of the element and the amount this amino acid contains
	 */
	@XmlElement(name = "element")
	private List<Name2Count> elementList;
	/**
	 * Stores the name of the isotope and the amount this amino acid contains
	 */
	@XmlElement(name = "isotope", required = false)
	private List<Name2Count> isotopeList;

	public AminoAcidComposition(){}

	public AminoAcidComposition(String symbol, String shortName, String name, List<Name2Count> elementList,
			List<Name2Count> isotopeList){
		if(symbol.length() != 1){
			throw new Error("Symbol attribute must be of length 1.");
		}
		this.symbol = symbol.toUpperCase();
		this.shortName = shortName;
		this.name = name;
		this.elementList = elementList;
		this.isotopeList = isotopeList;
	}

	@Override
	public String toString(){
		return symbol + ", " + shortName + ", " + name;
	}

	public String getSymbol() {
		return symbol;
	}

	public void setSymbol(String symbol) {
		this.symbol = symbol;
	}

	public void setShortName(String shortName) {
		this.shortName = shortName;
	}

	public String getShorName(){
		return this.shortName;
	}

	public void setName(String name) {
		this.name = name;
	}

	public String getName(){
		return this.name;
	}

	public List<Name2Count> getElementList(){
		return this.elementList;
	}

	public List<Name2Count> getIsotopeList(){
		return this.isotopeList;
	}
}
