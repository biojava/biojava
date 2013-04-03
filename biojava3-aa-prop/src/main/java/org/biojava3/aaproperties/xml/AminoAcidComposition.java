package org.biojava3.aaproperties.xml;

import java.util.List;

import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlAccessType;

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
