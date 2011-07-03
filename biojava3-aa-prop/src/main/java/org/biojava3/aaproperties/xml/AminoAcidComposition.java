package org.biojava3.aaproperties.xml;

import java.util.Map;

import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlAccessType;

@XmlRootElement(name = "AminoAcidComposition", namespace ="http://biojava.org")
@XmlAccessorType(XmlAccessType.FIELD)
public class AminoAcidComposition {
	/**
	 * Amino acid symbol - single character
	 */
	@XmlAttribute(name = "symbol", required = true)
	private String symbol;
	/**
	 * Amino acid short name - three characters
	 */
	@XmlAttribute(name = "shortname", required = true)
	private String shortName;
	/**
	 * Amino acid full name
	 */
	@XmlAttribute(name = "name", required = true)
	private String name;
	/**
	 * Stores the name of the element and the amount this amino acid contains
	 */
	//@XmlElement(type=HashMap.class, nillable = false, defaultValue = "1", name = "elements")
	private Map<String, Integer> elementName2Count;
	/**
	 * Stores the name of the isotope and the amount this amino acid contains
	 */
	//@XmlElement(type=HashMap.class, name = "isotopes")
	private Map<String, Integer> isotopeName2Count;
	
	public AminoAcidComposition(){}
	
	public AminoAcidComposition(String symbol, String shortName, String name, Map<String, Integer> elementName2Count, 
			Map<String, Integer> isotopeName2Count){
		this.symbol = symbol;
		this.shortName = shortName;
		this.name = name;
		this.elementName2Count = elementName2Count;
		this.isotopeName2Count = isotopeName2Count;
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
	
	public void setName(String name) {
		this.name = name;
	}
	
	public Map<String, Integer> getIsotope2Count() {
		return isotopeName2Count;
	}
	
	public int getIsotopeCount(String isotopeName){
		return this.isotopeName2Count.get(isotopeName);
	}
	
	public void setIsotope2Count(Map<String, Integer> isotopeName2Count) {
		this.isotopeName2Count = isotopeName2Count;
	}
	
	public void setElement2Count(Map<String, Integer> elementName2Count) {
		this.elementName2Count = elementName2Count;
	}
	
	public Map<String, Integer> getElement2Count() {
		return elementName2Count;
	}
	
	public int getElementCount(String elementName){
		return this.elementName2Count.get(elementName);
	}
}
