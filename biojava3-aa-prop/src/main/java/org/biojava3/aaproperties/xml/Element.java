package org.biojava3.aaproperties.xml;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlElement;

/**
 * One way to model the elements
 * @author pvtroshin
 *
 */
@XmlAccessorType(XmlAccessType.NONE)
public class Element {

	/**
	 * Element name as per periodic table e.g. "Hydrogen"
	 */
	@XmlAttribute(name = "name", required = true)
	private String name; 
	/**
	 * Element short name as in periodic table e.g. "H"
	 */
	@XmlAttribute(name = "symbol")
	private String symbol; 
	/**
	 * The atomic number of the element = number of protons. 
	 */
	@XmlAttribute(name = "atomicnumber")
	private int atomicNumber; 
	/**
	 * The computed mass based on isotopes and their abundances
	 */
	@XmlAttribute(name = "mass", required = true)
	private double mass;
	/**
	 * List of common isotopes of the element
	 */
	@XmlElement
	private List<Isotope> isotope;
	
	/**
	 * To enable quick retrieval of Isotope from its name
	 */
	private Map<String, Isotope> name2Isotope;
	
	public Element(){}
	
	public Element(String name, String symbol, int atomicNumber, List<Isotope> isotopes, double mass){
		if(atomicNumber <= 0){
			throw new Error("Atomic number of Elements must be > 0.");
		}
		if(mass <= 0){
			throw new Error("Mass of Elements must be > 0.");
		}
		this.setName(name);
		this.setSymbol(symbol);
		this.setAtomicNumber(atomicNumber);
		this.setIsotopes(isotopes);
		this.setMass(mass);
	}
	
	public String toString(){
		return symbol + ", " + name + ", " + atomicNumber;
	}
	
	public void setMass(double mass){
		this.mass = mass;
	}
	
	public double getMass(){
		return this.mass;
	}

	public void setName(String name) {
		this.name = name;
	}
	
	public String getName(){
		return this.name;
	}

	public void setSymbol(String symbol) {
		this.symbol = symbol;
	}

	public void setAtomicNumber(int atomicNumber) {
		this.atomicNumber = atomicNumber;
	}


	public List<Isotope> getIsotopes() {
		return isotope;
	}


	public void setIsotopes(List<Isotope> isotopes) {
		this.isotope = isotopes;
		this.name2Isotope = new HashMap<String, Isotope>();
		if(isotopes != null){
			for(Isotope i:isotopes){
				name2Isotope.put(i.getName(), i);
			}
		}
	}
}
