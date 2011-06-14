package org.biojava3.aaproperties.xml;

import java.util.List;

import org.biojava3.aaproperties.Utils;


/**
 * One way to model the elements
 * @author pvtroshin
 *
 */
public class Element {

	/**
	 * Element name as per periodic table e.g. "Hydrogen"
	 */
	String name; 
	/**
	 * Element short name as in periodic table e.g. "H"
	 */
	String shortName; 
	/**
	 * The atomic number of the element = number of protons. 
	 */
	int atomicNumber; 
	/**
	 * List of common isotopes of the element
	 */
	List<Isotope> isotopes; 
	
	public Element(){}
	
	public Element(String name, String shortName, int atomicNumber, List<Isotope> isotopes){
		this.setName(name);
		this.setShortName(shortName);
		this.setAtomicNumber(atomicNumber);
		this.setIsotopes(isotopes);
	}
	
	/**
	 * Returns the standard atomic weight. It is computed by the relative abundance of the isotope multiply by its atomic weight
	 * @return the standard atomic weight 
	 */
	public double getMass() { 
		return getMass(-1);
	}
	
	public double getMass(int decimalPlace){
		double total = 0.0;
		for(Isotope i:isotopes){
			total += i.getWeight() * i.getAbundance();
		}
		return Utils.roundToDecimals(total, decimalPlace);
	}
	
	public String getName() {
		return name;
	}


	public void setName(String name) {
		this.name = name;
	}


	public String getShortName() {
		return shortName;
	}


	public void setShortName(String shortName) {
		this.shortName = shortName;
	}


	public int getAtomicNumber() {
		return atomicNumber;
	}


	public void setAtomicNumber(int atomicNumber) {
		this.atomicNumber = atomicNumber;
	}


	public List<Isotope> getIsotopes() {
		return isotopes;
	}


	public void setIsotopes(List<Isotope> isotopes) {
		this.isotopes = isotopes;
	}
}
