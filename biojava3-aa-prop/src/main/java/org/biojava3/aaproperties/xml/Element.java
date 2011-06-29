package org.biojava3.aaproperties.xml;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * One way to model the elements
 * @author pvtroshin
 *
 */
public class Element {

	/**
	 * Element name as per periodic table e.g. "Hydrogen"
	 */
	private String name; 
	/**
	 * Element short name as in periodic table e.g. "H"
	 */
	private String symbol; 
	/**
	 * The atomic number of the element = number of protons. 
	 */
	private int atomicNumber; 
	/**
	 * The computed mass based on isotopes and their abundances
	 */
	private double mass;
	/**
	 * List of common isotopes of the element
	 */
	private List<Isotope> isotope;
	
	/**
	 * To enable quick retrieval of Isotope from its name
	 */
	private Map<String, Isotope> name2Isotope;
	
	public Element(){}
	
	public Element(String name, String symbol, int atomicNumber, List<Isotope> isotopes, double mass){
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
	
	/**
	 * Returns the standard atomic weight. It is computed by the relative abundance of the isotope multiply by its atomic weight
	 * @return the standard atomic weight 
	 */
	public double getMass() { 
		/*double total = 0.0;
		for(Isotope i:isotope){
			total += (i.getWeight() * i.getAbundance());
		}
		this.mass = total;*/
		return mass;
	}
	
	public String getName() {
		return name;
	}


	public void setName(String name) {
		this.name = name;
	}


	public String getSymbol() {
		return symbol;
	}


	public void setSymbol(String symbol) {
		this.symbol = symbol;
	}


	public int getAtomicNumber() {
		return atomicNumber;
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
