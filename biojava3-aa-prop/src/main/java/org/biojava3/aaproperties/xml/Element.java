package org.biojava3.aaproperties.xml;

import java.util.List;

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
	int atomic_number; 
	/**
	 * List of common isotopes of the element
	 */
	List<Isotope> isotopes; 

	/**
	 * Standard Atomic Weight
	 *   
	 */
	double getMass() { 
		throw new UnsupportedOperationException();
	}
	
	
	class Isotope { 
		/**
		 * The isotope name i.e. T (tritium)
		 */
		String name; 
		/**
		 * Number of neutrons 
		 */
		int neutrons_num;
		/**
		 * Relative Atomic Mass of the isotope
		 */
		double weight; 
		/**
		 * Natural abundance 
		 */
		double abandunce; 
	}
	
}
