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

import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlElement;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

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

	@Override
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
