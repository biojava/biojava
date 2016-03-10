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
import javax.xml.bind.annotation.XmlType;


@XmlType(name = "Iostope", propOrder = {"name","neutronsNum","mass"})
@XmlAccessorType(XmlAccessType.NONE)
public class Isotope {
	/**
	 * The isotope name i.e. T (tritium)
	 */
	@XmlAttribute(name = "name", required = true)
	private String name;
	/**
	 * Number of neutrons
	 */
	@XmlAttribute(name = "neutronsnum", required = true)
	private int neutronsNum;
	/**
	 * Relative Atomic Mass of the isotope
	 */
	@XmlAttribute(name = "mass", required = true)
	private double mass;

	public Isotope(){}

	public Isotope(String name, int neutronsNum, double mass){
		if(neutronsNum <= 0){
			throw new Error("Neutrons number of Isotopes must be > 0.");
		}
		if(mass <= 0){
			throw new Error("Mass of Isotopes must be > 0.");
		}
		this.setName(name);
		this.setNeutronsNum(neutronsNum);
		this.setMass(mass);
	}

	public String getName(){
		return this.name;
	}

	public void setName(String name) {
		this.name = name;
	}

	public int getNeutronsNum() {
		return neutronsNum;
	}

	public void setNeutronsNum(int neutronsNum) {
		this.neutronsNum = neutronsNum;
	}

	public double getMass() {
		return mass;
	}

	public void setMass(double weight) {
		this.mass = weight;
	}
}
