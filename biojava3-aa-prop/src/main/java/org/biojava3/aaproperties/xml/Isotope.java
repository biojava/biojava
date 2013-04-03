package org.biojava3.aaproperties.xml;

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
