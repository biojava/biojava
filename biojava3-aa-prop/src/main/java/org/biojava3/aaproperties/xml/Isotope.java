package org.biojava3.aaproperties.xml;

public class Isotope {
	/**
	 * The isotope name i.e. T (tritium)
	 */
	String name; 
	/**
	 * Number of neutrons 
	 */
	int neutronsNum;
	/**
	 * Relative Atomic Mass of the isotope
	 */
	double weight; 
	/**
	 * Natural abundance 
	 */
	double abundance;

	public Isotope(){}

	public Isotope(String name, int neutronsNum, double weight, double abundance){
		this.setName(name);
		this.setNeutronsNum(neutronsNum);
		this.setWeight(weight);
		this.setAbundance(abundance);
	}

	public String getName() {
		return name;
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
	public double getWeight() {
		return weight;
	}
	public void setWeight(double weight) {
		this.weight = weight;
	}
	public double getAbundance() {
		return abundance;
	}
	public void setAbundance(double abundance) {
		this.abundance = abundance;
	} 
}
