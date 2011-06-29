package org.biojava3.aaproperties.xml;

public class Isotope {
	/**
	 * The isotope name i.e. T (tritium)
	 */
	private String name; 
	/**
	 * Number of neutrons 
	 */
	private int neutronsNum;
	/**
	 * Relative Atomic Mass of the isotope
	 */
	private double weight; 

	public Isotope(){}

	public Isotope(String name, int neutronsNum, double weight){
		this.setName(name);
		this.setNeutronsNum(neutronsNum);
		this.setWeight(weight);
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
}
