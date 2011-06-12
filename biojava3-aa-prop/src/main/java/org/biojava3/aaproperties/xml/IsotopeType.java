package org.biojava3.aaproperties.xml;

public class IsotopeType {
    private int numOfProton;
    private String shortName;
    private double atomWeight;
    
    public IsotopeType(){}
    
    public IsotopeType(String shortName, int numOfProton, double atomWeight){
    	this.setShortName(shortName);
    	this.setNumOfProton(numOfProton);
    	this.setAtomWeight(atomWeight);
    }
	
    public int getNumOfProton() {
        return this.numOfProton;
    }
    public void setNumOfProton(int value) {
        this.numOfProton = value;
    }
    
    public String getShortName(){
    	return this.shortName;
    }
    
    public void setShortName(String value){
    	this.shortName = value;
    }
    
    public double getAtomWeight() {
    	return this.atomWeight;
    }
    
    public void setAtomWeight(double value){
    	this.atomWeight = value;
    }
}
