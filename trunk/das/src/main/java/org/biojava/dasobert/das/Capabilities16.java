package org.biojava.dasobert.das;

public class Capabilities16 {
	
	enum Capabilities{
		ENTRY_POINTS("entry_points"), FEATURES("features"),STYLESHEET("styesheet"),SEQUENCE("sequence"),ALIGNMENT("alignment"),STRUCTURE("structure"),TYPES("types"),INTERACTION("interaction"),SOURCES("sources");
	
	
	private String name;
	
	Capabilities(String name){
		this.name=name;
	}
	
	public String getName(){
		return this.name;
	}
	}
}

	

