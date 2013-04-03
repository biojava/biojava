package org.biojava3.aaproperties.xml;

import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlAttribute;

@XmlAccessorType(XmlAccessType.FIELD)
public class Name2Count{
	@XmlAttribute(name = "name", required = true)
	private String name;
	@XmlAttribute(name = "count", required = true)
	private int count;
	
	public Name2Count(){}
	
	public Name2Count(String n, int c){
		if(c <= 0){
			throw new Error("Count must be > 0.");
		}
		this.name = n;
		this.count = c;
	}
	
	public String getName(){
		return this.name;
	}
	
	public int getCount(){
		return this.count;
	}
}