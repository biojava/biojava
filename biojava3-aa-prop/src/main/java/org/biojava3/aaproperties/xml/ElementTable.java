package org.biojava3.aaproperties.xml;

import java.util.List;

import javax.xml.bind.annotation.XmlRootElement;

@XmlRootElement(name="elements", namespace="http://biojava.org")
public class ElementTable {
	
	private List<Element> element;
	
	public ElementTable(){}

	public ElementTable(List<Element> eList){
		this.setElement(eList);
	}
	
	public void setElement(List<Element> eList){
		this.element = eList;
	}

	public List<Element> getElement(){
		return this.element;
	}
}
