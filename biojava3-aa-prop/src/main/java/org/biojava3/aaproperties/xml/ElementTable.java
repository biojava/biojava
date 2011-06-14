package org.biojava3.aaproperties.xml;

import java.util.List;

import javax.xml.bind.annotation.XmlRootElement;

@XmlRootElement(name="atom", namespace="http://biojava.org")
public class ElementTable {
	
	private List<Element> elementList;
	
	public ElementTable(){}

	public ElementTable(List<Element> eList){
		this.setElementList(eList);
	}
	
	public void setElementList(List<Element> eList){
		this.elementList = eList;
	}

	public List<Element> getElementList(){
		return this.elementList;
	}
}
