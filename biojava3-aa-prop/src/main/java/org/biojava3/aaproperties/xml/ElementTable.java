package org.biojava3.aaproperties.xml;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.xml.bind.annotation.XmlRootElement;

@XmlRootElement(name="elements", namespace="http://biojava.org")
public class ElementTable {
	
	/**
	 * Contain the full list of elements
	 */
	private List<Element> element;
	
	/**
	 * Enable quick retrieval of Element from its name
	 */
	private Map<String, Element> elementName2Element;
	
	/**
	 * Enable quick retrieval of Isotop from its name
	 */
	private Map<String, Isotope> isotopeName2Isotope;
	
	public ElementTable(){}

	public ElementTable(List<Element> eList){
		this.setElement(eList);
	}
	
	public void setElement(List<Element> eList){
		this.element = eList;
		populateMaps();
	}
	
	/**
	 * Populate the Maps for quick retrieval
	 */
	public void populateMaps(){
		this.elementName2Element = new HashMap<String, Element>();
		this.isotopeName2Isotope = new HashMap<String, Isotope>();
		if(this.element != null){
			for(Element e:this.element){
				this.elementName2Element.put(e.getName(), e);
				if(e.getIsotopes() != null){
					for(Isotope i:e.getIsotopes()){
						this.isotopeName2Isotope.put(i.getName(), i);
					}
				}
			}
		}
	}

	public List<Element> getElement(){
		return this.element;
	}
	
	public Element getElement(String name){
		return this.elementName2Element.get(name);
	}
	
	public Isotope getIsotope(String name){
		return this.isotopeName2Isotope.get(name);
	}
}
