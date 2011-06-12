package org.biojava3.aaproperties.xml;

import java.io.StringWriter;
import java.util.ArrayList;
import java.util.List;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.JAXBException;
import javax.xml.bind.Marshaller;
import javax.xml.bind.annotation.XmlRootElement;

@XmlRootElement(name="atom", namespace="http://biojava.org")
public class IsotopeTable {
	private List<IsotopeType> isotopeList;
	
	public IsotopeTable(){}

	public IsotopeTable(List<IsotopeType> iList){
		this.setIsotopeList(iList);
	}
	
	public void setIsotopeList(List<IsotopeType> iList){
		this.isotopeList = iList;
	}

	public List<IsotopeType> getIsotopeList(){
		return this.isotopeList;
	}

	public static void generateXml() throws JAXBException {
		List<IsotopeType> iList = new ArrayList<IsotopeType>();
		iList.add(new IsotopeType("H", 1, 1.00782503207));
		iList.add(new IsotopeType("D", 2, 2.0141017778));
		iList.add(new IsotopeType("T", 3, 3.0160492777));
		iList.add(new IsotopeType("He", 3, 3.0160293191));
		iList.add(new IsotopeType("He", 4, 4.00260325415));
		
		IsotopeTable iTable = new IsotopeTable(iList);

		// Get a JAXB Context for the object we created above
		JAXBContext context = JAXBContext.newInstance(iTable.getClass());

		// To convert ex to XML, I need a JAXB Marshaller
		Marshaller marshaller = context.createMarshaller();

		// Make the output pretty
		marshaller.setProperty(Marshaller.JAXB_FORMATTED_OUTPUT, true);
		StringWriter sw = new StringWriter();

		// marshall the object to XML
		marshaller.marshal(iTable, sw);

		// print it out for this example
		System.out.println(sw.toString());
	}

	public static void main(String[] args){
		try{
			generateXml();
		}catch(JAXBException e){e.printStackTrace();}
	}
}
