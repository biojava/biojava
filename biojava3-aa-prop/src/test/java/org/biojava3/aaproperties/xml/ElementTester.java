package org.biojava3.aaproperties.xml;

import static junit.framework.Assert.assertEquals;

import java.io.StringWriter;
import java.util.ArrayList;
import java.util.List;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.JAXBException;
import javax.xml.bind.Marshaller;

import org.junit.Test;

public class ElementTester {
	@Test
	public void generateXml() throws JAXBException {
		List<Isotope> iList = new ArrayList<Isotope>();
		iList.add(new Isotope("Hydrogen", 1, 1.00782503207, 0.999885));
		iList.add(new Isotope("Deuterium", 2, 2.0141017778, 0.000115));
		iList.add(new Isotope("Tritium", 3, 3.0160492777, 0.0));
		Element hydrogen = new Element("Hydrogen", "H", 1, iList);
		assertEquals(1.00794, hydrogen.getMass(5));
		
		iList.clear();
		iList.add(new Isotope("Helium-3", 3, 3.0160293191, 0.00000134));
		iList.add(new Isotope("Helium-4", 4, 4.00260325415, 0.99999866));
		Element helium = new Element("Helium", "He", 2, iList);
		assertEquals(4.002602, helium.getMass(6));
		
		List<Element> eList = new ArrayList<Element>();
		eList.add(hydrogen);
		eList.add(helium);
		
		ElementTable iTable = new ElementTable(eList);

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

}
