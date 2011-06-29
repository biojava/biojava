package org.biojava3.aaproperties.xml;

import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.JAXBException;
import javax.xml.bind.Marshaller;
import javax.xml.bind.Unmarshaller;

import org.junit.Test;

public class AminoAcidTester {
	
	@Test
	public void generateSchema() throws JAXBException, IOException{
		JAXBContext context = JAXBContext.newInstance(AminoAcidCompositionTable.class);
		context.generateSchema(new SchemaGenerator("./src/main/resources/AminoAcidComposition.xsd"));
	}
	
	@Test
	public void readXml() throws JAXBException, IOException{
		ElementTable iTable = new ElementTable();
		// Get a JAXB Context for the object we created above
		JAXBContext jc = JAXBContext.newInstance(iTable.getClass());
		Unmarshaller u = jc.createUnmarshaller();
		iTable = (ElementTable)u.unmarshal(new FileInputStream("./src/main/resources/ElementMass.xml" ) );
		iTable.populateMaps();
		
		AminoAcidCompositionTable aTable = new AminoAcidCompositionTable();
		// Get a JAXB Context for the object we created above
		JAXBContext jc2 = JAXBContext.newInstance(aTable.getClass());
		Unmarshaller u2 = jc2.createUnmarshaller();
		
		aTable = (AminoAcidCompositionTable)u2.unmarshal(new FileInputStream("./src/main/resources/AminoAcidComposition.xml" ) );
		aTable.computeMolecularWeight(iTable);
		for(AminoAcidComposition a:aTable.getAminoacid()){
			System.out.println(a + ", " + aTable.getMolecularWeight(a.getSymbol()));
		}
	}
	
	@Test
	public void generateXml() throws JAXBException, IOException{
		Map<String, Integer> elementName2Count = new HashMap<String, Integer>();
		Map<String, Integer> isotopeName2Count = new HashMap<String, Integer>();
		elementName2Count.put("Carbon", 3);
		elementName2Count.put("Hydrogen", 5);
		elementName2Count.put("Nitrogen", 1);
		elementName2Count.put("Oxygen", 1);
		AminoAcidComposition alanine = new AminoAcidComposition("A", "Ala", "Alanine", elementName2Count, isotopeName2Count);
		
		elementName2Count = new HashMap<String, Integer>();
		isotopeName2Count = new HashMap<String, Integer>();
		elementName2Count.put("Carbon", 6);
		elementName2Count.put("Hydrogen", 12);
		elementName2Count.put("Nitrogen", 4);
		elementName2Count.put("Oxygen", 1);
		AminoAcidComposition arginine = new AminoAcidComposition("R", "Arg", "Arginine", elementName2Count, isotopeName2Count);
		
		elementName2Count = new HashMap<String, Integer>();
		isotopeName2Count = new HashMap<String, Integer>();
		elementName2Count.put("Carbon", 4);
		elementName2Count.put("Hydrogen", 6);
		elementName2Count.put("Nitrogen", 2);
		elementName2Count.put("Oxygen", 2);
		AminoAcidComposition asparagine = new AminoAcidComposition("N", "Asn", "Asparagine", elementName2Count, isotopeName2Count);
		
		elementName2Count = new HashMap<String, Integer>();
		isotopeName2Count = new HashMap<String, Integer>();
		elementName2Count.put("Carbon", 4);
		elementName2Count.put("Hydrogen", 5);
		elementName2Count.put("Nitrogen", 1);
		elementName2Count.put("Oxygen", 3);
		AminoAcidComposition asparticAcid = new AminoAcidComposition("D", "Asp", "Aspartic Acid", elementName2Count, isotopeName2Count);
		
		elementName2Count = new HashMap<String, Integer>();
		isotopeName2Count = new HashMap<String, Integer>();
		elementName2Count.put("Carbon", 3);
		elementName2Count.put("Hydrogen", 5);
		elementName2Count.put("Nitrogen", 1);
		elementName2Count.put("Oxygen", 1);
		elementName2Count.put("Sulfur", 1);
		AminoAcidComposition cysteine = new AminoAcidComposition("C", "Cys", "Cysteine", elementName2Count, isotopeName2Count);
		
		elementName2Count = new HashMap<String, Integer>();
		isotopeName2Count = new HashMap<String, Integer>();
		elementName2Count.put("Carbon", 5);
		elementName2Count.put("Hydrogen", 7);
		elementName2Count.put("Nitrogen", 1);
		elementName2Count.put("Oxygen", 3);
		AminoAcidComposition glutamicAcid = new AminoAcidComposition("E", "Glu", "Glutamic Acid", elementName2Count, isotopeName2Count);
		
		elementName2Count = new HashMap<String, Integer>();
		isotopeName2Count = new HashMap<String, Integer>();
		elementName2Count.put("Carbon", 5);
		elementName2Count.put("Hydrogen", 8);
		elementName2Count.put("Nitrogen", 2);
		elementName2Count.put("Oxygen", 2);
		AminoAcidComposition glutamine = new AminoAcidComposition("Q", "Gln", "Glutamine", elementName2Count, isotopeName2Count);
		
		elementName2Count = new HashMap<String, Integer>();
		isotopeName2Count = new HashMap<String, Integer>();
		elementName2Count.put("Carbon", 2);
		elementName2Count.put("Hydrogen", 3);
		elementName2Count.put("Nitrogen", 1);
		elementName2Count.put("Oxygen", 1);
		AminoAcidComposition glycine = new AminoAcidComposition("G", "Gly", "Glycine", elementName2Count, isotopeName2Count);
		
		elementName2Count = new HashMap<String, Integer>();
		isotopeName2Count = new HashMap<String, Integer>();
		elementName2Count.put("Carbon", 6);
		elementName2Count.put("Hydrogen", 7);
		elementName2Count.put("Nitrogen", 3);
		elementName2Count.put("Oxygen", 1);
		AminoAcidComposition histidine = new AminoAcidComposition("H", "His", "Histidine", elementName2Count, isotopeName2Count);
		
		elementName2Count = new HashMap<String, Integer>();
		isotopeName2Count = new HashMap<String, Integer>();
		elementName2Count.put("Carbon", 6);
		elementName2Count.put("Hydrogen", 11);
		elementName2Count.put("Nitrogen", 1);
		elementName2Count.put("Oxygen", 1);
		AminoAcidComposition isoleucine = new AminoAcidComposition("I", "Ile", "Isoleucine", elementName2Count, isotopeName2Count);
		
		elementName2Count = new HashMap<String, Integer>();
		isotopeName2Count = new HashMap<String, Integer>();
		elementName2Count.put("Carbon", 6);
		elementName2Count.put("Hydrogen", 11);
		elementName2Count.put("Nitrogen", 1);
		elementName2Count.put("Oxygen", 1);
		AminoAcidComposition leucine = new AminoAcidComposition("L", "Leu", "Leucine", elementName2Count, isotopeName2Count);
		
		elementName2Count = new HashMap<String, Integer>();
		isotopeName2Count = new HashMap<String, Integer>();
		elementName2Count.put("Carbon", 6);
		elementName2Count.put("Hydrogen", 12);
		elementName2Count.put("Nitrogen", 2);
		elementName2Count.put("Oxygen", 1);
		AminoAcidComposition lysine = new AminoAcidComposition("K", "Lys", "Lysine", elementName2Count, isotopeName2Count);
		
		elementName2Count = new HashMap<String, Integer>();
		isotopeName2Count = new HashMap<String, Integer>();
		elementName2Count.put("Carbon", 5);
		elementName2Count.put("Hydrogen", 9);
		elementName2Count.put("Nitrogen", 1);
		elementName2Count.put("Oxygen", 1);
		elementName2Count.put("Sulfur", 1);
		AminoAcidComposition methionine = new AminoAcidComposition("M", "Met", "Methionine", elementName2Count, isotopeName2Count);
		
		elementName2Count = new HashMap<String, Integer>();
		isotopeName2Count = new HashMap<String, Integer>();
		elementName2Count.put("Carbon", 9);
		elementName2Count.put("Hydrogen", 9);
		elementName2Count.put("Nitrogen", 1);
		elementName2Count.put("Oxygen", 1);
		AminoAcidComposition phenylalanine = new AminoAcidComposition("F", "Phe", "Phenylalanine", elementName2Count, isotopeName2Count);
		
		elementName2Count = new HashMap<String, Integer>();
		isotopeName2Count = new HashMap<String, Integer>();
		elementName2Count.put("Carbon", 5);
		elementName2Count.put("Hydrogen", 7);
		elementName2Count.put("Nitrogen", 1);
		elementName2Count.put("Oxygen", 1);
		AminoAcidComposition proline = new AminoAcidComposition("P", "Pro", "Proline", elementName2Count, isotopeName2Count);
		
		elementName2Count = new HashMap<String, Integer>();
		isotopeName2Count = new HashMap<String, Integer>();
		elementName2Count.put("Carbon", 3);
		elementName2Count.put("Hydrogen", 5);
		elementName2Count.put("Nitrogen", 1);
		elementName2Count.put("Oxygen", 2);
		AminoAcidComposition serine = new AminoAcidComposition("S", "Ser", "Serine", elementName2Count, isotopeName2Count);
		
		elementName2Count = new HashMap<String, Integer>();
		isotopeName2Count = new HashMap<String, Integer>();
		elementName2Count.put("Carbon", 4);
		elementName2Count.put("Hydrogen", 7);
		elementName2Count.put("Nitrogen", 1);
		elementName2Count.put("Oxygen", 2);
		AminoAcidComposition threonine = new AminoAcidComposition("T", "Thr", "Threonine", elementName2Count, isotopeName2Count);
		
		elementName2Count = new HashMap<String, Integer>();
		isotopeName2Count = new HashMap<String, Integer>();
		elementName2Count.put("Carbon", 11);
		elementName2Count.put("Hydrogen", 10);
		elementName2Count.put("Nitrogen", 2);
		elementName2Count.put("Oxygen", 1);
		AminoAcidComposition tryptophan = new AminoAcidComposition("W", "Trp", "Tryptophan", elementName2Count, isotopeName2Count);
		
		elementName2Count = new HashMap<String, Integer>();
		isotopeName2Count = new HashMap<String, Integer>();
		elementName2Count.put("Carbon", 9);
		elementName2Count.put("Hydrogen", 9);
		elementName2Count.put("Nitrogen", 1);
		elementName2Count.put("Oxygen", 2);
		AminoAcidComposition tyrosine = new AminoAcidComposition("Y", "Tyr", "Tyrosine", elementName2Count, isotopeName2Count);
		
		elementName2Count = new HashMap<String, Integer>();
		isotopeName2Count = new HashMap<String, Integer>();
		elementName2Count.put("Carbon", 5);
		elementName2Count.put("Hydrogen", 9);
		elementName2Count.put("Nitrogen", 1);
		elementName2Count.put("Oxygen", 1);
		AminoAcidComposition valine = new AminoAcidComposition("V", "Val", "Valine", elementName2Count, isotopeName2Count);
		
		List<AminoAcidComposition> aList = new ArrayList<AminoAcidComposition>();
		aList.add(alanine);
		aList.add(arginine);
		aList.add(asparagine);
		aList.add(asparticAcid);
		aList.add(cysteine);
		aList.add(glutamicAcid);
		aList.add(glutamine);
		aList.add(glycine);
		aList.add(histidine);
		aList.add(isoleucine);
		aList.add(leucine);
		aList.add(lysine);
		aList.add(methionine);
		aList.add(phenylalanine);
		aList.add(proline);
		aList.add(serine);
		aList.add(threonine);
		aList.add(tryptophan);
		aList.add(tyrosine);
		aList.add(valine);
		
		AminoAcidCompositionTable aTable = new AminoAcidCompositionTable(aList);
		// Get a JAXB Context for the object we created above
		JAXBContext context = JAXBContext.newInstance(aTable.getClass());
		// To convert ex to XML, I need a JAXB Marshaller
		Marshaller marshaller = context.createMarshaller();
		// Make the output pretty
		marshaller.setProperty(Marshaller.JAXB_FORMATTED_OUTPUT, true);
		StringWriter sw = new StringWriter();
		// marshall the object to XML
		marshaller.marshal(aTable, sw);
		// print it out for this example
		BufferedWriter output = new BufferedWriter(new FileWriter("./src/main/resources/AminoAcidComposition.xml"));
		output.write(sw.toString());
		output.close();
	}
}
