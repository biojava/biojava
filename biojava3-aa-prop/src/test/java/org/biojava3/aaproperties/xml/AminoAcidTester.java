package org.biojava3.aaproperties.xml;

import static junit.framework.Assert.assertEquals;

import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.List;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.JAXBException;
import javax.xml.bind.Marshaller;
import javax.xml.bind.Unmarshaller;

import org.biojava3.aaproperties.PeptideProperties;
import org.junit.Test;

public class AminoAcidTester {
	@Test
	public void generateSchema() throws JAXBException, IOException{
		JAXBContext context = JAXBContext.newInstance(AminoAcidCompositionTable.class);
		context.generateSchema(new SchemaGenerator("./src/main/resources/AminoAcidComposition.xsd"));
	}
	
	@Test
	public void readAdvancedXml() throws JAXBException, IOException{
		ElementTable iTable = new ElementTable();
		// Get a JAXB Context for the object we created above
		JAXBContext jc = JAXBContext.newInstance(iTable.getClass());
		Unmarshaller u = jc.createUnmarshaller();
		iTable = (ElementTable)u.unmarshal(new FileInputStream("./src/main/resources/AdvancedElementMass.xml" ) );
		iTable.populateMaps();
		
		AminoAcidCompositionTable aTable = new AminoAcidCompositionTable();
		// Get a JAXB Context for the object we created above
		JAXBContext jc2 = JAXBContext.newInstance(aTable.getClass());
		Unmarshaller u2 = jc2.createUnmarshaller();
		
		aTable = (AminoAcidCompositionTable)u2.unmarshal(new FileInputStream("./src/main/resources/AdvancedAminoAcidComposition.xml" ) );
		aTable.computeMolecularWeight(iTable);
		//Assert the weight of the radioactives
		String sequence = "00000";
		assertEquals(398.558744445, PeptideProperties.getMolecularWeightBasedOnXML(sequence, aTable));
		sequence = "1111";
		assertEquals(702.335483556, PeptideProperties.getMolecularWeightBasedOnXML(sequence, aTable));
		sequence = "JJJJ";
		assertEquals(0.0, PeptideProperties.getMolecularWeightBasedOnXML(sequence, aTable));
	}
	
	@Test
	public void readWithIDXml() throws JAXBException, IOException{
		ElementTable iTable = new ElementTable();
		// Get a JAXB Context for the object we created above
		JAXBContext jc = JAXBContext.newInstance(iTable.getClass());
		Unmarshaller u = jc.createUnmarshaller();
		iTable = (ElementTable)u.unmarshal(new FileInputStream("./src/main/resources/MinElementMass.xml" ) );
		iTable.populateMaps();
		
		AminoAcidCompositionTable aTable = new AminoAcidCompositionTable();
		// Get a JAXB Context for the object we created above
		JAXBContext jc2 = JAXBContext.newInstance(aTable.getClass());
		Unmarshaller u2 = jc2.createUnmarshaller();
		
		aTable = (AminoAcidCompositionTable)u2.unmarshal(new FileInputStream("./src/main/resources/AminoAcidCompositionWithID.xml" ) );
		aTable.computeMolecularWeight(iTable);
	}
	
	
	@Test
	public void readMinXml() throws JAXBException, IOException{
		ElementTable iTable = new ElementTable();
		// Get a JAXB Context for the object we created above
		JAXBContext jc = JAXBContext.newInstance(iTable.getClass());
		Unmarshaller u = jc.createUnmarshaller();
		iTable = (ElementTable)u.unmarshal(new FileInputStream("./src/main/resources/MinElementMass.xml" ) );
		iTable.populateMaps();
		
		AminoAcidCompositionTable aTable = new AminoAcidCompositionTable();
		// Get a JAXB Context for the object we created above
		JAXBContext jc2 = JAXBContext.newInstance(aTable.getClass());
		Unmarshaller u2 = jc2.createUnmarshaller();
		
		aTable = (AminoAcidCompositionTable)u2.unmarshal(new FileInputStream("./src/main/resources/MinAminoAcidComposition.xml" ) );
		aTable.computeMolecularWeight(iTable);
		
		String sequence = "AAAAA";
		assertEquals(373.4047, PeptideProperties.getMolecularWeightBasedOnXML(sequence, aTable));
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
		
		String sequence = "AAAAA";
		assertEquals(373.4047, PeptideProperties.getMolecularWeightBasedOnXML(sequence, aTable));
	}
	
	@Test
	public void generateXml() throws JAXBException, IOException{
		List<Name2Count> elementList = new ArrayList<Name2Count>();
		elementList.add(new Name2Count("Carbon", 3));
		elementList.add(new Name2Count("Hydrogen", 5));
		elementList.add(new Name2Count("Nitrogen", 1));
		elementList.add(new Name2Count("Oxygen", 1));
		AminoAcidComposition alanine = new AminoAcidComposition("A", "Ala", "Alanine", elementList, null);
		
		elementList = new ArrayList<Name2Count>();
		elementList.add(new Name2Count("Carbon", 6));
		elementList.add(new Name2Count("Hydrogen", 12));
		elementList.add(new Name2Count("Nitrogen", 4));
		elementList.add(new Name2Count("Oxygen", 1));
		AminoAcidComposition arginine = new AminoAcidComposition("R", "Arg", "Arginine", elementList, null);
		
		elementList = new ArrayList<Name2Count>();
		elementList.add(new Name2Count("Carbon", 4));
		elementList.add(new Name2Count("Hydrogen", 6));
		elementList.add(new Name2Count("Nitrogen", 2));
		elementList.add(new Name2Count("Oxygen", 2));
		AminoAcidComposition asparagine = new AminoAcidComposition("N", "Asn", "Asparagine", elementList, null);
	
		elementList = new ArrayList<Name2Count>();
		elementList.add(new Name2Count("Carbon", 4));
		elementList.add(new Name2Count("Hydrogen", 5));
		elementList.add(new Name2Count("Nitrogen", 1));
		elementList.add(new Name2Count("Oxygen", 3));
		AminoAcidComposition asparticAcid = new AminoAcidComposition("D", "Asp", "Aspartic Acid", elementList, null);
		
		elementList = new ArrayList<Name2Count>();
		
		elementList.add(new Name2Count("Carbon", 3));
		elementList.add(new Name2Count("Hydrogen", 5));
		elementList.add(new Name2Count("Nitrogen", 1));
		elementList.add(new Name2Count("Oxygen", 1));
		elementList.add(new Name2Count("Sulfur", 1));
		AminoAcidComposition cysteine = new AminoAcidComposition("C", "Cys", "Cysteine", elementList, null);
		
		elementList = new ArrayList<Name2Count>();
		
		elementList.add(new Name2Count("Carbon", 5));
		elementList.add(new Name2Count("Hydrogen", 7));
		elementList.add(new Name2Count("Nitrogen", 1));
		elementList.add(new Name2Count("Oxygen", 3));
		AminoAcidComposition glutamicAcid = new AminoAcidComposition("E", "Glu", "Glutamic Acid", elementList, null);
		
		elementList = new ArrayList<Name2Count>();
		
		elementList.add(new Name2Count("Carbon", 5));
		elementList.add(new Name2Count("Hydrogen", 8));
		elementList.add(new Name2Count("Nitrogen", 2));
		elementList.add(new Name2Count("Oxygen", 2));
		AminoAcidComposition glutamine = new AminoAcidComposition("Q", "Gln", "Glutamine", elementList, null);
		
		elementList = new ArrayList<Name2Count>();
		
		elementList.add(new Name2Count("Carbon", 2));
		elementList.add(new Name2Count("Hydrogen", 3));
		elementList.add(new Name2Count("Nitrogen", 1));
		elementList.add(new Name2Count("Oxygen", 1));
		AminoAcidComposition glycine = new AminoAcidComposition("G", "Gly", "Glycine", elementList, null);
		
		elementList = new ArrayList<Name2Count>();
		
		elementList.add(new Name2Count("Carbon", 6));
		elementList.add(new Name2Count("Hydrogen", 7));
		elementList.add(new Name2Count("Nitrogen", 3));
		elementList.add(new Name2Count("Oxygen", 1));
		AminoAcidComposition histidine = new AminoAcidComposition("H", "His", "Histidine", elementList, null);
		
		elementList = new ArrayList<Name2Count>();
		
		elementList.add(new Name2Count("Carbon", 6));
		elementList.add(new Name2Count("Hydrogen", 11));
		elementList.add(new Name2Count("Nitrogen", 1));
		elementList.add(new Name2Count("Oxygen", 1));
		AminoAcidComposition isoleucine = new AminoAcidComposition("I", "Ile", "Isoleucine", elementList, null);
		
		elementList = new ArrayList<Name2Count>();
		
		elementList.add(new Name2Count("Carbon", 6));
		elementList.add(new Name2Count("Hydrogen", 11));
		elementList.add(new Name2Count("Nitrogen", 1));
		elementList.add(new Name2Count("Oxygen", 1));
		AminoAcidComposition leucine = new AminoAcidComposition("L", "Leu", "Leucine", elementList, null);
		
		elementList = new ArrayList<Name2Count>();
		
		elementList.add(new Name2Count("Carbon", 6));
		elementList.add(new Name2Count("Hydrogen", 12));
		elementList.add(new Name2Count("Nitrogen", 2));
		elementList.add(new Name2Count("Oxygen", 1));
		AminoAcidComposition lysine = new AminoAcidComposition("K", "Lys", "Lysine", elementList, null);
		
		elementList = new ArrayList<Name2Count>();
		
		elementList.add(new Name2Count("Carbon", 5));
		elementList.add(new Name2Count("Hydrogen", 9));
		elementList.add(new Name2Count("Nitrogen", 1));
		elementList.add(new Name2Count("Oxygen", 1));
		elementList.add(new Name2Count("Sulfur", 1));
		AminoAcidComposition methionine = new AminoAcidComposition("M", "Met", "Methionine", elementList, null);
		
		elementList = new ArrayList<Name2Count>();
		
		elementList.add(new Name2Count("Carbon", 9));
		elementList.add(new Name2Count("Hydrogen", 9));
		elementList.add(new Name2Count("Nitrogen", 1));
		elementList.add(new Name2Count("Oxygen", 1));
		AminoAcidComposition phenylalanine = new AminoAcidComposition("F", "Phe", "Phenylalanine", elementList, null);
		
		elementList = new ArrayList<Name2Count>();
		
		elementList.add(new Name2Count("Carbon", 5));
		elementList.add(new Name2Count("Hydrogen", 7));
		elementList.add(new Name2Count("Nitrogen", 1));
		elementList.add(new Name2Count("Oxygen", 1));
		AminoAcidComposition proline = new AminoAcidComposition("P", "Pro", "Proline", elementList, null);
		
		elementList = new ArrayList<Name2Count>();
		
		elementList.add(new Name2Count("Carbon", 3));
		elementList.add(new Name2Count("Hydrogen", 5));
		elementList.add(new Name2Count("Nitrogen", 1));
		elementList.add(new Name2Count("Oxygen", 2));
		AminoAcidComposition serine = new AminoAcidComposition("S", "Ser", "Serine", elementList, null);
		
		elementList = new ArrayList<Name2Count>();
		
		elementList.add(new Name2Count("Carbon", 4));
		elementList.add(new Name2Count("Hydrogen", 7));
		elementList.add(new Name2Count("Nitrogen", 1));
		elementList.add(new Name2Count("Oxygen", 2));
		AminoAcidComposition threonine = new AminoAcidComposition("T", "Thr", "Threonine", elementList, null);
		
		elementList = new ArrayList<Name2Count>();
		
		elementList.add(new Name2Count("Carbon", 11));
		elementList.add(new Name2Count("Hydrogen", 10));
		elementList.add(new Name2Count("Nitrogen", 2));
		elementList.add(new Name2Count("Oxygen", 1));
		AminoAcidComposition tryptophan = new AminoAcidComposition("W", "Trp", "Tryptophan", elementList, null);
		
		elementList = new ArrayList<Name2Count>();
		
		elementList.add(new Name2Count("Carbon", 9));
		elementList.add(new Name2Count("Hydrogen", 9));
		elementList.add(new Name2Count("Nitrogen", 1));
		elementList.add(new Name2Count("Oxygen", 2));
		AminoAcidComposition tyrosine = new AminoAcidComposition("Y", "Tyr", "Tyrosine", elementList, null);
		
		elementList = new ArrayList<Name2Count>();
		
		elementList.add(new Name2Count("Carbon", 5));
		elementList.add(new Name2Count("Hydrogen", 9));
		elementList.add(new Name2Count("Nitrogen", 1));
		elementList.add(new Name2Count("Oxygen", 1));
		AminoAcidComposition valine = new AminoAcidComposition("V", "Val", "Valine", elementList, null);
		
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
