package org.biojava3.protmod.structure;

import java.io.IOException;
import java.util.Set;

import org.biojava.bio.structure.Structure;
import org.biojava3.protmod.ProteinModification;
import org.biojava3.protmod.io.ProteinModificationXMLConverter;
import org.biojava3.protmod.structure.ModifiedCompound;
import org.biojava3.protmod.structure.ProteinModificationIdentifier;

import junit.framework.TestCase;

public class TestProteinModificationSerialilzation extends TestCase {

	public void test1CAD(){
		String pdbId = "1CAD";
		
		testXMLSerialization(pdbId);
	}
	
	
	public void testXMLSerialization(String pdbId){

		//ProteinModification.init();
		//ProteinModification fromEmpty = ProteinModification.getById("0252.CROSS_LINK_4");
		//assertNotNull("Could not get ProteinModification by ID!", fromEmpty);
		
		
		try {
			
			Structure struc = TmpAtomCache.cache.getStructure(pdbId);

			ProteinModificationIdentifier parser = new ProteinModificationIdentifier();
			parser.identify(struc, ProteinModification.allModifications());
			Set<ModifiedCompound> mcs = parser.getIdentifiedModifiedCompound();

			for (ModifiedCompound mc : mcs){
				ProteinModification modifcation = mc.getModification();
				String xml =  doXMLSerialization(modifcation) ;
				//System.out.println("got XML: " + String.format("%n") + xml + String.format("%n") + " for object " + modifcation);
				ProteinModification newMod = getProteinModificationFromXML(xml);
				assertEquals("The two objects are not equal before and after XML serialization" , modifcation, newMod);
				
				String xml2 = doXMLSerialization(newMod);
				assertEquals(xml,xml2);
				//assertEquals(modifcation, fromEmpty);
				
			}
		} catch (Exception e){
			e.printStackTrace();
			fail(e.getMessage());
		}
	}

	private ProteinModification getProteinModificationFromXML(String xml) {
		return ProteinModificationXMLConverter.fromXML(xml);
	}

	private String doXMLSerialization(ProteinModification modification) throws IOException {
		
		return ProteinModificationXMLConverter.toXML(modification);
	}
	
	
}
