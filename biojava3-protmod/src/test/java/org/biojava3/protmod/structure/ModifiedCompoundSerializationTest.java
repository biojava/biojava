/*
 *                    BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 *
 *      http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *      http://www.biojava.org/
 *
 * Created on Jul 31, 2010
 * Author: Jianjiong Gao 
 *
 */

package org.biojava3.protmod.structure;


import java.io.File;
import java.io.IOException;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;

import java.util.Set;

import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;

import org.biojava3.protmod.ProteinModification;
import org.biojava3.protmod.io.ModifiedCompoundXMLConverter;
import org.biojava3.protmod.structure.ProteinModificationIdentifier;

import junit.framework.TestCase;

public class ModifiedCompoundSerializationTest extends TestCase {
	boolean allOK = false;

	private String[][] strucs;
	
	public void setUp() {
		//strucs = ProteinModificationParserTest.setUpShortTest();
		strucs = ProteinModificationParserTest.setUpLongTest();
	}
	
	public void testMulti() {
		int count = 0;
		for ( String[] name : strucs){
			try {
				
				testXMLSerialization(name[0]);
				count++;
			} catch (Exception e){
				System.err.println("failed after running " + count + " serializations at PDB ID " + name[0]);
				e.printStackTrace();
				fail(e.getMessage());
			}
		}
	}
	
	@SuppressWarnings("unchecked") 
	public void testSerialization() throws StructureException, IOException, ClassNotFoundException {
		String pdbId = "1CAD";
		Structure struc = TmpAtomCache.cache.getStructure(pdbId);

		ProteinModificationIdentifier parser = new ProteinModificationIdentifier();
		parser.identify(struc, ProteinModification.allModifications());
		Set<ModifiedCompound> mcs = parser.getIdentifiedModifiedCompound();

		String file = System.getProperty("java.io.tmpdir") + File.separatorChar + pdbId;
		FileOutputStream fos = new FileOutputStream(file);
		ObjectOutputStream oos = new ObjectOutputStream(fos);
		oos.writeObject(mcs);
		oos.close();

		FileInputStream fin = new FileInputStream(file);
		ObjectInputStream ois = new ObjectInputStream(fin); 
		mcs = (Set<ModifiedCompound>) ois.readObject();
		ois.close();

		//		System.out.println(mcs);

		try {
			new File(file).delete();
		} catch (Exception e) {
			fail(e.getMessage());
		}

		//System.out.println(mcs);
	}

	public void test1CAD(){
		String pdbId = "1CAD";
		testXMLSerialization(pdbId);
	}
	
	public void test1UIS(){
		String pdbId = "1UIS";
		testXMLSerialization(pdbId);
	}
	
	public void testXMLSerialization(String pdbId){
		String xml = null;
		ModifiedCompound currentMC = null;
		try {
			
			Structure struc = TmpAtomCache.cache.getStructure(pdbId);

			ProteinModificationIdentifier parser = new ProteinModificationIdentifier();
			parser.identify(struc, ProteinModification.allModifications());
			Set<ModifiedCompound> mcs = parser.getIdentifiedModifiedCompound();
			
			for (ModifiedCompound mc : mcs){
				currentMC = mc;
				 xml =  doXMLSerialization(mc) ;
				//System.out.println("got XML: " + String.format("%n") + xml);
				ModifiedCompound newMC = getModifiedCompoundFromXML(xml);
				String xml2 = doXMLSerialization(newMC);
				assertEquals(xml,xml2);
				//System.out.println(xml2);
				//assertEquals("The two objects are not equal before and after XML serialization" , mc, newMC);
				//System.out.println(mc.getDescription());
				//System.out.println(newMC.getDescription());
			}
		} catch (Exception e){
			System.out.println("error when serializing " + pdbId);
			System.out.println(currentMC.getDescription());
			System.out.println(xml);
			e.printStackTrace();
			fail(e.getMessage());
		}
		xml = null;
		currentMC =null;
	}

	private ModifiedCompound getModifiedCompoundFromXML(String xml) {
		return ModifiedCompoundXMLConverter.fromXML(xml);

	}

	private String doXMLSerialization(ModifiedCompound mc) throws IOException{
		return ModifiedCompoundXMLConverter.toXML(mc);
	}
}
