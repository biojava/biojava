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


import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.Reader;
import java.io.StringWriter;
import java.io.Writer;

import java.util.Set;

import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.Structure;

import org.biojava3.protmod.ProteinModification;
import org.biojava3.protmod.io.ModifiedCompoundXMLConverter;
import org.biojava3.protmod.structure.ProteinModificationIdentifier;

import junit.framework.TestCase;

public class ModifiedCompoundSerializationTest extends TestCase {
	boolean allOK = false;

	private String[][] strucs;

	public void setUp() {
		strucs = new String[1][1];

		strucs[0][0]="2TMD";
		//strucs = ProteinModificationParserTest.setUpShortTest();
		//strucs = ProteinModificationParserTest.setUpLongTest();
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

	public void test1CAD(){
		String pdbId = "1CAD";
		testXMLSerialization(pdbId);
	}

	public void test1UIS(){
		String pdbId = "1UIS";
		testXMLSerialization(pdbId);
	}
	
//	public void test2TMD(){
//		String pdbId = "2TMD";
//		testXMLSerialization(pdbId);
//	}
	
	public void test1CDG(){
		String pdbId = "1CDG";
		testXMLSerialization(pdbId);
	}

	public void testXMLSerialization(String pdbId){
		String xml = null;
		ModifiedCompound currentMC = null;
		try {

			Structure struc = TmpAtomCache.cache.getStructure(pdbId);

			ProteinModificationIdentifier parser = new ProteinModificationIdentifier();

			for (Chain c : struc.getChains()) {

				parser.identify(c, ProteinModification.allModifications());
				Set<ModifiedCompound> mcs = parser.getIdentifiedModifiedCompound();

				for (ModifiedCompound mc : mcs){
					currentMC = mc;
					xml =  doXMLSerialization(mc) ;
					//System.out.println( pdbId + " got XML: " + String.format("%n") + xml);
					ModifiedCompound newMC = getModifiedCompoundFromXML(xml);
					String xml2 = doXMLSerialization(newMC);
					assertEquals(xml,xml2);
					//System.out.println(xml2);
					//assertEquals("The two objects are not equal before and after XML serialization" , mc, newMC);
					//System.out.println(mc.getDescription());
					//System.out.println(newMC.getDescription());
				}
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

	public void testFlatFileParsing(){
		InputStream inStream = this.getClass().getResourceAsStream("/org/biojava3/protmod/parser/modifiedCompound.xml");
		assertNotNull(inStream);
		try {
			String xml = convertStreamToString(inStream);
			//System.out.println(xml);
			ModifiedCompound newMC = getModifiedCompoundFromXML(xml);

			assertNotNull(newMC);
		} catch (Exception e){
			fail(e.getMessage());
		}

	}

	public String convertStreamToString(InputStream is)
	throws IOException {
		/*
		 * To convert the InputStream to String we use the
		 * Reader.read(char[] buffer) method. We iterate until the
		 * Reader return -1 which means there's no more data to
		 * read. We use the StringWriter class to produce the string.
		 */
		if (is != null) {
			Writer writer = new StringWriter(); 
			char[] buffer = new char[1024];
			try {
				Reader reader = new BufferedReader(
						new InputStreamReader(is, "UTF-8"));
				int n;
				while ((n = reader.read(buffer)) != -1) {
					writer.write(buffer, 0, n);
				}
			} finally {
				is.close();
			}
			return writer.toString();
		} else {       
			return "";
		}
	}
}
