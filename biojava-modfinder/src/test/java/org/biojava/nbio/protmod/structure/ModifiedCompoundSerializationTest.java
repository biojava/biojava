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

package org.biojava.nbio.protmod.structure;

import junit.framework.TestCase;
import org.biojava.nbio.protmod.ModificationCategory;
import org.biojava.nbio.protmod.ProteinModificationRegistry;
import org.biojava.nbio.protmod.io.ModifiedCompoundXMLConverter;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Structure;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;

public class ModifiedCompoundSerializationTest extends TestCase {
	
	private static final Logger logger = LoggerFactory.getLogger(ModifiedCompoundSerializationTest.class);

	boolean allOK = false;

	private String[][] strucs;

	@Override
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
				logger.error("Failed after running {} serializations at PDB ID: {}", count, name[0], e);
				fail(e.getMessage());
			}
		}
	}

	public void test1CAD(){
		String pdbId = "1CAD";
		testXMLSerialization(pdbId);
	}

	public void test1a4w(){
		String pdbId = "1a4w";
		List<ModifiedCompound> all = testXMLSerialization(pdbId);

		try {
		mainLoop:
		for ( ModifiedCompound mc : all){
			Set<StructureGroup> groups = mc.getGroups();
			
			for (StructureGroup g: groups){
				if (! g.getChainId().equals("H"))
					continue mainLoop;
			}
			ModificationCategory cat = mc.getModification().getCategory();
			
			// all modifications on chain H should be crosslink 2
			
//			if (  groups.size() != 2  ) {
//				logger.info(ModifiedCompoundXMLConverter.toXML(mc));
//				logger.info(cat);
//				logger.info(mc);
//			
//				logger.info(mc.getAtomLinkages());
//				
//				for (StructureGroup structureGroup : groups) {
//					logger.info(structureGroup);
//				}
//			}
//			assertEquals("Not the right number of groups! should be 2, but got " + groups.size() + " in: " + ModifiedCompoundXMLConverter.toXML(mc),2,groups.size());
//			
			if (!cat.equals(ModificationCategory.CROSS_LINK_2)) {
				logger.info(ModifiedCompoundXMLConverter.toXML(mc));
				logger.info(cat.toString());
				logger.info(mc.toString());
			}
			
			assertEquals(ModificationCategory.CROSS_LINK_2, cat);
		}
		} catch (Exception e){
			logger.error("Exception: ", e);
			fail(e.getMessage());
		}

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

	public List<ModifiedCompound> testXMLSerialization(String pdbId){
		String xml = null;
		ModifiedCompound currentMC = null;
		List<ModifiedCompound> all = new ArrayList<ModifiedCompound>();		
		try {

			Structure struc = TmpAtomCache.cache.getStructure(pdbId);

			ProteinModificationIdentifier parser = new ProteinModificationIdentifier();

			for (Chain c : struc.getChains()) {

				parser.identify(c, ProteinModificationRegistry.allModifications());
				Set<ModifiedCompound> mcs = parser.getIdentifiedModifiedCompound();

				for (ModifiedCompound mc : mcs){
					currentMC = mc;
					xml =  doXMLSerialization(mc) ;
					//logger.info( pdbId + " got XML: " + String.format("%n") + xml);
					ModifiedCompound newMC = getModifiedCompoundFromXML(xml);
					String xml2 = doXMLSerialization(newMC);
					assertEquals(xml,xml2);
					//logger.info(xml2);
					//assertEquals("The two objects are not equal before and after XML serialization" , mc, newMC);
					//logger.info(mc.getDescription());
					//logger.info(newMC.getDescription());
					all.add(mc);
				}
			}
		} catch (Exception e){
			logger.error("Error when serializing {}", pdbId);
			logger.error(currentMC.getDescription());
			logger.error(xml, e);
			fail(e.getMessage());
		}
		xml = null;
		currentMC =null;
		return all;
	}

	private ModifiedCompound getModifiedCompoundFromXML(String xml) {
		return ModifiedCompoundXMLConverter.fromXML(xml);

	}

	private String doXMLSerialization(ModifiedCompound mc) throws IOException{
		return ModifiedCompoundXMLConverter.toXML(mc);
	}

	public void testFlatFileParsing(){
		InputStream inStream = this.getClass().getResourceAsStream("/org/biojava/nbio/protmod/parser/modifiedCompound.xml");
		assertNotNull(inStream);
		try {
			String xml = convertStreamToString(inStream);
			//logger.info(xml);
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
