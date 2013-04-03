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
 * Created on May 18, 2010
 * Author: Andreas Prlic 
 *
 */

package org.biojava.bio.structure.align.ce;

import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringWriter;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.align.StructureAlignment;
import org.biojava.bio.structure.align.StructureAlignmentFactory;

import org.biojava.bio.structure.align.client.JFatCatClient;

import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.align.xml.AFPChainFlipper;
import org.biojava.bio.structure.align.xml.AFPChainXMLConverter;
import org.biojava.bio.structure.align.xml.AFPChainXMLParser;
import org.biojava3.core.util.PrettyXMLWriter;


import junit.framework.TestCase;


public class TestWebStartClient extends TestCase
{

	@SuppressWarnings("unused")
	public void testCPAlignment(){

		//String name1="1cdg.A";
		//String name2="1tim.A";
		String name1="1VHR.A";
		String name2="2IHB.A";
		
		try {
			//StructureAlignment algorithm = StructureAlignmentFactory.getAlgorithm(CeCPMain.algorithmName);
			for (StructureAlignment algorithm : StructureAlignmentFactory.getAllAlgorithms()){
			   // disable for now
				//align(name1,name2,algorithm);
			}
		} catch (Exception e){
			e.printStackTrace();
			fail (e.getMessage());
		}
	}

	@SuppressWarnings("unused")
	private void align(String name1, String name2, StructureAlignment algorithm) 
	throws StructureException, IOException {
		if ( algorithm.getAlgorithmName().startsWith("Smith")) {
			System.err.println("not testing SW, no need to run that on server...");
			return;
		}
			
		//System.out.println("testing " + name1 + " " + name2 + " " + algorithm.getAlgorithmName());
		AtomCache cache = new AtomCache();


		Atom[] ca1 = cache.getAtoms(name1);
		Atom[] ca2 = cache.getAtoms(name2);

		
		AFPChain afpChain = algorithm.align(ca1,ca2);
		afpChain.setName1(name1);  
		afpChain.setName2(name2);

		assertNotNull(afpChain);
		assertNotNull(afpChain.getAlgorithmName());
		assertTrue(afpChain.getAlgorithmName().equals(algorithm.getAlgorithmName()));

		String xml = AFPChainXMLConverter.toXML(afpChain,ca1,ca2);

		/// SERVER part
		String serverLocation = "http://beta.rcsb.org/pdb/rest/";
		AFPChain afpServer = JFatCatClient.getAFPChainFromServer(serverLocation,algorithm.getAlgorithmName(), name1, name2, ca1, ca2, 5000);
		assertNotNull(afpServer); 

		assertTrue("Algorithm names don't match!", afpServer.getAlgorithmName().equals(algorithm.getAlgorithmName()));
		assertTrue("Alignment blockNum < 1" , afpServer.getBlockNum() >= 1);

		String xml2 = AFPChainXMLConverter.toXML(afpServer, ca1, ca2);
		//System.err.println(" tmp disabled comparison of server and client XML, a minor rounding diff...");
		assertEquals("The server and the locally calculated XML representations don;t match!", xml,xml2);

		AFPChain afpFlip = AFPChainFlipper.flipChain(afpChain);
		String xmlFlipped = AFPChainXMLConverter.toXML(afpFlip, ca2, ca1);
		//System.out.println(xmlFlipped);
		AFPChain fromXmlFlipped = AFPChainXMLParser.fromXML(xmlFlipped, ca2, ca1);
		assertEquals("The alignment lengths don't match", afpFlip.getNrEQR(), fromXmlFlipped.getNrEQR());

		String xmlFromFlippled = AFPChainXMLConverter.toXML(fromXmlFlipped,ca2,ca1);
		assertEquals("The XML of the flipped and the recreated from that XML don't match!", xmlFlipped, xmlFromFlippled);

		AFPChain afpBackToOrig = AFPChainFlipper.flipChain(fromXmlFlipped);
		//String xml5 = AFPChainXMLConverter.toXML(afpBackToOrig, ca1, ca2);
		// ok in the double flipping there are some minor after comma mismatches.

		String xmlShortOrig = getShortXML(afpChain,ca1, ca2);
		String xmlShortFinal = getShortXML(afpBackToOrig, ca1, ca2);
		assertEquals("The 2 x flipped alignment does not match the original", xmlShortOrig,xmlShortFinal);




	}

	private String getShortXML(AFPChain afpChain, Atom[] ca1, Atom[] ca2) throws IOException {

		StringWriter result = new StringWriter();
		PrintWriter writer = new PrintWriter(result);
		PrettyXMLWriter xml = new PrettyXMLWriter(writer);
		xml.openTag("AFPChain");
		//AFPChainXMLConverter.printXMLHeader(xml, afpChain);
		int blockNum = afpChain.getBlockNum();
		for(int bk = 0; bk < blockNum; bk ++) {
			
			xml.openTag("block");
			AFPChainXMLConverter.printXMLEQRInferPositions(xml, afpChain, bk, ca1, ca2);
			xml.closeTag("block");
		}
		xml.closeTag("AFPChain");
		return result.toString();
	}
}
