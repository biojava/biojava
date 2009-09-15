/*
 *                  BioJava development code
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
 * Created on Jun 13, 2007
 * 
 */

package org.biojava.dasobert;

import java.io.InputStream;
import javax.xml.parsers.SAXParser;
import javax.xml.parsers.SAXParserFactory;

import org.biojava.dasobert.das.DAS_Sequence_Handler;
import org.xml.sax.InputSource;
import org.xml.sax.XMLReader;

import junit.framework.TestCase;

public class TestSequenceParse extends TestCase {

	public void testParseSequenceResponse(){

		InputStream inStream = this.getClass().getResourceAsStream("/sequence.xml");

		assertNotNull(inStream);


		boolean errorHappend = false;
		try {

			SAXParserFactory spfactory =
				SAXParserFactory.newInstance();



			SAXParser saxParser = null ;


			saxParser =
				spfactory.newSAXParser();

		
			XMLReader xmlreader = saxParser.getXMLReader();

			xmlreader.setFeature("http://xml.org/sax/features/validation", false);
			xmlreader.setFeature("http://apache.org/xml/features/nonvalidating/load-external-dtd",false);
		       
			
			DAS_Sequence_Handler cont_handle = new DAS_Sequence_Handler() ;
			xmlreader.setContentHandler(cont_handle);
			xmlreader.setErrorHandler(new org.xml.sax.helpers.DefaultHandler());
			InputSource insource = new InputSource() ;
			insource.setByteStream(inStream);

			xmlreader.parse(insource);
			String sequence = cont_handle.get_sequence();		        
			
			assertNotNull(sequence);

			String version = cont_handle.getVersion();
			assertNotNull(version);
			
			assertEquals("the sequence did not match the expected length.",149,sequence.length());

			String shouldSequence = "MLAKATLAIVLSAASLPVLAAQCEATIESNDAMQYNLKEMVVDKSCKQFTVHLKHVGKMAKVAMGHNWVLTKEADKQGVATDGMNAGLAQDYVKAGDTRVIAHTKVIGGGESDSVTFDVSKLTPGEAYAYFCSFPGHWAMMKGTLKLSN";

			assertEquals("the sequence is not the right one." , shouldSequence,sequence);
			
		} catch (Exception ex) {
			ex.printStackTrace();
			errorHappend = true;
		}

		assertTrue(errorHappend != true);

	}
}