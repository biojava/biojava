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
 * Created on Aug 21, 2007
 * 
 */

package org.biojava.dasobert;

import java.io.InputStream;
import java.util.List;

import javax.xml.parsers.ParserConfigurationException;
import javax.xml.parsers.SAXParser;
import javax.xml.parsers.SAXParserFactory;

import org.biojava.dasobert.das.DASInteractionXMLParser;
import org.xml.sax.InputSource;
import org.xml.sax.XMLReader;

import de.mpg.mpiinf.ag3.dasmi.model.Interaction;

import junit.framework.TestCase;

public class TestInteractionParse extends TestCase {

	
	public void testParseInteractionFile(){


		InputStream inStream = this.getClass().getResourceAsStream("/interaction.xml");
		
		assertNotNull(inStream);


		boolean errorHappend = false;
		try {


			SAXParserFactory spfactory =
				SAXParserFactory.newInstance();

			spfactory.setValidating(false);

			SAXParser saxParser = null ;

			try{
				saxParser =
					spfactory.newSAXParser();
			} catch (ParserConfigurationException e) {
				e.printStackTrace();
			}

			String vali = System.getProperty("XMLVALIDATION");

			boolean validation = false ;
			if ( vali != null )
				if ( vali.equals("true") ) 
					validation = true ;


			XMLReader xmlreader = saxParser.getXMLReader();

			//XMLReader xmlreader = XMLReaderFactory.createXMLReader();

			xmlreader.setFeature("http://xml.org/sax/features/validation", validation);	            	            	            
			xmlreader.setFeature("http://apache.org/xml/features/nonvalidating/load-external-dtd",validation);


			DASInteractionXMLParser cont_handle = new DASInteractionXMLParser();
			
			xmlreader.setContentHandler(cont_handle);
			xmlreader.setErrorHandler(new org.xml.sax.helpers.DefaultHandler());
			InputSource insource = new InputSource() ;
			insource.setByteStream(inStream);


			xmlreader.parse(insource);			
			Interaction[] interactions = cont_handle.getInteractions();
			assertEquals("did not find the right number of interactions.", 1,interactions.length);
			
			Interaction i1 = interactions[0];
			List partners = i1.getParticipants();
			assertEquals("did not find the right number of interaction partners.", 2,partners.size());
			

					


		}
		catch (Exception ex) {
			ex.printStackTrace();
			errorHappend = true;
		}
		
		assertTrue(errorHappend != true);



	
	}
}
