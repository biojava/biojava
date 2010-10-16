

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
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import javax.xml.parsers.ParserConfigurationException;
import javax.xml.parsers.SAXParser;
import javax.xml.parsers.SAXParserFactory;

import org.biojava.dasobert.das.DAS_Feature_Handler;
import org.xml.sax.InputSource;
import org.xml.sax.XMLReader;

import junit.framework.TestCase;
public class TestUCSCFeatures extends TestCase {


	protected void setUp(){


	}

	public void testParseFeaturesResponse(){

		InputStream inStream = this.getClass().getResourceAsStream("/ucscfeatures.xml");

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


			DAS_Feature_Handler cont_handle = new DAS_Feature_Handler() ;
			cont_handle.setDASCommand("");
			xmlreader.setContentHandler(cont_handle);
			xmlreader.setErrorHandler(new org.xml.sax.helpers.DefaultHandler());
			InputSource insource = new InputSource() ;
			insource.setByteStream(inStream);


			xmlreader.parse(insource);			
			List features = cont_handle.get_features();


			assertTrue("did not find the right number of features.", 35 < features.size());

			Iterator iter = features.iterator();
			while (iter.hasNext()){
				Map f = (Map) iter.next();

				String type = (String)f.get("TYPE");

				assertNotNull(type);

			}


		}
		catch (Exception ex) {
			ex.printStackTrace();
			errorHappend = true;
		}

		assertTrue(errorHappend != true);



	}

}


