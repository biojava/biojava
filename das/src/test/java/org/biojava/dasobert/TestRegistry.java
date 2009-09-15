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

import java.net.MalformedURLException;
import java.net.URL;

import org.biojava.dasobert.das2.io.DasSourceReaderImpl;
import org.biojava.dasobert.dasregistry.DasSource;

import junit.framework.TestCase;

public class TestRegistry extends TestCase {

	public static final String REGISTRY_LOCATION =  "http://www.dasregistry.org/das1/sources";
	//public static final String REGISTRY_LOCATION =  "http://localhost:8088/dasregistry/das1/sources";
	protected void setUp(){
		// if you are behind a proxy, please change the following lines
		// for your proxy setup ...
		//System.setProperty("proxySet", "true");
		//System.setProperty("proxyHost", "wwwcache.sanger.ac.uk");
		//System.setProperty("proxyPort", "3128");


//		make sure we use the Xerces XML parser..
		System.setProperty("javax.xml.parsers.DocumentBuilderFactory",
		"org.apache.xerces.jaxp.DocumentBuilderFactoryImpl");
		System.setProperty("javax.xml.parsers.SAXParserFactory",
		"org.apache.xerces.jaxp.SAXParserFactoryImpl");


	}

	
	public void testPublicDASRegistry(){

		DasSourceReaderImpl reader = new DasSourceReaderImpl();

		
		URL url = null ;
		try {
			url = new URL(REGISTRY_LOCATION);
		} catch (MalformedURLException e) {
			
		}
		assertNotNull(url);
		
		DasSource[] sources = reader.readDasSource(url);

		assertTrue(sources.length > 200);





	}

}
