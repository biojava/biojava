
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
 * Created on Jan 16, 2009
 *
 */
package org.biojava.dasobert.das2.io;

import java.io.InputStream;
import java.net.MalformedURLException;
import java.net.URL;
import javax.xml.parsers.ParserConfigurationException;
import javax.xml.parsers.SAXParser;
import javax.xml.parsers.SAXParserFactory;

import org.biojava.dasobert.dasregistry.DasCoordinateSystem;
import org.biojava.dasobert.dasregistry.DasSource;
import org.biojava.dasobert.util.HttpConnectionTools;
import org.xml.sax.InputSource;
import org.xml.sax.SAXException;
import org.xml.sax.SAXNotRecognizedException;
import org.xml.sax.XMLReader;


public class DASRegistryCoordinatesReaderXML implements DASRegistryCoordinatesReader {

	Exception loggedException;
	private static String PUBLIC_REGISTRY="http://www.dasregistry.org/das1/coordinatesystem";
	private static String INTERNAL_REGISTRY="http://deskpro20727.dynamic.sanger.ac.uk:8080/dasregistryOID/das1/coordinatesystem";
	public DASRegistryCoordinatesReaderXML() {
		super();
		loggedException = null;
	}

	/**
	 * Use this method if you want to read from the public registry instance www.dasregistry.org
	 * @return
	 */
	public DasCoordinateSystem[] readRegistryDas1CoorinateSystems()  {
		URL url=null;
		try {
			url=new URL(PUBLIC_REGISTRY);
		} catch (MalformedURLException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		DasCoordinateSystem[] regCoords=this.readRegistryDas1CoordinateSystems(url);
		return regCoords;
	}
	
	public DasCoordinateSystem[] readRegistryDas1CoordinateSystems(URL url){
		DasCoordinateSystem[] coords = new DasCoordinateSystem[0];
		System.setProperty("proxySet", "true");
		System.setProperty("proxyHost", "wwwcache.sanger.ac.uk");
		System.setProperty("proxyPort", "3128");
		
	
		System.setProperty("javax.xml.parsers.DocumentBuilderFactory", "org.apache.xerces.jaxp.DocumentBuilderFactoryImpl");
		System.setProperty("javax.xml.parsers.SAXParserFactory", "org.apache.xerces.jaxp.SAXParserFactoryImpl");
		InputStream stream=null;
		try {
			System.out.println("getting url for validation"+url);
			stream = HttpConnectionTools.getInputStream(url);

			
		} catch (Exception e){
			e.printStackTrace();
			loggedException = e;
		}
		coords = readRegistryCoordinates(stream);
		return coords;
	}

	/** read a DAS2 coordinates response and return a list of coordinate systems.
	 * 
	 */
	public DasCoordinateSystem[] readRegistryCoordinates(InputStream stream)  {

		DasCoordinateSystem[] regCoords = new DasCoordinateSystem[0];

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
				loggedException = e;
			}

			String vali = System.getProperty("XMLVALIDATION");

			boolean validation = vali != null && vali.equals("true") ? true : false;

			XMLReader xmlreader = saxParser.getXMLReader();

			//XMLReader xmlreader = XMLReaderFactory.createXMLReader();
			try {
				xmlreader.setFeature("http://xml.org/sax/features/validation", validation);
			} catch (SAXException e) {
				//logger.log(Level.FINE,"Cannot set validation " + validation); 
				e.printStackTrace();
			}

			try {

				xmlreader.setFeature("http://apache.org/xml/features/nonvalidating/load-external-dtd",validation);

			} catch (SAXNotRecognizedException e){
				e.printStackTrace();
				//logger.log(Level.FINE,"Cannot set load-external-dtd "+validation); 

			}

			Das1RegistryCoordinatesHandler cont_handle = new Das1RegistryCoordinatesHandler();

			xmlreader.setContentHandler(cont_handle);
			xmlreader.setErrorHandler(new org.xml.sax.helpers.DefaultHandler());

			InputSource insource = new InputSource() ;
			insource.setByteStream(stream);


			xmlreader.parse(insource);    

			regCoords = cont_handle.getRegistryCoordinates();

		} catch (Exception e) {
			e.printStackTrace();
			loggedException = e;
		}
		return regCoords;
	}

	public Exception getLoggedException(){
		return loggedException;
	}

	


	public static void main (String[] args){
		String url = "http://das.sanger.ac.uk/registry/das1/sources/";
		DasSourceReaderImpl reader = new DasSourceReaderImpl();
		try {
			URL u = new URL(url);
			DasSource[] sources = reader.readDasSource(u);
			for (int i=0; i< sources.length;i++){
				DasSource ds = sources[i];
				System.out.println(ds.toString());
			}

		} catch (Exception e){
			e.printStackTrace();
		}

	}



}
