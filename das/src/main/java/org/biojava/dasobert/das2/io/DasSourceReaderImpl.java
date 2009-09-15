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
 * Created on Feb 24, 2006
 *
 */
package org.biojava.dasobert.das2.io;

import java.io.IOException;
import java.io.InputStream;
import java.lang.reflect.Method;
import java.net.ConnectException;
import java.net.HttpURLConnection;
import java.net.URL;
import javax.xml.parsers.ParserConfigurationException;
import javax.xml.parsers.SAXParser;
import javax.xml.parsers.SAXParserFactory;
import org.biojava.dasobert.dasregistry.DasSource;
import org.biojava.dasobert.util.DasobertDefaults;
import org.biojava.dasobert.util.HttpConnectionTools;
import org.xml.sax.InputSource;
import org.xml.sax.SAXException;
import org.xml.sax.SAXNotRecognizedException;
import org.xml.sax.XMLReader;


public class DasSourceReaderImpl implements DasSourceReader {

	Exception loggedException;

	public DasSourceReaderImpl() {
		super();
		loggedException = null;
	}


	public DasSource[] readDasSource(URL url){
		DasSource[] sources = new DasSource[0];

		try {
			InputStream stream = HttpConnectionTools.getInputStream(url);

			sources = readDasSource(stream);
		} catch (Exception e){
			e.printStackTrace();
			loggedException = e;
		}
		return sources;
	}

	/** read a DAS2 sources response and return a list of DAS sources.
	 * 
	 */
	public DasSource[] readDasSource(InputStream stream)  {

		DasSource[] sources = new DasSource[0];

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

			boolean validation = false ;
			if ( vali != null )
				if ( vali.equals("true") ) 
					validation = true ;


			XMLReader xmlreader = saxParser.getXMLReader();

			//XMLReader xmlreader = XMLReaderFactory.createXMLReader();
			try {
				xmlreader.setFeature("http://xml.org/sax/features/validation", validation);
			} catch (SAXException e) {
				//logger.log(Level.FINE,"Cannot set validation " + validation); 
			}

			try {

				xmlreader.setFeature("http://apache.org/xml/features/nonvalidating/load-external-dtd",validation);

			} catch (SAXNotRecognizedException e){
				e.printStackTrace();
				//logger.log(Level.FINE,"Cannot set load-external-dtd "+validation); 

			}

			Das2SourceHandler cont_handle = new Das2SourceHandler() ;

			xmlreader.setContentHandler(cont_handle);
			xmlreader.setErrorHandler(new org.xml.sax.helpers.DefaultHandler());

			InputSource insource = new InputSource() ;
			insource.setByteStream(stream);


			xmlreader.parse(insource);    

			sources = cont_handle.getSources();

		} catch (Exception e) {
			e.printStackTrace();
			loggedException = e;
		}
		return sources;
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
