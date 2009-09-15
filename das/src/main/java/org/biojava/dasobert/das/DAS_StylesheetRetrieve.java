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
 * Created on Aug 3, 2005
 *
 */
package org.biojava.dasobert.das;

import java.io.InputStream;
import java.net.*;
import java.util.*;
import java.util.logging.Logger;

import javax.xml.parsers.ParserConfigurationException;
import javax.xml.parsers.SAXParser;
import javax.xml.parsers.SAXParserFactory;

import org.biojava.dasobert.util.HttpConnectionTools;
import org.xml.sax.InputSource;
import org.xml.sax.XMLReader;

/** this stores the stylesheet config for a DAS source.
 * 
 * @author Andreas Prlic
 *
 */
public class DAS_StylesheetRetrieve {
    static Logger logger = Logger.getLogger("org.biojava.spice");
    /**
     * 
     */
    Map[] t3DMap;
    public DAS_StylesheetRetrieve() {
        super();
       
    }
    
    /** retrieve a StyleSheet from a URL
     * The style sheet is represented as a Map[],
     *  where each Map contains the description of how to draw a features of a particular type.
     *  @param dasStylesheetRequest url to get stylesheet from
     *  @return Map[] containing the stylesheet
     *  
     *  */
    public Map[] retrieve(URL dasStylesheetRequest){
        
        logger.finest("requesting stylesheet from " + dasStylesheetRequest);
        
        InputStream inStream = null;

		try {
		    HttpURLConnection huc = HttpConnectionTools.openHttpURLConnection(dasStylesheetRequest);

		    logger.finest("got connection: "+huc.getResponseMessage());
		    //String contentEncoding = huc.getContentEncoding();
		    inStream = huc.getInputStream();
		    
		    
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
		    
			DAS_Stylesheet_Handler cont_handle = new DAS_Stylesheet_Handler() ;
			XMLReader xmlreader = saxParser.getXMLReader();
			
			xmlreader.setContentHandler(cont_handle);
			xmlreader.setErrorHandler(new org.xml.sax.helpers.DefaultHandler());
			InputSource insource = new InputSource() ;
			insource.setByteStream(inStream);
			
		
		    xmlreader.parse(insource);			
			Map[] typeMap = cont_handle.getTypeStyles();
			 
			t3DMap = cont_handle.get3DStyles();
			return typeMap;
		    
		} catch (Exception e) {
		    logger.finest(e.getMessage());
		    return null;
		}        
    }

    
    public Map[] get3DStyle(){
        return t3DMap;
    }
}
