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
 * Created on Oct 26, 2006
 * 
 */

package org.biojava.dasobert.das2;

import java.io.IOException;
import java.io.InputStream;
import java.net.MalformedURLException;
import java.net.URL;
import java.net.URLEncoder;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import javax.xml.parsers.ParserConfigurationException;
import javax.xml.parsers.SAXParser;
import javax.xml.parsers.SAXParserFactory;

import org.biojava.bio.program.das.dasalignment.DASException;
import org.biojava.dasobert.das2.io.Das2ValidationHandler;
import org.biojava.dasobert.util.HttpConnectionTools;
import org.xml.sax.InputSource;
import org.xml.sax.SAXException;
import org.xml.sax.SAXNotRecognizedException;
import org.xml.sax.XMLReader;

public class Das2ValidatorCall {

	
	static final String dasypus="http://cgi.biodas.org:8080/validate_url";
	
	public static void main(String[] args){
		
		try {
	      
			//String check = "http://www.dasregistry.org/registry/das1/sources";
			//String check = "http://das.biopackages.net/das/genome";
			//String check = "http://das.biopackages.net/das/genome/human/17/type";
			String check = "http://das.biopackages.net/das/genome/human/17/feature?segment=http://www.ncbi.nlm.nih.gov/genome/H_sapiens/B36.1/dna/chr1:overlaps=1:1000";
			URL u = new URL(check);
		
			
			Das2ValidatorCall das2Validator = new Das2ValidatorCall();
			das2Validator.validate(u);
			
			
	        
		} catch (Exception e){
			e.printStackTrace();
		}
	}
	
	
	public Das2ValidatorCall() {
		
	}
	
	public Map[] validate(URL checkMe) 
	throws MalformedURLException,
	IOException,
	SAXException,DASException{
		
		URL u = new URL(dasypus + "?url=" +  URLEncoder.encode(checkMe.toString(),"UTF-8" ));
		
		InputStream stream = HttpConnectionTools.getInputStream(u);
		
		XMLReader xmlreader = getXMLReader();
		
     
        Das2ValidationHandler contentHandler = new Das2ValidationHandler();
                     
        xmlreader.setContentHandler(contentHandler);
        xmlreader.setErrorHandler(new org.xml.sax.helpers.DefaultHandler());
                
        InputSource insource = new InputSource() ;        
        insource.setByteStream(stream);
        
        xmlreader.parse(insource);
        
        List messages = contentHandler.getMessages();
        Iterator iter = messages.iterator();
        while (iter.hasNext()){
        	Map m = (Map) iter.next();
        	
        	String text = (String)m.get("text");
        	String severity = (String)m.get("severity");
        	
        	System.out.println("severity:" + severity + " : " + text);
        }
        return (Map[]) messages.toArray(new Map[messages.size()]); 
	}
	
	private XMLReader getXMLReader() 
	throws SAXException{
		SAXParserFactory spfactory =
            SAXParserFactory.newInstance();
        
		// never do this 
      
        String vali = "false";
        boolean validate = false ;
        if ((vali != null) && ( vali.equals("true")) ) 
            validate = true ;
        spfactory.setValidating(validate);
        
        SAXParser saxParser = null ;
        
        try{
            saxParser =
                spfactory.newSAXParser();
        } catch (ParserConfigurationException e) {
            //e.printStackTrace();
         	System.err.println("got exception " + e.getMessage() );
        }
        
        XMLReader xmlreader = saxParser.getXMLReader();
        
        try {
            xmlreader.setFeature("http://xml.org/sax/features/validation", validate);
        } catch (SAXException e) {
        	System.err.println("Cannot set validation to " + validate);
            System.err.println("got exception " + e.getMessage());            
        }
        
        try {
            xmlreader.setFeature("http://apache.org/xml/features/nonvalidating/load-external-dtd",validate);
        } catch (SAXNotRecognizedException e){
            //e.printStackTrace();
            System.err.println("Cannot set load-external-dtd to " + validate); 
            System.err.println("got exception" + e.getMessage()); 
        }
        
        return xmlreader;
	}
	
}
