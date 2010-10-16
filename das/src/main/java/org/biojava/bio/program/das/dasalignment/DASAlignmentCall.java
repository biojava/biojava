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
 * Created on 15.05.2004
 * @author Andreas Prlic
 *
 */



 

package org.biojava.bio.program.das.dasalignment ;

import java.io.IOException;
import java.io.InputStream;
import java.net.HttpURLConnection;
import java.net.URL;
import java.util.zip.GZIPInputStream;

import javax.xml.parsers.ParserConfigurationException;
import javax.xml.parsers.SAXParser;
import javax.xml.parsers.SAXParserFactory;

import org.xml.sax.InputSource;
import org.xml.sax.SAXException;
import org.xml.sax.SAXNotRecognizedException;
import org.xml.sax.XMLReader;

/** takes care of the communication with a DAS Alignment service.
 *
 * @author Andreas Prlic
 * @version %I% %G%
 * @since 1.4
 */
public class DASAlignmentCall {
    
    String serverurl;

    /**
     * Constructs a DASAlignmentCall object.
     */
    public DASAlignmentCall() {
	serverurl = "" ;
    }
    
    /**
     * Constructs a DASAlignmentCall object.
     *
     * @param url  a String ...
     */
    public DASAlignmentCall(String url){
	serverurl = url;
    }
    
    /** set url of aligmnent service.
     *
     * @param s  a String specifying the serverurl value
     *
     * @see #getServerurl
     */
    public void   setServerurl(String s) { serverurl = s;     }

    /** get url of alignment service.
     *
     * @return a String representing the serverurl value
     *
     * @see #setServerurl
     */
    public String getServerurl(        ) { return serverurl;}
    
    /** connect to a DAS alignment service and retreive alignments.
     * return Alignment objects.
     * uses the serverurl specified in the constructore to create http request
    
     * @return an array of Alignment objects
     * @throws IOException ...
     */
    
    public Alignment[] getAlignments()
	throws IOException
    {
	URL dasUrl = null ;
	try {
	    dasUrl = new URL(serverurl);
	} catch (Exception e) {
	    throw new IOException("error during creation of URL " + e.getMessage());
	}
	System.out.println("connecting to "+serverurl);
	InputStream inStream = connectDASServer(dasUrl);
	

	Alignment[] ali = null;
	try{
	    ali =  parseDASResponse(inStream) ;
	} catch (Exception e) {
	    throw new IOException("error during creation of URL " + e.getMessage());
	}
	return ali;	
    }
    /** connect to a DAS alignment  service and retreive data.
     * return a biojava Alignment object.
     *
     * @param query  a String
     * @return an array of Alignment objects
     * @throws IOException ...
     */
    
    public Alignment[] getAlignments(String query)
	throws IOException
    {
	/* now connect to DAS server */
	String connstr = serverurl + query ;
	URL dasUrl = null ;
	try {
	    dasUrl = new URL(connstr);
	} catch (Exception e) {
	    throw new IOException("error during creation of URL " + e.getMessage());
	}
	System.out.println("connecting to "+connstr);
	InputStream inStream = connectDASServer(dasUrl);
	

	Alignment[] ali = null;
	try{
	    ali =  parseDASResponse(inStream) ;
	} catch (Exception e) {
	    throw new IOException("error during creation of URL " + e.getMessage());
	}
	return ali;
	
    }


    /** connect to DAS server and return result as an InputStream.
     *
     */    
    private InputStream connectDASServer(URL url) 
	throws IOException
    {
	InputStream inStream = null ;
				
	System.out.println("opening connection to "+url);
	HttpURLConnection huc = null;
	huc = (HttpURLConnection) url.openConnection();	    
	 

	//System.out.println("temporarily disabled: accepting gzip encoding ");
	// should make communication much faster!
	huc.setRequestProperty("Accept-Encoding", "gzip");
	
	System.out.println("response code " +huc.getResponseCode());
	String contentEncoding = huc.getContentEncoding();
	System.out.println("getting InputStream");
	inStream = huc.getInputStream();
	if (contentEncoding != null) {
	    if (contentEncoding.indexOf("gzip") != -1) {
		// we have gzip encoding
		inStream = new GZIPInputStream(inStream);
		System.out.println("using gzip encoding!");
	    }
	}
	System.out.println("got InputStream from  DAS Alignment server");
	System.out.println("encoding: " + contentEncoding);

	return inStream;
	
    }

    /** parse the Response of a DAS ALignment service and return a
     * biojava Alignment object.
     *
     */
    public Alignment[] parseDASResponse(InputStream inStream) 
	throws IOException, SAXException
    {
	
	
	
	SAXParserFactory spfactory =
	    SAXParserFactory.newInstance();
	
	spfactory.setValidating(true);

	SAXParser saxParser = null ;

	try{
	    saxParser =
		spfactory.newSAXParser();
	} catch (ParserConfigurationException e) {
	    e.printStackTrace();
	}
	
	XMLReader xmlreader = saxParser.getXMLReader();

	
	try {
	    xmlreader.setFeature("http://xml.org/sax/features/validation", true);
	} catch (SAXException e) {
	    System.err.println("Cannot activate validation."); 
	}
       	
	try {
	    xmlreader.setFeature("http://apache.org/xml/features/nonvalidating/load-external-dtd",true);
	} catch (SAXNotRecognizedException e){
	    e.printStackTrace();
	}
	

	DASAlignmentXMLResponseParser cont_handle = new DASAlignmentXMLResponseParser() ;
	xmlreader.setContentHandler(cont_handle);
	xmlreader.setErrorHandler(new org.xml.sax.helpers.DefaultHandler());
	InputSource insource = new InputSource() ;
	insource.setByteStream(inStream);

	//System.out.println("DASAlignmentCall parse XML response ...");
	xmlreader.parse(insource);
	//System.out.println("DASAlignmentCall parse XML response done.");

	return cont_handle.getAlignments();
	
    }
    
}
