

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
 * Created on 06.05.2004
 * @author Andreas Prlic
 *
 */

/** takes care of the communication with a DAS Structure service
 */

package org.biojava.bio.program.das.dasstructure ;

import java.io.IOException;
import java.io.InputStream;
import java.net.HttpURLConnection;
import java.net.URL;
import java.util.Calendar;
import java.util.zip.GZIPInputStream;

import javax.xml.parsers.ParserConfigurationException;
import javax.xml.parsers.SAXParser;
import javax.xml.parsers.SAXParserFactory;

import org.biojava.bio.structure.Structure;
import org.xml.sax.InputSource;
import org.xml.sax.SAXException;
import org.xml.sax.XMLReader;

/**
 *  Calls a DAS - Structure server.
 *
 * @author Andreas Prlic
 * @version %I% %G%
 * @since 1.4
 */
public class DASStructureCall {
    
    String serverurl;

    /**
     * Constructs a DASStructureCall object.
     */

    public DASStructureCall() {
	serverurl = "" ;
    }
    
    /**
     * Constructs a DASStructureCall object.
     *
     * @param url  a String ...
     */
    public DASStructureCall(String url){
	serverurl = url;
    }

    /**
     * Returns a time stamp .
     * @return a String representing the current system time 
     */
    protected String getTimeStamp(){

	Calendar cal = Calendar.getInstance() ;
	// Get the components of the time
	int hour24 = cal.get(Calendar.HOUR_OF_DAY);     // 0..23
	int min = cal.get(Calendar.MINUTE);             // 0..59
	int sec = cal.get(Calendar.SECOND);             // 0..59
	String s = "time: "+hour24+" "+min+" "+sec;
	return s ;
    }

    
    /** set url of structure service. 
     * @param s  a String specifying the serverurl value
     * @see #getServerurl
     */
    public void   setServerurl(String s) { serverurl=s;     }

    /** get url of structure service.

     * @return a Structure object
     * @throws IOException ...
     * @see #setServerurl
    */
    public String getServerurl(        ) { return serverurl;}
    

    /** connect to a DAS structure service and retreive 3D data.
	return a biojava Structure object
	*
	* @param pdb_code  a String
	* @return a Structure object
	* @throws IOException ...
	*/
    
    public Structure getStructure(String pdb_code)
	throws IOException
    {
	/* now connect to DAS server */
	String connstr = serverurl + pdb_code ;
	System.out.println("DASStructureCall: connstr" + connstr);
	URL dasUrl = null ;
	try {
	    dasUrl = new URL(connstr);
	} catch (Exception e) {
	    throw new IOException("error during creation of URL " + e.getMessage());
	    //e.printStackTrace();
	    //return null;
	}
	//System.out.println("connecting to "+connstr);
	InputStream inStream = connectDASServer(dasUrl);
	

	Structure structure = null;
	//System.out.println(getTimeStamp());
	//System.out.println("starting to parse DAS structure response");
	try{
	    structure = parseDASResponse(inStream) ;
	} catch (Exception e) {
	    e.printStackTrace() ;
	    throw new IOException("error during parsing of DAS response " + e.getMessage());
	    
	}
	//System.out.println("finished parsing DAS structure response");
	//System.out.println(getTimeStamp());
	//System.out.println(structure.toPDB());
	return structure;
	
    }

    /** connect to DAS server and return result as an InputStream */
    private InputStream connectDASServer(URL url) 
	throws IOException
    {
	
				
	System.out.println(getTimeStamp() );
	System.out.println("opening connection to DAS Structure server");
	HttpURLConnection huc = null;	
	
	//huc = (HttpURLConnection) dasUrl.openConnection();	    
	//huc = proxyUrl.openConnection();	    
	//System.out.println("opening "+url);
	
	huc = (HttpURLConnection) url.openConnection();
	System.out.println(huc);
	//System.out.println("deacactivated encoding temporarily");
	// should make communication much faster!
	huc.setRequestProperty("Accept-Encoding", "gzip");
	
	//System.out.println(huc.getResponseMessage());
	
	

	
	System.out.println("getContentEncoding");
	String contentEncoding = huc.getContentEncoding();
	System.out.println("getInputStream");
	InputStream inStream = huc.getInputStream();	

	if (contentEncoding != null) {
	    if (contentEncoding.indexOf("gzip") != -1) {
		// we have gzip encoding
		inStream = new GZIPInputStream(inStream);
		System.out.println("using gzip encoding!");
	    }
	}
	
	
	System.out.println(getTimeStamp() );
	System.out.println("got InputStream from  DAS Structure server");
	System.out.println("encoding: " + contentEncoding);
	System.out.println("code:" + huc.getResponseCode());
	//System.out.println("message:" + huc.getResponseMessage());
	//inStream = huc.getInputStream();
	
	
	return inStream;
	
    }

    
    /** parse the Response of a DAS Structure service and return a
     * biojava Structure */
    private Structure parseDASResponse(InputStream inStream) 
	throws IOException, SAXException
    {
	
	
	// System.setProperty("org.xml.sax.driver", 
	// "org.apache.crimson.parser.XMLReaderImpl");


	SAXParserFactory spfactory =  SAXParserFactory.newInstance();	
	//spfactory.setValidating(false);
	//spfactory.setNamespaceAware(false);
	
	/*
	  try {
	  spfactory.setFeature("http://apache.org/xml/features/nonvalidating/load-external-dtd",false);
	  } catch (ParserConfigurationException e){
	  e.printStackTrace();
	  }
	*/
	SAXParser saxParser = null ;
	try{
	    saxParser =
		spfactory.newSAXParser();
	} catch (ParserConfigurationException e) {
	    e.printStackTrace();
	}
	
       
	XMLReader xmlreader = saxParser.getXMLReader();
	
	//http://apache.org/xml/features/nonvalidating/load-external-dtd

	//XMLReader xmlreader = XMLReaderFactory.createXMLReader();	
	//xmlreader.setValidating(false);

	// try to deactivate validation
	/*
	try {
	    xmlreader.setFeature("http://xml.org/sax/features/validation", false);
	    xmlreader.setFeature("http://apache.org/xml/features/nonvalidating/load-external-dtd",false); 
	} catch (SAXException e) {
	    System.err.println("Cannot deactivate validation."); 
	}
	*/
	/*
	try {
	    //System.out.println("deactivating validation");
	    	    //System.out.println("done...");
	} catch (SAXException e) {
	    System.err.println("Cannot deactivate validation."); 
	}
       	
	try {
	    xmlreader.setLoadDTDGrammar(false);
	    // otherwise :java.io.FileNotFoundException: http://www.sanger.ac.uk/Users/ap3/DAS/XML/protodas.dtd
	    //	    xmlreader.setFeature("http://apache.org/xml/features/nonvalidating/load-external-dtd",false);
	} catch (SAXNotRecognizedException e){
	    e.printStackTrace();
	    //System.out.println("continuing ...");
	}
	*/
	//System.out.println("DASStructureCall setting DASStructureXMLResponseParser");

	DASStructureXMLResponseParser cont_handle = new DASStructureXMLResponseParser() ;
	xmlreader.setContentHandler(cont_handle);
	xmlreader.setErrorHandler(new org.xml.sax.helpers.DefaultHandler());
	InputSource insource = new InputSource() ;
	insource.setByteStream(inStream);

	//System.out.println("DASStructureCall parse XML response ...");
	xmlreader.parse(insource);

	/*
	BufferedReader buf = new BufferedReader (new InputStreamReader (inStream));
	String line = buf.readLine ();
	while (line != null) {
	    line = buf.readLine ();
	    System.out.println(line);
	}
	*/
	return cont_handle.get_structure();
	
    }
    

}

