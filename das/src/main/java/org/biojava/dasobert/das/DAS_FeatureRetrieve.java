/**
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
 * Created on 19.03.2004
 * @author Andreas Prlic
 *
 */
package org.biojava.dasobert.das;


import java.net.URL                         ;
import java.io.InputStream                  ;

import org.biojava.dasobert.util.HttpConnectionTools;
import org.xml.sax.InputSource              ;
import org.xml.sax.XMLReader                ;
import javax.xml.parsers.*                  ;
import org.xml.sax.*                        ;
import java.util.ArrayList                  ;
import java.util.List;
import java.util.Map;
import java.util.logging.*                  ;
import java.net.HttpURLConnection           ;


/**
 * A class to perform a DAS features request
 * 
 * @author Andreas Prlic
 *
 */
public class DAS_FeatureRetrieve {
    
	String version;
	
    List<Map<String,String>> features ;
    Logger logger     ;
    int comeBackLater;
    URL url;
    /**
     * @param url the URL the features should be downloaded from
     *  
     */
    public DAS_FeatureRetrieve(URL url) {
        super();
        
        logger = Logger.getLogger("org.biojava.spice");
        features = new ArrayList<Map<String,String>>() ;
        comeBackLater = -1;
        version = "";
        this.url=url;
        reload();
    }
    
    
    /** contact the DAS-feature server again. Usually
     * it is not necessary to call this again, because the constructor already does, but
     * if comeBackLater > -1 this should be called again.
     *
     */
    public void reload(){
        
        try {
            
            InputStream dasInStream = null;
            try {
                dasInStream	= open(url); 
            } catch (Exception e ){
                comeBackLater = -1;
                if (logger.isLoggable(Level.FINE)){
                	logger.log(Level.FINE,"could not open connection to " + url,e);
                }
                return ;
            }
            
            
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
            try {
                xmlreader.setFeature("http://xml.org/sax/features/validation", validation);
            } catch (SAXException e) {
            	   if (logger.isLoggable(Level.FINE)){
            		   logger.log(Level.FINE,"Cannot set validation " + validation);
            	   }
            }
            
            try {
                xmlreader.setFeature("http://apache.org/xml/features/nonvalidating/load-external-dtd",validation);
            } catch (SAXNotRecognizedException e){
                e.printStackTrace();
                if (logger.isLoggable(Level.FINE)){
                	logger.log(Level.FINE,"Cannot set load-external-dtd "+validation);
                }
                
            }
            
            DAS_Feature_Handler cont_handle = new DAS_Feature_Handler() ;
            cont_handle.setDASCommand(url.toString());
            xmlreader.setContentHandler(cont_handle);
            xmlreader.setErrorHandler(new org.xml.sax.helpers.DefaultHandler());
            InputSource insource = new InputSource() ;
            insource.setByteStream(dasInStream);
            
            try {
                xmlreader.parse(insource);	
                
                features = cont_handle.get_features();
                version  = cont_handle.getVersion();
                
                comeBackLater = cont_handle.getComBackLater();
            } 
            catch ( Exception e){
                e.printStackTrace();
                if (logger.isLoggable(Level.FINE)){
                	logger.log(Level.FINE,"error while parsing response from "+ url);
                }
                comeBackLater = -1;
                features = new ArrayList<Map<String,String>>();
            }
        }
        catch (Exception ex) {
            ex.printStackTrace();
            comeBackLater = -1;
        }
    }
    
    /** open a HttpURLConnection to a URL and return an InputStream
     * 
     * @param url
     * @return an open stream
     * @throws java.io.IOException
     * @throws java.net.ConnectException
     */
    private InputStream open(URL url)
    throws java.io.IOException, java.net.ConnectException
    {
        InputStream inStream = null;
                
        HttpURLConnection huc = HttpConnectionTools.openHttpURLConnection(url);
        
        inStream = huc.getInputStream();		
        
        return inStream;
        
    }
    
    /** returns a List of Features 
     * @return a List of Maps containing the features*/
    public List<Map<String,String>> get_features() {
      
        return features;
    }
    
    
    /** Get the version string of the reference object.
     * If it does not match the version string that is obtained from the 
     * reference server there is a version problem!
     *  
     * @return version string. (e.g. a MD5 digest of the reference sequence)
     */
    public String getVersion() {
		return version;
	}


	public void setVersion(String version) {
		this.version = version;
	}


	/** returns the comeBackLater value - if a server returned suchh - 
     * 
     * @return comeBackLater in seconds, or -1 if not provided by server 
     */
    public int getComeBackLater(){
        
        return comeBackLater;
        
    }
    
    
}
