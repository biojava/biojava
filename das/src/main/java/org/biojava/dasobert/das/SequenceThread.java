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
 * Created on Nov 20, 2005
 *
 */
package org.biojava.dasobert.das;

import java.io.InputStream;
import java.net.HttpURLConnection;
import java.net.URL;
import java.util.Iterator;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.xml.parsers.ParserConfigurationException;
import javax.xml.parsers.SAXParser;
import javax.xml.parsers.SAXParserFactory;

import org.biojava.dasobert.dasregistry.Das1Source;
import org.biojava.dasobert.eventmodel.SequenceEvent;
import org.biojava.dasobert.eventmodel.SequenceListener;
import org.biojava.dasobert.util.HttpConnectionTools;
import org.xml.sax.InputSource;
import org.xml.sax.SAXException;
import org.xml.sax.SAXNotRecognizedException;
import org.xml.sax.XMLReader;
import java.util.*;

/** a thread that gets the sequence from a DAS server
 * 
 * @author Andreas Prlic
 *
 */
public class SequenceThread  
extends Thread {
    
    Das1Source[] sequenceServers;
    String sp_accession;
    List seqListeners;
    String version ;
    
    static Logger logger = Logger.getLogger("org.biojava.spice");
    
     public SequenceThread(String sp_accession,Das1Source ds ) {
        super();
        Das1Source[] dss =new Das1Source[1];
	dss[0] = ds;
        this.sp_accession = sp_accession;
        this.sequenceServers =dss ;
        clearSequenceListeners();
        version = "";
    }
    public SequenceThread(String sp_accession,Das1Source[] ds ) {
        super();
        
        this.sp_accession = sp_accession;
        this.sequenceServers =ds ;
        clearSequenceListeners();
    }
    
    public void clearSequenceListeners(){
        seqListeners = new ArrayList();
    }
    
    public void addSequenceListener(SequenceListener lis){
        seqListeners.add(lis);
    }
    
    public void run() {
        getSequence();
    }
    
    public void getSequence( ){
        
        boolean gotSequence = false ;
        
        for ( int i = 0 ; i< sequenceServers.length; i++){
            
            if ( gotSequence ) break ;
            
            Das1Source ds = sequenceServers[i];
            String url = ds.getUrl() ;
            char lastChar = url.charAt(url.length()-1);      
            if ( ! (lastChar == '/') ) 
                url +="/" ;
            String dascmd = url + "sequence?segment=";
            String connstr = dascmd + sp_accession ;
            
            try {
                version = "";
                
                String seq= retrieveSequence(connstr);
                // bug in aristotle das source?
                seq.replaceAll(" ","");
                gotSequence = true ;
                // set the sequence ...
                
                triggerNewSequence(sp_accession,seq,ds,version);
                
                
                return;
            }
            catch (Exception ex) {
                ex.printStackTrace();     
                logger.warning(ex.getMessage());
                
                //triggerException(ex);
                
            }       
        }
        
        logger.log(Level.WARNING,"could not retreive UniProt sequence from any available DAS sequence server");
        
        triggerNoSequence(sp_accession);
        
    }

   
    
//    private void triggerException(Exception e){
//        Iterator iter = seqListeners.iterator();
//        while (iter.hasNext()){
//           SequenceListener li = (SequenceListener)iter.next();
//           li.exceptionOccured(e);
//        }
//    }
    
    private void triggerNewSequence(String sp_accession,String sequence,Das1Source source,String version){

        Iterator iter = seqListeners.iterator();
        while (iter.hasNext()){
           SequenceListener li = (SequenceListener)iter.next();
            //SequenceEvent event = new SequenceEvent(sequence);
           SequenceEvent event = new SequenceEvent(sp_accession,sequence,version); 
           event.setSource(source);
           li.newSequence(event);
        }
    }
    
    private void triggerNoSequence(String ac){

        Iterator iter = seqListeners.iterator();
        while (iter.hasNext()){
            SequenceListener li = (SequenceListener)iter.next();
            li.noObjectFound(ac);
        }
    
    }
    
    /** retrieve the Sequence from a DAS server.  
     * 
     * @param connstr - the DAS - request string. e.g. http://www.ebi.ac.uk/das-srv/uniprot/das/aristotle/sequence?segment=P00280
     * @return the requested Sequence
     * @throws Exception
     */
    public String retrieveSequence( String connstr) 
    throws Exception 
    {
        
        //logger.finest("trying: " + connstr) ;
        URL dasUrl = new URL(connstr);
        //DAS_httpConnector dhtp = new DAS_httpConnector() ;
        logger.info("requesting sequence from " + connstr);
        InputStream dasInStream =open(dasUrl); 
        
        
        SAXParserFactory spfactory =
            SAXParserFactory.newInstance();
        
	// never do this 
        //String vali = System.getProperty("XMLVALIDATION");
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
            logger.log(Level.FINER,"Uncaught exception", e);
        }
        
        XMLReader xmlreader = saxParser.getXMLReader();
        
        try {
            xmlreader.setFeature("http://xml.org/sax/features/validation", validate);
        } catch (SAXException e) {
            logger.finer("Cannot set validation to " + validate); 
            logger.log(Level.FINER,"Uncaught exception", e);
        }
        
        try {
            xmlreader.setFeature("http://apache.org/xml/features/nonvalidating/load-external-dtd",validate);
        } catch (SAXNotRecognizedException e){
            //e.printStackTrace();
            logger.finer("Cannot set load-external-dtd to" + validate); 
            logger.log(Level.FINER,"Uncaught exception", e);
            //System.err.println("Cannot set load-external-dtd to" + validate); 
        }
        
        
        //DAS_DNA_Handler cont_handle = new DAS_DNA_Handler() ;
        DAS_Sequence_Handler cont_handle = new DAS_Sequence_Handler() ;
        xmlreader.setContentHandler(cont_handle);
        xmlreader.setErrorHandler(new org.xml.sax.helpers.DefaultHandler());
        InputSource insource = new InputSource() ;
        insource.setByteStream(dasInStream);
        
        xmlreader.parse(insource);
        String sequence = cont_handle.get_sequence();
        version = cont_handle.getVersion();
        //logger.finest("Got sequence from DAS: " +sequence);
        
        logger.exiting(this.getClass().getName(), "retreiveSequence",  sequence);
        return sequence ;
    }
    
    private InputStream open(URL url) {
        {
            
            InputStream inStream = null;
            try{
                
             
               
                
                HttpURLConnection huc = null;
             
                huc = HttpConnectionTools.openHttpURLConnection(url);
                
                
                logger.finest(huc.getResponseMessage());
                
                inStream = huc.getInputStream();
                
                
            }
            catch ( Exception ex){
                ex.printStackTrace();
                logger.log(Level.WARNING,"exception occured", ex);
            }
            
            return inStream;
        }
        
    }
    
}
