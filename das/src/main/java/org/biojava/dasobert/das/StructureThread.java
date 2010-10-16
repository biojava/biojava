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
 * Created on Nov 7, 2005
 *
 */
package org.biojava.dasobert.das;

import java.io.IOException;
import java.util.Iterator;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.List;
import java.util.ArrayList;

import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.io.DASStructureClient;
import org.biojava.bio.structure.io.PDBFileReader;
import org.biojava.dasobert.dasregistry.Das1Source;
import org.biojava.dasobert.eventmodel.StructureEvent;
import org.biojava.dasobert.eventmodel.StructureListener;


/** a thread that gets the protein structure from a das server.
 * If several servers are provided it takes the structure from the first server.
 * 
 * <pre>
 *      Das1Source dasSource = new Das1Source();
        
        dasSource.setUrl("http://das.sanger.ac.uk/das/structure/");
 * 
 * 
 *  String pdbCode = "1boi";
        
        // now we create the thread that will fetch the structure
        StructureThread thread = new StructureThread(pdbCode,dasSource);
        
        // add a structureListener that simply prints the PDB code
        StructureListener listener = new MyStructureListener();
        thread.addStructureListener(listener);

        // and now start the DAS request
        thread.start();
     }   
        class MyStructureListener 
    implements StructureListener {


    public synchronized void newStructure(StructureEvent event){
        Structure s = event.getStructure();
        System.out.println(s.toPDB()); 
        //System.out.println(s.toPDB()); 
        System.exit(0);
    }
    


    // the methods below are required by the interface but not needed here
    public void selectedChain(StructureEvent event){}
    public void newObjectRequested(String name){}
    public void noObjectFound(String accessionCode){}

    }
        
 * </pre>
 * 
 * @author Andreas Prlic
 *
 */
public class StructureThread
extends Thread{
    
    
    Das1Source[] dasSources;
    String accessionCode;
    static Logger logger = Logger.getLogger("org.biojava.spice");
    List structureListeners;
     
    public StructureThread(String accessionCode, Das1Source ds) {
	Das1Source[] dss = new Das1Source[1];
	dss[0] = ds;
	dasSources = dss;
	this.accessionCode = accessionCode;
	structureListeners = new ArrayList();
    }


    public StructureThread(String accessionCode, Das1Source[] dss) {
        super();
        dasSources = dss;
        
        // if the accessioncode has a chain, remove it ...
        int pos = accessionCode.indexOf(".");
        System.out.println("ac " +  accessionCode +" " + pos);
        if ( pos > 0)                 
            this.accessionCode = accessionCode.substring(0,pos);
        else
            this.accessionCode = accessionCode;
        
        structureListeners = new ArrayList();
    }
    
    public void addStructureListener(StructureListener li){
        if ( ! structureListeners.contains(li))
            structureListeners.add(li);
    }
    
    public void triggerNoStructure(String ac){
        
        Iterator iter = structureListeners.iterator();
        while (iter.hasNext()){
           //StructureListener li = (StructureListener) iter.next();
           //li.newStructure(event);
            StructureListener li = (StructureListener) iter.next();
            li.noObjectFound(ac);
           
        }
}
    
    public void triggerNewStructure(StructureEvent event){
        
            Iterator iter = structureListeners.iterator();
            while (iter.hasNext()){
               //StructureListener li = (StructureListener) iter.next();
               //li.newStructure(event);
                StructureListener li = (StructureListener) iter.next();
                li.newStructure(event);
            }
    }
    
    public void run() {
    	
        Structure structure = null ;
        
        for (int i=0 ; i < dasSources.length; i++){
        	
            Das1Source ds = dasSources[i];
            
            String url = ds.getUrl();
            logger.finest(url);
            
            if ( url.substring(0,7).equals("file://") ) {
                // load local PDB file
                String dir  = url.substring(7);
                structure = getLocalPDB(dir,accessionCode);
            } else {
                char lastChar = url.charAt(url.length()-1);      
                if ( ! (lastChar == '/') ) 
                    url +="/" ;
                
                String dasstructurecommand = url + "structure?model=1&query=";
                
                
                DASStructureClient dasc= new DASStructureClient(dasstructurecommand);
                logger.info("requesting structure from "+dasstructurecommand  +accessionCode);     
                try {
                    structure = dasc.getStructureById(accessionCode); 
                    
                    // return as soon as we found a structure
                    if ( structure != null ){
                    	logger.info("found structure " + structure.getPDBCode());
                        StructureEvent event = new StructureEvent(structure);
                        event.setSource(ds);
                        triggerNewStructure(event);
                        return;
                    }
                }
                catch (Exception e) {
                    logger.log(Level.WARNING,"could not retreive structure from "+dasstructurecommand ,e);
                    triggerNoStructure(accessionCode);
                }
            }
        }
        
        if ( structure != null ){
        	System.out.println("StrucutureThread got structure, but not sure from which DAS source...");
            StructureEvent event = new StructureEvent(structure);            
            triggerNewStructure(event);
        } else {
        	 triggerNoStructure(accessionCode);
        }
       
        
    }
    
    private Structure getLocalPDB(String dir,String pdbcode) {
        PDBFileReader parser = new PDBFileReader() ;
        //TODO make extensions configurable
        String[] extensions = {".ent",".pdb"};
        for (int i =0; i< extensions.length ; i++){
            String s = extensions[i];
            parser.addExtension(s);
        }
        parser.setPath(dir+ java.io.File.separator);
        Structure struc = null ;
        try {
            struc = parser.getStructureById(pdbcode);
        }  catch ( IOException e) {
            logger.log(Level.INFO,"local structure "+pdbcode+" not found, trying somewhere else");
            //e.printStackTrace();
            return null ;
        }
        return struc ;
    }
}
