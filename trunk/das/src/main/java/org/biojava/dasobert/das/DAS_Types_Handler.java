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
 * Created on 19.03.2004
 * @author Andreas Prlic
 *
 */
package org.biojava.dasobert.das ;

import org.xml.sax.helpers.DefaultHandler;
import org.xml.sax.Attributes;

import java.util.ArrayList ;
import java.util.List;


/** a class to parse the reponse of a DAS - types request 
 */
public class DAS_Types_Handler extends DefaultHandler {
    List types;
    boolean dastypesPresent;
    boolean gffPresent;
    boolean segmentPresent;
	private int maxNrFeaturesOntolgy;
    
    public DAS_Types_Handler() {
	super();
	types = new ArrayList();
	dastypesPresent = false;
	gffPresent=false;
	segmentPresent=false;
    }

    public void startElement (String uri, String name, String qName, Attributes atts){
	if ( qName.equals("DASTYPES")) {
	    dastypesPresent = true;
	    
	} else if ( qName.equals("GFF")) {
	    gffPresent = true;
	    
	} else if ( qName.equals("SEGMENT")) {
	    segmentPresent = true;	
	 
	    String id = atts.getValue("id");
	    // id is optional here
	    if ( id != null ) {
		types.add(id);
	    }
	} else if ( qName.equals("TYPE")){
	    String type = atts.getValue("id");
	    // id is mandatory ...	    
	    types.add(type);
	    
	}
	
    }

    public String[] getTypes(){
	return (String[])types.toArray(new String[types.size()]);
    }

	public void setMaxFeatures(int maxNrFeaturesOntology) {
		// TODO Auto-generated method stub
		this.maxNrFeaturesOntolgy=maxNrFeaturesOntology;
	}
}

