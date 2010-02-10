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


/** a class to parse the reponse of a DAS - types request 
 */
public class DAS_Entry_Points_Handler extends DefaultHandler {

	String version ;

	// This is to prevent problems with excessive DAS servers...
	// We had one case where one DAS server hat ~200.000 entry points.
	// the returned file was 18 M big and took the DAS server 1 minute to create.
	
	public static final int MAX_NUMBER_ENTRY_POINTS=1000;
	
	
	int entry_points_counter;
	
	public DAS_Entry_Points_Handler() {
		super();

		version = null;
		
		entry_points_counter = 0;
	}

	public void startElement (String uri, String name, String qName, Attributes atts){
		if ( qName.equals("DASEP")) {

		}  else if ( qName.equals("ENTRY_POINTS")) {

			String v = atts.getValue("version");
			version = v;	    
		} 	
	}

	/** returns a String if the server returns an entry points 
	 * @return a String containing the version
	 * */
	public String getVersion() {
		return version;
	}

}

