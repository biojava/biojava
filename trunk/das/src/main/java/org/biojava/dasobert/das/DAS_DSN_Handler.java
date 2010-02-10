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
 * Created on Dec 7, 2006
 * 
 */

package org.biojava.dasobert.das;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.xml.sax.Attributes;
import org.xml.sax.helpers.DefaultHandler;

/** XML parser for a DAS - DSN response
 * 
 * @author Andreas Prlic
 *
 */
public class DAS_DSN_Handler 
	extends DefaultHandler
	{
	
	List dsnSources;
	Map currentDSN;
	StringBuffer characterData;
	
	public DAS_DSN_Handler(){
		
		dsnSources = new ArrayList();
		currentDSN = new HashMap();
		characterData = new StringBuffer();
	}
	
	public List getDsnSources(){
		return dsnSources;
	}



	public void endElement(String uri, String name, String qName)  {
		
		if (qName.equals("MAPMASTER"))
			currentDSN.put(qName,characterData.toString());
		
		if ( qName.equals("DESCRIPTION"))
			currentDSN.put(qName,characterData.toString());
		
		if ( qName.equals("DSN"))
			dsnSources.add(currentDSN);
		
	}

	public void startDocument()  {
	
		dsnSources = new ArrayList();
		
	}

	public void startElement(String uri, String name, String qName, Attributes atts) {
		
		
		if ( qName.equals("DSN") ){
			currentDSN = new HashMap();
		} else if (qName.equals("SOURCE")){
			String id = atts.getValue("id");
			if (! (id == null)){
				currentDSN.put("id",id);
			}
		}
		
		characterData = new StringBuffer();
		
	}
	
	public void characters (char ch[], int start, int length){
		for (int i = start; i < start + length; i++) {

			characterData.append(ch[i]);
		}
	}
	
	
	
	
}
