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

package org.biojava.dasobert.das2.io;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.xml.sax.Attributes;
import org.xml.sax.helpers.DefaultHandler;

public class Das2ValidationHandler extends DefaultHandler{
	
	List messages ;
	
	public void startDocument(){
		messages = new ArrayList();

	}

	public List getMessages(){
		return messages;
	}
	
	private void addMessage(String uri, String name, String qName, Attributes atts){
		String text = atts.getValue("text");
		String severity = atts.getValue("severity");
		
		Map m = new HashMap();
		m.put("text", text);
		m.put("severity", severity);
		
		messages.add(m);
		
	}
	
	public void startElement (String uri, String name, String qName, Attributes atts){
		if (qName.equals("MESSAGE")){
			addMessage(uri,name,qName,atts);
		}
	}
}
