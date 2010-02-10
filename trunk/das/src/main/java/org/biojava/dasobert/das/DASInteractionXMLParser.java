/*
 * 
 * This is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2.1 of
 * the License, or (at your option) any later version.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this software; if not, write to the Free
 * Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
 * 02110-1301 USA, or see the FSF site: http://www.fsf.org.
 * 
 */
package org.biojava.dasobert.das;

import org.xml.sax.helpers.DefaultHandler;
import org.xml.sax.Attributes;

import java.util.ArrayList;
import java.util.HashMap;

import de.mpg.mpiinf.ag3.dasmi.model.Interaction;
import de.mpg.mpiinf.ag3.dasmi.model.Interactor;
import de.mpg.mpiinf.ag3.dasmi.model.Detail;
import de.mpg.mpiinf.ag3.dasmi.model.Range;
import de.mpg.mpiinf.ag3.dasmi.model.Sequence;
import de.mpg.mpiinf.ag3.dasmi.model.Participant;

/**
 * Class parsing the interaction DAS XML response. Converts the XML content into
 * java objects based on the DASMI model
 * 
 * @author Hagen Blankenburg, Max Planck Institute for Informatics
 * 
 */
public class DASInteractionXMLParser extends DefaultHandler {

	private ArrayList <Interaction>interactions = null;
	private HashMap <String,Interactor>interactors = null;
	private Interaction interaction = null;
	private Interactor interactor = null;
	private Detail detail = null;
	private Range range = null;
	private Sequence sequence = null;
	private Participant participant = null;
	private String elementContent = null;
	private String detailParent = null;
	
	private boolean unknownSegment = false;

	
	/**
	 * Constructor. Performs some basic initialisations
	 * 
	 */
	public DASInteractionXMLParser() {
		super();
		interactors = new HashMap<String,Interactor>();
		interactions = new ArrayList<Interaction>();
		//System.out.println("constructed DASInteractionXMLParser");
	}

	/**
	 * Returns all parsed interactions.
	 * 
	 * @return An array of interaction objects
	 */
	public Interaction[] getInteractions() {
		return (Interaction[]) interactions.toArray(new Interaction[interactions.size()]);
	}

	
	/**
	 * Returns the interaction at a specific position within all interactions
	 * 
	 * @param position the psotion in the interactiolist
	 * @return the interaction
	 */
	public Interaction getInteraction(int position) {
		Interaction interaction = (Interaction) interactions.get(position);
		return interaction;
	}
	

	/**
	 * Called whenever an opening element is read from the input. Distributes
	 * the tasks to the proper private handling methods
	 * 
	 * @param uri
	 * @param name
	 * @param qName
	 * @param atts
	 */
	public void startElement(String uri, String name, String qName,	Attributes atts) {
		//System.out.println("<" + qName + ">");
		if (qName.equalsIgnoreCase("INTERACTION")) {
			interactionStart(atts);
		} else if (qName.equalsIgnoreCase("INTERACTOR")) {
			interactorStart(atts);
		} else if (qName.equalsIgnoreCase("PARTICIPANT")) {
			participantStart(atts);
		} else if (qName.equalsIgnoreCase("DETAIL")) {
			detailStart(atts);
		} else if (qName.equalsIgnoreCase("RANGE")) {
			rangeStart(atts);
		} else if (qName.equalsIgnoreCase("SEQUENCE")) {
			sequenceStart(atts);
		} else if (qName.equalsIgnoreCase("UNKNOWNSEGMENT")) {
			
		}
	}

	/**
	 * Called whenever an enclosing element is read from the input. Distributes
	 * the tasks to the proper private handling methods
	 * 
	 */
	public void endElement(String uri, String name, String qName) {
		//System.out.println("</" + qName + ">");
		if (qName.equalsIgnoreCase("INTERACTION")) {
			interactionEnd();
		} else if (qName.equalsIgnoreCase("INTERACTOR")) {
			interactorEnd();
		} else if (qName.equalsIgnoreCase("PARTICIPANT")) {
			participantEnd();
		} else if (qName.equalsIgnoreCase("DETAIL")) {
			detailEnd();
		} else if (qName.equalsIgnoreCase("RANGE")) {
			rangeEnd();
		} else if (qName.equalsIgnoreCase("SEQUENCE")) {
			sequenceEnd();
		} else if (qName.equalsIgnoreCase("UNKNOWNSEGMENT")) {
			unknownSegment = true;
		} 

	}

	/**
	 * Interaction start. Reads the interaction parameters from the input stream 
	 * @param atts
	 */
	private void interactionStart(Attributes atts) {
		interaction = new Interaction();
		interaction.setDbAccessionId(atts.getValue("dbAccessionId"));
		interaction.setName(atts.getValue("name"));
		interaction.setDbSource(atts.getValue("dbSource"));
		try {
			interaction.setDbSourceCvId(atts.getValue("dbSourceCvId"));
		} catch (Exception e) {
		}
		try {
			interaction.setDbVersion(atts.getValue("dbVersion"));
		} catch (Exception e) {
		}
		detailParent = "interaction";
	}

	/**
	 * Interaction end. Adds the interaction to all interactions. 
	 * 
	 */
	private void interactionEnd() {
		interactions.add(interaction);
		interaction = null;
	}

	/**
	 * Interactor start. Reads the interactor attributes from th input stream.
	 * @param atts
	 */
	private void interactorStart(Attributes atts) {
		interactor = new Interactor();
		if (atts.getValue("intId") == null){
			interactor.setId(atts.getValue("id"));
			//System.out.println("setting id by getting id");
		}else{
			interactor.setId(atts.getValue("intId"));
			//System.out.println("setting id by getting intId");
		}
		interactor.setName(atts.getValue("shortLabel"));
		interactor.setDbAccessionId(atts.getValue("dbAccessionId"));
		interactor.setDbCoordSys(atts.getValue("dbCoordSys"));
		try{
			interactor.setDbSource(atts.getValue("dbSource"));
		} catch (Exception e) {}
		interactor.setDbSourceCvId(atts.getValue("dbSourceCvId"));
		try {
			interactor.setDbVersion(atts.getValue("dbVersion"));
		} catch (Exception e) {}
		detailParent = "interactor";
	}

	/**
	 * Interactor end. Binds the sequence, if present and puts the 
	 * interactor into a hash for later assignment
	 * 
	 */
	private void interactorEnd() {
		if (sequence != null) {
			interactor.setSequence(sequence);
		}
		interactors.put(interactor.getId(), interactor);
		interactor = null;
	}

	/**
	 * Participant start. Reads the participant id.
	 * @param atts
	 */
	private void participantStart(Attributes atts) {
		participant = new Participant();
		//System.out.println("intId= "+atts.getValue("intId"));
		//System.out.println("id= "+atts.getValue("id"));
		if (atts.getValue("intId") == null){
			participant.setId(atts.getValue("id"));
			//System.out.println("participant set id from id");
			
		}else{
			participant.setId(atts.getValue("intId"));
			//System.out.println("participant set id using intId");
		}
		detailParent = "participant";
	}

	/**
	 * Participant end. Assigns the proper interactor and adds the
	 *  participant to the associated interaction
	 * 
	 */
	private void participantEnd() {
		Interactor inter = (Interactor) interactors.get(participant.getId());
		participant.setInteractor(inter);
		interaction.addParticipant(participant);
		//System.out.println(participant);
		participant = null;
	}

	/**
	 * Sequence start. Reads the sequence attributes from the input stream
	 * @param atts
	 */
	private void sequenceStart(Attributes atts) {
		sequence = new Sequence();
		try {
			sequence.setStart(Integer.valueOf(atts.getValue("start")).intValue());
		} catch (Exception e) {
		}
		try {
			sequence.setEnd(Integer.valueOf(atts.getValue("end")).intValue());
		} catch (Exception e) {
		}
	}

	/**
	 * Sequence end. Adds the sequence to the object
	 * 
	 */
	private void sequenceEnd() {
		sequence.setSequence(elementContent);
		elementContent = "";
	}

	/**
	 * Range start. Reads the range attributs form the input stream
	 * @param atts
	 */
	private void rangeStart(Attributes atts) {
		range = new Range();
		range.setStart(Integer.valueOf(atts.getValue("start")).intValue());
		range.setEnd(Integer.valueOf(atts.getValue("end")).intValue());
		try {
			range.setStartStatus(atts.getValue("startStatus"));
		} catch (Exception e) {
		}
		try {
			range.setStartStatusCvId(atts.getValue("startStatusCvId"));
		} catch (Exception e) {
		}
		try {
			range.setEndStatus(atts.getValue("endStatus"));
		} catch (Exception e) {
		}
		try {
			range.setEndStatusCvId(atts.getValue("endStatusCvId"));
		} catch (Exception e) {
		}

	}

	/**
	 * Range end. Does nothing.
	 *
	 */
	private void rangeEnd() {
	}

	/**
	 * Detail start. Sets the detail attributes.
	 * @param atts
	 */
	private void detailStart(Attributes atts) {
		detail = new Detail();
		detail.setProperty(atts.getValue("property"));
		try {
			detail.setPropertyCvId(atts.getValue("propertyCvId"));
		} catch (Exception e) {
		}
		detail.setValue(atts.getValue("value"));
		try {
			detail.setValueCvId(atts.getValue("valueCvId"));
		} catch (Exception e) {
		}
	}

	/**
	 * Detail end. Adds a range (if present) and the detail to the current parent element
	 *
	 */
	private void detailEnd() {
		if (range != null) {
			detail.setRange(range);
			range = null;
		}
		if (detailParent.equals("interactor")) {
			interactor.addDetail(detail);
		} else if (detailParent.equals("interaction")) {
			interaction.addDetail(detail);
		} else if (detailParent.equals("participant")) {
			participant.addDetail(detail);
		}
		detail = null;
	}

	/**
	 * Called when the beginning of the document is reached
	 */
	public void startDocument() {
		//System.out.println("start document");
	}

	/**
	 * Called when the end of the docuemnt is reached 
	 */
	public void endDocument() {
		//System.out.println("end document");
	}

	/**
	 * Called whenever element text is read
	 */
	public void characters(char ch[], int start, int length) {
		elementContent = "";
		for (int i = start; i < start + length; i++) {
			elementContent += ch[i];
		}
		//System.out.println("content"+elementContent);
	}

}
