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
 */
package org.biojava.nbio.aaproperties.xml;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import javax.xml.bind.ValidationEvent;
import javax.xml.bind.ValidationEventHandler;
import javax.xml.bind.ValidationEventLocator;

public class MyValidationEventHandler implements ValidationEventHandler{

	private final static Logger logger = LoggerFactory.getLogger(MyValidationEventHandler.class);

	@Override
	public boolean handleEvent(ValidationEvent ve) {
		if (ve.getSeverity() == ValidationEvent.FATAL_ERROR ||  ve.getSeverity() == ValidationEvent.ERROR){
			ValidationEventLocator  locator = ve.getLocator();
			//print message from valdation event
			logger.info("Message is {}", ve.getMessage());
			//output line and column number
			logger.info("Column is {} at line number {}", locator.getColumnNumber(), locator.getLineNumber());

			return false;
		}
		return true;
	}
}
