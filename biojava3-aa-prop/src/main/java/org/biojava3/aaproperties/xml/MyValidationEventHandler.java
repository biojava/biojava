package org.biojava3.aaproperties.xml;

import javax.xml.bind.ValidationEvent;
import javax.xml.bind.ValidationEventHandler;
import javax.xml.bind.ValidationEventLocator;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

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