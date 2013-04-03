package org.biojava3.aaproperties.xml;

import javax.xml.bind.ValidationEvent;
import javax.xml.bind.ValidationEventHandler;
import javax.xml.bind.ValidationEventLocator;

public class MyValidationEventHandler implements ValidationEventHandler{
	@Override
	public boolean handleEvent(ValidationEvent ve) {            
		if (ve.getSeverity() == ValidationEvent.FATAL_ERROR ||  ve.getSeverity() == ValidationEvent.ERROR){
			ValidationEventLocator  locator = ve.getLocator();
			//print message from valdation event
			System.out.println("Message is " + ve.getMessage());
			//output line and column number
			System.out.println("Column is " + locator.getColumnNumber() + " at line number " + locator.getLineNumber());
			return false;
		}
		return true;
	}
}