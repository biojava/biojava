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
package org.biojava.utils.bytecode;

class ExceptionMemento {
    SimpleReference startHandled;
    SimpleReference endHandled;
    CodeClass eClass;
    SimpleReference handler;

    public ExceptionMemento(SimpleReference startHandled,
			    SimpleReference endHandled,
			    CodeClass eClass,
			    SimpleReference handler)
    {
	this.startHandled = startHandled;
	this.endHandled = endHandled;
	this.eClass = eClass;
	this.handler = handler;
    }

    public boolean isFullyResolved() {
	return startHandled.isResolved() &&
	       endHandled.isResolved() &&
	       handler.isResolved();
    }
}
