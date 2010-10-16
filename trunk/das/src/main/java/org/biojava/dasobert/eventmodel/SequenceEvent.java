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
 * Created on Nov 20, 2005
 *
 */
package org.biojava.dasobert.eventmodel;


public class SequenceEvent
extends AbstractDasEvent{

    String sequence;
    String accessionCode;
    String version;
    
    public SequenceEvent(String accessionCode, String seq, String version) {
        super();
        sequence = seq;
        this.accessionCode = accessionCode;        
    }
    
    public String getAccessionCode(){
        return accessionCode;
    }
    
    public String getSequence(){
        return sequence;
    }

	public String getVersion() {
		return version;
	}

	public void setVersion(String version) {
		this.version = version;
	}

    
    
}
