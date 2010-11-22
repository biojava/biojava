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
 * created at Sep 7, 2007
 */
package org.biojava.bio.structure.server;

import org.biojava.bio.structure.Structure;

public class StructureEventImpl 
implements StructureEvent{

	Structure s;
	String pdbCode;
	
	public StructureEventImpl(Structure s){
		this.s = s;
	}
	
	public Structure getStructure() {
	
		return s;
	}

	public String getPDBCode() {
		return pdbCode;
	}

	public void setPDBCode(String pdbCode) {
		this.pdbCode = pdbCode;
		
	}
	
	

}
