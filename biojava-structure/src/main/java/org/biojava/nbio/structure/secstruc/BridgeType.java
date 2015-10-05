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
package org.biojava.nbio.structure.secstruc;

import java.io.Serializable;

public enum BridgeType implements Serializable{
	
	parallel("parallel",'p'),
	antiparallel("antiparallel",'a');
	
	public final Character type;
	public final String name;


	private BridgeType(String name,Character stype){
		this.name = name;
		this.type = stype;
	}

	public static BridgeType fromCharacter(Character stype){

		for ( BridgeType c : BridgeType.values()){
			if ( c.type.equals(stype)){
				return c;
			}
		}
		return null;
	}

	@Override
	public String toString(){
		return type.toString();
	}

}
