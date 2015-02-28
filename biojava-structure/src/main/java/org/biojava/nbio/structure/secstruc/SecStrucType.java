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

public enum SecStrucType implements Serializable{

	coil("Coil",' '),
	helix4("alpha Helix",'H'),
	helix3("3-10 Helix",'G'),
	helix5("pi helix",'I'),
	turn("Turn",'T'),	
	bend("Bend",'S'),
	extended("Extended",'E'),
	bridge("Bridge",'B');
	;

	public final Character type;
	public final String name;


	private SecStrucType(String name,Character stype){
		this.name = name;
		this.type = stype;
	}

	public static SecStrucType fromCharacter(Character stype){

		for ( SecStrucType c : SecStrucType.values()){
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

	public boolean isHelixType() {
		if ( type.equals(helix4.type) || type.equals(helix3.type) || type.equals(helix5.type))
			return true;
		return false;
	}
	

}
