package org.biojava.bio.structure.secstruc;

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

	public String toString(){
		return type.toString();
	}

	public boolean isHelixType() {
		if ( type.equals(helix4.type) || type.equals(helix3.type) || type.equals(helix5.type))
			return true;
		return false;
	}
	

}
