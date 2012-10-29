package org.biojava.bio.structure.secstruc;

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

	public String toString(){
		return type.toString();
	}

}
