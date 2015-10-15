package org.biojava.nbio.structure.secstruc;

/**
 * A bridge is formed by two non-overlapping stretches of three 
 * residues each (i-1,i,i+1) and (j-1,j,j+1), where i<j.
 * <p>
 * Depending on two basic patterns, a Bridge can be either of 
 * type parallel 
 * (H bonds in {(i-1,j) and (j,i+1)} OR {(j-1,i) and (i,j-1)})
 * or antiparallel 
 * (H bonds in {(i1,j) and (j,i} OR {(i-1,j+1) and (j-1,i+1)})
 * 
 * @author Andreas Prlic
 * @author Aleix Lafita
 *
 */
public enum BridgeType {
	
	parallel("parallel",'p'),
	antiparallel("antiparallel",'a');
	
	public final Character type;
	public final String name;

	private BridgeType(String name, Character stype){
		this.name = name;
		this.type = stype;
	}

	public static BridgeType fromCharacter(Character stype){

		for (BridgeType c : BridgeType.values()){
			if (c.type.equals(stype)) return c;
		}
		return null;
	}

	@Override
	public String toString(){
		return type.toString();
	}

}
