package org.biojava.bio.structure.align.client;


import java.io.Serializable;
import java.io.StringWriter;


/** A utility class that makes working with names of strucutres, domains and ranges easier
 * 
 * @param name the name. e.g. 4hhb, 4hhb.A, d4hhba_, etc.
 */
public class StructureName implements Comparable<StructureName>, Serializable{

	/**
	 * 
	 */
	private static final long serialVersionUID = 4021229518711762954L;
	String name;

	
	public StructureName(String name){
		if ( name.length() <  4)
			throw new IllegalArgumentException("This is not a valid StructureName:" + name);
		this.name = name;
	}

	/** PDB IDs are always returned as upper case
	 * 
	 * @return upper case PDB ID
	 */
	public String getPdbId(){
		if ( isScopName() ) {
			return name.substring(1,5).toUpperCase();
		}
		else  {
			// all other names start with PDB id
			return name.substring(0,4).toUpperCase();
		}

	}
	
	public String getName(){
		return name;
	}
	
	public String toString(){
		StringWriter s = new StringWriter();
		
		s.append(name);
		
		s.append(" PDB ID: ");
		s.append(getPdbId());
		
		if ( isScopName()) {
			s.append(" is a SCOP name");
		}
		
		if ( hasChainID()) {
			s.append(" has chain ID: ");
			s.append(getChainId());
					
		}
		
		return s.toString();
		
	}

	public boolean isScopName() {
		if (name.startsWith("d") && name.length() >6)		
			return true;
		return false;
	}
	public boolean hasChainID(){
		return name.contains(".");
	}

	public String getChainId() {
		// looks like PDB.chainID syntax
		String[] spl = name.split("\\.");
		String chain = null;

		if ( spl.length == 2) {

			chain = spl[1];

			if ( chain.length()>1)
				chain = chain.substring(0,1);

		} else if ( isScopName()){
			// chain ID is the 
			chain = name.substring(5,6);
		}

		return chain;
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((name == null) ? 0 : name.hashCode());
		return result;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		StructureName other = (StructureName) obj;
		if (name == null) {
			if (other.name != null)
				return false;
		} else if (!name.equals(other.name))
			return false;
		return true;
	}

	public int compareTo(StructureName o) {
		
		return name.compareTo(o.getName());
		
	}


	
	
}
