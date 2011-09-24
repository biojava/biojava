package org.biojava.bio.structure.align.client;


import java.io.Serializable;
import java.io.StringWriter;
import java.util.regex.Matcher;

import org.biojava.bio.structure.align.util.AtomCache;


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
		
		String chainID= getChainId();
		if ( chainID != null) {
			s.append(" has chain ID: ");
			s.append(chainID);
					
		}
		
		if ( isPDPDomain())
			s.append(" is a PDP domain");
		
		return s.toString();
		
	}

	public boolean isScopName() {
		if (name.startsWith("d") && name.length() >6)		
			return true;
		return false;
	}
	
	public boolean hasChainID(){
		//return name.contains(AtomCache.CHAIN_SPLIT_SYMBOL);
		
		String id= getChainId();
		if ( id != null)
			return true;
		return false;
	}
	
	public String getChainId(){
		if (name.length() == 6){
			// name is PDB.CHAINID style (e.g. 4hhb.A)
			
		
			if ( name.substring(4,5).equals(AtomCache.CHAIN_SPLIT_SYMBOL)) {
				return name.substring(5,6);
			}
		} else  if ( name.startsWith("d")){
			

			Matcher scopMatch = AtomCache.scopIDregex.matcher(name);
			if( scopMatch.matches() ) {
				//String pdbID = scopMatch.group(1);
				String chainID = scopMatch.group(2);
				//String domainID = scopMatch.group(3);
				return chainID;
			}
			
			
		} else if ( name.startsWith(AtomCache.PDP_DOMAIN_IDENTIFIER)){
			// eg. PDP:4HHBAa
			String chainID = name.substring(8,9);
			//System.out.println("chain " + chainID + " for " + name);
			return chainID;
		}
		
		return null;
	}

	
	
	public boolean isPDPDomain(){
		return name.startsWith(AtomCache.PDP_DOMAIN_IDENTIFIER);
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
