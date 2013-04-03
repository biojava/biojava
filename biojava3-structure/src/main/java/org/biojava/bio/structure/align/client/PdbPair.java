package org.biojava.bio.structure.align.client;

/** A pair for structure alignment
 * 
 * @author Andreas Prlic
 * 
 * name1 is always < name2
 *
 */
public class PdbPair implements Comparable<PdbPair> {

	
	StructureName name1;
	StructureName name2;
	public PdbPair(String name1, String name2) {
		super();
		
	
		// always make sure the first 4 chars are upper case...
			
		this.name1= getCheckName(name1);
		this.name2= getCheckName(name2);
				
	}
	private StructureName getCheckName(String name) {
		
		StructureName rname = new StructureName(name);
		
		return rname;
		
		
	}
	public String getName1() {
		return name1.getName();
	}
	public void setName1(String name1) {
		this.name1 = new StructureName(name1);
	}
	public String getName2() {
		return name2.getName();
	}
	public void setName2(String name2) {
		this.name2 = new StructureName(name2);
	}
	
	public String toString() {
		return "PdbPair [name1=" + name1 + ", name2=" + name2 + "]";
	}
	
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((name1 == null) ? 0 : name1.hashCode());
		result = prime * result + ((name2 == null) ? 0 : name2.hashCode());
		return result;
	}
	
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		PdbPair other = (PdbPair) obj;
		if (name1 == null) {
			if (other.name1 != null)
				return false;
		} else if (!name1.equals(other.name1))
			return false;
		if (name2 == null) {
			if (other.name2 != null)
				return false;
		} else if (!name2.equals(other.name2))
			return false;
		return true;
	}
	
	public int compareTo(PdbPair o) {
		if ( this.equals(o))
			return 0;
		
		int c = name1.getName().compareTo(o.getName1()); 
		if ( c == 0 )
			return name2.getName().compareTo(o.getName2());
		 else
			return c;
	}
	
	public String getPDBCode1() {
		return name1.getPdbId();
	
	}
	

	public String getPDBCode2(){
		return name2.getPdbId();
		
	}
	public String getChainId1(){
		return  name1.getChainId();
		
	}
	public String getChainId2(){
		return name2.getChainId();
		
	}
	
	public PdbPair getReverse() {
		PdbPair newPair = new PdbPair(name2.getName(), name1.getName());
		return newPair;
	}
	
	
	
}
