package org.biojava.bio.structure.align.client;

/** A pair for structure alignment
 * 
 * @author Andreas Prlic
 * 
 * name1 is always < name2
 *
 */
public class PdbPair implements Comparable<PdbPair> {

	
	String name1;
	String name2;
	public PdbPair(String name1, String name2) {
		super();
		
		/*if ( name1.compareTo(name2) < 0){
			String tmp = name2;
			name2= name1;
			name1 = tmp;
		}*/
		// always make sure the first 4 chars are upper case...
		
		String pdb1 = name1.substring(0,4);
		String pdb2 = name2.substring(0,4);
		
		this.name1 = pdb1.toUpperCase() + name1.substring(4,name1.length());
		this.name2 = pdb2.toUpperCase() + name2.substring(4,name2.length());
	}
	public String getName1() {
		return name1;
	}
	public void setName1(String name1) {
		this.name1 = name1;
	}
	public String getName2() {
		return name2;
	}
	public void setName2(String name2) {
		this.name2 = name2;
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
		
		int c = name1.compareTo(o.getName1()); 
		if ( c == 0 )
			return name2.compareTo(o.getName2());
		 else
			return c;
	}
	
	public String getPDBCode1() {
		return getPDBFromCode(name1);
	
	}
	
	private String getPDBFromCode(String name) {
		if ( name.length() > 4)
			return name.substring(0,4);
		return null;
	}
	public String getPDBCode2(){
		return getPDBFromCode(name2);
		
	}
	public String getChainId1(){
		return getChainFromName(name1);
		
	}
	public String getChainId2(){
		return getChainFromName(name2);
		
	}
	private String getChainFromName(String name) {

		if ( name.length() == 6)
			return name.substring(5,6);
		
		return null;
	}
	
}
