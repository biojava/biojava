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
package org.biojava.nbio.structure.align.client;

import org.biojava.nbio.structure.StructureException;

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
		this(new StructureName(name1),new StructureName(name2));
	}
	public PdbPair(StructureName name1, StructureName name2) {
		super();
		this.name1 = name1;
		this.name2 = name2;
	}

	public String getName1() {
		return name1.getIdentifier();
	}
	public void setName1(String name1) {
		this.name1 = new StructureName(name1);
	}
	public String getName2() {
		return name2.getIdentifier();
	}
	public void setName2(String name2) {
		this.name2 = new StructureName(name2);
	}

	@Override
	public String toString() {
		return "PdbPair [name1=" + name1 + ", name2=" + name2 + "]";
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((name1 == null) ? 0 : name1.hashCode());
		result = prime * result + ((name2 == null) ? 0 : name2.hashCode());
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

	@Override
	public int compareTo(PdbPair o) {
		if ( this.equals(o))
			return 0;
		// Use StructureName's compareTo method
		int c = name1.compareTo(o.name1);
		if ( c != 0 )
			return c;
		return name2.compareTo(o.name2);
	}

	public String getPDBCode1() throws StructureException {
		return name1.getPdbId();
	}
	public String getPDBCode2() throws StructureException{
		return name2.getPdbId();
	}

	public String getChainId1(){
		return  name1.getChainId();
	}
	public String getChainId2(){
		return name2.getChainId();
	}

	public PdbPair getReverse() {
		PdbPair newPair = new PdbPair(name2, name1);
		return newPair;
	}
}
