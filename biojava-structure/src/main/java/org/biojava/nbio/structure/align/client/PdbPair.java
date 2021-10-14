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

import java.util.Objects;

import org.biojava.nbio.structure.PdbId;
import org.biojava.nbio.structure.StructureException;

/**
 * A pair for structure alignment.
 * a pair is considered equal to another pair if their two respective tuple poles are equal either in their original or reversed order.
 * i.e. both <code>new PdbPair("1abc", "2abc").equals(new PdbPair("1abc", "2abc"))</code> and 
 * <code>new PdbPair("1abc", "2abc").equals(new PdbPair("2abc", "1abc"))</code> are <code>true</code>.
 * @author Andreas Prlic
 *
 */
public class PdbPair implements Comparable<PdbPair> {

	private StructureName name1;
	private StructureName name2;

	public PdbPair(String name1, String name2) {
		this(new StructureName(Objects.requireNonNull(name1)),
				new StructureName(Objects.requireNonNull(name1)));
	}

	public PdbPair(StructureName name1, StructureName name2) {
		this.name1 = Objects.requireNonNull(name1);
		this.name2 = Objects.requireNonNull(name2);
	}

	public String getName1() {
		return name1.getIdentifier();
	}
	public void setName1(String name1) {
		this.name1 = new StructureName(Objects.requireNonNull(name1));
	}
	public String getName2() {
		return name2.getIdentifier();
	}
	public void setName2(String name2) {
		this.name2 = new StructureName(Objects.requireNonNull(name2));
	}

	@Override
	public String toString() {
		return "PdbPair [name1=" + name1 + ", name2=" + name2 + "]";
	}

	@Override
	public int hashCode() {
		return Objects.hashCode(name1) + Objects.hashCode(name2);
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
		return (this.name1.equals(other.name1) && this.name2.equals(other.name2)) ||
				(this.name1.equals(other.name2) && this.name2.equals(other.name1));
	}

	@Override
	public int compareTo(PdbPair o) {
		//make sure they are not just reverse.
		if ( this.equals(o))
			return 0;

		// Use StructureName's compareTo method
		int c = name1.compareTo(o.name1);
		if ( c != 0 )
			return c;
		return name2.compareTo(o.name2);
	}

	/**
	 * @deprecated use {@link #getPDBCode1()} instead
	 * @return
	 * @throws StructureException
	 */
	public String getPDBCode1() throws StructureException {
		PdbId pdbId = name1.getPdbId();
		return pdbId == null? null: pdbId.getId();
	}

	/**
	 * @deprecated use {@link #getPDBCode2()} instead
	 * @return
	 * @throws StructureException
	 */
	@Deprecated
	public String getPDBCode2() throws StructureException{
		PdbId pdbId = name2.getPdbId();
		return pdbId == null? null: pdbId.getId();
	}
	
	/**
	 * @since 6.0.0
	 * @return
	 * @throws StructureException
	 */
	public PdbId getPdbId1() throws StructureException{
		return name1.getPdbId();
	}

	/**
	 * @since 6.0.0
	 * @return
	 * @throws StructureException
	 */
	public PdbId getPdbId2() throws StructureException{
		return name2.getPdbId();
	}
	
	public String getChainId1(){
		return  name1.getChainId();
	}

	public String getChainId2(){
		return name2.getChainId();
	}

	public PdbPair getReverse() {
		return new PdbPair(name2, name1);
	}
}
