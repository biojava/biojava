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
 * Created on Aug 2, 2010
 * Author: Jianjiong Gao 
 *
 */

package org.biojava3.protmod.structure;

public class StructureAtomLinkage {
	
	private final StructureAtom atom1;
	private final StructureAtom atom2;
	private final double distance;
	
	public StructureAtomLinkage(final StructureAtom atom1,
			final StructureAtom atom2, final double distance) {
		if (atom1 == null || atom2 == null)
			throw new IllegalArgumentException("Null atom(s)");
		this.atom1 = atom1;
		this.atom2 = atom2;
		this.distance = distance;
	}

	public StructureAtom getAtom1() {
		return atom1;
	}
	
	public StructureAtom getAtom2() {
		return atom2;
	}
	
	public double getDistance() {
		return distance;
	}
	
	public boolean equals(Object obj) {
		if (obj == this)
			return true;
		
		if (!(obj instanceof StructureAtomLinkage))
			return false;
		
		StructureAtomLinkage aLink = (StructureAtomLinkage) obj;
		if (aLink.atom1.equals(atom1) && aLink.atom2.equals(atom2))
			return true;
		
		if (aLink.atom1.equals(atom2) && aLink.atom2.equals(atom1))
			return true;
		
		return false;
	}
	
	public int hashCode() {
		int result = 17;
		result = result * 31 + atom1.hashCode() + atom2.hashCode();
		return result;
	}
	
	public String toString() {
		String dat =  atom1.toString() + "-" + atom2.toString() + " distance: " + String.format("%.2f",distance);
		dat = dat.replaceAll("\t"," ");
		return dat;
	}
}
