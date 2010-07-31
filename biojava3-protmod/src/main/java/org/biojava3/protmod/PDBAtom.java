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
 * Created on Jul 29, 2010
 * Author: Jianjiong Gao 
 *
 */

package org.biojava3.protmod;

import java.io.Serializable;

import org.biojava.bio.structure.PDBResidueNumber;

/**
 * Everything that is needed to uniquely describe a atom.
 * @author Jianjiong Gao
 * @since 3.0
 */
public class PDBAtom
implements Serializable {
	private static final long serialVersionUID = -3586772436145093984L;
	
	private final PDBResidueNumber group;
	private final String atomName;
	
	public PDBAtom(final PDBResidueNumber group, final String atomName) {
		if (group==null || atomName==null) {
			throw new IllegalArgumentException("Null argument(s).");
		}
		this.group = group;
		this.atomName = atomName;
	}
	
	public PDBResidueNumber getGroup() {
		return group;
	}
	
	public String getAtomName() {
		return atomName;
	}
	
	public boolean equals(Object obj) {
		if (!(obj instanceof PDBAtom))
			return false;
		
		if (obj == this)
			return true;
		
		PDBAtom anAtom = (PDBAtom)obj;
		
		if (!anAtom.getGroup().equals(group))
			return false;
		
		if (!anAtom.getAtomName().equals(atomName))
			return false;
		
		return true;
	}
	
	public int hashCode() {
		int result = 17;
		result = result * 32 + group.hashCode();
		result = result * 32 + atomName.hashCode();
		return result;
	}
	
	public String toString() {
		return group.toString() + ", Atom:" + atomName; 
	}
}
