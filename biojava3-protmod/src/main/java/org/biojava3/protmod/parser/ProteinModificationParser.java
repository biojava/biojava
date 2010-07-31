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
 * Created on Jun 6, 2010
 * Author: Jianjiong Gao 
 *
 */

package org.biojava3.protmod.parser;

import java.util.List;
import java.util.Set;

import org.biojava.bio.structure.PDBResidueNumber;
import org.biojava.bio.structure.Structure;

import org.biojava3.protmod.ModifiedCompound;
import org.biojava3.protmod.PDBAtom;
import org.biojava3.protmod.ProteinModification;

/**
 * Identify protein modifications in a 3-d structure.
 * 
 * @author Jianjiong Gao
 * @since 3.0
 */
public interface ProteinModificationParser {
	/**
	 * Parse modifications in a structure.
	 * @param structure query {@link Structure}.
	 * @param potentialModifications query {@link ProteinModification}s.
	 */
	public void parse(Structure structure, Set<ProteinModification> potentialModifications);
	
	/**
	 * 
	 * @return a list of identified {@link ModifiedCompound}s from
	 *  the last parse result.
	 */
	public List<ModifiedCompound> getIdentifiedModifiedCompound();
	
	/**
	 * 
	 * @param recordUnidentifiableAtomLinkages true if choosing to record unidentifiable
	 *  atoms; false, otherwise.
	 * @see #getUnidentifiableModifiedResidues
	 * @see #getUnidentifiableAtomLinkages
	 */
	public void setRecordUnidentifiableCompounds(boolean recordUnidentifiableModifiedCompounds);
	
	/**
	 * 
	 * @return a list of modified residues that were not covered by
	 *  the identified ModifiedCompounds from the last parse 
	 *  result.
	 *  @see #setRecordUnidentifiableCompounds
	 *  @see #getIdentifiedModifiedCompound
	 */
	public List<PDBResidueNumber> getUnidentifiableModifiedResidues();
	
	/**
	 * 
	 * @return a list of atom pairs, which represent the 
	 *  atom bonds that were not covered by the identified 
	 *  {@link ModifiedCompound}s from the last parse result.
	 *  Each element of the list is a array containing two atoms.
	 * @see PDBAtom
	 * @see #setRecordUnidentifiableCompounds
	 */
	public List<PDBAtom[]> getUnidentifiableAtomLinkages();
}
