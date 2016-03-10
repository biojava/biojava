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
package org.biojava.nbio.aaproperties.xml;

import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.core.sequence.template.CompoundSet;
import org.biojava.nbio.core.sequence.template.Sequence;

import java.util.*;

public class ModifiedAminoAcidCompoundSet implements CompoundSet<AminoAcidCompound> {

	private final Map<String, AminoAcidCompound> aminoAcidCompoundCache = new HashMap<String, AminoAcidCompound>();

	public ModifiedAminoAcidCompoundSet(List<AminoAcidComposition> aaList, Map<Character, Double> aaSymbol2MolecularWeight) {
		this.aminoAcidCompoundCache.put("-", new AminoAcidCompound(null, "-", "", "", 0.0f));
		for (AminoAcidComposition aa : aaList) {
			this.aminoAcidCompoundCache.put(aa.getSymbol(),
					new AminoAcidCompound(null, aa.getSymbol(), aa.getShorName(), aa.getName(),
							aaSymbol2MolecularWeight.get(aa.getSymbol().charAt(0)).floatValue()));
		}
	}

	@Override
	public int getMaxSingleCompoundStringLength() {
		return 1;
	}

	@Override
	public boolean isCompoundStringLengthEqual() {
		return true;
	}

	@Override
	public AminoAcidCompound getCompoundForString(String string) {
		if (string.length() == 0) {
			return null;
		}
		if (string.length() > this.getMaxSingleCompoundStringLength()) {
			throw new IllegalArgumentException("String supplied (" + string + ") is too long. Max is " + getMaxSingleCompoundStringLength());
		}
		return this.aminoAcidCompoundCache.get(string);
	}

	@Override
	public String getStringForCompound(AminoAcidCompound compound) {
		return compound.toString();
	}

	@Override
	public boolean compoundsEquivalent(AminoAcidCompound compoundOne, AminoAcidCompound compoundTwo) {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public Set<AminoAcidCompound> getEquivalentCompounds(AminoAcidCompound compound) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public boolean hasCompound(AminoAcidCompound compound) {
		return aminoAcidCompoundCache.containsValue(compound);
	}

	@Override
	public List<AminoAcidCompound> getAllCompounds() {
		return new ArrayList<AminoAcidCompound>(aminoAcidCompoundCache.values());
	}

	@Override
	public boolean isComplementable() {
		return false;
	}

	@Override
	public boolean isValidSequence(Sequence<AminoAcidCompound> sequence) {
		for (AminoAcidCompound c : sequence) {
			if (!hasCompound(c)) {
				return false;
			}
		}
		return true;
	}

}
