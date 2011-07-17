package org.biojava3.aaproperties.xml;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.biojava3.core.exceptions.CompoundNotFoundError;
import org.biojava3.core.sequence.compound.AminoAcidCompound;
import org.biojava3.core.sequence.template.CompoundSet;
import org.biojava3.core.sequence.template.Sequence;

public class ModifiedAminoAcidCompoundSet implements CompoundSet<AminoAcidCompound>{
	private final Map<String, AminoAcidCompound> aminoAcidCompoundCache = new HashMap<String, AminoAcidCompound>();
	
	public ModifiedAminoAcidCompoundSet(List<AminoAcidComposition> aaList, Map<Character, Double> aaSymbol2MolecularWeight){
		this.aminoAcidCompoundCache.put("-", new  AminoAcidCompound(null, "-", "", "", 0.0f));
		for(AminoAcidComposition aa:aaList){
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
            throw new IllegalArgumentException("String supplied ("+string+") is too long. Max is "+getMaxSingleCompoundStringLength());
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
	public void verifySequence(Sequence<AminoAcidCompound> sequence) throws CompoundNotFoundError {
		for (AminoAcidCompound compound : sequence) {
            if (!hasCompound(compound)) {
                throw new CompoundNotFoundError("Compound (" + compound + ") not found in AminoAcidCompoundSet.");
            }
        }		
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

}
