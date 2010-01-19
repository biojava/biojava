package org.biojava3.core.sequence.impl;

import org.biojava3.core.sequence.Compound;

public class DNACompound implements Compound {
	private char base;
	private DNACompoundSet compoundSet;
	private CharSequence complementStr;
	
	public DNACompound(char base, DNACompoundSet compoundSet, CharSequence complementStr) {
		this.base = base;
		this.compoundSet = compoundSet;
		this.complementStr = complementStr;
	}
	
	public DNACompound getComplement() {
		return this.compoundSet.getCompoundForCharSequence(this.complementStr);
	}
	
	public String toString() {
		return ""+base;
	}
	
	public boolean equals(Object obj) {
		if (obj==null) {
			return false;
		}
		if (!(obj instanceof DNACompound)) {
			return false;
		}
		DNACompound them = (DNACompound)obj;
		return this.base == them.base;
	}
	
	public int hashCode() {
		return this.base;
	}
	
	public boolean equalsIgnoreCase(Compound compound) {
		if (compound==null) {
			return false;
		}
		if (!(compound instanceof DNACompound)) {
			return false;
		}
		DNACompound them = (DNACompound)compound;
		return this.toString().equalsIgnoreCase(them.toString());
	}
}
