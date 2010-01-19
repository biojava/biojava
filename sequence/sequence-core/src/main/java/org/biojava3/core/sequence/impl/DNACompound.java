package org.biojava3.core.sequence.impl;

import org.biojava3.core.sequence.Compound;

public class DNACompound implements Compound {
	private char base;
	private char complement;
	
	public DNACompound(char base, char complement) {
		this.base = base;
	}
	
	public char getComplement() {
		return this.complement;
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
