package org.biojava3.core.sequence.impl;

import java.util.HashMap;
import java.util.Map;

import org.biojava3.core.sequence.CompoundSet;

public class DNACompoundSet implements CompoundSet<DNACompound> {

	private final Map<CharSequence,DNACompound> dnaCompoundCache = new HashMap<CharSequence,DNACompound>();
	
	public DNACompoundSet() {
		dnaCompoundCache.put("A", new DNACompound('A',this,"T"));
		dnaCompoundCache.put("T", new DNACompound('T',this,"A"));
		dnaCompoundCache.put("C", new DNACompound('C',this,"G"));
		dnaCompoundCache.put("G", new DNACompound('G',this,"C"));
		dnaCompoundCache.put("a", new DNACompound('a',this,"t"));
		dnaCompoundCache.put("t", new DNACompound('t',this,"a"));
		dnaCompoundCache.put("c", new DNACompound('c',this,"g"));
		dnaCompoundCache.put("g", new DNACompound('g',this,"c"));
	}
	
	public CharSequence getCharSequenceForCompound(DNACompound compound) {
		return compound.toString();
	}

	public DNACompound getCompoundForCharSequence(CharSequence string) {
		if (string.length()==0) {
			return null;
		}
		if (string.length()>this.getMaxSingleCompoundCharSequenceLength()) {
			throw new IllegalArgumentException("CharSequence supplied is too long.");
		}
		return this.dnaCompoundCache.get(string);
	}

	public int getMaxSingleCompoundCharSequenceLength() {
		return 1;
	}

}
