package org.biojava3.core.sequence.impl;

import java.util.HashMap;
import java.util.Map;

import org.biojava3.core.sequence.CompoundSet;

public class DNACompoundSet implements CompoundSet<DNACompound> {

	private final Map<CharSequence,DNACompound> dnaCompoundCache = new HashMap<CharSequence,DNACompound>();
	
	public DNACompoundSet() {
		dnaCompoundCache.put("A", new DNACompound('A','T'));
		dnaCompoundCache.put("T", new DNACompound('T','A'));
		dnaCompoundCache.put("C", new DNACompound('C','G'));
		dnaCompoundCache.put("G", new DNACompound('G','C'));
		dnaCompoundCache.put("a", new DNACompound('a','t'));
		dnaCompoundCache.put("t", new DNACompound('t','a'));
		dnaCompoundCache.put("c", new DNACompound('c','g'));
		dnaCompoundCache.put("g", new DNACompound('g','c'));
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
