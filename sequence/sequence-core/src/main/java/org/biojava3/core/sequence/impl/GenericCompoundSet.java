package org.biojava3.core.sequence.impl;

import java.util.HashMap;
import java.util.Map;

import org.biojava3.core.sequence.CompoundSet;

public class GenericCompoundSet implements CompoundSet<GenericCompound> {

	private final Map<CharSequence,GenericCompound> genericCompoundCache = new HashMap<CharSequence,GenericCompound>();
	
	public CharSequence getCharSequenceForCompound(GenericCompound compound) {
		return compound.toString();
	}

	public GenericCompound getCompoundForCharSequence(CharSequence string) {
		if (string.length()==0) {
			return null;
		}
		if (string.length()>this.getMaxSingleCompoundCharSequenceLength()) {
			throw new IllegalArgumentException("CharSequence supplied is too long.");
		}
		if (this.genericCompoundCache.containsKey(string)) {
			return this.genericCompoundCache.get(string);
		}
		GenericCompound compound = new GenericCompound(string);
		this.genericCompoundCache.put(string, compound);
		return compound;
	}

	public int getMaxSingleCompoundCharSequenceLength() {
		return 1;
	}

}
