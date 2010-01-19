package org.biojava3.core.sequence.impl;

import org.biojava3.core.sequence.Compound;

public class GenericCompound implements Compound {
	private CharSequence string;
	
	public GenericCompound(CharSequence string) {
		this.string = string;
	}
	
	public String toString() {
		return string.toString();
	}
	
	public boolean equals(Object obj) {
		if (obj==null) {
			return false;
		}
		if (!(obj instanceof GenericCompound)) {
			return false;
		}
		GenericCompound them = (GenericCompound)obj;
		return this.string.equals(them.string);
	}
	
	public int hashCode() {
		return this.string.hashCode();
	}
	
	public boolean equalsIgnoreCase(Compound compound) {
		if (compound==null) {
			return false;
		}
		if (!(compound instanceof GenericCompound)) {
			return false;
		}
		GenericCompound them = (GenericCompound)compound;
		return this.string.toString().equalsIgnoreCase(them.string.toString());
	}
}
