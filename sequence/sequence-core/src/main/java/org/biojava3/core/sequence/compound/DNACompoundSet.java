package org.biojava3.core.sequence.compound;

import java.util.HashMap;
import java.util.Map;

import org.biojava3.core.sequence.template.CompoundSet;


public class DNACompoundSet implements CompoundSet<NucleotideCompound>{

	private final Map<String,NucleotideCompound> dnaCompoundCache = new HashMap<String,NucleotideCompound>();
	
	public DNACompoundSet() {
		dnaCompoundCache.put("A", new NucleotideCompound("A",this,"T"));
		dnaCompoundCache.put("T", new NucleotideCompound("T",this,"A"));
		dnaCompoundCache.put("C", new NucleotideCompound("C",this,"G"));
		dnaCompoundCache.put("G", new NucleotideCompound("G",this,"C"));
		dnaCompoundCache.put("a", new NucleotideCompound("a",this,"t"));
		dnaCompoundCache.put("t", new NucleotideCompound("t",this,"a"));
		dnaCompoundCache.put("c", new NucleotideCompound("c",this,"g"));
		dnaCompoundCache.put("g", new NucleotideCompound("g",this,"c"));
	}
	
	public String getStringForCompound(NucleotideCompound compound) {
		return compound.toString();
	}

	public NucleotideCompound getCompoundForString(String string) {
		if (string.length()==0) {
			return null;
		}
		if (string.length()>this.getMaxSingleCompoundStringLength()) {
			throw new IllegalArgumentException("String supplied is too long.");
		}
		return this.dnaCompoundCache.get(string);
	}

	public int getMaxSingleCompoundStringLength() {
		return 1;
	}

 
        private final static DNACompoundSet dnaCompoundSet = new DNACompoundSet();

        static public DNACompoundSet getDNACompoundSet(){
            return dnaCompoundSet;
        }

}
