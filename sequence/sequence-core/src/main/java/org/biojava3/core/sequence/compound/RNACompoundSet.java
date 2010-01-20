package org.biojava3.core.sequence.compound;

import java.util.HashMap;
import java.util.Map;

import org.biojava3.core.sequence.template.CompoundSet;

public class RNACompoundSet implements CompoundSet<NucleotideCompound> {

    private final Map<String, NucleotideCompound> rnaCompoundCache = new HashMap<String, NucleotideCompound>();

    public RNACompoundSet() {
        rnaCompoundCache.put("U", new NucleotideCompound("U", this, "U"));
        rnaCompoundCache.put("A", new NucleotideCompound("A", this, "A"));
        rnaCompoundCache.put("G", new NucleotideCompound("G", this, "G"));
        rnaCompoundCache.put("C", new NucleotideCompound("C", this, "C"));
        rnaCompoundCache.put("u", new NucleotideCompound("u", this, "u"));
        rnaCompoundCache.put("a", new NucleotideCompound("a", this, "a"));
        rnaCompoundCache.put("g", new NucleotideCompound("g", this, "g"));
        rnaCompoundCache.put("c", new NucleotideCompound("c", this, "c"));
    }

    public String getStringForCompound(NucleotideCompound compound) {
        return compound.toString();
    }

    public NucleotideCompound getCompoundForString(String string) {
        if (string.length() == 0) {
            return null;
        }
        if (string.length() > this.getMaxSingleCompoundStringLength()) {
            throw new IllegalArgumentException("Sequence supplied is too long.");
        }
        return this.rnaCompoundCache.get(string);
    }

    public int getMaxSingleCompoundStringLength() {
        return 1;
    }
    private final static RNACompoundSet rnaCompoundSet = new RNACompoundSet();

    static public RNACompoundSet getRNACompoundSet() {
        return rnaCompoundSet;
    }
}
