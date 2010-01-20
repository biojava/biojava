package org.biojava3.core.sequence.compound;

import java.util.HashMap;
import java.util.Map;

import org.biojava3.core.sequence.template.CompoundSet;

public class AminoAcidCompoundSet implements CompoundSet<AminoAcidCompound> {

    private final Map<String, AminoAcidCompound> aminoAcidCompoundCache = new HashMap<String, AminoAcidCompound>();

    public AminoAcidCompoundSet() {
        aminoAcidCompoundCache.put("A", new AminoAcidCompound(this, "A", "Ala", "Alanine", 0.0f));
        aminoAcidCompoundCache.put("R", new AminoAcidCompound(this, "R", "Arg", "Arginine", 0.0f));
        aminoAcidCompoundCache.put("N", new AminoAcidCompound(this, "N", "Asn", "Asparagine", 0.0f));
        aminoAcidCompoundCache.put("D", new AminoAcidCompound(this, "D", "Asp", "Aspartic acid", 0.0f));
        aminoAcidCompoundCache.put("C", new AminoAcidCompound(this, "C", "Cys", "Cysteine", 0.0f));
        aminoAcidCompoundCache.put("E", new AminoAcidCompound(this, "E", "Glu", "Glutamic acid", 0.0f));
        aminoAcidCompoundCache.put("Q", new AminoAcidCompound(this, "Q", "Gln", "Glutamine", 0.0f));
        aminoAcidCompoundCache.put("G", new AminoAcidCompound(this, "G", "Gly", "Glycine", 0.0f));
        aminoAcidCompoundCache.put("H", new AminoAcidCompound(this, "H", "His", "Histidine", 0.0f));
        aminoAcidCompoundCache.put("I", new AminoAcidCompound(this, "I", "Ile", "Isoleucine", 0.0f));
        aminoAcidCompoundCache.put("L", new AminoAcidCompound(this, "L", "Leu", "Leucine", 0.0f));
        aminoAcidCompoundCache.put("K", new AminoAcidCompound(this, "K", "Lys", "Lysine", 0.0f));
        aminoAcidCompoundCache.put("M", new AminoAcidCompound(this, "M", "Met", "Methionine", 0.0f));
        aminoAcidCompoundCache.put("F", new AminoAcidCompound(this, "F", "Phe", "Phenylalanine", 0.0f));
        aminoAcidCompoundCache.put("P", new AminoAcidCompound(this, "P", "Pro", "Proline", 0.0f));
        aminoAcidCompoundCache.put("S", new AminoAcidCompound(this, "S", "Ser", "Serine", 0.0f));
        aminoAcidCompoundCache.put("T", new AminoAcidCompound(this, "T", "Thr", "Threonine", 0.0f));
        aminoAcidCompoundCache.put("W", new AminoAcidCompound(this, "W", "Trp", "Tryptophan", 0.0f));
        aminoAcidCompoundCache.put("Y", new AminoAcidCompound(this, "Y", "Tyr", "Tyrosine", 0.0f));
        aminoAcidCompoundCache.put("V", new AminoAcidCompound(this, "V", "Val", "Valine", 0.0f));
        aminoAcidCompoundCache.put("B", new AminoAcidCompound(this, "B", "Asx", "Asparagine or Aspartic acid", null));
        aminoAcidCompoundCache.put("Z", new AminoAcidCompound(this, "Z", "Glx", "Glutamine or glutamic acid", null));
        aminoAcidCompoundCache.put("J", new AminoAcidCompound(this, "J", "Xle", "Leucine or Isoleucine", null));
        aminoAcidCompoundCache.put("X", new AminoAcidCompound(this, "Z", "Xaa", "Unspecified", null));

    }

    public String getStringForCompound(AminoAcidCompound compound) {
        return compound.toString();
    }

    public AminoAcidCompound getCompoundForString(String string) {
        if (string.length() == 0) {
            return null;
        }
        if (string.length() > this.getMaxSingleCompoundStringLength()) {
            throw new IllegalArgumentException("String supplied is too long.");
        }
        return this.aminoAcidCompoundCache.get(string);
    }

    public int getMaxSingleCompoundStringLength() {
        return 1;
    }
    private final static AminoAcidCompoundSet aminoAcidCompoundSet = new AminoAcidCompoundSet();

    static public AminoAcidCompoundSet getAminoAcidCompoundSet() {
        return aminoAcidCompoundSet;
    }
}
