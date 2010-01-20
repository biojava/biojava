package org.biojava3.core.sequence.compound;

import org.biojava3.core.sequence.template.Compound;
import org.biojava3.core.sequence.template.CompoundSet;
import org.biojava3.core.sequence.template.NucleotideCompoundInterface;

public class NucleotideCompound<C extends Compound> implements NucleotideCompoundInterface {

    private String base;
    private CompoundSet<NucleotideCompound> compoundSet;
    private String complementStr;

    public NucleotideCompound(String base, CompoundSet<NucleotideCompound> compoundSet, String complementStr) {
        this.base = base;
        this.compoundSet = compoundSet;
        this.complementStr = complementStr;
    }

    public NucleotideCompound getComplement() {
        return this.compoundSet.getCompoundForString(this.complementStr);
    }

    public String toString() {
        return base.toString();
    }

    public boolean equals(Object obj) {
        if (obj == null) {
            return false;
        }
        if (!(obj instanceof NucleotideCompound)) {
            return false;
        }
        NucleotideCompound them = (NucleotideCompound) obj;
        return this.base.equals(them.base);
    }

    public int hashCode() {
        return this.base.hashCode();
    }

    public boolean equalsIgnoreCase(Compound compound) {
        if (compound == null) {
            return false;
        }
        if (!(compound instanceof NucleotideCompound)) {
            return false;
        }
        NucleotideCompound them = (NucleotideCompound) compound;
        return this.base.toString().equalsIgnoreCase(them.base.toString());
    }
    private String shortName = null;
    private String longName = null;
    private String description = null;
    private Float molecularWeight = null;

    public String getDescription() {
        return description;
    }

    public void setDescription(String _description) {
        description = _description;
    }

    public String getShortName() {
        return shortName;
    }

    public void setShortName(String _shortName) {
        shortName = _shortName;
    }

    public String getLongName() {
        return longName;
    }

    public void setLongName(String _longName) {
        longName = _longName;
    }

    public Float getMolecularWeight() {
        return molecularWeight;
    }

    public void setMolecularWeight(Float _molecularWeight) {
        molecularWeight = _molecularWeight;
    }

    public void setComplement(NucleotideCompoundInterface nucleotideCompound) {
        throw new UnsupportedOperationException("Not supported yet.");
    }
}
