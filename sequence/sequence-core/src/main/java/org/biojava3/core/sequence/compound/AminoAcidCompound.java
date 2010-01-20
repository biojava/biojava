package org.biojava3.core.sequence.compound;

import org.biojava3.core.sequence.template.Compound;

public class AminoAcidCompound implements Compound {

    
    private AminoAcidCompoundSet compoundSet;
    

    public AminoAcidCompound(AminoAcidCompoundSet compoundSet,String shortName,String longName,String description,Float molecularWeight ) {
        setShortName(shortName);
        setLongName(longName);
        setDescription(description);
        setMolecularWeight(molecularWeight);
        this.compoundSet = compoundSet;
        
    }



    public String toString() {
        return shortName.toString();
    }
//TODO need to allow for modified name
    public boolean equals(Object obj) {
        if (obj == null) {
            return false;
        }
        if (!(obj instanceof AminoAcidCompound)) {
            return false;
        }
        AminoAcidCompound them = (AminoAcidCompound) obj;
        if(shortName.equals(them.shortName)){
            return true;
        }
        return longName.equals(them.longName);

    }

//    public int hashCode() {
//        return this.base.hashCode();
//    }

    public boolean equalsIgnoreCase(Compound compound) {
        if (compound == null) {
            return false;
        }
        if (!(compound instanceof AminoAcidCompound)) {
            return false;
        }
        AminoAcidCompound them = (AminoAcidCompound) compound;
        if(shortName.toString().equalsIgnoreCase(them.shortName.toString())){
            return true;
        }
        return longName.toString().equalsIgnoreCase(them.shortName.toString());
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
}
