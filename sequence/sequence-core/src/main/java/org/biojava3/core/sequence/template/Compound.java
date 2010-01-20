package org.biojava3.core.sequence.template;

public interface Compound {

    public boolean equalsIgnoreCase(Compound compound);

    public String getDescription();

    public void setDescription(String description);

    public String getShortName();

    public void setShortName(String shortName);

    public String getLongName();

    public void setLongName(String longName);

    public Float getMolecularWeight();

    public void setMolecularWeight(Float molecularWeight);
}
