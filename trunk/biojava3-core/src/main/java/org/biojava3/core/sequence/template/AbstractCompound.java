package org.biojava3.core.sequence.template;

/**
 *
 * @author Andy Yates
 */
public abstract class AbstractCompound implements Compound {

  private String base;
  private String shortName = null;
  private String longName = null;
  private String description = null;
  private Float molecularWeight = null;

  public AbstractCompound(String base) {
    this.base = base;
  }

  public String getBase() {
    return base;
  }

  public String getDescription() {
    return description;
  }

  public void setDescription(String description) {
    this.description = description;
  }

  public String getShortName() {
    return shortName;
  }

  public void setShortName(String shortName) {
    this.shortName = shortName;
  }

  public String getLongName() {
    return longName;
  }

  public void setLongName(String longName) {
    this.longName = longName;
  }

  public Float getMolecularWeight() {
    return molecularWeight;
  }

  public void setMolecularWeight(Float molecularWeight) {
    this.molecularWeight = molecularWeight;
  }

  public String toString() {
    return base;
  }

//  public boolean equals(Object obj) {
//    if (obj == null) {
//      return false;
//    }
//    if (!(obj instanceof NucleotideCompound)) {
//      return false;
//    }
//    NucleotideCompound them = (NucleotideCompound) obj;
//    return this.base.equals(them.base);
//  }

//  public int hashCode() {
//    return this.base.hashCode();
//  }

//  public boolean equalsIgnoreCase(Compound compound) {
//    if (compound == null) {
//      return false;
//    }
//    if (!(compound instanceof NucleotideCompound)) {
//      return false;
//    }
//    NucleotideCompound them = (NucleotideCompound) compound;
//    return this.base.toString().equalsIgnoreCase(them.base.toString());
//  }
}
