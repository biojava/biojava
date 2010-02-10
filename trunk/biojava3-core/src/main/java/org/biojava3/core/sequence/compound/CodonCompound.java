package org.biojava3.core.sequence.compound;

import org.biojava3.core.sequence.template.AbstractCompound;
import org.biojava3.core.sequence.template.Compound;

public class CodonCompound extends AbstractCompound {

  private final NucleotideCompound one;
  private final NucleotideCompound two;
  private final NucleotideCompound three;
  private final boolean start;

  public CodonCompound(NucleotideCompound one, NucleotideCompound two,
      NucleotideCompound three, boolean start) {
    super(one.toString()+two.toString()+three.toString());
    this.one = one;
    this.two = two;
    this.three = three;
    this.start = start;
  }

  public boolean equalsIgnoreCase(Compound compound) {
    if (compound == null) {
      return false;
    }
    if (!(compound instanceof CodonCompound)) {
      return false;
    }
    CodonCompound them = (CodonCompound) compound;
    return toString().equalsIgnoreCase(them.toString());
  }

  public boolean equals(Object obj) {
    if (obj == null) {
      return false;
    }
    if (!(obj instanceof CodonCompound)) {
      return false;
    }
    CodonCompound them = (CodonCompound)obj;
    return toString().equals(them.toString());
  }

  public int hashCode() {
    return toString().hashCode();
  }

  public NucleotideCompound getOne() {
    return one;
  }

  public NucleotideCompound getTwo() {
    return two;
  }

  public NucleotideCompound getThree() {
    return three;
  }

  public boolean isStart() {
    return start;
  }

  public String getDescription() {
    // TODO Auto-generated method stub
    return null;
  }

  public String getLongName() {
    // TODO Auto-generated method stub
    return null;
  }

  public Float getMolecularWeight() {
    // TODO Auto-generated method stub
    return null;
  }

  public String getShortName() {
    // TODO Auto-generated method stub
    return null;
  }

  public void setDescription(String description) {
    // TODO Auto-generated method stub

  }

  public void setLongName(String longName) {
    // TODO Auto-generated method stub

  }

  public void setMolecularWeight(Float molecularWeight) {
    // TODO Auto-generated method stub

  }

  public void setShortName(String shortName) {
    // TODO Auto-generated method stub

  }

}
