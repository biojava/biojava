package org.biojava3.core.sequence.compound;

import org.biojava3.core.sequence.template.AbstractNucleotideCompoundSet;

/**
 * @author Andy Yates
 */
public class DNACompoundSet extends AbstractNucleotideCompoundSet<NucleotideCompound> {

  private static class InitaliseOnDemand {
    public static final DNACompoundSet INSTANCE = new DNACompoundSet();
  }

  public static DNACompoundSet getDNACompoundSet() {
    return InitaliseOnDemand.INSTANCE;
  }

  public DNACompoundSet() {
    addNucleotideCompound("A", "T");
    addNucleotideCompound("T", "A");
    addNucleotideCompound("G", "C");
    addNucleotideCompound("C", "G");
  }

  protected NucleotideCompound newNucleotideCompound(String base, String complement) {
    return new NucleotideCompound(base, this, complement);
  }
}
