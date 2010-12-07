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
    addNucleotideCompound("N", "N");
    addNucleotideCompound("-", "-");
  }

  public NucleotideCompound newNucleotideCompound(String base, String complement, String... equivalents) {
    if(equivalents.length == 0) {
      return new NucleotideCompound(base, this, complement);
    }
    else {
      NucleotideCompound[] compounds = new NucleotideCompound[equivalents.length];
      for(int i=0; i<compounds.length; i++) {
        compounds[i] = getCompoundForString(equivalents[i]);
      }
      return new NucleotideCompound(base, this, complement, compounds);
    }
  }
}
