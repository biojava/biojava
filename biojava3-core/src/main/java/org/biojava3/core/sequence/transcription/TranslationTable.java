package org.biojava3.core.sequence.transcription;

import org.biojava3.core.sequence.compound.AminoAcidCompound;
import org.biojava3.core.sequence.compound.NucleotideCompound;
import org.biojava3.core.sequence.template.CompoundSet;

public class TranslationTable {

  private final CompoundSet<NucleotideCompound> nucleotides;
  private final CompoundSet<AminoAcidCompound> aminoAcids;

  public TranslationTable(CompoundSet<NucleotideCompound> nucleotides,
      CompoundSet<AminoAcidCompound> aminoAcids) {
    this.nucleotides = nucleotides;
    this.aminoAcids = aminoAcids;
  }

  public static class Codon {

    private final NucleotideCompound baseOne;
    private final NucleotideCompound baseTwo;
    private final NucleotideCompound baseThree;
    private boolean isStartCodon;
    private final AminoAcidCompound aminoAcid;

    public Codon(NucleotideCompound baseOne, NucleotideCompound baseTwo,
        NucleotideCompound baseThree, boolean isStartCodon, AminoAcidCompound aminoAcid) {
      super();
      this.baseOne = baseOne;
      this.baseTwo = baseTwo;
      this.baseThree = baseThree;
      this.isStartCodon = isStartCodon;
      this.aminoAcid = aminoAcid;
    }

    public NucleotideCompound getBaseOne() {
      return baseOne;
    }

    public NucleotideCompound getBaseTwo() {
      return baseTwo;
    }

    public NucleotideCompound getBaseThree() {
      return baseThree;
    }

    public boolean isStartCodon() {
      return isStartCodon;
    }

    public AminoAcidCompound getAminoAcid() {
      return aminoAcid;
    }
  }



}
