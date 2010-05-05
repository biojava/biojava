package org.biojava3.core.sequence.transcription;

import java.util.List;

import org.biojava3.core.sequence.compound.NucleotideCompound;
import org.biojava3.core.sequence.io.template.SequenceCreatorInterface;
import org.biojava3.core.sequence.template.AbstractCompoundTranslator;
import org.biojava3.core.sequence.template.CompoundSet;
import org.biojava3.core.sequence.template.Sequence;

public class DNAToRNATranslator extends AbstractCompoundTranslator<NucleotideCompound, NucleotideCompound>{

  public DNAToRNATranslator(SequenceCreatorInterface<NucleotideCompound> rnaCreator,
      CompoundSet<NucleotideCompound> dna, CompoundSet<NucleotideCompound> rna) {
    super(rnaCreator, dna, rna);
    defaultMappings();
    thyamineToUracil();
  }

  private void defaultMappings() {
    NucleotideCompound thymine = getFromCompoundSet().getCompoundForString("T");
    for(NucleotideCompound dnaBase: getFromCompoundSet().getAllCompounds()) {
      if(dnaBase.equalsIgnoreCase(thymine)) {
        continue;
      }
      NucleotideCompound rnaBase = getToCompoundSet().getCompoundForString(
          dnaBase.toString());
      addCompounds(dnaBase, rnaBase);
    }

  }

  private void thyamineToUracil() {
    addCompounds(getFromCompoundSet().getCompoundForString("T"),
        getToCompoundSet().getCompoundForString("U"));
    addCompounds(getFromCompoundSet().getCompoundForString("t"),
        getToCompoundSet().getCompoundForString("u"));
  }

  public Sequence<NucleotideCompound> createSequence(Sequence<NucleotideCompound> originalSequence, Frame frame) {
    Sequence<NucleotideCompound> wrapped = frame.wrap(originalSequence);
    return super.createSequence(wrapped);
  }

  @Override
  public Sequence<NucleotideCompound> createSequence(Sequence<NucleotideCompound> originalSequence) {
    return createSequence(originalSequence, Frame.getDefaultFrame());
  }

  @Override
  protected void postProcessCompoundLists(
      List<List<NucleotideCompound>> compoundLists) {
    //No post processing needed
  }
}
