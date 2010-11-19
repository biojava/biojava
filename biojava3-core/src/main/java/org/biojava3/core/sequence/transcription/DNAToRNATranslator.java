package org.biojava3.core.sequence.transcription;

import java.util.ArrayList;
import java.util.List;
import org.biojava3.core.sequence.RNASequence;

import org.biojava3.core.sequence.compound.NucleotideCompound;
import org.biojava3.core.sequence.io.template.SequenceCreatorInterface;
import org.biojava3.core.sequence.template.AbstractCompoundTranslator;
import org.biojava3.core.sequence.template.CompoundSet;
import org.biojava3.core.sequence.template.ProxySequenceReader;
import org.biojava3.core.sequence.template.Sequence;
import org.biojava3.core.sequence.views.RnaSequenceView;

/**
 * Performs the first stage of transcription by going from DNA to RNA. This
 * class will first delegate to {@link Frame} in order to be in the correctly
 * specified translation frame and then translates T to U. The other
 * translation carried out is to convert an equivalent compound in DNA to RNA
 * i.e. for the base A in DNA fetching the equivalent A base in the RNA
 * {@link CompoundSet}.
 *
 * @author ayates
 */
public class DNAToRNATranslator extends AbstractCompoundTranslator<NucleotideCompound, NucleotideCompound> {

    private final boolean shortCutTranslation;

  public DNAToRNATranslator(SequenceCreatorInterface<NucleotideCompound> rnaCreator,
      CompoundSet<NucleotideCompound> dna, CompoundSet<NucleotideCompound> rna,
      boolean shortCutTranslation) {
    super(rnaCreator, dna, rna);
    this.shortCutTranslation = shortCutTranslation;
    defaultMappings();
    thyamineToUracil();
  }

    /**
     * Overloaded local version which delegates to an optional translator
     * when told to (specified during construction).
     *
     * @param originalSequence The DNA sequence to translate
     * @return The translated single sequence
     */
    @Override
    public List<Sequence<NucleotideCompound>> createSequences(Sequence<NucleotideCompound> originalSequence) {
        if(shortCutTranslation) {
            List<Sequence<NucleotideCompound>> result = new ArrayList<Sequence<NucleotideCompound>>(1);
            result.add(wrapToRna(originalSequence));
            return result;
        }
        else {
            return super.createSequences(originalSequence);
        }
    }

    /**
     * Takes in the given DNA Sequence and returns an instance of RNASequence
     * which is using {@link RnaSequenceView} as a
     * {@link ProxySequenceReader}.
     */
    protected RNASequence wrapToRna(Sequence<NucleotideCompound> dna) {
        ProxySequenceReader<NucleotideCompound> rnaView = new RnaSequenceView(dna);
        return new RNASequence(rnaView);
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
