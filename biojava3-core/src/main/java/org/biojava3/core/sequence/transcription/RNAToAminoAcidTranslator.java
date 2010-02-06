package org.biojava3.core.sequence.transcription;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.biojava3.core.sequence.compound.AminoAcidCompound;
import org.biojava3.core.sequence.compound.NucleotideCompound;
import org.biojava3.core.sequence.io.template.SequenceCreatorInterface;
import org.biojava3.core.sequence.template.AbstractCompoundTranslator;
import org.biojava3.core.sequence.template.CompoundSet;
import org.biojava3.core.sequence.template.Sequence;
import org.biojava3.core.sequence.transcription.Table.Codon;
import org.biojava3.core.sequence.views.WindowedSequence;

public class RNAToAminoAcidTranslator extends
    AbstractCompoundTranslator<NucleotideCompound, AminoAcidCompound> {

  private final boolean                              trimStops;
  private final Map<List<NucleotideCompound>, Codon> quickLookup;

  public RNAToAminoAcidTranslator(
      SequenceCreatorInterface<AminoAcidCompound> creator,
      CompoundSet<NucleotideCompound> nucleotides, CompoundSet<Codon> codons,
      CompoundSet<AminoAcidCompound> aminoAcids, Table table, boolean trimStops) {

    super(creator, nucleotides, aminoAcids);
    this.trimStops = trimStops;

    quickLookup = new HashMap<List<NucleotideCompound>, Codon>(codons
        .getAllCompounds().size());
    for (Codon codon : table.getCodons(nucleotides, aminoAcids)) {
      quickLookup.put(codon.getAsList(), codon);
    }
  }

  @Override
  protected void addCompoundToLists(List<List<AminoAcidCompound>> list,
      AminoAcidCompound compound) {
    if (trimStops && compound.getShortName().equals("*")) {
      return;
    }
    super.addCompoundToLists(list, compound);
  }

  @Override
  public List<Sequence<AminoAcidCompound>> createSequences(
      Sequence<NucleotideCompound> originalSequence) {

    List<List<AminoAcidCompound>> workingList = new ArrayList<List<AminoAcidCompound>>();
    Iterable<List<NucleotideCompound>> iter = new WindowedSequence<NucleotideCompound>(
        originalSequence, 3);
    for (List<NucleotideCompound> element : iter) {
      Codon target = quickLookup.get(element);
      addCompoundsToList(Arrays.asList(target.getAminoAcid()), workingList);
    }

    return workingListToSequences(workingList);
  }
}
