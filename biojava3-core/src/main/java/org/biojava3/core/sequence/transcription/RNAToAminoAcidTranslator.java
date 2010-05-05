package org.biojava3.core.sequence.transcription;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.biojava3.core.sequence.RNASequence;
import org.biojava3.core.sequence.compound.AminoAcidCompound;
import org.biojava3.core.sequence.compound.NucleotideCompound;
import org.biojava3.core.sequence.io.template.SequenceCreatorInterface;
import org.biojava3.core.sequence.template.AbstractCompoundTranslator;
import org.biojava3.core.sequence.template.CompoundSet;
import org.biojava3.core.sequence.template.Sequence;
import org.biojava3.core.sequence.transcription.Table.Codon;
import org.biojava3.core.sequence.views.WindowedSequence;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Multimap;

/**
 * Takes a {@link Sequence} of {@link NucleotideCompound} which should
 * represent an RNA sequence ({@link RNASequence} is good for this) and returns
 * a list of {@link Sequence} which hold {@link AminoAcidCompound}. The
 * translator can also trim stop codons as well as changing any valid
 * start codon to an initiating met.
 *
 * @author ayates
 */
public class RNAToAminoAcidTranslator extends
    AbstractCompoundTranslator<NucleotideCompound, AminoAcidCompound> {

  private final boolean                                trimStops;
  private final boolean                                initMetOnly;
  private final Map<List<NucleotideCompound>, Codon>   quickLookup;
  private final Multimap<AminoAcidCompound, Codon> aminoAcidToCodon;

  public RNAToAminoAcidTranslator(
      SequenceCreatorInterface<AminoAcidCompound> creator,
      CompoundSet<NucleotideCompound> nucleotides, CompoundSet<Codon> codons,
      CompoundSet<AminoAcidCompound> aminoAcids, Table table,
      boolean trimStops, boolean initMetOnly) {

    super(creator, nucleotides, aminoAcids);
    this.trimStops = trimStops;
    this.initMetOnly = initMetOnly;

    quickLookup = new HashMap<List<NucleotideCompound>, Codon>(codons
        .getAllCompounds().size());
    aminoAcidToCodon = ArrayListMultimap.create();

    List<Codon> codonList = table.getCodons(nucleotides, aminoAcids);
    for (Codon codon : codonList) {
      quickLookup.put(codon.getAsList(), codon);
      aminoAcidToCodon.put(codon.getAminoAcid(), codon);
    }
  }

  /**
   * Refuses to add the * compound (stop) to the available compounds if
   * trimStops is on
   */
  @Override
  protected void addCompoundToLists(List<List<AminoAcidCompound>> list,
      AminoAcidCompound compound) {
    if (trimStops && compound.getShortName().equals("*")) {
      return;
    }
    super.addCompoundToLists(list, compound);
  }

  /**
   * Performs the core conversion of RNA to Peptide. It does this by walking
   * a windowed version of the given sequence. Any trailing DNA base pairs
   * are ignored according to the specification of {@link WindowedSequence}.
   */
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

  /**
   * Performs the trimming of stop codons and the conversion of a valid start
   * amino acid to M
   */
  @Override
  protected void postProcessCompoundLists(
      List<List<AminoAcidCompound>> compoundLists) {
    for (List<AminoAcidCompound> compounds : compoundLists) {
      if (initMetOnly) {
        initMet(compounds);
      }
      if (trimStops) {
        trimStop(compounds);
      }
    }
  }

  private void initMet(List<AminoAcidCompound> sequence) {
    AminoAcidCompound initMet = getToCompoundSet().getCompoundForString("M");
    AminoAcidCompound start = sequence.get(0);
    boolean isStart = false;
    for(Codon c: aminoAcidToCodon.get(start)) {
      if(c.isStart()) {
        isStart = true;
        break;
      }
    }

    if(isStart) {
      sequence.set(0, initMet);
    }
  }

  private void trimStop(List<AminoAcidCompound> sequence) {
    AminoAcidCompound stop = sequence.get(sequence.size() - 1);
    boolean isStop = false;
    for(Codon c: aminoAcidToCodon.get(stop)) {
      if(c.isStop()) {
        isStop = true;
        break;
      }
    }

    if(isStop) {
      sequence.remove(sequence.size()-1);
    }
  }
}
