package org.biojava3.core.sequence.transcription;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.biojava3.core.sequence.RNASequence;
import org.biojava3.core.sequence.compound.AminoAcidCompound;
import org.biojava3.core.sequence.transcription.CaseInsensitiveCompound;
import org.biojava3.core.sequence.compound.NucleotideCompound;
import org.biojava3.core.sequence.io.template.SequenceCreatorInterface;
import org.biojava3.core.sequence.template.AbstractCompoundTranslator;
import org.biojava3.core.sequence.template.CompoundSet;
import org.biojava3.core.sequence.template.Sequence;
import org.biojava3.core.sequence.transcription.Table.Codon;
import org.biojava3.core.sequence.views.WindowedSequence;

/**
 * Takes a {@link Sequence} of {@link NucleotideCompound} which should
 * represent an RNA sequence ({@link RNASequence} is good for this) and returns
 * a list of {@link Sequence} which hold {@link AminoAcidCompound}. The
 * translator can also trim stop codons as well as changing any valid
 * start codon to an initiating met.
 *
 * @author ayates
 */
public class RNAToAminoAcidTranslator extends AbstractCompoundTranslator<NucleotideCompound, AminoAcidCompound> {

    private final boolean trimStops;
    private final boolean initMetOnly;
    private final Map<List<CaseInsensitiveCompound>, Codon> quickLookup;
    private final Map<AminoAcidCompound, List<Codon>> aminoAcidToCodon;
    private final NucleotideCompound nCompound;
    private final AminoAcidCompound unknownAminoAcidCompound;
    private final boolean translateNCodons;

    public RNAToAminoAcidTranslator(
            SequenceCreatorInterface<AminoAcidCompound> creator,
            CompoundSet<NucleotideCompound> nucleotides, CompoundSet<Codon> codons,
            CompoundSet<AminoAcidCompound> aminoAcids, Table table,
            boolean trimStops, boolean initMetOnly, boolean translateNCodons) {

        super(creator, nucleotides, aminoAcids);
        this.trimStops = trimStops;
        this.initMetOnly = initMetOnly;
        this.translateNCodons = translateNCodons;

        quickLookup = new HashMap<List<CaseInsensitiveCompound>, Codon>(codons.getAllCompounds().size());
        aminoAcidToCodon = new HashMap<AminoAcidCompound, List<Codon>>();

        List<Codon> codonList = table.getCodons(nucleotides, aminoAcids);
        for (Codon codon : codonList) {
            quickLookup.put(codon.getAsList(), codon);
            
            List<Codon> codonL = aminoAcidToCodon.get(codon.getAminoAcid());
            if ( codonL == null){
            	codonL = new ArrayList<Codon>();
            	aminoAcidToCodon.put(codon.getAminoAcid(), codonL);
            }
            codonL.add(codon);
            
        }

        nCompound = nucleotides.getCompoundForString("N");
        unknownAminoAcidCompound = aminoAcids.getCompoundForString("X");
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
        
        Iterable<List<NucleotideCompound>> iter = 
            new WindowedSequence<NucleotideCompound>(originalSequence, 3);
                
        for (List<NucleotideCompound> element : iter) {
            AminoAcidCompound aminoAcid;
            if(hasN(element)) {
                aminoAcid = unknownAminoAcidCompound;
            }
            else {
              List<CaseInsensitiveCompound> c = wrap(element);
              Codon target = quickLookup.get(c);
              aminoAcid = target.getAminoAcid();
            }
            addCompoundsToList(Arrays.asList(aminoAcid), workingList);
        }

        return workingListToSequences(workingList);
    }

    protected boolean hasN(List<NucleotideCompound> compounds) {
        if(! translateNCodons()) {
            return false;
        }
        for(NucleotideCompound c: compounds) {
            if(c.equalsIgnoreCase(nCompound)) {
                return true;
            }
        }
        return false;
    }
    
    protected List<CaseInsensitiveCompound> wrap(List<NucleotideCompound> list) {
      List<CaseInsensitiveCompound> output = new ArrayList<CaseInsensitiveCompound>(list.size());
      for(NucleotideCompound c: list) {
        output.add(new CaseInsensitiveCompound(c));
      }
      return output;
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
        for (Codon c : aminoAcidToCodon.get(start)) {
            if (c.isStart()) {
                isStart = true;
                break;
            }
        }

        if (isStart) {
            sequence.set(0, initMet);
        }
    }

    private void trimStop(List<AminoAcidCompound> sequence) {
        AminoAcidCompound stop = sequence.get(sequence.size() - 1);
        boolean isStop = false;
        for (Codon c : aminoAcidToCodon.get(stop)) {
            if (c.isStop()) {
                isStop = true;
                break;
            }
        }

        if (isStop) {
            sequence.remove(sequence.size() - 1);
        }
    }

    /**
     * Indicates if we want to force exact translation of compounds or not i.e.
     * those with internal N RNA bases. This will cause a translation to an
     * X amino acid
     */
    public boolean translateNCodons() {
        return translateNCodons;
    }
}
