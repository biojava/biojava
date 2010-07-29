package org.biojava3.alignment;

import java.net.URL;
import java.util.ArrayList;
import java.util.List;

import org.biojava3.alignment.Alignments.PairwiseSequenceAlignerType;
import org.biojava3.alignment.template.SequencePair;
import org.biojava3.alignment.template.SubstitutionMatrix;
import org.biojava3.core.sequence.ProteinSequence;
import org.biojava3.core.sequence.compound.AminoAcidCompound;
import org.biojava3.core.sequence.io.FastaReaderHelper;
import org.biojava3.core.util.ConcurrencyTools;

public class CookbookAlignAllGlobal {

    public static void main(String[] args) {
        String[] ids = new String[] {"Q21691", "Q21495", "O48771"};
        try {
            alignAllGlobal(ids);
        } catch (Exception e){
            e.printStackTrace();
        }
    }

    private static void alignAllGlobal(String[] ids) throws Exception {
        List<ProteinSequence> lst = new ArrayList<ProteinSequence>();
        for (String id : ids) {
            lst.add(getSequenceForId(id));
        }
        SubstitutionMatrix<AminoAcidCompound> matrix = new SimpleSubstitutionMatrix<AminoAcidCompound>();
        List<SequencePair<ProteinSequence, AminoAcidCompound>> alig = Alignments.getAllPairsAlignments(lst,
                PairwiseSequenceAlignerType.GLOBAL, new SimpleGapPenalty(), matrix);
        for (SequencePair<ProteinSequence, AminoAcidCompound> pair : alig) {
            System.out.printf("%n%s vs %s%n%s", pair.getQuery().getAccession(), pair.getTarget().getAccession(), pair);
        }
        ConcurrencyTools.shutdown();
    }

    private static ProteinSequence getSequenceForId(String uniProtId) throws Exception {
        URL uniprotFasta = new URL(String.format("http://www.uniprot.org/uniprot/%s.fasta", uniProtId));
        ProteinSequence seq = FastaReaderHelper.readFastaProteinSequence(uniprotFasta.openStream()).get(uniProtId);
        System.out.printf("id : %s %s%n%s%n", uniProtId, seq, seq.getOriginalHeader());
        return seq;
    }

}
