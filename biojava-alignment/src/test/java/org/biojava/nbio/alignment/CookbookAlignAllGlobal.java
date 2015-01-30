/*
 *                    BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 *
 *      http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *      http://www.biojava.org/
 *
 */
package org.biojava.nbio.alignment;

import org.biojava.nbio.alignment.Alignments.PairwiseSequenceAlignerType;
import org.biojava.nbio.alignment.template.SequencePair;
import org.biojava.nbio.alignment.template.SubstitutionMatrix;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.core.sequence.io.FastaReaderHelper;
import org.biojava.nbio.core.util.ConcurrencyTools;

import java.net.URL;
import java.util.ArrayList;
import java.util.List;

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
        SubstitutionMatrix<AminoAcidCompound> matrix = SimpleSubstitutionMatrix.getBlosum62();
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
