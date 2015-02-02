/**
 *                    BioJava development code
 *
 * This code may be freely distributed and modified under the terms of the GNU
 * Lesser General Public Licence. This should be distributed with the code. If
 * you do not have a copy, see:
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
 * Created on Oct 5, 2011 
 * Created by Andreas Prlic
 *
 * @since 3.0.2
 */
package org.biojava.nbio.alignment;

import junit.framework.TestCase;
import org.biojava.nbio.alignment.Alignments.PairwiseSequenceAlignerType;
import org.biojava.nbio.alignment.template.PairwiseSequenceAligner;
import org.biojava.nbio.alignment.template.Profile;
import org.biojava.nbio.alignment.template.SequencePair;
import org.biojava.nbio.alignment.template.SubstitutionMatrix;
import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.compound.AmbiguityDNACompoundSet;
import org.biojava.nbio.core.sequence.compound.DNACompoundSet;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;
import org.biojava.nbio.core.sequence.io.FastaReaderHelper;
import org.biojava.nbio.core.util.ConcurrencyTools;

import java.io.InputStream;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;

public class TestDNAAlignment extends TestCase {

	private static final double PRECISION = 0.00000001;
	

    public void testDNAAlignment() {

        try {
            List<DNASequence> lst = getDNAFASTAFile();

            Profile<DNASequence, NucleotideCompound> profile = Alignments.getMultipleSequenceAlignment(lst);

            assertTrue(profile.getSize() == 10);

            assertTrue(profile.getAlignedSequence(1).getSequenceAsString().length() > 50);


            // here how to print the MSA:

            //System.out.printf("MSA:%n%s%n", profile);
        } catch (Exception e) {
            e.printStackTrace();
            fail(e.getMessage());
        }
        ConcurrencyTools.shutdown();
    }

    private static List<DNASequence> getDNAFASTAFile() throws Exception {

        InputStream inStream = TestDNAAlignment.class.getResourceAsStream(String.format("/dna-fasta.txt"));
        LinkedHashMap<String, DNASequence> fastas = FastaReaderHelper.readFastaDNASequence(inStream);

        List<DNASequence> sequences = new ArrayList<DNASequence>();

        for (String key : fastas.keySet()) {
            DNASequence seq = fastas.get(key);
            sequences.add(seq);
        }

        return sequences;
    }

    /**
     * @author brandstaetter 
     */
    public void testDNAMultipleAlignmentWithMixedCompoundSets() throws CompoundNotFoundException {

        DNASequence target = new DNASequence("ACTGACGTGTAGCTGACTGA", DNACompoundSet.getDNACompoundSet());
        DNASequence query = new DNASequence("ACTGACGTGTAGCTGACTGTA", AmbiguityDNACompoundSet.getDNACompoundSet());

        List<DNASequence> lst = new ArrayList<DNASequence>();
        lst.add(target);
        lst.add(query);

        try {
        	@SuppressWarnings("unused")
			Profile<DNASequence, NucleotideCompound> profile = Alignments.getMultipleSequenceAlignment(lst);
        	fail("Alignments.getMultipleSequenceAlignment(lst) expected exception with differing compound sets");
        } catch (IllegalArgumentException ex) {
        	// expected exception
        }
    }

    /**
     * @author brandstaetter
     */
    public void testDNAPairwiseAlignmentWithMixedCompoundSets() throws CompoundNotFoundException {
        DNASequence target = new DNASequence("ACTGACGTGTAGCTGACTGA", DNACompoundSet.getDNACompoundSet());
        DNASequence query = new DNASequence("ACTGACGTGTAGCTGACTGT", AmbiguityDNACompoundSet.getDNACompoundSet());
        SubstitutionMatrix<NucleotideCompound> matrix = SubstitutionMatrixHelper.getNuc4_4();
        SimpleGapPenalty gapP = new SimpleGapPenalty();
        gapP.setOpenPenalty((short) 5);
        gapP.setExtensionPenalty((short) 2);
        
        try {
        	@SuppressWarnings("unused")
			SequencePair<DNASequence, NucleotideCompound> psa = Alignments.getPairwiseAlignment(query, target, PairwiseSequenceAlignerType.LOCAL, gapP, matrix);
        	fail("Alignments.getPairwiseAlignment() expected exception with differing compound sets");
        } catch (IllegalArgumentException ex) {
        	// expected exception
        }
    }
    /**
     * @author Daniel Cameron
     */
    public void testMixedCaseInputStringsMatchUnderlyingBases() throws CompoundNotFoundException {
        DNASequence target = new DNASequence("AAAAAAAAGTC", DNACompoundSet.getDNACompoundSet());
        DNASequence query = new DNASequence("aaaaaaaagtc", DNACompoundSet.getDNACompoundSet());
        SubstitutionMatrix<NucleotideCompound> matrix = SubstitutionMatrixHelper.getNuc4_4();
        SimpleGapPenalty gapP = new SimpleGapPenalty((short)5, (short)2);
        // should be a full match with +5 per match
        assertEquals(5.0 * query.getLength(), Alignments.getPairwiseAligner(query, target, PairwiseSequenceAlignerType.LOCAL, gapP, matrix).getScore(), PRECISION);
    }
    /**
     * @author Daniel Cameron
     */
    public void testNoAlignedBases() throws CompoundNotFoundException {
        DNASequence target = new DNASequence("A", DNACompoundSet.getDNACompoundSet());
        DNASequence query = new DNASequence("T", DNACompoundSet.getDNACompoundSet());
        SubstitutionMatrix<NucleotideCompound> matrix = SubstitutionMatrixHelper.getNuc4_4();
        SimpleGapPenalty gapP = new SimpleGapPenalty((short)0, (short)1);
        PairwiseSequenceAligner<DNASequence, NucleotideCompound> aligner = Alignments.getPairwiseAligner(query, target, PairwiseSequenceAlignerType.GLOBAL, gapP, matrix);
        assertEquals(2, aligner.getPair().getLength());
    }
    /**
    * @author Daniel Cameron
    */
	public void testLinearAlignment() throws CompoundNotFoundException {
		DNASequence query = new DNASequence("GTAAAAG", DNACompoundSet.getDNACompoundSet());
		DNASequence target = new DNASequence("GAAAACGTTTTTTTTTT", DNACompoundSet.getDNACompoundSet());
		SubstitutionMatrix<NucleotideCompound> matrix = SubstitutionMatrixHelper.getNuc4_4();
		SimpleGapPenalty gapP = new SimpleGapPenalty((short)0, (short)3);
		PairwiseSequenceAligner<DNASequence, NucleotideCompound> aligner = Alignments.getPairwiseAligner(query, target, PairwiseSequenceAlignerType.GLOBAL, gapP, matrix);
		assertEquals(String.format("GTAAAA-G----------%nG-AAAACGTTTTTTTTTT%n"), aligner.getPair().toString());;
	}
}
