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
 * Created on June 17, 2010
 * Author: Mark Chapman
 */

package org.biojava3.alignment;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.biojava3.alignment.template.GapPenalty;
import org.biojava3.alignment.template.SequencePair;
import org.biojava3.alignment.template.SubstitutionMatrix;
import org.biojava3.core.exceptions.CompoundNotFoundException;
import org.biojava3.core.sequence.DNASequence;
import org.biojava3.core.sequence.ProteinSequence;
import org.biojava3.core.sequence.compound.AmbiguityDNACompoundSet;
import org.biojava3.core.sequence.compound.AminoAcidCompound;
import org.biojava3.core.sequence.compound.DNACompoundSet;
import org.biojava3.core.sequence.compound.NucleotideCompound;
import org.junit.Before;
import org.junit.Test;

public class NeedlemanWunschTest {

	private static final double PRECISION = 0.00000001;
	
    private ProteinSequence query, target;
    private GapPenalty gaps;
    private SubstitutionMatrix<AminoAcidCompound> blosum62;
    private NeedlemanWunsch<ProteinSequence, AminoAcidCompound> alignment, self;

    @Before
    public void setup() throws CompoundNotFoundException { 
        query = new ProteinSequence("ARND");
        target = new ProteinSequence("RDG");
        gaps = new SimpleGapPenalty(10, 1);
        blosum62 = SubstitutionMatrixHelper.getBlosum62();
        alignment = new NeedlemanWunsch<ProteinSequence, AminoAcidCompound>(query, target, gaps, blosum62);
        self = new NeedlemanWunsch<ProteinSequence, AminoAcidCompound>(query, query, gaps, blosum62);
    }

    @Test
    public void testNeedlemanWunsch() {
        NeedlemanWunsch<ProteinSequence, AminoAcidCompound> nw =
                new NeedlemanWunsch<ProteinSequence, AminoAcidCompound>();
        nw.setQuery(query);
        nw.setTarget(target);
        nw.setGapPenalty(gaps);
        nw.setSubstitutionMatrix(blosum62);
        assertEquals(nw.getScore(), alignment.getScore(), PRECISION);
    }

    @Test
    public void testGetQuery() {
        assertEquals(alignment.getQuery(), query);
        assertEquals(self.getQuery(), query);
    }

    @Test
    public void testGetTarget() {
        assertEquals(alignment.getTarget(), target);
        assertEquals(self.getTarget(), query);
    }

    @Test
    public void testGetGapPenalty() {
        assertEquals(alignment.getGapPenalty(), gaps);
        assertEquals(self.getGapPenalty(), gaps);
    }

    @Test
    public void testGetSubstitutionMatrix() {
        assertEquals(alignment.getSubstitutionMatrix(), blosum62);
        assertEquals(self.getSubstitutionMatrix(), blosum62);
    }

    @Test
    public void testIsStoringScoreMatrix() {
        assertFalse(alignment.isStoringScoreMatrix());
        assertFalse(self.isStoringScoreMatrix());
    }

    @Test
    public void testGetScoreMatrix() {
        int[][][] scores = alignment.getScoreMatrix();
        assertEquals(4, scores[2][1][0]);
        scores = self.getScoreMatrix();
        assertEquals(4, scores[3][4][2]);
    }

    @Test
    public void testGetScoreMatrixAsString() {
    	//System.out.println(alignment.getScoreMatrixAsString());
    	//System.out.println(self.getScoreMatrixAsString());
        assertEquals(String.format(
                "Substitution%n" +
                "        R   D   G%n" +
                "    0  -\u221E  -\u221E  -\u221E%n" +
                "A  -\u221E  -1  -3  -2%n" +
                "R  -\u221E   4  -3  -5%n" +
                "N  -\u221E  -2   5  -3%n" +
                "D  -\u221E  -5   4   4%n" +
                "%nDeletion%n" +
                "        R   D   G%n" +
                "    0  -\u221E  -\u221E  -\u221E%n" +
                "A  -1  -\u221E  -\u221E  -\u221E%n" +
                "R  -2 -12 -14 -13%n" +
                "N  -3  -7 -14 -14%n" +
                "D  -4  -8  -6 -14%n" +
                "%nInsertion%n" +
                "        R   D   G%n" +
                "    0  -1  -2  -3%n" +
                "A  -\u221E  -\u221E -12 -13%n" +
                "R  -\u221E  -\u221E  -7  -8%n" +
                "N  -\u221E  -\u221E -13  -6%n" +
                "D  -\u221E  -\u221E -16  -7%n"),
                alignment.getScoreMatrixAsString());
        assertEquals(String.format(
                "Substitution%n" +
                "        A   R   N   D%n" +
                "    0  -\u221E  -\u221E  -\u221E  -\u221E%n" +
                "A  -\u221E   4  -2  -4  -5%n" +
                "R  -\u221E  -2   9  -2  -6%n" +
                "N  -\u221E  -4  -2  15  -1%n" +
                "D  -\u221E  -5  -6  -1  21%n" +
                "%nDeletion%n" +
                "        A   R   N   D%n" +
                "    0  -\u221E  -\u221E  -\u221E  -\u221E%n" +
                "A  -1  -\u221E  -\u221E  -\u221E  -\u221E%n" +
                "R  -2  -7 -13 -15 -16%n" +
                "N  -3  -8  -2 -13 -17%n" +
                "D  -4  -9  -3   4 -12%n" +
                "%nInsertion%n" +
                "        A   R   N   D%n" +
                "    0  -1  -2  -3  -4%n" +
                "A  -\u221E  -\u221E  -7  -8  -9%n" +
                "R  -\u221E  -\u221E -13  -2  -3%n" +
                "N  -\u221E  -\u221E -15 -13   4%n" +
                "D  -\u221E  -\u221E -16 -17 -12%n"),
                self.getScoreMatrixAsString());
    }

    @Test
    public void testGetComputationTime() {
        assertTrue(alignment.getComputationTime() > 0);
        assertTrue(self.getComputationTime() > 0);
    }

    @Test
    public void testGetProfile() {
        assertEquals(String.format("ARND%n-RDG%n"), alignment.getProfile().toString());
        assertEquals(String.format("ARND%nARND%n"), self.getProfile().toString());
    }

    @Test
    public void testGetMaxScore() {
        assertEquals(21, alignment.getMaxScore(), PRECISION);
        assertEquals(21, self.getMaxScore(), PRECISION);
    }

    @Test
    public void testGetMinScore() {
        assertEquals(-27, alignment.getMinScore(), PRECISION);
        assertEquals(-28, self.getMinScore(), PRECISION);
    }

    @Test
    public void testGetScore() {
        assertEquals(4, alignment.getScore(), PRECISION);
        assertEquals(21, self.getScore(), PRECISION);
    }

    @Test
    public void testGetPair() {
        assertEquals(String.format("ARND%n-RDG%n"), alignment.getPair().toString());
        assertEquals(String.format("ARND%nARND%n"), self.getPair().toString());
    }
    
    /**
     * @author Daniel Cameron 
     */
    @Test
	public void should_align_all_anchored() throws CompoundNotFoundException {
    	DNASequence query = new DNASequence("ACG", AmbiguityDNACompoundSet.getDNACompoundSet());
		DNASequence target = new DNASequence("CGT", AmbiguityDNACompoundSet.getDNACompoundSet());
		NeedlemanWunsch<DNASequence, NucleotideCompound> aligner = new NeedlemanWunsch<DNASequence, NucleotideCompound>(query, target, new SimpleGapPenalty((short)0, (short)10), SubstitutionMatrixHelper.getNuc4_4());
		aligner.setAnchors(new int[] { 0, 1, 2} );
		assertEquals(String.format("ACG%nCGT%n"), aligner.getPair().toString());
    }
    
    /**
     * @author Daniel Cameron 
     */
    @Test
	public void should_align_starting_anchor() throws CompoundNotFoundException {
    	DNASequence query = new DNASequence("AAT", AmbiguityDNACompoundSet.getDNACompoundSet());
		DNASequence target = new DNASequence("AATT", AmbiguityDNACompoundSet.getDNACompoundSet());
		NeedlemanWunsch<DNASequence, NucleotideCompound> aligner = new NeedlemanWunsch<DNASequence, NucleotideCompound>(query, target, new SimpleGapPenalty((short)0, (short)10), SubstitutionMatrixHelper.getNuc4_4());
		aligner.setAnchors(new int[] { 1, -1, -1} );
		assertEquals(String.format("-AAT%nAATT%n"), aligner.getPair().toString());
    }
    
    /**
     * @author Daniel Cameron 
     */
    @Test
	public void should_align_ending_anchor() throws CompoundNotFoundException {
    	DNASequence query = new DNASequence("AAG", AmbiguityDNACompoundSet.getDNACompoundSet());
		DNASequence target = new DNASequence("AATT", AmbiguityDNACompoundSet.getDNACompoundSet());
		NeedlemanWunsch<DNASequence, NucleotideCompound> aligner = new NeedlemanWunsch<DNASequence, NucleotideCompound>(query, target, new SimpleGapPenalty((short)0, (short)10), SubstitutionMatrixHelper.getNuc4_4());
		aligner.addAnchor(2, 3);
		assertEquals(String.format("AA-G%nAATT%n"), aligner.getPair().toString());
    }
    
    /**
     * @author Daniel Cameron 
     */
    @Test
	public void should_align_middle_anchor() throws CompoundNotFoundException {
    	DNASequence query = new DNASequence("ACTTT", AmbiguityDNACompoundSet.getDNACompoundSet());
		DNASequence target = new DNASequence("ACGTTT", AmbiguityDNACompoundSet.getDNACompoundSet());
		NeedlemanWunsch<DNASequence, NucleotideCompound> aligner = new NeedlemanWunsch<DNASequence, NucleotideCompound>(query, target, new SimpleGapPenalty((short)0, (short)10), SubstitutionMatrixHelper.getNuc4_4());
		aligner.setAnchors(new int[] { -1, 2, -1} );
		assertEquals(String.format("A-CTTT%nACGTTT%n"), aligner.getPair().toString());
    }
    
    /**
     * @author Daniel Cameron 
     */
    @Test
	public void should_align_multiple_anchors() throws CompoundNotFoundException {
    	DNASequence query = new DNASequence("ACGT", AmbiguityDNACompoundSet.getDNACompoundSet());
		DNASequence target = new DNASequence("ATACGT", AmbiguityDNACompoundSet.getDNACompoundSet());
		NeedlemanWunsch<DNASequence, NucleotideCompound> aligner = new NeedlemanWunsch<DNASequence, NucleotideCompound>(query, target, new SimpleGapPenalty((short)0, (short)10), SubstitutionMatrixHelper.getNuc4_4());
		aligner.addAnchor(0, 0);
		aligner.addAnchor(1, 1);
		aligner.addAnchor(2, 2);
		aligner.addAnchor(3, 5);
		assertEquals(String.format("ACG--T%nATACGT%n"), aligner.getPair().toString());
    }
    
    /**
     * @author Daniel Cameron 
     */
    @Test
	public void anchors_should_not_change_score() throws CompoundNotFoundException {
    	DNASequence query = new DNASequence("ACGT", AmbiguityDNACompoundSet.getDNACompoundSet());
		DNASequence target = new DNASequence("ACGT", AmbiguityDNACompoundSet.getDNACompoundSet());
		NeedlemanWunsch<DNASequence, NucleotideCompound> aligner = new NeedlemanWunsch<DNASequence, NucleotideCompound>(query, target, new SimpleGapPenalty((short)5, (short)10), SubstitutionMatrixHelper.getNuc4_4());
		NeedlemanWunsch<DNASequence, NucleotideCompound> anchored = new NeedlemanWunsch<DNASequence, NucleotideCompound>(query, target, new SimpleGapPenalty((short)5, (short)10), SubstitutionMatrixHelper.getNuc4_4());
		anchored.addAnchor(0, 0);
		anchored.addAnchor(1, 1);
		anchored.addAnchor(2, 2);
		anchored.addAnchor(3, 3);
		assertEquals(aligner.getScore(), anchored.getScore(), PRECISION);
    }
    
    /**
     * @author Daniel Cameron 
     */
	@Test
	public void testAnchoredDNAAlignment() throws CompoundNotFoundException {
		DNASequence query = new DNASequence(  "ACGTACCGGTTTT", DNACompoundSet.getDNACompoundSet());
		DNASequence target = new DNASequence("TACGTCCGGTTACGTACGTT", DNACompoundSet.getDNACompoundSet());
		NeedlemanWunsch<DNASequence, NucleotideCompound> aligner = new NeedlemanWunsch<DNASequence, NucleotideCompound>(query, target, new SimpleGapPenalty((short)5, (short)2), SubstitutionMatrixHelper.getNuc4_4());
		assertEquals(String.format("-ACGTACCGGTT-------TT%nTACGT-CCGGTTACGTACGTT%n"), aligner.getPair().toString());
	}
	
	/**
	 * See issue #202 in github 
	 * @author Jose M Duarte
	 */
	@Test
	public void testIntOverflowBug() throws CompoundNotFoundException {

		SubstitutionMatrix<NucleotideCompound> matrix = SubstitutionMatrixHelper.getNuc4_4();
		SimpleGapPenalty gap = new SimpleGapPenalty();

		String str1 =
				"AGATATATCTGAAGCTTAAAGGGCAGTGACAATGGCTGGCTCGGTTAACGGGAATCATAGTGCTGTAGGACCTGGTATAAATTATGAGACGGTGTCTCAAGTGGATGAGTTCTGTAAAGCACTTAGAGGGAAAAGGCCGATCCATAGTATTTTGATAGCTAACAATGGAATGGCGGCTGTGAAGTTTATACGTAGTGTCAGAACATGGGCTTATGAAACATTTGGTACGGAAAAAGCCATATTGTTGGTGGGGATGGCAACCCCTGAAGACATGCGGATCAATGCGGAGCATATCAGAATCGCTGATCAGTTTGTTGAGGTTCCCGGAGGAACCAACAATAACAATTATGCTAACGTTCAGCTGATTGTGGAGATGGCTGAAGTAACACGCGTGGATGCAGTTTGGCCTGGTTGGGGTCATGCATCTGAAAACCCCGAATTACCTGATGCCCTAGATGCAAAAGGAATCATATTTCTTGGTCCTCCAGCATCTTCAATGGCAGCACTGGGAGATAAGATTGGTTCTTCGTTGATTGCACAAGCTGCTGATGTACCCACTCTGCCATGGAGTGGTTCCCATGTTAAAATACCTCCTAATAGCAACTTGGTAACCATCCCAGAGGAGATCTACCGGCAAGCATGTGTCTACACAACTGAAGAAGCGATTGCTAGCTGTCAAGTTGTCGGTTACCCAGCAATGATCAAAGCATCGTGGGGTGGTGGTGGTAAAGGAATCAGGAAGGTTCATAATGATGATGAGGTTAGGGCTCTATTCAAGCAAGTTCAGGGTGAGGTCCCAGGCTCACCAATATTCATAATGAAGGTTGCGTCACAGAGTCGGCATCTAGAGGTCCAGCTGCTCTGTGACAAGCATGGAAATGTTTCAGCTCTGCATAGCCGTGATTGTAGCGTCCAGAGAAGACATCAAAAGATCATAGAGGAGGGTCCAATTACTGTGGCTCCGCCAGAAACTGTCAAGAAACTTGAACAAGCAGCTAGAAGGTTGGCTAAGAGTGTTAACTATGTTGGAGCTGCTACTGTTGAGTATCTCTACAGTATGGACACTGGGGAGTACTACTTCTTAGAGCTTAACCCTCGCTTACAGGTTGAGCATCCTGTCACTGAGTGGATTGCCGAGATAAATCTTCCTGCTGCCCAAGTTGCTGTGGGGATGGGAATTCCTCTCTGGCAAATCCCTGAGATAAGACGGTTCTATGGAATAGAACATGGTGGAGGTTATGATTCTTGGCGAAAAACATCTGTTGTAGCCTTCCCTTTTGATTTTGATAAAGCTCAATCTATAAGGCCAAAAGGTCATTGTGTGGCTGTACGTGTGACAAGTGAGGATCCTGATGACGGGTTCAAACCAACCAGCGGTAGAGTTCAGGAGTTGAGTTTTAAGAGCAAGCCAAATGTGTGGGCGTACTTCTCTGTCAAGTCTGGTGGAGGCATCCACGAGTTCTCGGATTCCCAGTTTGGACATGTTTTTGCATTTGGGGAATCCAGAGCCCTGGCGATAGCGAATATGGTTCTTGGGCTAAAAGAAATTCAGATCCGTGGAGAAATTAGGACTAACGTTGACTACACGATCGACCTTTTACATGCTTCTGATTACCGTGATAACAAAATTCACACTGGTTGGTTGGATAGTAGGATTGCTATGCGGGTCAGAGCTGAGAGGCCTCCATGGTATCTCTCTGTTGTCGGCGGAGCTCTCTATAAAGCATCAGCGACCAGTGCTGCTGTGGTTTCAGATTACGTTGGTTATCTGGAGAAGGGGCAAATCCCTCCAAAGCATATATCTCTTGTACATTCTCAAGTGTCTCTGAATATTGAAGGAAGTAAATATACGATTGATGTAGTCCGGGGTGGATCAGGAACCTACAGGCTAAGAATGAACAAGTCAGAAGTGGTAGCAGAAATACACACTCTACGTGATGGAGGTCTGTTGATGCAGTTGGATGGCAAAAGCCATGTGATATATGCAGAGGAAGAAGCTGCAGGAACTCGTCTTCTCATTGATGGAAGAACTTGTTTGCTACAGAATGACCACGATCCATCAAAGTTAATGGCTGAGACACCGTGCAAGTTGATGAGGTATTTGATTTCCGACAACAGCAATATTGACGCTGATACGCCTTATGCCGAAGTTGAGGTCATGAAGATGTGCATGCCACTTCTTTCACCTGCTTCAGGAGTTATCCATTTTAAAATGTCTGAAGGACAAGCCATGCAGGCTGGTGAACTTATAGCCAATCTTGATCTTGATGATCCTTCTGCTGTAAGAAAGGCCGAACCCTTCCATGGAAGTTTCCCAAGATTAGGGCTTCCAACTGCAATATCCGGTAGAGTTCATCAGAGATGTGCCGCAACATTAAATGCTGCACGCATGATTCTTGCTGGCTATGAGCATAAAGTAGATGAGGTTGTTCAAGACTTACTTAATTGCCTTGATAGCCCTGAACTCCCATTTCTTCAGTGGCAAGAGTGCTTTGCAGTTCTGGCGACACGACTACCTAAAAATCTCAGGAACATGCTAGAATCAAAGTATAGGGAATTTGAGAGTATTTCCAGAAACTCTTTGACCACCGATTTCCCTGCCAAACTTTTAAAAGGCATTCTTGAGGCACATTTATCTTCTTGTGATGAGAAAGAGAGAGGTGCCCTTGAAAGGCTCATTGAACCATTGATGAGCCTTGCAAAATCTTATGAAGGTGGTAGAGAAAGTCATGCCCGTGTTATTGTTCATTCTCTCTTTGAAGAATATCTATCAGTAGAAGAATTATTCAATGATAACATGCTGGCTGATGTTATAGAACGCATGCGTCAGCTATACAAGAAAGATCTGTTGAAAATTGTGGATATAGTGCTCTCACACCAGGGCATAAAAAACAAAAACAAACTCGTTCTCCGGCTCATGGAGCAGCTTGTTTACCCTAATCCTGCTGCTTACAGAGATAAACTTATTCGATTCTCAACACTTAACCATACTAACTACTCTGAGTTGGCGCTCAAGGCGAGTCAATTACTTGAACAGACCAAACTAAGTGAGCTTCGTTCAAACATTGCTAGAAGCCTTTCAGAGTTAGAAATGTTTACAGAGGACGGAGAAAATATGGATACTCCCAAGAGGAAAAGTGCCATTAATGAAAGAATAGAAGATCTTGTAAGCGCATCTTTAGCTGTTGAAGACGCTCTCGTGGGACTATTTGACCATAGCGATCACACACTTCAAAGACGGGTTGTTGAGACTTATATTCGCAGATTATACCAGCCCTACGTCGTTAAAGATAGCGTGAGGATGCAGTGGCACCGTTCTGGTCTTCTTGCTTCCTGGGAGTTCCTAGAGGAGCATATGGAAAGAAAAAACATTGGCTTAGACGATCCCGACACATCTGAAAAAGGATTGGTTGAGAAGCGTAGTAAGAGAAAATGGGGGGCTATGGTTATAATCAAATCTTTGCAGTTTCTTCCAAGTATAATAAGTGCAGCATTGAGAGAAACAAAGCACAACGACTATGAAACTGCCGGAGCTCCTTTATCTGGCAATATGATGCACATTGCTATTGTGGGCATCAACAACCAGATGAGTCTGCTTCAGGACAGTGGGGATGAAGACCAAGCTCAGGAAAGAGTAAACAAGTTGGCCAAAATTCTTAAAGAGGAAGAAGTGAGTTCAAGCCTCTGTTCTGCCGGTGTTGGTGTAATCAGCTGTATAATTCAGCGAGATGAAGGACGAACACCCATGAGACATTCTTTCCATTGGTCGTTGGAGAAACAGTATTATGTAGAAGAGCCGTTGCTGCGTCATCTTGAACCTCCTCTGTCCATTTACCTTGAGTTGGATAAGCTGAAAGGATACTCAAATATACAATATACGCCTTCTCGAGATCGTCAATGGCATCTGTATACTGTTACAGACAAGCCAGTGCCAATCAAGAGGATGTTCCTGAGATCTCTTGTTCGACAGGCTACAATGAACGATGGATTTATATTGCAGCAAGGGCAGGATAAGCAGCTTAGCCAAACACTGATCTCCATGGCGTTTACGTCGAAATGTGTTCTGAGGTCTTTGATGGATGCCATGGAGGAACTGGAACTGAATGCCCATAATGCTGCAATGAAACCAGATCACGCACATATGTTTCTTTGCATATTGCGTGAGCAGCAGATAGATGATCTTGTGCCTTTCCCCAGGAGAGTTGAAGTGAATGCGGAGGATGAAGAAACTACAGTTGAAATGATCTTAGAAGAAGCAGCACGAGAGATACATAGATCTGTTGGAGTGAGAATGCATAGGTTGGGCGTGTGCGAGTGGGAAGTGCGGCTGTGGTTGGTGTCCTCTGGACTGGCATGTGGTGCTTGGAGGGTTGTGGTTGCAAACGTGACAGGCCGTACATGCACTGTCCACATATACCGAGAAGTTGAAACTCCTGGAAGAAACAGTTTAATCTACCACTCAATAACCAAGAAGGGACCTTTGCATGAAACACCAATCAGTGATCAATATAAGCCCCTGGGATATCTCGACAGGCAACGTTTAGCAGCAAGGAGGAGTAACACTACTTATTGCTATGACTTCCCGTTGGCATTTGGGACAGCCTTGGAACTGTTGTGGGCATCACAACACCCAGGAGTTAAGAAACCATATAAGGATACTCTGATCAATGTTAAAGAGCTTGTATTCTCAAAACCAGAAGGTTCTTCGGGTACATCTCTAGATCTGGTTGAAAGACCACCCGGTCTCAACGACTTTGGAATGGTTGCCTGGTGCCTAGATATGTCGACCCCAGAGTTTCCTATGGGGCGGAAACTTCTCGTGATTGCGAATGATGTCACCTTCAAAGCTGGTTCTTTTGGTCCTAGAGAGGACGCGTTTTTCCTTGCTGTTACTGAACTCGCTTGTGCCAAGAAGCTTCCCTTGATTTACTTGGCAGCAAATTCTGGTGCCCGACTTGGGGTTGCTGAAGAAGTCAAAGCCTGCTTCAAAGTTGGATGGTCGGATGAAATTTCCCCTGAGAATGGTTTTCAGTATATATACCTAAGCCCTGAAGACCACGAAAGGATTGGATCATCTGTCATTGCCCATGAAGTAAAGCTCTCTAGTGGGGAAACTAGGTGGGTGATTGATACGATCGTTGGCAAAGAAGATGGTATTGGTGTAGAGAACTTAACAGGAAGTGGGGCCATAGCGGGTGCTTACTCAAAGGCATACAATGAAACTTTTACTTTAACCTTTGTTAGTGGAAGAACGGTTGGAATTGGTGCTTATCTTGCCCGCCTAGGTATGCGGTGCATACAGAGACTTGATCAGCCGATCATCTTGACTGGCTTCTCTACACTCAACAAGTTACTTGGGCGTGAGGTCTATAGCTCTCACATGCAACTGGGTGGCCCGAAAATCATGGGCACAAATGGTGTTGTTCATCTTACAGTCTCAGATGATCTTGAAGGCGTATCAGCAATTCTCAACTGGCTCAGCTACATTCCTGCTTACGTGGGTGGTCCTCTTCCTGTTCTTGCCCCTTTAGATCCACCGGAGAGAATTGTGGAGTATGTCCCAGAGAACTCTTGCGACCCACGAGCGGCTATAGCTGGGGTCAAAGACAATACCGGTAAATGGCTTGGAGGTATCTTTGATAAAAATAGTTTCATTGAGACTCTTGAAGGCTGGGCAAGGACGGTAGTGACTGGTAGAGCCAAGCTCGGGGGAATACCCGTTGGAGTTGTTGCAGTTGAGACACAGACTGTCATGCAGATCATCCCAGCCGATCCTGGACAGCTTGACTCTCATGAAAGAGTGGTTCCGCAAGCAGGGCAAGTCTGGTTTCCTGATTCAGCGGCCAAGACTGCTCAAGCGCTTATGGATTTCAACCGGGAAGAGCTTCCATTGTTTATCCTAGCGAACTGGAGAGGGTTTTCAGGTGGGCAGAGAGATCTTTTCGAAGGAATACTTCAGGCAGGTTCAACTATAGTAGAAAATCTGAGAACCTATCGTCAGCCAGTGTTTGTGTACATCCCAATGATGGGAGAGCTGCGCGGTGGAGCGTGGGTTGTTGTTGACAGCCAGATAAATTCGGATTATGTTGAAATGTATGCTGATGAAACAGCTCGTGGAAATGTGCTTGAGCCAGAAGGGACAATAGAGATAAAATTTAGAACAAAAGAGCTATTAGAGTGCATGGGAAGGTTGGACCAGAAGCTAATCAGTCTGAAAGCAAAACTGCAAGATGCCAAGCAAAGCGAGGCCTATGCAAACATCGAGCTTCTCCAGCAACAGATTAAAGCCCGAGAGAAACAGCTTTTACCAGTTTATATCCAAATCGCCACCAAATTTGCAGAACTTCATGACACTTCCATGAGAATGGCTGCAAAGGGAGTGATCAAAAGTGTTGTGGAATGGAGCGGCTCGCGGTCCTTCTTCTACAAAAAGCTCAATAGGAGAATCGCTGAGAGCTCTCTTGTGAAAAACGTAAGAGAAGCATCTGGAGACAACTTAGCATATAAATCTTCAATGCGTCTGATTCAGGATTGGTTCTGCAACTCTGATATTGCAAAGGGGAAAGAAGAAGCTTGGACAGACGACCAAGTGTTCTTTACATGGAAGGACAATGTTAGTAACTACGAGTTGAAGCTGAGCGAGTTGAGAGCGCAGAAACTACTGAACCAACTTGCAGAGATTGGGAATTCCTCAGATTTGCAAGCTCTGCCACAAGGACTTGCTAATCTTCTAAACAAGGTGGAGCCGTCGAAAAGAGAAGAGCTGGTGGCTGCTATTCGAAAGGTCTTGGGTTGACTGA";
		String str2 =
				"TAAAGTCTTCGATATCAGTCAACCCAAGACCTTTCGAATAGCAGCCACCAGCTCTTCTCTTTTCGACGGCTCCACCTTGTTTAGAAGATTAGCA";
		//System.out.println("Lengths: "+str1.length()+" "+str2.length());

		DNASequence target = new DNASequence(str1,
				AmbiguityDNACompoundSet.getDNACompoundSet());

		DNASequence query = new DNASequence(str2,
				AmbiguityDNACompoundSet.getDNACompoundSet());

		NeedlemanWunsch<DNASequence, NucleotideCompound> aligner = 
				new NeedlemanWunsch<DNASequence, NucleotideCompound>(query, target, gap, matrix);

		
		//System.out.println("getScore: " + aligner.getScore());
		//System.out.println("getMaxScore: " + aligner.getMaxScore());
		//System.out.println("getMinScore: " + aligner.getMinScore());
		//System.out.println("getSimilarity: " + aligner.getSimilarity());
		
		assertTrue("Similarity must be positive, this must be an integer overflow bug!", aligner.getSimilarity()>0);
	} 
	
	/**
	 * The starting gap in a global alignment should not be penalised as a gap opening,
	 * other implementations of Needleman-Wunsch behave in that way (see for instance needle from 
	 * EMBOSS, http://www.ebi.ac.uk/Tools/psa/emboss_needle/)
	 * @author Jose M Duarte
	 * @throws CompoundNotFoundException
	 */
	@Test
	public void testGapAtStartIssue() throws CompoundNotFoundException {
		
		String str1 = "HHHAAAA";
		String str2 = "MAAAA";

		
		SubstitutionMatrix<AminoAcidCompound> matrix = SubstitutionMatrixHelper.getBlosum50();
		//System.out.println(matrix.toString());
		
		GapPenalty penalty = new SimpleGapPenalty(10, 1);

		ProteinSequence s1 = new ProteinSequence(str1);
		ProteinSequence s2 = new ProteinSequence(str2);

		NeedlemanWunsch<ProteinSequence,AminoAcidCompound> nw = 
				new NeedlemanWunsch<ProteinSequence,AminoAcidCompound>(s1,s2, penalty, matrix);
		
		//System.out.println(nw.getScoreMatrixAsString());
		//System.out.println("Score: "+nw.getScore());
		
		SequencePair<ProteinSequence, AminoAcidCompound> pair = nw.getPair();
		//System.out.println(pair.toString(100));
		assertEquals("The first position in target aligned seq should be a gap","-",pair.getCompoundInTargetAt(1).toString());

		assertEquals(String.format("HHHAAAA%n--MAAAA%n"), pair.toString());
	}
}
