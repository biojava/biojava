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
 * Created on August 11, 2010
 * Author: Mark Chapman
 */

package org.biojava.nbio.alignment.routines;

import org.biojava.nbio.alignment.SimpleGapPenalty;
import org.biojava.nbio.core.alignment.matrices.SubstitutionMatrixHelper;
import org.biojava.nbio.alignment.template.GapPenalty;
import org.biojava.nbio.core.alignment.template.SubstitutionMatrix;
import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.compound.AmbiguityDNACompoundSet;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;
import org.junit.Before;
import org.junit.Test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

public class GuanUberbacherTest {

	private static final double PRECISION = 0.00000001;
	
    private ProteinSequence query, target;
    private GapPenalty gaps;
    private SubstitutionMatrix<AminoAcidCompound> blosum62;
    private GuanUberbacher<ProteinSequence, AminoAcidCompound> alignment, self;

    @Before
    public void setup() throws CompoundNotFoundException { 
        query = new ProteinSequence("ARND");
        target = new ProteinSequence("RDG");
        gaps = new SimpleGapPenalty(10, 1);
        blosum62 = SubstitutionMatrixHelper.getBlosum62();
        alignment = new GuanUberbacher<ProteinSequence, AminoAcidCompound>(query, target, gaps, blosum62);
        self = new GuanUberbacher<ProteinSequence, AminoAcidCompound>(query, query, gaps, blosum62);
    }

    @Test
    public void testGuanUberbacher() {
        GuanUberbacher<ProteinSequence, AminoAcidCompound> gu =
                new GuanUberbacher<ProteinSequence, AminoAcidCompound>();
        gu.setQuery(query);
        gu.setTarget(target);
        gu.setGapPenalty(gaps);
        gu.setSubstitutionMatrix(blosum62);
        assertEquals(gu.getScore(), alignment.getScore(), PRECISION);
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
        assertEquals(-6, alignment.getScore(), PRECISION);
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
	public void should_align_shorter_query() throws CompoundNotFoundException {
    	DNASequence query = new DNASequence("A", AmbiguityDNACompoundSet.getDNACompoundSet());
		DNASequence target = new DNASequence("AT", AmbiguityDNACompoundSet.getDNACompoundSet());
		GuanUberbacher<DNASequence, NucleotideCompound> aligner = new GuanUberbacher<DNASequence, NucleotideCompound>(query, target, new SimpleGapPenalty((short)5, (short)2), SubstitutionMatrixHelper.getNuc4_4());
		assertEquals(String.format("A-%nAT%n"), aligner.getPair().toString());
    }
    /**
     * @author Daniel Cameron 
     */
    @Test
	public void should_align_shorter_target() throws CompoundNotFoundException {
    	DNASequence query = new DNASequence("AT", AmbiguityDNACompoundSet.getDNACompoundSet());
		DNASequence target = new DNASequence("A", AmbiguityDNACompoundSet.getDNACompoundSet());
		GuanUberbacher<DNASequence, NucleotideCompound> aligner = new GuanUberbacher<DNASequence, NucleotideCompound>(query, target, new SimpleGapPenalty((short)5, (short)2), SubstitutionMatrixHelper.getNuc4_4());
		assertEquals(String.format("AT%nA-%n"), aligner.getPair().toString());
    }
    /**
     * @author Daniel Cameron 
     */
    @Test
	public void should_align_multiple_cuts() throws CompoundNotFoundException {
    	DNASequence query = new DNASequence("AAT", AmbiguityDNACompoundSet.getDNACompoundSet());
		DNASequence target = new DNASequence("AATG", AmbiguityDNACompoundSet.getDNACompoundSet());
		GuanUberbacher<DNASequence, NucleotideCompound> aligner = new GuanUberbacher<DNASequence, NucleotideCompound>(query, target, new SimpleGapPenalty((short)0, (short)2), SubstitutionMatrixHelper.getNuc4_4());
		aligner.setCutsPerSection(2); // 3 bases with 2 cuts
		assertEquals(String.format("AAT-%nAATG%n"), aligner.getPair().toString());
    }
}
