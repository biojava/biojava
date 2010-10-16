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

package org.biojava3.alignment.routines;

import static org.junit.Assert.*;

import org.biojava3.alignment.SimpleGapPenalty;
import org.biojava3.alignment.SubstitutionMatrixHelper;
import org.biojava3.alignment.template.GapPenalty;
import org.biojava3.alignment.template.SubstitutionMatrix;
import org.biojava3.core.sequence.ProteinSequence;
import org.biojava3.core.sequence.compound.AminoAcidCompound;
import org.junit.Before;
import org.junit.Test;

public class GuanUberbacherTest {

    private ProteinSequence query, target;
    private GapPenalty gaps;
    private SubstitutionMatrix<AminoAcidCompound> blosum62;
    private GuanUberbacher<ProteinSequence, AminoAcidCompound> alignment, self;

    @Before
    public void setup() {
        query = new ProteinSequence("ARND");
        target = new ProteinSequence("RDG");
        gaps = new SimpleGapPenalty((short) 10, (short) 1);
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
        assertEquals(gu.getScore(), alignment.getScore());
    }

    @Test
    public void testGetComputationTime() {
        assertTrue(alignment.getComputationTime() > 0);
        assertTrue(self.getComputationTime() > 0);
    }

    @Test
    public void testGetProfile() {
        assertEquals(alignment.getProfile().toString(), String.format("ARND%n-RDG%n"));
        assertEquals(self.getProfile().toString(), String.format("ARND%nARND%n"));
    }

    @Test
    public void testGetMaxScore() {
        assertEquals(alignment.getMaxScore(), 21);
        assertEquals(self.getMaxScore(), 21);
    }

    @Test
    public void testGetMinScore() {
        assertEquals(alignment.getMinScore(), -27);
        assertEquals(self.getMinScore(), -28);
    }

    @Test
    public void testGetScore() {
        assertEquals(alignment.getScore(), -6);
        assertEquals(self.getScore(), 21);
    }

    @Test
    public void testGetPair() {
        assertEquals(alignment.getPair().toString(), String.format("ARND%n-RDG%n"));
        assertEquals(self.getPair().toString(), String.format("ARND%nARND%n"));
    }

}
