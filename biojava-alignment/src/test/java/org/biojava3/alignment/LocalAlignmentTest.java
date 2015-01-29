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
 * Created on December 10, 2013
 * Author: Daniel Cameron
 */
package org.biojava3.alignment;

import org.biojava3.alignment.Alignments.PairwiseSequenceAlignerType;
import org.biojava3.alignment.template.PairwiseSequenceAligner;
import org.biojava3.alignment.template.SubstitutionMatrix;
import org.biojava3.core.exceptions.CompoundNotFoundException;
import org.biojava3.core.sequence.DNASequence;
import org.biojava3.core.sequence.compound.DNACompoundSet;
import org.biojava3.core.sequence.compound.NucleotideCompound;
import org.junit.Test;

import static org.junit.Assert.assertEquals;

/**
 * @author Daniel Cameron
 *
 */
public class LocalAlignmentTest {

	private static final double PRECISION = 0.00000001;
	
	@Test
	public void shouldAllowZeroLengthMatches() throws CompoundNotFoundException { 
        DNASequence query = new DNASequence("C", DNACompoundSet.getDNACompoundSet());
        DNASequence target = new DNASequence("A", DNACompoundSet.getDNACompoundSet());
        SubstitutionMatrix<NucleotideCompound> matrix = SubstitutionMatrixHelper.getNuc4_4();
        SimpleGapPenalty gapP = new SimpleGapPenalty((short)5, (short)2);
        PairwiseSequenceAligner<DNASequence, NucleotideCompound> result = Alignments.getPairwiseAligner(query, target, PairwiseSequenceAlignerType.LOCAL, gapP, matrix);
        assertEquals(0, result.getScore(), PRECISION);
        assertEquals(0, result.getProfile().getLength());
	}
}
