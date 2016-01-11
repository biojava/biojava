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
 * Created on June 10, 2010
 * Author: Mark Chapman
 */

package org.biojava.nbio.alignment;

import org.biojava.nbio.core.alignment.matrices.SimpleSubstitutionMatrix;
import org.biojava.nbio.core.alignment.template.SubstitutionMatrix;
import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava.nbio.core.sequence.compound.DNACompoundSet;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;
import org.junit.Test;

import java.io.File;
import java.io.FileNotFoundException;

import static org.junit.Assert.assertEquals;

public class SimpleSubstitutionMatrixTest {

    @Test(expected=FileNotFoundException.class)
    public void testSimpleSubstitutionMatrixNotFound() throws FileNotFoundException {
        new SimpleSubstitutionMatrix<AminoAcidCompound>(AminoAcidCompoundSet.getAminoAcidCompoundSet(),
                new File("blosum63.txt"));
    }

    @Test
    public void test() throws CompoundNotFoundException {
        NucleotideCompound A = new DNASequence("A").getCompoundAt(1);
        NucleotideCompound a = new DNASequence("a").getCompoundAt(1);
        NucleotideCompound c = new DNASequence("c").getCompoundAt(1);
        SubstitutionMatrix<NucleotideCompound> matrix = new SimpleSubstitutionMatrix<NucleotideCompound>(DNACompoundSet.getDNACompoundSet(), (short)1, (short)0);
        assertEquals(1, (matrix.getValue(A, A)));
        assertEquals(1, (matrix.getValue(a, a)));
        assertEquals(1, (matrix.getValue(A, a)));
        assertEquals(0, (matrix.getValue(a, c)));
    }

    @Test()
    public void testSimpleSubstitutionMatrix() {
        SubstitutionMatrix<AminoAcidCompound> matrix = SimpleSubstitutionMatrix.getBlosum62();
        assertEquals(matrix.getCompoundSet(), AminoAcidCompoundSet.getAminoAcidCompoundSet());
        assertEquals(matrix.getName(), "blosum62");
        assertEquals(matrix.getMaxValue(), 11);
        assertEquals(matrix.getMinValue(), -4);
    }

    @Test
    public void testSimpleSubstitutionMatrixCompoundSetOfCStringString() {
        DNACompoundSet dnacs = DNACompoundSet.getDNACompoundSet();
        SubstitutionMatrix<NucleotideCompound> dnaTest = new SimpleSubstitutionMatrix<NucleotideCompound>(dnacs,
                "# Test\nA C G T\nA 5 0 0 0\nC 0 5 0 0\nG 0 0 5 0\nT 0 0 0 1\n", "DNA Test");
        short[][] matrix = dnaTest.getMatrix();
        assertEquals(matrix[1][1], 5);
        assertEquals(matrix[3][3], 1);
        assertEquals(matrix[3][1], 0);
        assertEquals(dnaTest.getMatrixAsString().substring(2,9), "A C G T");
        assertEquals(dnaTest.getValue(dnacs.getCompoundForString("G"), dnacs.getCompoundForString("G")), 5);
        assertEquals(dnaTest.getValue(dnacs.getCompoundForString("A"), dnacs.getCompoundForString("G")), 0);
    }

    @Test
    public void testSimpleSubstitutionMatrixCompoundSetOfCShortShort() {
        SubstitutionMatrix<AminoAcidCompound> matrix = new SimpleSubstitutionMatrix<AminoAcidCompound>(
                AminoAcidCompoundSet.getAminoAcidCompoundSet(), (short) 5, (short) 1);
        assertEquals(matrix.getName(), "IDENTITY_5_1");
    }

    @Test
    public void testSetDescription() {
        SubstitutionMatrix<AminoAcidCompound> matrix = SimpleSubstitutionMatrix.getBlosum62();
        assertEquals(matrix.getDescription().substring(0, 2), "# ");
        matrix.setDescription("blah");
        assertEquals(matrix.getDescription().substring(0, 2), "bl");
    }

    @Test
    public void testSetName() {
        SubstitutionMatrix<AminoAcidCompound> matrix = SimpleSubstitutionMatrix.getBlosum62();
        assertEquals(matrix.getName(), "blosum62");
        matrix.setName("blah");
        assertEquals(matrix.getName(), "blah");
    }

    @Test
    public void testToString() {
        SubstitutionMatrix<NucleotideCompound> matrix = new SimpleSubstitutionMatrix<NucleotideCompound>(
                DNACompoundSet.getDNACompoundSet(),
                "# Test\nA C  G T\nA 5  0 0 0\nC 0 5 0 0 \nG 0 0 5 0\n    T 0     0 0 1\n", "DNAtest");
        assertEquals(matrix.toString().substring(0,6), "# Test");
        assertEquals(matrix.toString(),
                String.format("# Test%n  A C G T%nA 5 0 0 0%nC 0 5 0 0%nG 0 0 5 0%nT 0 0 0 1%n"));
    }
    /*
     * Author: Daniel Cameron
     */
    @Test
    public void testCaseEquivalence() {
    	DNACompoundSet dnacs = DNACompoundSet.getDNACompoundSet();
        SubstitutionMatrix<NucleotideCompound> dnaTest = new SimpleSubstitutionMatrix<NucleotideCompound>(dnacs,
                "# Test\nA C G T\nA 5 0 0 0\nC 0 5 0 0\nG 0 0 5 0\nT 0 0 0 1\n", "DNA Test");
        @SuppressWarnings("unused")
		short[][] matrix = dnaTest.getMatrix();
        assertEquals(dnaTest.getValue(dnacs.getCompoundForString("G"), dnacs.getCompoundForString("g")), 5);
        assertEquals(dnaTest.getValue(dnacs.getCompoundForString("A"), dnacs.getCompoundForString("g")), 0);
        assertEquals(dnaTest.getValue(dnacs.getCompoundForString("g"), dnacs.getCompoundForString("G")), 5);
        assertEquals(dnaTest.getValue(dnacs.getCompoundForString("g"), dnacs.getCompoundForString("A")), 0);
    }

}
