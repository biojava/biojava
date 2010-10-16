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

package org.biojava3.alignment;

import static org.junit.Assert.*;

import java.io.File;
import java.io.FileNotFoundException;

import org.biojava3.alignment.template.SubstitutionMatrix;
import org.biojava3.core.sequence.compound.AminoAcidCompound;
import org.biojava3.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava3.core.sequence.compound.DNACompoundSet;
import org.biojava3.core.sequence.compound.NucleotideCompound;
import org.junit.Ignore;
import org.junit.Test;

public class SimpleSubstitutionMatrixTest {

    @Test(expected=FileNotFoundException.class)
    public void testSimpleSubstitutionMatrixNotFound() throws FileNotFoundException {
        new SimpleSubstitutionMatrix<AminoAcidCompound>(AminoAcidCompoundSet.getAminoAcidCompoundSet(),
                new File("blosum63.txt"));
    }

    @Ignore // TODO why does this not cause ClassCastException? loses typing at runtime?
    @Test(expected=ClassCastException.class)
    public void testSimpleSubstitutionMatrixWrong() {
        new SimpleSubstitutionMatrix<NucleotideCompound>();
    }

    @Test()
    public void testSimpleSubstitutionMatrix() {
        SubstitutionMatrix<AminoAcidCompound> matrix = new SimpleSubstitutionMatrix<AminoAcidCompound>();
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
        SubstitutionMatrix<AminoAcidCompound> matrix = new SimpleSubstitutionMatrix<AminoAcidCompound>();
        assertEquals(matrix.getDescription().substring(0, 2), "# ");
        matrix.setDescription("blah");
        assertEquals(matrix.getDescription().substring(0, 2), "bl");
    }

    @Test
    public void testSetName() {
        SubstitutionMatrix<AminoAcidCompound> matrix = new SimpleSubstitutionMatrix<AminoAcidCompound>();
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

}
