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
 * Created on Jun 10, 2010
 * Author: Mark 
 *
 */

package org.biojava3.alignment;

import static org.junit.Assert.*;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.InputStreamReader;

import org.biojava3.core.sequence.compound.AminoAcidCompound;
import org.biojava3.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava3.core.sequence.compound.CodonCompound;
import org.biojava3.core.sequence.compound.DNACompoundSet;
import org.biojava3.core.sequence.compound.NucleotideCompound;
import org.biojava3.core.sequence.template.CompoundSet;
import org.junit.Ignore;
import org.junit.Test;

public class SimpleSubstitutionMatrixTest {

    @Test(expected=FileNotFoundException.class)
    public void testSimpleSubstitutionMatrixNotFound() throws FileNotFoundException {
        new SimpleSubstitutionMatrix<AminoAcidCompound>(new AminoAcidCompoundSet(), new File("blosum63.txt"));
    }

    @Test(expected=IllegalStateException.class)
    public void testSimpleSubstitutionMatrixNull() {
        new SimpleSubstitutionMatrix<AminoAcidCompound>();
    }

    @Ignore
    @Test(expected=ClassCastException.class)
    public void testSimpleSubstitutionMatrixWrong() {
        SimpleSubstitutionMatrix.Default.set(new AminoAcidCompoundSet(), new InputStreamReader(
                SimpleSubstitutionMatrix.class.getResourceAsStream("/blosum62.txt") ), "blosum62");
        new SimpleSubstitutionMatrix<CodonCompound>();
        // TODO why does this not cause ClassCastException? loses typing at runtime?
    }

    @Test()
    public void testSimpleSubstitutionMatrix() throws FileNotFoundException {
        CompoundSet<AminoAcidCompound> compoundSet = new AminoAcidCompoundSet();
        SimpleSubstitutionMatrix.Default.set(compoundSet, new InputStreamReader(
                SimpleSubstitutionMatrix.class.getResourceAsStream("/blosum62.txt") ), "blosum62");
        SimpleSubstitutionMatrix<AminoAcidCompound> matrix = new SimpleSubstitutionMatrix<AminoAcidCompound>();
        assertEquals(matrix.getCompoundSet(), compoundSet);
        assertEquals(matrix.getName(), "blosum62");
        assertEquals(matrix.getMaxValue(), 11);
        assertEquals(matrix.getMinValue(), -4);
    }

    @Test
    public void testSimpleSubstitutionMatrixSubstitutionMatrixOfC() {
        SimpleSubstitutionMatrix<AminoAcidCompound> matrix1 =
                new SimpleSubstitutionMatrix<AminoAcidCompound>(new AminoAcidCompoundSet(), new InputStreamReader(
                SimpleSubstitutionMatrix.class.getResourceAsStream("/blosum62.txt") ), "blosum62"),
                matrix2 = new SimpleSubstitutionMatrix<AminoAcidCompound>(matrix1);
        assertEquals(matrix2.getCompoundSet(), matrix1.getCompoundSet());
        assertEquals(matrix2.getName(), "blosum62");
        assertEquals(matrix2.getMaxValue(), 11);
        assertEquals(matrix2.getMinValue(), -4);
    }

    @Test
    public void testSimpleSubstitutionMatrixCompoundSetOfCStringString() {
        SimpleSubstitutionMatrix<NucleotideCompound> matrix = new SimpleSubstitutionMatrix<NucleotideCompound>(
                new DNACompoundSet(), "# Test\nA C G T\nA 5 0 0 0\nC 0 5 0 0\nG 0 0 5 0\nT 0 0 0 1\n", "SimpleDNA");
        assertEquals(matrix.getMatrixAsString().substring(2,9), "A C G T");
    }

    @Test
    public void testSimpleSubstitutionMatrixCompoundSetOfCShortShort() {
        SimpleSubstitutionMatrix<AminoAcidCompound> matrix =
                new SimpleSubstitutionMatrix<AminoAcidCompound>(new AminoAcidCompoundSet(), (short) 5, (short) 1);
        assertEquals(matrix.getName(), "IDENTITY_5_1");
    }

    @Test
    public void testSetDescription() {
        SimpleSubstitutionMatrix<AminoAcidCompound> matrix =
                new SimpleSubstitutionMatrix<AminoAcidCompound>(new AminoAcidCompoundSet(), new InputStreamReader(
                SimpleSubstitutionMatrix.class.getResourceAsStream("/blosum62.txt") ), "blosum62");
        assertEquals(matrix.getDescription().substring(0, 2), "# ");
        matrix.setDescription("blah");
        assertEquals(matrix.getDescription().substring(0, 2), "bl");
    }

    @Test
    public void testSetName() {
        SimpleSubstitutionMatrix<AminoAcidCompound> matrix =
                new SimpleSubstitutionMatrix<AminoAcidCompound>(new AminoAcidCompoundSet(), new InputStreamReader(
                SimpleSubstitutionMatrix.class.getResourceAsStream("/blosum62.txt") ), "blosum62");
        assertEquals(matrix.getName(), "blosum62");
        matrix.setName("blah");
        assertEquals(matrix.getName(), "blah");
    }

    @Test
    public void testToString() {
        SimpleSubstitutionMatrix<NucleotideCompound> matrix = new SimpleSubstitutionMatrix<NucleotideCompound>(
                new DNACompoundSet(), "# Test\nA C  G T\nA 5  0 0 0\nC 0 5 0 0 \nG 0 0 5 0\n T 0 0 0 1\n", "DNAtest");
        assertEquals(matrix.toString().substring(0,6), "# Test");
        String nl = System.getProperty("line.separator");
        assertEquals(matrix.toString(), "# Test" + nl + "  A C G T" + nl + "A 5 0 0 0" + nl + "C 0 5 0 0" + nl +
                "G 0 0 5 0" + nl + "T 0 0 0 1" + nl);
    }

}
