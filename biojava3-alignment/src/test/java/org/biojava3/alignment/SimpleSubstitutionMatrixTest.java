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

import java.io.FileNotFoundException;
import java.io.InputStreamReader;

import org.biojava3.core.sequence.compound.AminoAcidCompound;
import org.biojava3.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava3.core.sequence.compound.CodonCompound;
import org.biojava3.core.sequence.template.CompoundSet;
import org.junit.Test;

public class SimpleSubstitutionMatrixTest {

    @Test(expected=IllegalStateException.class)
    public void testSimpleSubstitutionMatrixNull() {
        new SimpleSubstitutionMatrix<AminoAcidCompound>();
    }

    @Test(expected=ClassCastException.class)
    public void testSimplSubstitutionMatrixWrong() {
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
        fail("Not yet implemented");
    }

    @Test
    public void testSimpleSubstitutionMatrixCompoundSetOfCFile() {
        fail("Not yet implemented");
    }

    @Test
    public void testSimpleSubstitutionMatrixCompoundSetOfCReaderString() {
        fail("Not yet implemented");
    }

    @Test
    public void testSimpleSubstitutionMatrixCompoundSetOfCStringString() {
        fail("Not yet implemented");
    }

    @Test
    public void testSimpleSubstitutionMatrixCompoundSetOfCShortShort() {
        fail("Not yet implemented");
    }

    @Test
    public void testSetDescription() {
        fail("Not yet implemented");
    }

    @Test
    public void testSetName() {
        fail("Not yet implemented");
    }

    @Test
    public void testToString() {
        fail("Not yet implemented");
    }

}
