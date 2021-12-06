package org.biojava.nbio.core.alignment.matrices;

import org.biojava.nbio.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava.nbio.core.sequence.compound.DNACompoundSet;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;


class SubstitutionMatrixTest {
    AminoAcidCompoundSet aaSet = new AminoAcidCompoundSet();
    DNACompoundSet dnaSet = DNACompoundSet.getDNACompoundSet();
    final short MATCH = 5;
    final short REPLACE = -10;
    SimpleSubstitutionMatrix<NucleotideCompound> sm = null;
    @BeforeEach
    void before (){
        sm = new SimpleSubstitutionMatrix<NucleotideCompound>(
            dnaSet, (short)MATCH, REPLACE);
    }

    @Test
    void createIdentityMatrix() {
        assertEquals(MATCH, sm.getMaxValue());
        assertEquals(REPLACE, sm.getMinValue());
        assertEquals(dnaSet, sm.getCompoundSet());
        short value = sm.getValue(dnaSet.getCompoundForString("T"),dnaSet.getCompoundForString("T"));
        assertEquals(MATCH, value);

        value = sm.getValue(dnaSet.getCompoundForString("T"),dnaSet.getCompoundForString("A"));
        assertEquals(REPLACE, value);
    }

    @Test
    void matrixDimensions(){
        int dnaSetSize = dnaSet.getAllCompounds().size();
        NucleotideCompound thy = dnaSet.getCompoundForString("T");
        assertEquals(dnaSetSize, sm.getColumn(thy).size());
        assertEquals(dnaSetSize, sm.getRow(thy).size());
    }

    @Test
    void getMatrixReturnsCopy(){
        
        short [][] matrix = sm.getMatrix();
        assertEquals(MATCH, matrix[0][0]);
        matrix [0][0]= 100; // new value doesn't affect internal matrix
        assertEquals(MATCH, sm.getMatrix()[0][0]);
    }

    @Test
    void matrixToString(){
        String asString = sm.toString();
        // description + 5*2 for ATCGNatcgn + 1 header + 1  for '-'
        assertEquals(13, asString.split("\\n").length);
        String header = asString.split("\\R")[1];
        assertTrue(header.replaceAll(" ","").matches("[ATCGatcgNn-]+"));
    }

    @Test
    void matrixAsString(){
        String asString = sm.getMatrixAsString();

        // 5*2 for ATCGNatcgn + 1 header + 1  for '-'
        assertEquals(12, asString.split("\\n").length);
        String header = asString.split("\\R")[0];
        assertTrue(header.replaceAll(" ","").matches("[ATCGatcgNn-]+"));
    }
}