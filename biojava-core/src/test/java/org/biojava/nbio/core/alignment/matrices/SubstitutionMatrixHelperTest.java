package org.biojava.nbio.core.alignment.matrices;

import org.biojava.nbio.core.alignment.template.SubstitutionMatrix;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompoundSet;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

class SubstitutionMatrixHelperTest {
    AminoAcidCompoundSet aaSet = new AminoAcidCompoundSet();


    @Test
    void getMatrixFromAAINDEX() {
        SubstitutionMatrix<AminoAcidCompound> aaIndex = SubstitutionMatrixHelper.getMatrixFromAAINDEX("ALTS910101");
        assertNotNull(aaIndex);
        assertEquals(-30, aaIndex.getValue(aaSet.getCompoundForString("R"),
                aaSet.getCompoundForString("A")));
    }

    @Test
    void getIdentity() {
        SubstitutionMatrix<AminoAcidCompound> identityMatrix = SubstitutionMatrixHelper.getIdentity();
        final String standard20 = "ARNDCQEGHILKMFPSTWYV";
        for (AminoAcidCompound from : aaSet.getAllCompounds()) {
            if (!standard20.contains(from.getShortName())) {
                continue;
            }
            for (AminoAcidCompound to : aaSet.getAllCompounds()) {
                if (!standard20.contains(to.getShortName())) {
                    continue;
                }
                if (from.equals(to)) {
                    assertEquals(1, identityMatrix.getValue(from, to));
                } else {
                    assertEquals(-10000, identityMatrix.getValue(from, to));
                }
            }
        }
    }

    @Test
    void getBlosum100() {
        assertNotNull(SubstitutionMatrixHelper.getBlosum100());
    }

    @Test
    void getBlosum30() {
        assertNotNull(SubstitutionMatrixHelper.getBlosum30());
    }

    @Test
    void getBlosum35() {
        assertNotNull(SubstitutionMatrixHelper.getBlosum35());
    }

    @Test
    void getBlosum40() {
        assertNotNull(SubstitutionMatrixHelper.getBlosum40());
    }

    @Test
    void getBlosum45() {
        assertNotNull(SubstitutionMatrixHelper.getBlosum45());
    }

    @Test
    void getBlosum50() {
        assertNotNull(SubstitutionMatrixHelper.getBlosum50());
    }

    @Test
    void getBlosum55() {
        assertNotNull(SubstitutionMatrixHelper.getBlosum55());
    }

    @Test
    void getBlosum60() {
        assertNotNull(SubstitutionMatrixHelper.getBlosum60());
    }

    @Test
    void getBlosum62() {
        SubstitutionMatrix<AminoAcidCompound> blosum62 = SubstitutionMatrixHelper.getBlosum62();
        assertNotNull(blosum62);
        AminoAcidCompound trypt = aaSet.getCompoundForString("W");
        assertEquals(11, blosum62.getValue(trypt, trypt));
    }

    @Test
    void getBlosum65() {
        assertNotNull(SubstitutionMatrixHelper.getBlosum65());
    }

    @Test
    void getBlosum70() {
        assertNotNull(SubstitutionMatrixHelper.getBlosum70());
    }

    @Test
    void getBlosum75() {
        assertNotNull(SubstitutionMatrixHelper.getBlosum75());
    }

    @Test
    void getBlosum80() {
        assertNotNull(SubstitutionMatrixHelper.getBlosum80());
    }

    @Test
    void getBlosum85() {
        assertNotNull(SubstitutionMatrixHelper.getBlosum85());
    }

    @Test
    void getBlosum90() {
        assertNotNull(SubstitutionMatrixHelper.getBlosum90());
    }

    @Test
    void getGonnet250() {
        assertNotNull(SubstitutionMatrixHelper.getGonnet250());
    }

    @Test
    void getNuc4_2() {
        assertNotNull(SubstitutionMatrixHelper.getNuc4_2());
    }

    @Test
    void getNuc4_4() {
        assertNotNull(SubstitutionMatrixHelper.getNuc4_4());
    }

    @Test
    void getPAM250() {
        assertNotNull(SubstitutionMatrixHelper.getPAM250());
    }

    @Test
    void getAminoAcidSubstitutionMatrix() {
        assertNotNull(SubstitutionMatrixHelper.getAminoAcidSubstitutionMatrix("blosum62"));
        assertNotNull(SubstitutionMatrixHelper.getAminoAcidSubstitutionMatrix("DAYM780301"));
    }

    @Test
    void unknownMatrixReturnsNull() {
        assertNull( SubstitutionMatrixHelper.getAminoAcidSubstitutionMatrix("?????"));
    }
}