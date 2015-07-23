/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.biojava.nbio.core.search.io;

import org.biojava.nbio.alignment.template.SequencePair;
import org.biojava.nbio.core.search.io.blast.BlastHspBuilder;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;
import org.junit.Ignore;

/**
 *
 * @author pavanpa
 */
public class HspTest {
    
    Hsp hspImpl = new BlastHspBuilder()
                .setHspNum(1)
                .setHspBitScore(377.211)
                .setHspEvalue(8.04143e-093)
                .setHspQueryFrom(1)
                .setHspQueryTo(224)
                .setHspHitFrom(1035)
                .setHspHitTo(811)
                .setHspQueryFrame(-1)
                .setHspIdentity(213)
                .setHspPositive(213)
                .setHspGaps(5)
                .setHspAlignLen(227)
                .setHspQseq("CTGACGACAGCCATGCACCACCTGTCTCGACTTTCCCCCGAAGGGCACCTAATGTATCTCTACCTCGTTAGTCGGATGTCAAGACCTGGTAAGGTTTTTTCGCGTATCTTCGAATTAAACCACATACTCCACTGCTTGTGCGG-CCCCCGTCAATTCCTTTGAGTTTCAACCTTGCGGCCGTACTCCC-AGGTGGA-TACTTATTGTGTTAACTCCGGCACGGAAGG")
                .setHspHseq("CTGACGACAACCATGCACCACCTGTCTCAACTTTCCCC-GAAGGGCACCTAATGTATCTCTACTTCGTTAGTTGGATGTCAAGACCTGGTAAGGTT-CTTCGCGTTGCTTCGAATTAAACCACATACTCCACTGCTTGTGCGGGCCCCCGTCAATTCCTTTGAGTTTCAACCTTGCGGTCGTACTCCCCAGGTGGATTACTTATTGTGTTAACTCCGGCACAGAAGG")
                .setHspIdentityString("||||||||| |||||||||||||||||| ||||||||| |||||||||||||||||||||||| |||||||| |||||||||||||||||||||||  |||||||  |||||||||||||||||||||||||||||||||||| |||||||||||||||||||||||||||||||||| ||||||||| ||||||| |||||||||||||||||||||||| |||||")
                .createBlastHsp();
    
    public HspTest() {
    }
    
    @BeforeClass
    public static void setUpClass() {
    }
    
    @AfterClass
    public static void tearDownClass() {
    }
    
    @Before
    public void setUp() {
    }
    
    @After
    public void tearDown() {
    }

    /**
     * Test of hashCode method, of class Hsp.
     */
    @Test
    @Ignore public void testHashCode() {
        System.out.println("hashCode");
        Hsp instance = null;
        int expResult = 0;
        int result = instance.hashCode();
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of equals method, of class Hsp.
     */
    @Test
    @Ignore public void testEquals() {
        System.out.println("equals");
        Object o = null;
        Hsp instance = null;
        boolean expResult = false;
        boolean result = instance.equals(o);
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of getAlignment method, of class Hsp.
     */
    @Test
    public void testGetAlignment() {
        System.out.println("getAlignment");
        
        SequencePair<DNASequence, NucleotideCompound> aln = hspImpl.getAlignment();
        
        StringBuilder s = new StringBuilder();
        s.append(hspImpl.getHspQseq());
        s.append(String.format("%n"));
        s.append(hspImpl.getHspHseq());
        s.append(String.format("%n"));
        
        String expResult = s.toString();
        
        String result = aln.toString();
        
        assertEquals(expResult, result);
    }
    
}
