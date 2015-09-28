package org.biojava.nbio.core.search.io;

import org.biojava.nbio.core.alignment.template.SequencePair;
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
 * Designed by Paolo Pavan.
 * You may want to find my contacts on Github and LinkedIn for code info 
 * or discuss major changes.
 * https://github.com/paolopavan
 * 
 * @author Paolo Pavan
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
    
    Hsp uncompleteHsp = new BlastHspBuilder()
                .setPercentageIdentity(100.00/100)
                .setHspAlignLen(48)
                .setMismatchCount(0)
                .setHspGaps(0)
                .setHspQueryFrom(1)
                .setHspQueryTo(48)
                .setHspHitFrom(344)
                .setHspHitTo(391)
                .setHspEvalue(4e-19)
                .setHspBitScore(95.6)
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
    public void testHashCode() {
        System.out.println("hashCode");
        Hsp instance;
        int expResult;
        int result;
        
        instance = hspImpl;
        expResult = 71782805;
        result = instance.hashCode();
        assertEquals(expResult, result);
        
        instance = uncompleteHsp;
        expResult = 679;
        result = instance.hashCode();
        assertEquals(expResult, result);
        
        Hsp uncompleteHsp2 = new BlastHspBuilder()
                .setPercentageIdentity(100.00/100)
                .setHspAlignLen(48)
                .setMismatchCount(0)
                .setHspGaps(0)
                .setHspQueryFrom(1)
                .setHspQueryTo(48)
                .setHspHitFrom(344)
                .setHspHitTo(391)
                .setHspEvalue(4e-19)
                .setHspBitScore(95.6)
                .createBlastHsp();
        
        assertFalse(uncompleteHsp.hashCode() == uncompleteHsp2.hashCode());
    }

    /**
     * Test of equals method, of class Hsp.
     */
    @Test
    public void testEquals() {
        System.out.println("equals");
        Object o;
        o = new BlastHspBuilder()
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
        Hsp instance = hspImpl;
        
        assertEquals(o, instance);
        
        // example of Hsp retrieved from uncomplete report. The result is null, hance false
        // (Those HSP may come from a tabular format, for example)
        o = new BlastHspBuilder()
                .setPercentageIdentity(100.00/100)
                .setHspAlignLen(48)
                .setMismatchCount(0)
                .setHspGaps(0)
                .setHspQueryFrom(1)
                .setHspQueryTo(48)
                .setHspHitFrom(344)
                .setHspHitTo(391)
                .setHspEvalue(4e-19)
                .setHspBitScore(95.6)
                .createBlastHsp();
        
        // At now, check will return false as it could be not determined
        assertFalse(uncompleteHsp.equals(o));
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
