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
 */
package org.biojava.nbio.core.sequence.io.embl;

import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import java.io.File;
import java.io.IOException;

/**
 * This class should test the parser of EmblReader class
 *
 * @author Noor Aldeen Al Mbaidin
 * @since 5.0.0
 */
public class EmblReaderTest {

    private String sequence;

    @Before
    public void initObjects() {
        sequence = "aaacaaaccaaatatggattttattgtagccatatttgctctgtttgttattagctcattcacaattacttcca" +
                "caaatgcagttgaagcttctactcttcttgacataggtaacctgagtcggagcagttttcctcgtggcttcatctttggtgctggatcttcagcatac" +
                "caatttgaaggtgcagtaaacgaaggcggtagaggaccaagtatttgggataccttcacccataaatatccagaaaaaataagggatggaagcaatgcaga" +
                "catcacggttgaccaatatcaccgctacaaggaagatgttgggattatgaaggatcaaaatatggattcgtatagattctcaatctcttggccaagaatactcc" +
                "caaagggaaagttgagcggaggcataaatcacgaaggaatcaaatattacaacaaccttatcaacgaactattggctaacggtatacaaccatttgtaactcttttt" +
                "cattgggatcttccccaagtcttagaagatgagtatggtggtttcttaaactccggtgtaataaatgattttcgagactatacggatctttgcttcaaggaatttgga" +
                "gatagagtgaggtattggagtactctaaatgagccatgggtgtttagcaattctggatatgcactaggaacaaatgcaccaggtcgatgttcggcctccaacgtggccaa" +
                "gcctggtgattctggaacaggaccttatatagttacacacaatcaaattcttgctcatgcagaagctgtacatgtgtataagactaaataccaggcatatcaaaagggaaa" +
                "gataggcataacgttggtatctaactggttaatgccacttgatgataatagcataccagatataaaggctgccgagagatcacttgacttccaatttggattgtttatggaac" +
                "aattaacaacaggagattattctaagagcatgcggcgtatagttaaaaaccgattacctaagttctcaaaattcgaatcaagcctagtgaatggttcatttgattttattggtat" +
                "aaactattactcttctagttatattagcaatgccccttcacatggcaatgccaaacccagttactcaacaaatcctatgaccaatatttcatttgaaaaacatgggatacc" +
                "cttaggtccaagggctgcttcaatttggatatatgtttatccatatatgtttatccaagaggacttcgagatcttttgttacatattaaaaataaatataacaatcctgcaatt" +
                "ttcaatcactgaaaatggtatgaatgaattcaacgatgcaacacttccagtagaagaagctcttttgaatacttacagaattgattactattaccgtcacttatactacattcgt" +
                "tctgcaatcagggctggctcaaatgtgaagggtttttacgcatggtcatttttggactgtaatgaatggtttgcaggctttactgttcgttttggattaaactttgtagattaga" +
                "aagatggattaaaaaggtaccctaagctttctgcccaatggtacaagaactttctcaaaagaaactagctagtattattaaaagaactttgtagtagattacagtacatcgtttg" +
                "aagttgagttggtgcacctaattaaataaaagaggttactcttaacatatttttaggccattcgttgtgaagttgttaggctgttatttctattatactatgttgtagtaataa" +
                "gtgcattgttgtaccagaagctatgatcataactataggttgatccttcatgtatcagtttgatgttgagaatactttgaattaaaagtctttttttatttttttaaaaaaaaaa" +
                "aaaaaaaaaaaaaaaaaaa";
    }

    @Test(expected = NullPointerException.class)
    public void givenNullFileParameterWhenProcessEmblFileThenThrowException() throws IOException {
        File file = new File(this.getClass().getResource(null).getFile());
        EmblReader.process(file);

    }

    @Test(expected = IllegalArgumentException.class)
    public void givenDirectoryWhenProcessEmblFileThenThrowException() throws IOException {
        File file = new File(this.getClass().getResource("/")
                .getPath());
        EmblReader.process(file);
    }

    @Test
    public void givenAnEmilFileWhenProcessEmilFileThanTheSequenceShouldReturnAsExpected() throws IOException {
        File file = new File(this.getClass().getResource("/test.embl").getFile());
        EmblRecord emblRecord = EmblReader.process(file);
        Assert.assertEquals(sequence, emblRecord.getSequence());
    }


}
