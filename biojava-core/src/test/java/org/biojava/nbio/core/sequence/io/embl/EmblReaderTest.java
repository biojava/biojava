package org.biojava.nbio.core.sequence.io.embl;

import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import java.io.File;
import java.io.IOException;

/**
 * This class should test the parser of EmblReader class
 * @since 5.0.0
 * @author Noor Aldeen Al Mbaidin
 */
public class EmblReaderTest {


    @Test(expected = NullPointerException.class)
    public void givenNullFileParameterWhenProcessEmblFileThenThrowException() throws IOException {
        File file = null;
        EmblReader.process(file);

    }

    @Test(expected = IllegalArgumentException.class)
    public void givenDirectoryWhenProcessEmblFileThenThrowException() throws IOException {
        File file = new File("./src/test/resources");
        EmblReader.process(file);
    }

    @Test
    public void givenAnEmilFileWhenProcessEmilFileThanTheSequenceShouldReturnAsExpected() throws IOException {
        File file = new File("./src/test/resources/test.embl");
        EmblRecord emblRecord = EmblReader.process(file);
        Assert.assertEquals("acaagatgccattgtcccccggcctcctgctgctg" +
                "ctgctctccggggccacggccaccgctgccctgcccctggagggtggccccaccggcc" +
                "gagacagcgagcatatgcaggaagcggcaggaataaggaaaagcagcctcctgactttcc" +
                "tcgcttggtggtttgagtggacctcccaggccagtgccgggcccctcataggagaggaagc" +
                "tcgggaggtggccaggcggcaggaaggcgcacccccccagcaatccgcgcgccgggacagaa" +
                "tgccctgcaggaacttcttctggaagaccttctcctcctgcaaataaaacctcacccatgaatgc" +
                "tcacgcaagtttaattacagacctgaa",emblRecord.getSequence());
    }




}
