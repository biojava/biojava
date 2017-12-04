package org.biojava.nbio.core.sequence.io.embl;

import org.junit.Assert;
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
        Assert.assertEquals("acaagatgccattgtcccccggcctcctgctgctg" +
                "ctgctctccggggccacggccaccgctgccctgcccctggagggtggccccaccggcc" +
                "gagacagcgagcatatgcaggaagcggcaggaataaggaaaagcagcctcctgactttcc" +
                "tcgcttggtggtttgagtggacctcccaggccagtgccgggcccctcataggagaggaagc" +
                "tcgggaggtggccaggcggcaggaaggcgcacccccccagcaatccgcgcgccgggacagaa" +
                "tgccctgcaggaacttcttctggaagaccttctcctcctgcaaataaaacctcacccatgaatgc" +
                "tcacgcaagtttaattacagacctgaa", emblRecord.getSequence());
    }


}
