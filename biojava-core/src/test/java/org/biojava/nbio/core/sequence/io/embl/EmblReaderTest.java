package org.biojava.nbio.core.sequence.io.embl;

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

    private EmblReader emblReader;

    @Before
    public void setUp() {
        emblReader = new EmblReader();
    }

    @Test(expected = NullPointerException.class)
    public void givenNullFileParameterWhenProcessEmblFileThenThrowException() throws IOException {
        File file = null;
        emblReader.process(file);

    }

    @Test(expected = IllegalArgumentException.class)
    public void givenDirectoryWhenProcessEmblFileThenThrowException() throws IOException {
        File file = new File("./src/test/resources");
        emblReader.process(file);
    }




}
