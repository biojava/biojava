package org.biojava.nbio.core.sequence.io.embl;

import org.junit.Test;

import java.io.File;
/**
 * This class should test the parser of EmblParser class
 * @author Noor Aldeen Al Mbaidin
 */
public class EmblParserTest {

    @Test(expected = NullPointerException.class)
    public void givenNullFileParameterWhenCreateEmblParserThenthrowException(){
        File file = null;
        EmblParser emblParser = new EmblParser(file);
    }

    @Test
    public void test(){
        File file = new File("/home/pslpt219/Desktop/Homo.dat");
        EmblParser emblParser = new EmblParser(file);
        emblParser.parse();
    }


}
