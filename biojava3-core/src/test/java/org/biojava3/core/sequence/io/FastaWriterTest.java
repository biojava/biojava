package org.biojava3.core.sequence.io;


import static org.hamcrest.CoreMatchers.is;
import static org.junit.Assert.assertThat;

import java.io.ByteArrayOutputStream;

import org.biojava3.core.sequence.AccessionID;
import org.biojava3.core.sequence.DNASequence;
import org.junit.Test;

public class FastaWriterTest {

  @Test
  public void writeBasicFasta() throws Exception {
    String id         = "Example";
    String dnaLineOne = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
    String dnaLineTwo = "T";

    DNASequence s = new DNASequence(dnaLineOne+dnaLineTwo);
    s.setAccession(new AccessionID(id));
    ByteArrayOutputStream baos = new ByteArrayOutputStream();
    FastaWriterHelper.writeSequence(baos, s);

    String actual = new String(baos.toByteArray());
    String expected = ">"+id+"\n"+dnaLineOne+"\n"+dnaLineTwo+"\n";

    assertThat("Writer not as expected", actual, is(expected));
  }

  @Test
  public void writeFastaEqualToLineLength() throws Exception {
    String id  = "Example";
    String dna = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAT";

    DNASequence s = new DNASequence(dna);
    s.setAccession(new AccessionID(id));
    ByteArrayOutputStream baos = new ByteArrayOutputStream();
    FastaWriterHelper.writeSequence(baos, s);

    String actual = new String(baos.toByteArray());
    String expected = ">"+id+"\n"+dna+"\n";

    assertThat("Writer not as expected", actual, is(expected));
  }

}
