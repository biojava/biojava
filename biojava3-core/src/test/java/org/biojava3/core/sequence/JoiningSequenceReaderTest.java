package org.biojava3.core.sequence;

import static org.junit.Assert.assertEquals;

import org.biojava3.core.sequence.compound.NucleotideCompound;
import org.biojava3.core.sequence.storage.JoiningSequenceReader;
import org.biojava3.core.sequence.template.SequenceMixin;
import org.junit.Test;

public class JoiningSequenceReaderTest {

  @SuppressWarnings("unchecked")
  @Test
  public void canScan() {
    JoiningSequenceReader<NucleotideCompound> seq =
      new JoiningSequenceReader<NucleotideCompound>(
          new DNASequence("AAAA"),
          new DNASequence("GGG"),
          new JoiningSequenceReader<NucleotideCompound>(new DNASequence("A"), new DNASequence("C")),
          new DNASequence("TT"),
          new DNASequence("C")
    );

    String expected = "AAAAGGGACTTC";

    StringBuilder builderByIndex = new StringBuilder();
    for(int i = 1; i <= seq.getLength(); i++) {
      builderByIndex.append(seq.getCompoundAt(i));
    }

    StringBuilder builderByIterator = SequenceMixin.toStringBuilder(seq);

    assertEquals("Index builder", expected, builderByIndex.toString());
    assertEquals("Iterator builder", expected, builderByIterator.toString());
  }

  @SuppressWarnings("unchecked")
  @Test
  public void empty() {
    JoiningSequenceReader<NucleotideCompound> seq =
      new JoiningSequenceReader<NucleotideCompound>(
          new DNASequence(""),
          new DNASequence(""),
          new DNASequence("A"),
          new DNASequence("")
      );
    assertEquals("Testing empty sequences", "A", seq.getSequenceAsString());
  }
}
