package org.biojava3.core.sequence;

import static org.junit.Assert.assertEquals;

import org.biojava3.core.sequence.compound.NucleotideCompound;
import org.biojava3.core.sequence.template.CompoundSet;
import org.biojava3.core.sequence.template.SequenceView;
import org.junit.Test;

public class SequenceViewTests {

  @Test
  public void testGetCompoundAt() {
    SequenceView<NucleotideCompound> s = new DNASequence("ATGC").getSubSequence(1, 4);
    assertEquals("Compound @ 1", s.getCompoundAt(1).toString(), "A");
    assertEquals("Compound @ 3", s.getCompoundAt(3).toString(), "G");
    assertEquals("Compound @ 3", s.getSubSequence(2,3).getCompoundAt(1).toString(), "T");
  }

  public void testLastIndexOf() {
    SequenceView<NucleotideCompound> s = new DNASequence("ATGC").getSubSequence(1, 4);
    CompoundSet<NucleotideCompound> cs = s.getCompoundSet();
    assertEquals("Last index of ", s.getLastIndexOf(cs.getCompoundForString("C")), 4);
    s = new DNASequence("GAAAAAAAAG").getSubSequence(0, 10);
    assertEquals("Last index of G is 10",
        s.getLastIndexOf(cs.getCompoundForString("G")), 10);
    assertEquals("Last index of G is 5",
        s.getSubSequence(5, 10).getLastIndexOf(cs.getCompoundForString("G")), 5);
  }
}
