package org.biojava3.core.sequence;

import static org.junit.Assert.assertEquals;

import org.biojava3.core.sequence.compound.NucleotideCompound;
import org.biojava3.core.sequence.template.CompoundSet;
import org.biojava3.core.sequence.template.SequenceView;
import org.junit.Test;

public class SequenceViewTest {

  @Test
  public void testGetCompoundAt() {
    SequenceView<NucleotideCompound> s = new DNASequence("ATGC").getSubSequence(1, 4);
    assertEquals("Compound @ 1", s.getCompoundAt(1).toString(), "A");
    assertEquals("Compound @ 3", s.getCompoundAt(3).toString(), "G");
    assertEquals("Compound @ 3", s.getSubSequence(2,3).getCompoundAt(1).toString(), "T");
  }

  @Test
  public void testLastIndexOf() {
    SequenceView<NucleotideCompound> s = new DNASequence("ATGC").getSubSequence(1, 4);
    CompoundSet<NucleotideCompound> cs = s.getCompoundSet();
    assertEquals("Last index of ", 4, s.getLastIndexOf(cs.getCompoundForString("C")));
    
    s = new DNASequence("GAAAAAAAAG").getSubSequence(1, 10);
    assertEquals("Last index of G is 10", 10,
        s.getLastIndexOf(cs.getCompoundForString("G")));
    assertEquals("Last index of G is 5", 6,
        s.getSubSequence(5, 10).getLastIndexOf(cs.getCompoundForString("G")));
  }

  @Test
  public void testInverse() {
    SequenceView<NucleotideCompound> s = new DNASequence("ATGC").getSubSequence(2, 3).getInverse();
    assertEquals("Reversed complementing view", s.getSequenceAsString(), "CA");
  }
}
