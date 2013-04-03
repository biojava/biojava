package org.biojava3.core.sequence.views;

import static org.hamcrest.CoreMatchers.is;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertThat;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import org.biojava3.core.sequence.RNASequence;
import org.biojava3.core.sequence.compound.NucleotideCompound;
import org.biojava3.core.sequence.template.SequenceView;
import org.junit.Test;

public class WindowViewTests {

  @Test
  public void basicWindow() {
    RNASequence rna = new RNASequence("AUGCCU");
    WindowedSequence<NucleotideCompound> window = new WindowedSequence<NucleotideCompound>(rna, 3);

    Iterator<SequenceView<NucleotideCompound>> iter = window.iterator();
    assertTrue("hasNext() returns true", iter.hasNext());

    int count = 0;
    for(SequenceView<NucleotideCompound> c: window) {
      count++;
      if(count == 0) {
        String extracted = c.getCompoundAt(0).toString() + c.getCompoundAt(1).toString() + c.getCompoundAt(2).toString();
        assertEquals("Checking codon string", "AUG", extracted);
      }
    }
    assertThat("Windowed iterator broken", count, is(2));
  }

  @Test
  public void reaminderWindow() {
    RNASequence rna = new RNASequence("AUGCC");
    WindowedSequence<NucleotideCompound> window = new WindowedSequence<NucleotideCompound>(rna, 3);
    List<SequenceView<NucleotideCompound>> list = new ArrayList<SequenceView<NucleotideCompound>>();
    for(SequenceView<NucleotideCompound> c: window) {
      list.add(c);
    }
    assertThat("First window is size 3", list.get(0).getLength(), is(3));
    assertThat("Only 1 window", list.size(), is(1));
  }
}
