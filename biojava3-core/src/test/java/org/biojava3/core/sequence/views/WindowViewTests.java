package org.biojava3.core.sequence.views;

import static org.hamcrest.CoreMatchers.is;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertThat;
import static org.junit.Assert.assertTrue;

import java.util.Iterator;
import java.util.List;

import org.biojava3.core.sequence.RNASequence;
import org.biojava3.core.sequence.compound.NucleotideCompound;
import org.junit.Test;

public class WindowViewTests {

  @Test
  public void basicWindow() {
    RNASequence rna = new RNASequence("AUGCCU");
    WindowedSequence<NucleotideCompound> window = new WindowedSequence<NucleotideCompound>(rna, 3);

    Iterator<List<NucleotideCompound>> iter = window.iterator();
    assertTrue("hasNext() returns true", iter.hasNext());

    int count = 0;
    for(List<NucleotideCompound> c: window) {
      count++;
      if(count == 0) {
        String extracted = c.get(0).getBase() + c.get(1).getBase() + c.get(2).getBase();
        assertEquals("Checking codon string", "AUG", extracted);
      }
    }
    assertThat("Windowed iterator broken", count, is(2));
  }

}
