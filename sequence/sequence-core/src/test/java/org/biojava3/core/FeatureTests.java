package org.biojava3.core;

import static org.junit.Assert.*;
import static org.hamcrest.CoreMatchers.*;

import org.biojava3.core.sequence.Compound;
import org.biojava3.core.sequence.Sequence;
import org.biojava3.core.sequence.impl.GenericCompound;
import org.biojava3.core.sequence.impl.GenericSequence;
import org.junit.Test;

public class FeatureTests {

  @Test
  public void testFeatureCostruction() {

  }

  public Sequence<? extends Compound> getSeq() {
    return new GenericSequence("ATGC");
  }

}
