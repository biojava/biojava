package org.biojava3.core.sequence.compound;

import static org.junit.Assert.assertEquals;

import java.util.ArrayList;
import java.util.List;

import org.junit.Test;

public class AmbiguityDNACompoundTest {

  private AmbiguityDNACompoundSet set = AmbiguityDNACompoundSet.getDNACompoundSet();

  @Test
  public void testAmbiguity() {
    NucleotideCompound actual = set.getAmbiguity(getCompounds("M","V"));
    assertEquals("Checking M & G = V", getCompounds("V")[0], actual);
  }

  @Test
  public void testBasicAmbiguity() {
    NucleotideCompound actual = set.getAmbiguity(getCompounds("A","C"));
    assertEquals("Checking A & C = M", getCompounds("M")[0], actual);
  }

  private NucleotideCompound[] getCompounds(String... compoundStrings) {
    List<NucleotideCompound> c = new ArrayList<NucleotideCompound>();
    for(String s: compoundStrings) {
      c.add(set.getCompoundForString(s));
    }
    return c.toArray(new NucleotideCompound[0]);
  }
}
