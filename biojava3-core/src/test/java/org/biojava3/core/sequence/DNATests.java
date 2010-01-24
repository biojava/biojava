package org.biojava3.core.sequence;

import static org.hamcrest.CoreMatchers.is;
import static org.junit.Assert.assertThat;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.biojava3.core.exceptions.CompoundNotFoundError;
import org.biojava3.core.sequence.compound.AmbiguityDNACompoundSet;
import org.biojava3.core.sequence.compound.DNACompoundSet;
import org.biojava3.core.sequence.compound.NucleotideCompound;
import org.biojava3.core.sequence.template.CompoundSet;
import org.biojava3.core.sequence.template.SequenceView;
import org.junit.Test;

public class DNATests {

  private DNACompoundSet set = new DNACompoundSet();
  private AmbiguityDNACompoundSet ambiguity = new AmbiguityDNACompoundSet();

  @Test
  public void reverseComplement() {
    String s = getSeq().getReverseComplement().getSequenceAsString();
    assertThat("Reversed Complemented sequence not as expected", s, is("GCAT"));
  }

  @Test
  public void complement() {
    String s = getSeq().getComplement().getSequenceAsString();
    assertThat("Complemented sequence not as expected", s, is("TACG"));
  }

  @Test
  public void reverse() {
    SequenceView<NucleotideCompound> r = getSeq().getReverse();
    assertThat("Reversed sequence not as expected", r.getSequenceAsString(), is("CGTA"));
    assertThat("Base at 2 not right", r.getCompoundAt(2).toString(), is("G"));

    List<String> actual = new ArrayList<String>();
    List<String> expected = Arrays.asList("C","G", "T", "A");
    for(NucleotideCompound c: r) {
      actual.add(c.toString());
    }
    assertThat("Iterator not behaving as expected", actual, is(expected));

    assertThat("Index of T not as expected", r.getIndexOf(set.getCompoundForString("T")), is(3));
  }

//
//  @Test
//  public void translateToRna() {
//    String s = getSeq("ATGGCGGCGCTGAGCGGT").toRNA().getSequenceAsString()();
//    assertThat("RNA as expected", s, is("AUGGCGGCGCUGAGCGGU"));
//  }
//

  @Test
  public void respectCase() {
    String s = "ATgc";
    DNASequence dna = getSeq(s);
    assertThat("Sequence does not remember casing", dna.getSequenceAsString(), is(s));

    StringBuilder reverse = new StringBuilder(s);
    assertThat("A reversed sequence does not remember case",
        dna.getReverse().getSequenceAsString(), is(reverse.reverse().toString()));

    assertThat("Reversed complement sequence forgets case",
        dna.getReverseComplement().getSequenceAsString(), is("gcAT"));
  }

  @Test(expected=CompoundNotFoundError.class)
  public void bogusSequence() {
    getSeq("ATGCx");
  }

  @Test
  public void basesEqual() {
    boolean equal = set.compoundsEqual(set.getCompoundForString("a"),
        set.getCompoundForString("A"));
    assertTrue("a & A should be equal bases", equal);
  }

  @Test
  public void basesEquivalent() {
    assertBaseEquivalence(ambiguity, "N", "A");
    assertBaseEquivalence(ambiguity, "G", "K");
    assertBaseEquivalence(ambiguity, "V", "C");
    assertBaseEquivalence(ambiguity, "W", "T");

    assertBaseEquivalence(ambiguity, "n", "A");
    assertBaseEquivalence(ambiguity, "g", "K");
    assertBaseEquivalence(ambiguity, "v", "C");
    assertBaseEquivalence(ambiguity, "w", "T");
  }

//  @Test
//  public void gc() {
//    assertThat("GC content not as expected", getSeq("GCGC").getGC(), is(1.0));
//    assertThat("GC content not as expected", getSeq("GAAC").getGC(), is(0.5));
//    assertThat("GC content not as expected",
//        getSeq("AATTTATATGAATTTATATGAATTTATATGAATTTATATGAATTTATATGAATTTATATGAATTTATATGAATTTATATGAATTTATATG").gc(),
//        is(0.1)
//    );
//  }

  private DNASequence getSeq() {
    return getSeq(null);
  }

  private DNASequence getSeq(final String seq) {
    String target = ( seq == null ) ? "ATGC" : seq;
    return new DNASequence(target);
  }

  private void assertBaseEquivalence(
      CompoundSet<NucleotideCompound> compoundSet, String one, String two) {
    boolean equal = compoundSet.compoundsEquivalent(
        compoundSet.getCompoundForString(one),
        compoundSet.getCompoundForString(two));
    assertTrue(one+" & "+two+" should be equivalent bases", equal);
  }

}
