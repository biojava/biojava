package org.biojava3.core.sequence;

import static org.junit.Assert.assertEquals;

import org.biojava3.core.sequence.compound.NucleotideCompound;
import org.biojava3.core.sequence.edits.Edit;
import org.biojava3.core.sequence.template.Compound;
import org.biojava3.core.sequence.template.Sequence;
import org.junit.Test;

public class EditSequenceTest {

  @Test
  public void substitute() {
    DNASequence seq = new DNASequence("ACGT");
    assertSeq(new Edit.Substitute<NucleotideCompound>("T", 2).edit(seq), "ATGT");
    assertSeq(new Edit.Substitute<NucleotideCompound>("TT", 2).edit(seq), "ATTT");
    assertSeq(new Edit.Substitute<NucleotideCompound>("T", 1).edit(seq), "TCGT");
    assertSeq(new Edit.Substitute<NucleotideCompound>("TTC", 2).edit(seq), "ATTC");
  }

  @Test(expected=IndexOutOfBoundsException.class)
  public void badSubstitute() {
    new Edit.Substitute<NucleotideCompound>("AAAA", 4).edit(new DNASequence("ACGT"));
  }

  @Test
  public void delete() {
    DNASequence seq = new DNASequence("ACGT");
    assertSeq(new Edit.Delete<NucleotideCompound>(1).edit(seq), "CGT");
    assertSeq(new Edit.Delete<NucleotideCompound>(4).edit(seq), "ACG");
    assertSeq(new Edit.Delete<NucleotideCompound>(2,3).edit(seq), "AT");

    //disabling this test, because we can't create a CompoundSet if we have no sequences...
    // assertSeq(new Edit.Delete<NucleotideCompound>(1,4).edit(seq), "");
  }

  @Test
  public void insert() {
    DNASequence seq = new DNASequence("ACGT");
    assertSeq(new Edit.Insert<NucleotideCompound>("TT", 1).edit(seq), "TTACGT");
    assertSeq(new Edit.Insert<NucleotideCompound>("TT", 2,3).edit(seq), "ACTTGT");
    assertSeq(new Edit.Insert<NucleotideCompound>("TT", 3,4).edit(seq), "ACGTTT");
    assertSeq(new Edit.Insert<NucleotideCompound>("A", 4).edit(seq), "ACGTA");

    //Original BioJava example
    assertSeq(
        new Edit.Insert<NucleotideCompound>("atgga", 3,4).edit(new DNASequence("gataca")),
        "gatatggaaca"
    );
  }

  private void assertSeq(Sequence<? extends Compound> seq, String expected) {
    assertEquals("Asserting sequence "+expected, expected, seq.getSequenceAsString());
  }
}
