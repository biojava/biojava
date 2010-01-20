package org.biojava3.core;

//import static org.hamcrest.CoreMatchers.is;
//import static org.junit.Assert.assertThat;

//import org.biojava3.core.sequence.template.Sequence;
//import org.junit.Test;

public class SequenceEditTests {

//  @Test
//  public void addSequenceToEnd() {
//    Sequence<DNA> seq = getSeq();
//    String s = seq.edit(new Add(seq.length(), "A")).toString();
//    assertThat(s, is("ATGCA"));
//  }
//
//  @Test
//  public void addSequenceToStart() {
//    Sequence<DNA> seq = getSeq();
//    String s = seq.edit(new Add(1, "A")).toString();
//    assertThat(s, is("AATGC"));
//  }
//
//  @Test
//  public void addSequenceToMiddle() {
//    Sequence<DNA> seq = getSeq();
//    String s = seq.edit(new Add(2, "A")).toString();
//    assertThat(s, is("ATAGC"));
//  }
//
//  @Test
//  public void replaceSequence() {
//    Sequence<DNA> seq = getSeq();
//    String s = seq.edit(new Replace(2, "A")).toString();
//    assertThat(s, is("AAGC"));
//
//    s = seq.edit(new Replace(2, "AAA")).toString();
//    assertThat(s, is("AAAA"));
//  }
//
//  @Test(expected=OutOfBoundsException.class)
//  public void replaceOutOfBounds() {
//    getSeq().edit(new Replace(10, "A"));
//  }
//
//  @Test(expected=OutOfBoundsException.class)
//  public void deleteOutOfBounds() {
//    getSeq().edit(new Delete(10));
//  }
//
//  @Test
//  public void deleteSequenceFromEnd() {
//    Sequence<DNA> seq = getSeq();
//    String s = seq.edit(new Delete(seq.getLength())).toString();
//    assertThat(s, is("ATG"));
//
//    int length = seq.getLength();
//    s = seq.edit(new Delete(length-1, length)).toString();
//    assertThat(s, is("AT"));
//  }
//
//  @Test
//  public void deleteSequenceRange() {
//    Sequence<DNA> seq = getSeq();
//    String s = seq.edit(new Delete(2, 4)).toString();
//    assertThat(s, is("A"));
//  }
//
//  public Sequence<DNA> getSeq(final String seq) {
//    String target = (seq == null) ? "ATGC" : seq;
//    return new DNASequence(target);
//  }

}
