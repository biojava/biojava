package org.biojava3.core.sequence;

import static org.hamcrest.CoreMatchers.is;
import static org.junit.Assert.assertThat;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;

import junit.framework.Assert;

import org.biojava3.core.exceptions.CompoundNotFoundError;
import org.biojava3.core.sequence.compound.AmbiguityDNACompoundSet;
import org.biojava3.core.sequence.compound.DNACompoundSet;
import org.biojava3.core.sequence.compound.NucleotideCompound;
import org.biojava3.core.sequence.storage.FourBitSequenceReader;
import org.biojava3.core.sequence.storage.SingleCompoundSequenceReader;
import org.biojava3.core.sequence.storage.TwoBitSequenceReader;
import org.biojava3.core.sequence.template.CompoundSet;
import org.biojava3.core.sequence.template.ProxySequenceReader;
import org.biojava3.core.sequence.template.SequenceMixin;
import org.biojava3.core.sequence.template.SequenceReader;
import org.biojava3.core.sequence.template.SequenceView;
import org.biojava3.core.sequence.transcription.Frame;
import org.biojava3.core.sequence.views.ComplementSequenceView;
import org.biojava3.core.sequence.views.ReversedSequenceView;
import org.junit.Test;

public class DNATest {

    private DNACompoundSet set = new DNACompoundSet();
    private AmbiguityDNACompoundSet ambiguity = new AmbiguityDNACompoundSet();

    @Test
    public void reverseComplement() {
        String s = getSeq().getInverse().getSequenceAsString();
        assertThat("Reversed Complemented sequence not as expected", s, is("GCAT"));
    }

    @Test
    public void complement() {
        String s = new ComplementSequenceView<NucleotideCompound>(getSeq()).getSequenceAsString();
        assertThat("Complemented sequence not as expected", s, is("TACG"));
    }

    @Test
    public void reverse() {
        SequenceView<NucleotideCompound> r = new ReversedSequenceView<NucleotideCompound>(getSeq());
        assertThat("Reversed sequence not as expected", r.getSequenceAsString(), is("CGTA"));
        assertThat("Base at 2 not right", r.getCompoundAt(2).toString(), is("G"));

        List<String> actual = new ArrayList<String>();
        List<String> expected = Arrays.asList("C", "G", "T", "A");
        for (NucleotideCompound c : r) {
            actual.add(c.toString());
        }
        assertThat("Iterator not behaving as expected", actual, is(expected));

        assertThat("Index of T not as expected", r.getIndexOf(set.getCompoundForString("T")), is(3));
    }

    @Test
    public void subSequence() {
        DNASequence seq = getSeq("ACGTGGC");
        SequenceView<NucleotideCompound> subSeq = seq.getSubSequence(2, 4);
        assertThat("Index 2 is the same as index 1 in the sub sequence",
                seq.getCompoundAt(2), is(subSeq.getCompoundAt(1)));
        assertThat("Length is equal to 3", subSeq.getLength(), is (3));
        assertThat("Index 4 is the same as index 3 in the sub sequence",
                seq.getCompoundAt(4), is(subSeq.getCompoundAt(3)));
        assertThat("Sub sequence contains only expected bases",
                subSeq.getSequenceAsString(), is("CGT"));
    }

    @Test
    public void translateToRna() {
        String s = getSeq("ATGGCGGCGCTGAGCGGT").getRNASequence().getSequenceAsString();
        assertThat("RNA as expected", s, is("AUGGCGGCGCUGAGCGGU"));
        String s2 = getSeq("ATGGCGGCGCTGAGCGGT").getRNASequence(Frame.TWO).getSequenceAsString();
        assertThat("RNA as expected", s2, is("UGGCGGCGCUGAGCGGU"));
    }

    @Test
    public void respectCase() {
        String s = "ATgc";
        DNASequence dna = getSeq(s);
        assertThat("Sequence does not remember casing", dna.getSequenceAsString(), is(s));
        assertThat("Reversed complement sequence forgets case",
                dna.getInverse().getSequenceAsString(), is("gcAT"));
    }

    @Test(expected = CompoundNotFoundError.class)
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

    @Test
    public void gc() {
        assertThat("GC content not as expected", SequenceMixin.countGC(getSeq("GCGC")), is(4));
        assertThat("GC content not as expected", getSeq("GCGC").getGCCount(), is(4));
        assertThat("GC content not as expected", SequenceMixin.countGC(getSeq("GAAC")), is(2));
        assertThat("GC content not as expected",
                SequenceMixin.countGC(getSeq("AATTTATATGAATTTATATGAATTTATATGAATTTATATGAATTTATATGAATTTATATGAATTTATATGAATTTATATGAATTTATATG")),
                is(9));
    }

    @Test
    public void at() {
        assertThat("AT content not as expected", SequenceMixin.countAT(getSeq("GCGC")), is(0));
        assertThat("AT content not as expected", SequenceMixin.countAT(getSeq("GCAT")), is(2));
        assertThat("AT content not as expected", SequenceMixin.countAT(getSeq("atAT")), is(4));
        assertThat("GC content not as expected",
                SequenceMixin.countAT(getSeq("AATTTATATGAATTTATATGAATTTATATGAATTTATATGAATTTATATGAATTTATATGAATTTATATGAATTTATATGAATTTATATG")),
                is(81));
    }

    @Test
    public void composition() {
        DNASequence seq = getSeq("ATTGGGCCCC");
        CompoundSet<NucleotideCompound> set = seq.getCompoundSet();
        Map<NucleotideCompound, Double> distribution = SequenceMixin.getDistribution(seq);
        assertThat("A distribution not as expected", distribution.get(set.getCompoundForString("A")), is(0.1));
        assertThat("T distribution not as expected", distribution.get(set.getCompoundForString("T")), is(0.2));
        assertThat("G distribution not as expected", distribution.get(set.getCompoundForString("G")), is(0.3));
        assertThat("C distribution not as expected", distribution.get(set.getCompoundForString("C")), is(0.4));
    }

    @Test
    public void twoBit() {
        String expected = "ATGCAACTGA";
        DNASequence seq = getSeq(expected);
        SequenceReader<NucleotideCompound> twoBitFromSeq =
                new TwoBitSequenceReader<NucleotideCompound>(seq);

        //being cheeky here & getting compound set from seq
        SequenceReader<NucleotideCompound> twoBitFromString =
                new TwoBitSequenceReader<NucleotideCompound>(expected, seq.getCompoundSet());

        assertThat("TwoBit from Sequence not as expected", twoBitFromSeq.getSequenceAsString(), is(expected));
        assertThat("TwoBit from String not as expected", twoBitFromString.getSequenceAsString(), is(expected));
    }

    @Test
    public void fourBit() {
        String expected = "ATGCAACTGA";
        DNASequence seq = getSeq(expected);
        SequenceReader<NucleotideCompound> bitFromSeq =
                new FourBitSequenceReader<NucleotideCompound>(seq);

        //being cheeky here & getting compound set from seq
        SequenceReader<NucleotideCompound> bitFromString =
                new FourBitSequenceReader<NucleotideCompound>(expected, seq.getCompoundSet());

        assertThat("FourBit from Sequence not as expected", bitFromSeq.getSequenceAsString(), is(expected));
        assertThat("FourBit from String not as expected", bitFromString.getSequenceAsString(), is(expected));
    }

    @Test(expected = IllegalStateException.class)
    public void badTwoBit() {
        DNASequence seq = getSeq();
        new TwoBitSequenceReader<NucleotideCompound>("ATNGC", seq.getCompoundSet());
    }

    @Test
    public void singleCompoundSequence() {
        CompoundSet<NucleotideCompound> cs = set;
        NucleotideCompound a = cs.getCompoundForString("A");
        NucleotideCompound n = cs.getCompoundForString("N");
        int length = 1000;

        ProxySequenceReader<NucleotideCompound> sr = new SingleCompoundSequenceReader<NucleotideCompound>(n, cs, length);
        DNASequence seq = new DNASequence(sr);

        int intCount = 0;
        int iteratorCount = 0;
        for (int i = 1; i <= seq.getLength(); i++) {
            if (seq.getCompoundAt(i).equals(n)) {
                intCount++;
            }
        }
        for (NucleotideCompound c : seq) {
            if (c.equals(n)) {
                iteratorCount++;
            }
        }

        assertThat("All positions from getCompoundAt() report N", intCount, is(length));
        assertThat("All positions from iterator report N", iteratorCount, is(length));
        assertThat("Non N compound reports -ve", seq.getIndexOf(a), is(-1));
        assertThat("Index of N compound reports 1", seq.getIndexOf(n), is(1));
        assertThat("Non N compound reports -ve", seq.getLastIndexOf(a), is(-1));
        assertThat("Last index of N compound reports length", seq.getLastIndexOf(n), is(length));
    }

    @Test
    public void kmerNonOverlap() {
        DNASequence d = new DNASequence("ATGTGCA");
        List<SequenceView<NucleotideCompound>> l =
                SequenceMixin.nonOverlappingKmers(d, 3);
        assertThat("Asserting we generate only 2 k-mers", l.size(), is(2));
        assertThat("Asserting first k-mer", l.get(0).getSequenceAsString(), is("ATG"));
        assertThat("Asserting second k-mer", l.get(1).getSequenceAsString(), is("TGC"));
    }

    @Test
    public void kmerOverlap() {
        DNASequence d = new DNASequence("ATGTT");
        List<SequenceView<NucleotideCompound>> l =
                SequenceMixin.overlappingKmers(d, 3);
        assertThat("Asserting we generate only 3 k-mers", l.size(), is(3));
        assertThat("Asserting first k-mer", l.get(0).getSequenceAsString(), is("ATG"));
        assertThat("Asserting second k-mer", l.get(1).getSequenceAsString(), is("TGT"));
        assertThat("Asserting second k-mer", l.get(2).getSequenceAsString(), is("GTT"));
    }

    @Test
    public void kmerOverlapExceedingSequenceLength() {
        DNASequence d = new DNASequence("ATGTT");
        List<SequenceView<NucleotideCompound>> l =
                SequenceMixin.overlappingKmers(d, 2);
        assertThat("Asserting we generate 4 k-mers", l.size(), is(4));
        assertThat("Asserting first k-mer", l.get(0).getSequenceAsString(), is("AT"));
        assertThat("Asserting second k-mer", l.get(2).getSequenceAsString(), is("GT"));
        assertThat("Asserting second k-mer", l.get(3).getSequenceAsString(), is("TT"));
    }

    @Test
    public void sequenceEquality() {
        DNASequence d = getSeq("ATGC");
        assertTrue("Asserting sequences are identical", SequenceMixin.sequenceEquality(d, d));
        Assert.assertFalse("Sequence identical but case different", SequenceMixin.sequenceEquality(d, getSeq("ATGc")));
        assertTrue("Asserting sequences are identical ignoring case", SequenceMixin.sequenceEqualityIgnoreCase(d, d));
        assertTrue("Asserting sequences are identical ignoring case & case different", SequenceMixin.sequenceEqualityIgnoreCase(d, getSeq("aTgC")));
        Assert.assertFalse("Sequence lengths differ", SequenceMixin.sequenceEquality(d, getSeq("ATG")));

        DNASequence bsr = new DNASequence(new TwoBitSequenceReader<NucleotideCompound>("ATGC", DNACompoundSet.getDNACompoundSet()));
        DNASequence bsrCI = new DNASequence(new TwoBitSequenceReader<NucleotideCompound>("ATGc", DNACompoundSet.getDNACompoundSet()));

        assertTrue("Asserting sequences are identical; backing stores differ", SequenceMixin.sequenceEquality(d, bsr));
        assertTrue("Asserting sequences are identical ignoring case; backing stores differ", SequenceMixin.sequenceEqualityIgnoreCase(d, bsrCI));
    }   

//  @Test
//  public void randomTwoBit() throws Exception {
//    int[] ar = new int[1000000];
//    Random r = new Random();
//    for(int i = 0; i < ar.length; i++) {
//      ar[i] = r.nextInt();
//    }
//
//    System.out.println(Runtime.getRuntime().freeMemory());
//    System.out.println(Runtime.getRuntime().totalMemory());
//    TwoBitArrayWorker<NucleotideCompound> worker =
//      new TwoBitArrayWorker<NucleotideCompound>(getSeq().getCompoundSet(), ar);
//    SequenceProxyLoader<NucleotideCompound> sbs =
//      new BitSequenceReader<NucleotideCompound>(worker, new AccessionID("barf"));
//
//    System.out.println(sbs.getLength());
//
//    System.out.println(Runtime.getRuntime().freeMemory());
//    System.out.println(Runtime.getRuntime().totalMemory());
//
//    List<NucleotideCompound> c = sbs.getAsList();
//
//    System.out.println(Runtime.getRuntime().freeMemory());
//    System.out.println(Runtime.getRuntime().totalMemory());
//
////    OutputStream os = new BufferedOutputStream(new FileOutputStream(new File("/tmp/random.fasta")));
////
////    List<DNASequence> seqs = Arrays.asList(new DNASequence(sbs, sbs.getCompoundSet()));
////    seqs.get(0).setAccession(sbs.getAccession());
////    FastaHeaderFormatInterface<DNASequence, NucleotideCompound> headerFormat =
////      new GenericFastaHeaderFormat<DNASequence, NucleotideCompound>();
////
////    FastaWriter<DNASequence, NucleotideCompound> writer =
////      new FastaWriter<DNASequence, NucleotideCompound>(os, seqs, headerFormat);
////
////    writer.process();
////
////    IOUtils.close(os);
//  }
    private DNASequence getSeq() {
        return getSeq(null);
    }

    private DNASequence getSeq(final String seq) {
        String target = (seq == null) ? "ATGC" : seq;
        return new DNASequence(target);
    }

    private void assertBaseEquivalence(
            CompoundSet<NucleotideCompound> compoundSet, String one, String two) {
        boolean equal = compoundSet.compoundsEquivalent(
                compoundSet.getCompoundForString(one),
                compoundSet.getCompoundForString(two));
        assertTrue(one + " & " + two + " should be equivalent bases", equal);
    }
}
