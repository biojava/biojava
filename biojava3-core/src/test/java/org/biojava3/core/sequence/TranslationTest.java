package org.biojava3.core.sequence;

import static org.biojava3.core.sequence.io.util.IOUtils.close;
import static org.hamcrest.CoreMatchers.is;
import static org.junit.Assert.assertThat;

import java.io.InputStream;

import org.biojava3.core.sequence.compound.AminoAcidCompound;
import org.biojava3.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava3.core.sequence.compound.DNACompoundSet;
import org.biojava3.core.sequence.io.DNASequenceCreator;
import org.biojava3.core.sequence.io.FastaReader;
import org.biojava3.core.sequence.io.GenericFastaHeaderParser;
import org.biojava3.core.sequence.io.IUPACParser;
import org.biojava3.core.sequence.io.ProteinSequenceCreator;
import org.biojava3.core.sequence.io.util.ClasspathResource;
import org.biojava3.core.sequence.template.Sequence;
import org.biojava3.core.sequence.transcription.RNAToAminoAcidTranslator;
import org.biojava3.core.sequence.transcription.TranscriptionEngine;
import org.junit.Assert;
import org.junit.BeforeClass;
import org.junit.Test;

public class TranslationTest {

  private static DNACompoundSet dnaCs = DNACompoundSet.getDNACompoundSet();
  private static AminoAcidCompoundSet aaCs = AminoAcidCompoundSet.getAminoAcidCompoundSet();

  private static DNASequence brca2Dna;
  private static Sequence<AminoAcidCompound> brca2Pep;

  @BeforeClass
  public static void parseSequences() {
    InputStream cdsIs = new ClasspathResource(
        "org/biojava3/core/sequence/BRCA2-cds.fasta").getInputStream();
    InputStream pepIs = new ClasspathResource(
        "org/biojava3/core/sequence/BRCA2-peptide.fasta").getInputStream();

    try {
      FastaReader<DNASequence> dnaReader = new FastaReader<DNASequence>(cdsIs,
          new GenericFastaHeaderParser(), new DNASequenceCreator(dnaCs));
      brca2Dna = dnaReader.process().iterator().next();
      FastaReader<ProteinSequence> pReader = new FastaReader<ProteinSequence>(
          pepIs, new GenericFastaHeaderParser(), new ProteinSequenceCreator(
              aaCs));
      brca2Pep = pReader.process().iterator().next();
    }
    catch (Exception e) {
      e.printStackTrace();
      Assert.fail("Encountered exception");
    }
    finally {
      close(cdsIs);
      close(pepIs);
    }
  }

  @Test
  public void getUniversal() {
    IUPACParser.getInstance().getTable(1);
    IUPACParser.getInstance().getTable("UNIVERSAL");
  }

//  @Test
//  public void basicTranslation() {
//    RNAToAminoAcidTranslator t = TranslationTables.getInstance()
//        .getTranslator("UNIVERSAL");
//    DNASequence dna = new DNASequence("ATG");
//    RNASequence rna = dna.getRNASequence();
//
//    WindowedSequenceView<CodonCompound, NucleotideCompound> window = new WindowedSequenceView<CodonCompound, NucleotideCompound>(
//        rna, 3, t.getFromCompoundSet());
//
//    CodonCompound one = window.getCompoundAt(1);
//
//    AminoAcidCompound initMet = t.translate(one);
//
//    assertThat("Initator methionine wrong", initMet.toString(), is("M"));
//  }

  @Test
  public void translateBrca2ExonOne() {
    TranscriptionEngine e = TranscriptionEngine.getDefault();
    DNASequence dna = new DNASequence(
        "ATGCCTATTGGATCCAAAGAGAGGCCAACATTTTTTGAAATTTTTAAGACACGCTGCAACAAAGCA");
    RNASequence rna = dna.getRNASequence(e);
    Sequence<AminoAcidCompound> peptide = rna.getProteinSequence(e);
    assertThat("Initator methionine wrong", peptide.getSequenceAsString(),
        is("MPIGSKERPTFFEIFKTRCNKA"));
  }

  @Test(timeout=1000)
  public void translateBrca2() {
    TranscriptionEngine e = TranscriptionEngine.getDefault();
    RNAToAminoAcidTranslator t = e.getRnaAminoAcidTranslator();
    for(int i =0; i < 100; i++) {
    RNASequence rna = brca2Dna.getRNASequence();
    Sequence<AminoAcidCompound> peptide = t.createSequence(rna);
    assertThat("BRCA2 does not translate", peptide.getSequenceAsString(),
        is(brca2Pep.getSequenceAsString()));
    }
  }

}
