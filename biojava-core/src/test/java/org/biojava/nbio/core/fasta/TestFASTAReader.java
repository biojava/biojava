package org.biojava.nbio.core.fasta;

import java.io.InputStream;
import java.util.LinkedHashMap;

import static org.junit.Assert.* ;
import static org.hamcrest.CoreMatchers.* ;

import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava.nbio.core.sequence.io.FastaReader;
import org.biojava.nbio.core.sequence.io.GenericFastaHeaderParser;
import org.biojava.nbio.core.sequence.io.ProteinSequenceCreator;
import org.biojava.nbio.core.sequence.io.util.ClasspathResource;
import org.junit.Test;


public class TestFASTAReader {

    @Test
    public void testProcessAll() throws Exception {
        ClasspathResource r = new ClasspathResource("org/biojava/nbio/core/fasta/P02768.fasta");
        FastaReader<ProteinSequence, AminoAcidCompound> fastaReader = null ;
        try( InputStream inStream = r.getInputStream() ) {
            fastaReader = new FastaReader<ProteinSequence, AminoAcidCompound>(
                    inStream,
                    new GenericFastaHeaderParser<ProteinSequence, AminoAcidCompound>(),
                    new ProteinSequenceCreator(AminoAcidCompoundSet.getAminoAcidCompoundSet()));
            LinkedHashMap<String, ProteinSequence> sequences = fastaReader.process();
            assertThat(sequences,is(notNullValue()));
            assertThat(sequences.size(),is(1));
        } finally {
            if(fastaReader != null) fastaReader.close();
        }
    }
    @Test
    public void testProcess1() throws Exception {
        ClasspathResource r = new ClasspathResource("org/biojava/nbio/core/fasta/P02768.fasta");
        FastaReader<ProteinSequence, AminoAcidCompound> fastaReader = null ;
        try( InputStream inStream = r.getInputStream() ) {
            fastaReader = new FastaReader<ProteinSequence, AminoAcidCompound>(
                    inStream,
                    new GenericFastaHeaderParser<ProteinSequence, AminoAcidCompound>(),
                    new ProteinSequenceCreator(AminoAcidCompoundSet.getAminoAcidCompoundSet()));
            assertThat(fastaReader.process(1),is(notNullValue()));
            assertThat(fastaReader.process(1),is(nullValue()));
        } finally {
            if(fastaReader != null) fastaReader.close();
        }
    }
    @Test
    public void testProcess1v2() throws Exception {
        ClasspathResource r = new ClasspathResource("org/biojava/nbio/core/fasta/TwoSequences.fasta");
        FastaReader<ProteinSequence, AminoAcidCompound> fastaReader = null ;
        try( InputStream inStream = r.getInputStream() ) {
            fastaReader = new FastaReader<ProteinSequence, AminoAcidCompound>(
                    inStream,
                    new GenericFastaHeaderParser<ProteinSequence, AminoAcidCompound>(),
                    new ProteinSequenceCreator(AminoAcidCompoundSet.getAminoAcidCompoundSet()));
            assertThat(fastaReader.process(1),is(notNullValue()));
            assertThat(fastaReader.process(1),is(notNullValue()));
            assertThat(fastaReader.process(1),is(nullValue()));
        } finally {
            if(fastaReader != null) fastaReader.close();
        }
    }

}
