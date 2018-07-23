/*
 *                    BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 *
 *      http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *      http://www.biojava.org/
 *
 */
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

	private void testProcessAll(String path) throws Exception {
        ClasspathResource r = new ClasspathResource(path);
        FastaReader<ProteinSequence, AminoAcidCompound> fastaReader = null ;
        try( InputStream inStream = r.getInputStream() ) {
            fastaReader = new FastaReader<ProteinSequence, AminoAcidCompound>(
                    inStream,
                    new GenericFastaHeaderParser<ProteinSequence, AminoAcidCompound>(),
                    new ProteinSequenceCreator(AminoAcidCompoundSet.getAminoAcidCompoundSet()));
            LinkedHashMap<String, ProteinSequence> sequences = fastaReader.process();
            assertThat(sequences,is(notNullValue()));
            assertThat(sequences.size(),is(1));
            assertThat(sequences.containsKey("P02768"),is(true));
            assertThat(sequences.get("P02768").getLength(),is(609));
        } finally {
            if(fastaReader != null) fastaReader.close();
        }
	}
	
	/**
	 * Test file contains one sequence (P02768 from swissprot). Read the whole
	 * file all at once by calling {@link FastaReader#process()} and verify that
	 * one sequence is read.
	 *
	 * @throws Exception
	 */
    @Test
    public void testProcessAll() throws Exception {
    	testProcessAll("org/biojava/nbio/core/fasta/P02768.fasta");
    }
    
    /**
     * Same as {@link #testProcessAll()} but input files contains blank lines
     * 
     * @throws Exception
     */
    @Test
    public void testProcessAllWithBlankLines() throws Exception {
    	testProcessAll("org/biojava/nbio/core/fasta/P02768_blank_lines.fasta");
    }
    
    private void testProcess1(String path) throws Exception {
        ClasspathResource r = new ClasspathResource(path);
        FastaReader<ProteinSequence, AminoAcidCompound> fastaReader = null ;
        try( InputStream inStream = r.getInputStream() ) {
            fastaReader = new FastaReader<ProteinSequence, AminoAcidCompound>(
                    inStream,
                    new GenericFastaHeaderParser<ProteinSequence, AminoAcidCompound>(),
                    new ProteinSequenceCreator(AminoAcidCompoundSet.getAminoAcidCompoundSet()));
            LinkedHashMap<String,ProteinSequence> out1 = fastaReader.process(1);
            assertThat(out1,is(notNullValue()));
            assertThat(out1.size(),is(1));
            assertThat(out1.containsKey("P02768"),is(true));
            assertThat(out1.get("P02768").getLength(),is(609));
            LinkedHashMap<String,ProteinSequence> out2 = fastaReader.process(1);
            assertThat(out2,is(nullValue()));
        } finally {
            if(fastaReader != null) fastaReader.close();
        }
    }
    
	/**
	 * Test file contains one sequence (P02768 from swissprot). Read one
	 * sequence at a time by calling {@link FastaReader#process(int)} and verify
	 * that the first call get one sequence and the second call get none.
	 * 
	 * @throws Exception
	 */
    @Test
    public void testProcess1() throws Exception {
    	testProcess1("org/biojava/nbio/core/fasta/P02768.fasta");
    }
    
    /**
     * Same as {@link #testProcess1()}, but input contains blank lines.
     * 
     * @throws Exception
     */
    @Test
    public void testProcess1WithBlankLines() throws Exception {
    	testProcess1("org/biojava/nbio/core/fasta/P02768_blank_lines.fasta");
    }
    
    private void testProcess2(String path) throws Exception {
        ClasspathResource r = new ClasspathResource(path);
        FastaReader<ProteinSequence, AminoAcidCompound> fastaReader = null ;
        try( InputStream inStream = r.getInputStream() ) {
            fastaReader = new FastaReader<ProteinSequence, AminoAcidCompound>(
                    inStream,
                    new GenericFastaHeaderParser<ProteinSequence, AminoAcidCompound>(),
                    new ProteinSequenceCreator(AminoAcidCompoundSet.getAminoAcidCompoundSet()));
            LinkedHashMap<String,ProteinSequence> out1 = fastaReader.process(1);
            assertThat(out1,is(notNullValue()));
            assertThat(out1.size(),is(1));
            assertThat(out1.containsKey("P02768"),is(true));
            assertThat(out1.get("P02768").getLength(),is(609));
            LinkedHashMap<String,ProteinSequence> out2 = fastaReader.process(1);
            assertThat(out2,is(notNullValue()));
            assertThat(out2.size(),is(1));
            assertThat(out2.containsKey("P00698"),is(true));
            assertThat(out2.get("P00698").getLength(),is(147));
            LinkedHashMap<String,ProteinSequence> out3 = fastaReader.process(1);
            assertThat(out3,is(nullValue()));
        } finally {
            if(fastaReader != null) fastaReader.close();
        }
    }
    
	/**
	 * Test file contains two sequences. Read one sequence at a time by calling
	 * {@link FastaReader#process(int)} and verify that the first and second
	 * call get one sequence each and the third call get none.
	 * 
	 * @throws Exception
	 */
    @Test
    public void testProcess2() throws Exception {
    	testProcess2("org/biojava/nbio/core/fasta/TwoSequences.fasta");
    }

    /**
     * Sane as {@link #testProcess2()} but input file contain blank lines
     * @throws Exception
     */
    @Test
    public void testProcess2WithBlankLines() throws Exception {
    	testProcess2("org/biojava/nbio/core/fasta/TwoSequences_blank_lines.fasta");
    }
}
