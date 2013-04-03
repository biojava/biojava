/* 
 * @(#)SequenceUtil.java 1.0 September 2009
 * 
 * Copyright (c) 2009 Peter Troshin
 *  
 *        BioJava development code
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

package org.biojava3.data.sequence;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Utility class for operations on sequences
 * 
 * @author Peter Troshin
 * @version 1.0
 * @since 3.0.2
 */
public final class SequenceUtil {

    /**
     * A whitespace character: [\t\n\x0B\f\r]
     */
    public static final Pattern WHITE_SPACE = Pattern.compile("\\s");

    /**
     * A digit
     */
    public static final Pattern DIGIT = Pattern.compile("\\d");

    /**
     * Non word
     */
    public static final Pattern NONWORD = Pattern.compile("\\W");

    /**
     * Valid Amino acids
     */
    public static final Pattern AA = Pattern.compile("[ARNDCQEGHILKMFPSTWYV]+",
	    Pattern.CASE_INSENSITIVE);

    /**
     * inversion of AA pattern
     */
    public static final Pattern NON_AA = Pattern.compile(
	    "[^ARNDCQEGHILKMFPSTWYVX]+", Pattern.CASE_INSENSITIVE);

    /**
     * Same as AA pattern but with two additional letters - XU
     */
    public static final Pattern AMBIGUOUS_AA = Pattern.compile(
	    "[ARNDCQEGHILKMFPSTWYVXU]+", Pattern.CASE_INSENSITIVE);

    /**
     * Nucleotides a, t, g, c, u
     */
    public static final Pattern NUCLEOTIDE = Pattern.compile("[AGTCU]+",
	    Pattern.CASE_INSENSITIVE);

    /**
     * Ambiguous nucleotide
     */
    public static final Pattern AMBIGUOUS_NUCLEOTIDE = Pattern.compile(
	    "[AGTCRYMKSWHBVDNU]+", Pattern.CASE_INSENSITIVE); // see IUPAC
    /**
     * Non nucleotide
     */
    public static final Pattern NON_NUCLEOTIDE = Pattern.compile("[^AGTCU]+",
	    Pattern.CASE_INSENSITIVE);

    private SequenceUtil() {
    } // utility class, no instantiation

    /*
     * public static void write_PirSeq(OutputStream os, FastaSequence seq)
     * throws IOException { BufferedWriter pir_out = new BufferedWriter(new
     * OutputStreamWriter(os)); pir_out.write(">P1;" + seq.getId() +
     * SysPrefs.newlinechar); pir_out.write(seq.getSequence() +
     * SysPrefs.newlinechar); pir_out.close(); }
     * 
     * public static void write_FastaSeq(OutputStream os, FastaSequence seq)
     * throws IOException { BufferedWriter fasta_out = new BufferedWriter( new
     * OutputStreamWriter(os)); fasta_out.write(">" + seq.getId() +
     * SysPrefs.newlinechar); fasta_out.write(seq.getSequence() +
     * SysPrefs.newlinechar); fasta_out.close(); }
     */

    /**
     * @return true is the sequence contains only letters a,c, t, g, u
     */
    public static boolean isNucleotideSequence(final FastaSequence s) {
	return SequenceUtil.isNonAmbNucleotideSequence(s.getSequence());
    }

    /**
     * Ambiguous DNA chars : AGTCRYMKSWHBVDN // differs from protein in only one
     * (!) - B char
     */
    public static boolean isNonAmbNucleotideSequence(String sequence) {
	sequence = SequenceUtil.cleanSequence(sequence);
	if (SequenceUtil.DIGIT.matcher(sequence).find()) {
	    return false;
	}
	if (SequenceUtil.NON_NUCLEOTIDE.matcher(sequence).find()) {
	    return false;
	    /*
	     * System.out.format("I found the text starting at " +
	     * "index %d and ending at index %d.%n", nonDNAmatcher .start(),
	     * nonDNAmatcher.end());
	     */
	}
	final Matcher DNAmatcher = SequenceUtil.NUCLEOTIDE.matcher(sequence);
	return DNAmatcher.find();
    }

    /**
     * Removes all whitespace chars in the sequence string
     * 
     * @param sequence
     * @return cleaned up sequence
     */
    public static String cleanSequence(String sequence) {
	assert sequence != null;
	final Matcher m = SequenceUtil.WHITE_SPACE.matcher(sequence);
	sequence = m.replaceAll("").toUpperCase();
	return sequence;
    }

    /**
     * Removes all special characters and digits as well as whitespace chars
     * from the sequence
     * 
     * @param sequence
     * @return cleaned up sequence
     */
    public static String deepCleanSequence(String sequence) {
	sequence = SequenceUtil.cleanSequence(sequence);
	sequence = SequenceUtil.DIGIT.matcher(sequence).replaceAll("");
	sequence = SequenceUtil.NONWORD.matcher(sequence).replaceAll("");
	final Pattern othernonSeqChars = Pattern.compile("[_-]+");
	sequence = othernonSeqChars.matcher(sequence).replaceAll("");
	return sequence;
    }

    /**
     * 
     * @param sequence
     * @return true is the sequence is a protein sequence, false overwise
     */
    public static boolean isProteinSequence(String sequence) {
	sequence = SequenceUtil.cleanSequence(sequence);
	if (SequenceUtil.isNonAmbNucleotideSequence(sequence)) {
	    return false;
	}
	if (SequenceUtil.DIGIT.matcher(sequence).find()) {
	    return false;
	}
	if (SequenceUtil.NON_AA.matcher(sequence).find()) {
		System.out.println("found non aa!");
	    return false;
	}
	final Matcher protmatcher = SequenceUtil.AA.matcher(sequence);
	return protmatcher.find();
    }

    /**
     * Check whether the sequence confirms to amboguous protein sequence
     * 
     * @param sequence
     * @return return true only if the sequence if ambiguous protein sequence
     *         Return false otherwise. e.g. if the sequence is non-ambiguous
     *         protein or DNA
     */
    public static boolean isAmbiguosProtein(String sequence) {
	sequence = SequenceUtil.cleanSequence(sequence);
	if (SequenceUtil.isNonAmbNucleotideSequence(sequence)) {
	    return false;
	}
	if (SequenceUtil.DIGIT.matcher(sequence).find()) {
	    return false;
	}
	if (SequenceUtil.NON_AA.matcher(sequence).find()) {
	    return false;
	}
	if (SequenceUtil.AA.matcher(sequence).find()) {
	    return false;
	}
	final Matcher amb_prot = SequenceUtil.AMBIGUOUS_AA.matcher(sequence);
	return amb_prot.find();
    }

    /**
     * Writes list of FastaSequeces into the outstream formatting the sequence
     * so that it contains width chars on each line
     * 
     * @param outstream
     * @param sequences
     * @param width
     *            - the maximum number of characters to write in one line
     * @throws IOException
     */
    public static void writeFasta(final OutputStream outstream,
	    final List<FastaSequence> sequences, final int width)
	    throws IOException {
	final OutputStreamWriter writer = new OutputStreamWriter(outstream);
	final BufferedWriter fastawriter = new BufferedWriter(writer);
	for (final FastaSequence fs : sequences) {
	    fastawriter.write(fs.getFormatedSequence(width));
	}
	outstream.flush();
	fastawriter.close();
	writer.close();
    }

    /**
     * Reads fasta sequences from inStream into the list of FastaSequence
     * objects
     * 
     * @param inStream
     *            from
     * @return list of FastaSequence objects
     * @throws IOException
     */
    public static List<FastaSequence> readFasta(final InputStream inStream)
	    throws IOException {
	final List<FastaSequence> seqs = new ArrayList<FastaSequence>();

	final BufferedReader infasta = new BufferedReader(
		new InputStreamReader(inStream, "UTF8"), 16000);
	final Pattern pattern = Pattern.compile("//s+");

	String line;
	String sname = "", seqstr = null;
	do {
	    line = infasta.readLine();
	    if ((line == null) || line.startsWith(">")) {
		if (seqstr != null) {
		    seqs.add(new FastaSequence(sname.substring(1), seqstr));
		}
		sname = line; // remove >
		seqstr = "";
	    } else {
		final String subseq = pattern.matcher(line).replaceAll("");
		seqstr += subseq;
	    }
	} while (line != null);

	infasta.close();
	return seqs;
    }

    /**
     * Writes FastaSequence in the file, each sequence will take one line only
     * 
     * @param os
     * @param sequences
     * @throws IOException
     */
    public static void writeFasta(final OutputStream os,
	    final List<FastaSequence> sequences) throws IOException {
	final OutputStreamWriter outWriter = new OutputStreamWriter(os);
	final BufferedWriter fasta_out = new BufferedWriter(outWriter);
	for (final FastaSequence fs : sequences) {
	    fasta_out.write(fs.getOnelineFasta());
	}
	fasta_out.close();
	outWriter.close();
    }

}
