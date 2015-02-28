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
 * Created on 01-21-2010
 */
package org.biojava.nbio.core.sequence.io.util;

import org.biojava.nbio.core.exceptions.ParserException;
import org.biojava.nbio.core.sequence.compound.AmbiguityDNACompoundSet;
import org.biojava.nbio.core.sequence.compound.AmbiguityRNACompoundSet;
import org.biojava.nbio.core.sequence.compound.DNACompoundSet;
import org.biojava.nbio.core.sequence.compound.RNACompoundSet;
import org.biojava.nbio.core.sequence.template.Compound;
import org.biojava.nbio.core.sequence.template.CompoundSet;
import org.biojava.nbio.core.sequence.template.Sequence;

import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.zip.GZIPInputStream;

public class IOUtils {

	private static final int BUFFER = 4096;

	/**
	 * Closes any Object which implements the interface {@link Closeable} and
	 * sending any error to the logger but not forcing any explicit catching of
	 * stream errors.
	 *
	 * @param c The stream to close
	 */
	public static void close(Closeable c) {
		try {
			if (c != null) {
				c.close();
			}
		} catch (IOException e) {
			Logger log = Logger.getLogger(IOUtils.class.getName());
			log.log(Level.WARNING, "Cannot close down the given Closeable object", e);
		}
	}

	/**
	 * Moves the bytes from input to output using a 4KB byte array.
	 *
	 * @param input Input stream of bytes
	 * @param output Output stream of bytes
	 * @throws IOException If anything occurs in the case of the reads and writes
	 */
	public static void copy(InputStream input, OutputStream output)
			throws IOException {
		byte[] buffer = new byte[BUFFER];
		int n = 0;
		while (-1 != (n = input.read(buffer))) {
			output.write(buffer, 0, n);
		}
	}

	/**
	 * Takes in a reader and a processor, reads every line from the given
	 * file and then invokes the processor. What you do with the lines is
	 * dependent on your processor.
	 *
	 * The code will automatically close the given BufferedReader.
	 *
	 * @param br The reader to process
	 * @param processor The processor to invoke on all lines
	 * @throws ParserException Can throw this if we cannot parse the given reader
	 */
	public static void processReader(BufferedReader br, ReaderProcessor processor) throws ParserException {
		String line;
		try {
			while( (line = br.readLine()) != null ) {
				processor.process(line);
			}
		}
		catch(IOException e) {
			throw new ParserException("Could not read from the given BufferedReader");
		}
		finally {
			close(br);
		}
	}

	/**
	 * Returns the contents of a buffered reader as a list of strings
	 *
	 * @param br BufferedReader to read from; <strong>will be closed</strong>
	 * @return List of Strings
	 * @throws ParserException Can throw this if we cannot parse the given reader
	 */
	public static List<String> getList(BufferedReader br) throws ParserException {
		final List<String> list = new ArrayList<String>();
		processReader(br, new ReaderProcessor() {
			@Override
			public void process(String line) {
				list.add(line);
			}
		});
		return list;
	}

	/**
	 * Delegates to {@link #getList(BufferedReader)} by wrapping the InputStream
	 * in a valid reader. No encoding is mentioned so if you need anything
	 * more advanced then use the other version of this method.
	 *
	 * @param is InputStream which is a text file
	 * @return List of Strings representing the lines of the files
	 * @throws ParserException Can throw this if the file is not a file or we
	 * cannot parse it
	 */
	public static List<String> getList(InputStream is) throws ParserException {
		return getList(new BufferedReader(new InputStreamReader(is)));
	}

	/**
	 * Delegates to {@link #getList(InputStream)} by wrapping the File
	 * in a valid stream. No encoding is mentioned so if you need anything
	 * more advanced then use the other version of this method. Since this
	 * uses {@link #openFile(File)} this code can support GZipped and plain
	 * files.
	 *
	 * @param file File which is a text file
	 * @return List of Strings representing the lines of the files
	 * @throws ParserException Can throw this if the file is not a file or we
	 * cannot parse it
	 */
	public static List<String> getList(File file) throws ParserException {
		return getList(openFile(file));
	}

	/**
	 * For a filename this code will check the extension of the file for a
	 * .gz extension. If it finds one then the InputStream given back
	 * is a {@link GZIPInputStream}. Otherwise we return a normal
	 * {@link FileInputStream}.
	 *
	 * @param file File which may or may not be GZipped
	 * @return The final stream
	 * @throws ParserException Can throw this if the file is not a file or we
	 * cannot open it for processing
	 */
	public static InputStream openFile(File file) throws ParserException {
		final InputStream is;
		if(!file.isFile()) {
			throw new ParserException("The file "+file+" is not a file.");
		}
		String name = file.getName();
		try {
			if(name.endsWith(".gz")) {
				is = new GZIPInputStream(new FileInputStream(file));
			}
			else {
				is = new FileInputStream(file);
			}
		}
		catch(IOException e) {
			throw new ParserException("Cannot open "+file+" for processing", e);
		}
		return is;
	}

	/**
	 * Closure interface used when working with
	 * {@link IOUtils#processReader(String)}. Each time a line is encountered
	 * the object that implements this interface will be invoked.
	 *
	 * @author ayates
	 */
	public static interface ReaderProcessor {
		void process(String line) throws IOException;
	}

	/**
	 * Calculates GCG checksum for entire list of sequences
	 *
	 * @param sequences list of sequences
	 * @return GCG checksum
	 */
	public static <S extends Sequence<C>, C extends Compound> int getGCGChecksum(List<S> sequences) {
		int check = 0;
		for (S as : sequences) {
			check += getGCGChecksum(as);
		}
		return check % 10000;
	}

	/**
	 * Calculates GCG checksum for a given sequence
	 *
	 * @param sequence given sequence
	 * @return GCG checksum
	 */
	public static <S extends Sequence<C>, C extends Compound> int getGCGChecksum(S sequence) {
		String s = sequence.toString().toUpperCase();
		int count = 0, check = 0;
		for (int i = 0; i < s.length(); i++) {
			count++;
			check += count * s.charAt(i);
			if (count == 57) {
				count = 0;
			}
		}
		return check % 10000;
	}

	/**
	 * Assembles a GCG file header
	 *
	 * @param sequences list of sequences
	 * @return GCG header
	 */
	public static <S extends Sequence<C>, C extends Compound> String getGCGHeader(List<S> sequences) {
		StringBuilder header = new StringBuilder();
		S s1 = sequences.get(0);
		header.append(String.format("MSA from BioJava%n%n MSF: %d  Type: %s  Check: %d ..%n%n",
				s1.getLength(), getGCGType(s1.getCompoundSet()), getGCGChecksum(sequences)));
		String format = " Name: " + getIDFormat(sequences) + " Len: " + s1.getLength() + "  Check: %4d  Weight: 1.0%n";
		for (S as : sequences) {
			header.append(String.format(format, as.getAccession(), getGCGChecksum(as)));
			// TODO show weights in MSF header
		}
		header.append(String.format("%n//%n%n"));
		// TODO? convert gap characters to '.'
		return header.toString();
	}

	/**
	 * Determines GCG type
	 * @param cs compound set of sequences
	 * @return GCG type
	 */
	public static <C extends Compound> String getGCGType(CompoundSet<C> cs) {
		return (cs == DNACompoundSet.getDNACompoundSet() || cs == AmbiguityDNACompoundSet.getDNACompoundSet()) ? "D" :
			(cs == RNACompoundSet.getRNACompoundSet() || cs == AmbiguityRNACompoundSet.getRNACompoundSet()) ? "R" : "P";
	}

	/**
	 * Creates format String for accession IDs
	 *
	 * @param sequences list of sequences
	 * @return format String for accession IDs
	 */
	public static <S extends Sequence<C>, C extends Compound> String getIDFormat(List<S> sequences) {
		int length = 0;
		for (S as : sequences) {
			length = Math.max(length, (as.getAccession() == null) ? 0 : as.getAccession().toString().length());
		}
		return (length == 0) ? null : "%-" + (length + 1) + "s";
	}

	/**
	 * Creates formatted String for a single character of PDB output
	 *
	 * @param web true for HTML display
	 * @param c1 character in first sequence
	 * @param c2 character in second sequence
	 * @param similar true if c1 and c2 are considered similar compounds
	 * @param c character to display
	 * @return formatted String
	 */
	public static String getPDBCharacter(boolean web, char c1, char c2, boolean similar, char c) {
		String s = c + "";
		return getPDBString(web, c1, c2, similar, s, s, s, s);
	}

	/**
	 * Creates formatted String for displaying conservation in PDB output
	 *
	 * @param web true for HTML display
	 * @param c1 character in first sequence
	 * @param c2 character in second sequence
	 * @param similar true if c1 and c2 are considered similar compounds
	 * @return formatted String
	 */
	public static String getPDBConservation(boolean web, char c1, char c2, boolean similar) {
		return getPDBString(web, c1, c2, similar, "|", ".", " ", web ? "&nbsp;" : " ");
	}

	// helper method for getPDBCharacter and getPDBConservation
	private static String getPDBString(boolean web, char c1, char c2, boolean similar, String m, String sm, String dm,
			String qg) {
		if (c1 == c2)
			return web ? "<span class=\"m\">" + m + "</span>" : m;                             
		else if (similar)
			return web ? "<span class=\"sm\">" + sm + "</span>" : sm;
		else if (c1 == '-' || c2 == '-')
			return web ? "<span class=\"dm\">" + dm + "</span>" : dm;
		else
			return web ? "<span class=\"qg\">" + qg + "</span>" : qg;
	}

	/**
	 * Creates formatted String for displaying conservation legend in PDB output
	 *
	 * @return legend String
	 */
	public static String getPDBLegend() {
		StringBuilder s = new StringBuilder();
		s.append("</pre></div>");
		s.append("          <div class=\"subText\">");
		s.append("          <b>Legend:</b>");
		s.append("          <span class=\"m\">Green</span> - identical residues |"); 
		s.append("          <span class=\"sm\">Pink</span> - similar residues | ");
		s.append("          <span class=\"qg\">Blue</span> - sequence mismatch |");
		s.append("          <span class=\"dm\">Brown</span> - insertion/deletion |");                  
		s.append("      </div>");
		s.append(String.format("%n"));
		return s.toString();
	}

	/**
	 * Prints {@code string} to {@code file}.
	 * @throws IOException If any I/O exception occurs while printing; this method does not catch any exceptions
	 */
	public static void print(String string, File file) throws IOException {
		PrintWriter out = null;
		try {
			out = new PrintWriter(new BufferedWriter(new FileWriter(file)));
			out.print(string);
			out.flush();
			out.close();
		} finally {
			if (out != null) out.close();
		}
	}

}
