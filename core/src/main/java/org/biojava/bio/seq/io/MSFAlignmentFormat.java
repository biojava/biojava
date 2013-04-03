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

package org.biojava.bio.seq.io;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.Vector;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.biojava.bio.BioException;
import org.biojava.bio.alignment.Alignment;
import org.biojava.bio.alignment.SimpleAlignment;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.ProteinTools;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.symbol.SymbolList;

/**
 * @author raemig
 * @author Thomas Down
 * @author Keith James
 * @author Nimesh Singh
 * @author Mark Schreiber
 * @author Matthew Pocock
 * @author Bradford Powell
 */

public class MSFAlignmentFormat implements AlignmentFormat {
	private static final boolean DEBUGPRINT = false;
	private static final int DNA = 1;
	private static final int PROTEIN = 2;

	public MSFAlignmentFormat() {
	}

	/**
	 * used to quick test the code
	 * 
	 * @param args
	 */
	public static void main(String[] args) {
		String filename;
		if (args.length < 1) {
			filename = "SimpleMSF.msf"; // change to your favorite
		} else {
			filename = args[0];
		}
		try {
			BufferedReader reader = new BufferedReader(new FileReader(filename));
			MSFAlignmentFormat MSFAlignmentFormat1 = new MSFAlignmentFormat();
			MSFAlignmentFormat1.read(reader);
		} catch (Exception E) {
		}
	}

	/**
	 * Reads an MSF Alignment File
	 * 
	 * @param reader
	 *            The file reader
	 * @return Alignment A SimpleAlignment consisting of the sequences in the
	 *         file.
	 */
	public Alignment read(BufferedReader reader) {
		Vector sequenceNames = new Vector();
		String sequenceName = null;
		StringBuffer sequenceData[] = null;
		int startOfData = 0; // the start of the sequence data in the line
		int currSeqCount = 0; // which sequence data you are currently trying to
		// get
		try {
			Pattern mtc = Pattern
					.compile("(Name:|NAME:)\\s+(.*?)\\s+(oo|OO|Len:|LEN:)");
			Pattern removewhitespace = Pattern.compile("\\s");
			// REMatch rem = null;
			String line = reader.readLine();
			// parse past header
			while (line.toUpperCase().indexOf("NAME:") == -1) {
				line = reader.readLine();
			}
			// read each name (between Name: and Len:
			while ((line.indexOf("//") == -1) && ((line.trim()).length() != 0)) {
				Matcher matcher = mtc.matcher(line);
				if (!matcher.find()) {
					break;
				} // end of sequence names
				// sequenceName = line.substring(rem.getSubStartIndex(1),
				// rem.getSubEndIndex(1));
				if ((line.trim()).length() == 0) {
					break;
				}
				sequenceName = matcher.group(2).trim();
				sequenceNames.add(sequenceName);

				line = reader.readLine();
			}
			sequenceData = new StringBuffer[sequenceNames.size()];
			for (int it = 0; it < sequenceNames.size(); it++) {
				sequenceData[it] = new StringBuffer();
			}
			// until you get a line that matches the first sequence
			while (line.indexOf((String) sequenceNames.get(0)) == -1) {
				line = reader.readLine();
			}
			// now you on the first line of the sequence data
			while (line != null) {
				for (currSeqCount = 0; currSeqCount < sequenceNames.size(); currSeqCount++) {// you
					// could
					// also
					// check
					// for
					// order
					// of
					// names
					if (line.indexOf((String) sequenceNames.get(currSeqCount)) == -1) {
						break;
					} // error

					startOfData = line.indexOf((String) sequenceNames
							.get(currSeqCount))
							+ ((String) sequenceNames.get(currSeqCount))
									.length();
					line = (line.substring(startOfData));
					line = removewhitespace.matcher(line).replaceAll("");
					sequenceData[currSeqCount].append(line); // make into string
					// buffer
					line = reader.readLine();
					if ((currSeqCount < sequenceNames.size() - 1)
							&& (line.trim().length() == 0)) {
						break;
					} // could be an error
				}
				// until you get a line that matches the first sequence
				while ((line != null)
						&& (line.indexOf((String) sequenceNames.get(0)) == -1)) // ||
				// (
				// (line.trim())
				// .length()>0
				// )
				// )
				{
					line = reader.readLine();
				}
			}
			// print them out for testing
			if (DEBUGPRINT) {
				for (currSeqCount = 0; currSeqCount < sequenceNames.size(); currSeqCount++) {
					System.out.println((String) sequenceNames.get(currSeqCount)
							+ ":" + sequenceData[currSeqCount]);
				}
			}
			// check DNA, RNA or Prot
			StringBuffer testString = new StringBuffer();
			for (currSeqCount = 0; currSeqCount < sequenceNames.size(); currSeqCount++) {
				testString.append(sequenceData[currSeqCount]);
			}
			String testStringUpper = testString.toString().toUpperCase();

			// now parse through them and create gapped symbol lists
			LinkedHashMap sequenceDataMap = new LinkedHashMap();
			FiniteAlphabet alph = null;

			for (int i = 0; i < testStringUpper.length(); i++) {
				char c = testStringUpper.charAt(i);
				if (c == 'F' || c == 'L' || c == 'I' || c == 'P' || c == 'Q'
						|| c == 'E') {
					alph = ProteinTools.getTAlphabet();
					break;
				}
			}
			if (alph == null) {
				alph = DNATools.getDNA();
			}
			for (currSeqCount = 0; currSeqCount < sequenceNames.size(); currSeqCount++) {
				String sd = sequenceData[currSeqCount].toString();
				// change stop codons to specified symbols
				sd = sd.replace('~', '-'); // sometimes this is a term signal
				// not a gap
				sd = sd.replace('.', '-'); // sometimes this is a term signal
				// not a gap
				sequenceDataMap.put((String) sequenceNames.get(currSeqCount),
						alph == ProteinTools.getTAlphabet() ? ProteinTools
								.createGappedProteinSequence(sd,
										(String) sequenceNames
												.get(currSeqCount)) : DNATools
								.createGappedDNASequence(sd,
										(String) sequenceNames
												.get(currSeqCount)));
			}
			SimpleAlignment sa = new SimpleAlignment(sequenceDataMap);
			return (sa);
		} catch (Exception e) {
			e.printStackTrace();
			System.err.println("MSFFormatReader " + e.getMessage());
			// throw (e);
		}
		return (null);
	} // end read it

	// This is where I am writing an alignment writer
	public void write(OutputStream os, Alignment align, int fileType)
			throws BioException, IllegalSymbolException {
		PrintStream out = new PrintStream(os);
		Object labels[] = align.getLabels().toArray();
		int numSeqs = labels.length;
		Iterator<?> seqIts[] = new Iterator<?>[numSeqs];
		int maxLabelLength = 0;
		for (int i = 0; i < numSeqs; i++) {
			seqIts[i] = align.symbolListForLabel(labels[i].toString())
					.iterator();
			if (((String) labels[i]).length() > maxLabelLength) {
				maxLabelLength = ((String) labels[i]).length();
			}
		}
		String nl = System.getProperty("line.separator");
		SymbolTokenization toke = null;

		// really should determine the filetype based on one of the seqeunces
		// alphabet

		if (align.symbolListForLabel(labels[0].toString()).getAlphabet() == DNATools
				.getDNA()) {
			fileType = DNA;

		} else if (align.symbolListForLabel(labels[0].toString()).getAlphabet() == ProteinTools
				.getAlphabet()
				|| align.symbolListForLabel(labels[0].toString()).getAlphabet() == ProteinTools
						.getTAlphabet()) {
			fileType = PROTEIN;
		}

		if (fileType == DNA) {
			out.print("PileUp" + nl);
			out.print(nl);
			out.print(" MSF: " + align.length() + "  Type: ");
			out.print("N");
			out.print("   Check: " + 0 + "   .." + nl);
			toke = DNATools.getDNA().getTokenization("token");
		} else if (fileType == PROTEIN) {
			out.print("PileUp" + nl);
			out.print(nl);
			out.print(" MSF: " + align.length() + "  Type: ");
			out.print("P");
			out.print("   Check: " + 0 + "   .." + nl);
			toke = ProteinTools.getTAlphabet().getTokenization("token");
		} else {
			System.out
					.println("MSFAlignment.write -- File type not recognized.");
			return;
		}
		out.print(nl);

		for (int i = 0; i < numSeqs; i++) {
			out.print(" Name: " + labels[i]);
			for (int j = 0; j < (maxLabelLength - ((String) labels[i]).length()); j++) {// padding
				out.print(" ");
			}
			out.print("  Len: " + align.length() + " 	Check: " + 0
					+ "	Weight: " + 0 + nl); // this really should be seq
			// length?
		}

		out.println(nl + "//" + nl + nl);
		// now should print the numbering line

		while (seqIts[0].hasNext()) {
			for (int i = 0; i < numSeqs; i++) {
				while (((String) labels[i]).length() < maxLabelLength + 1) {
					labels[i] = " " + labels[i];
				}
				out.print(labels[i] + " ");
				theLabel: for (int j = 0; j < 5; j++) {
					out.print(" ");
					for (int k = 0; k < 10; k++) {
						if (seqIts[i].hasNext()) {
							out.print(toke.tokenizeSymbol((Symbol) seqIts[i]
									.next()));
						} else {
							break theLabel;
						}
					}
				}
				out.print(nl);
			}
			out.print(nl);
		}

	} // end write

	public void writeDna(OutputStream os, Alignment align) throws BioException,
			IllegalSymbolException {
		write(os, align, DNA);
	}

	public void writeProtein(OutputStream os, Alignment align)
			throws BioException, IllegalSymbolException {
		write(os, align, PROTEIN);
	}

} // end class

