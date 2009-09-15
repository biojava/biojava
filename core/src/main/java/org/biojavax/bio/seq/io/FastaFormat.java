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

package org.biojavax.bio.seq.io;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.io.ParseException;
import org.biojava.bio.seq.io.SeqIOListener;
import org.biojava.bio.seq.io.SymbolTokenization;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.SimpleSymbolList;
import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.symbol.SymbolList;
import org.biojava.utils.ChangeVetoException;
import org.biojavax.Namespace;
import org.biojavax.RichObjectFactory;
import org.biojavax.SimpleNamespace;
import org.biojavax.bio.seq.RichSequence;


/**
 * Format object representing FASTA files. These files are almost pure
 * sequence data.
 * @author Thomas Down
 * @author Matthew Pocock
 * @author Greg Cox
 * @author Lukas Kall
 * @author Richard Holland
 * @author Mark Schreiber
 * @since 1.5
 */

public class FastaFormat extends RichSequenceFormat.HeaderlessFormat {

	// Register this format with the format auto-guesser.
	static {
		RichSequence.IOTools.registerFormat(FastaFormat.class);
	}

	/**
	 * The name of this format
	 */
	public static final String FASTA_FORMAT = "FASTA";

	// header line
	protected static final Pattern hp = Pattern.compile(">\\s*(\\S+)(\\s+(.*))?");
	// description chunk
	protected static final Pattern dp = Pattern.compile( "^(gi\\|(\\d+)\\|)?(\\w+)\\|(\\w+?)(\\.(\\d+))?\\|(\\w+)?$");

	protected static final Pattern readableFiles = Pattern.compile(".*(fa|fas)$");
	protected static final Pattern aminoAcids = Pattern.compile(".*[FLIPQE].*");

	private FastaHeader header = new FastaHeader();

	/**
	 * {@inheritDoc}
	 * A file is in FASTA format if the name ends with fa or fas, or the file starts with ">".
	 */
	@Override
	public boolean canRead(File file) throws IOException {
		if (readableFiles.matcher(file.getName()).matches()) return true;
		BufferedReader br = new BufferedReader(new FileReader(file));
		String firstLine = br.readLine();
		boolean readable = firstLine!=null && firstLine.startsWith(">");
		br.close();
		return readable;
	}

	/**
	 * {@inheritDoc}
	 * Returns an protein parser if the first line of sequence contains any of F/L/I/P/Q/E, 
	 * otherwise returns a DNA tokenizer.
	 */
	@Override
	public SymbolTokenization guessSymbolTokenization(File file) throws IOException {
		BufferedReader br = new BufferedReader(new FileReader(file));
		br.readLine(); // discard first line
		boolean aa = aminoAcids.matcher(br.readLine()).matches();
		br.close();
		if (aa) return RichSequence.IOTools.getProteinParser();
		else return RichSequence.IOTools.getDNAParser();
	}

	/**
	 * {@inheritDoc}
	 * A stream is in FASTA format if the stream starts with ">".
	 */
	public boolean canRead(BufferedInputStream stream) throws IOException {
		stream.mark(2000); // some streams may not support this
		BufferedReader br = new BufferedReader(new InputStreamReader(stream));
		String firstLine = br.readLine();
		boolean readable = firstLine!=null && firstLine.startsWith(">");
		// don't close the reader as it'll close the stream too.
		// br.close();
		stream.reset();
		return readable;
	}

	/**
	 * {@inheritDoc}
	 * Returns an protein parser if the first line of sequence contains any of F/L/I/P/Q/E, 
	 * otherwise returns a DNA tokenizer.
	 */
	public SymbolTokenization guessSymbolTokenization(BufferedInputStream stream) throws IOException {
		stream.mark(2000); // some streams may not support this
		BufferedReader br = new BufferedReader(new InputStreamReader(stream));
		br.readLine(); // discard first line
		boolean aa = aminoAcids.matcher(br.readLine()).matches();
		// don't close the reader as it'll close the stream too.
		// br.close();
		stream.reset();
		if (aa) return RichSequence.IOTools.getProteinParser();
		else return RichSequence.IOTools.getDNAParser();
	}

	/**
	 * {@inheritDoc}
	 */
	public boolean readSequence(
			BufferedReader reader,
			SymbolTokenization symParser,
			SeqIOListener listener
	)	throws
	IllegalSymbolException,
	IOException,
	ParseException {
		if (!(listener instanceof RichSeqIOListener)) throw new IllegalArgumentException("Only accepting RichSeqIOListeners today");
		return this.readRichSequence(reader,symParser,(RichSeqIOListener)listener,null);
	}

	/**
	 * {@inheritDoc}
	 * If namespace is null, then the namespace of the sequence in the fasta is used.
	 * If the namespace is null and so is the namespace of the sequence in the fasta,
	 * then the default namespace is used.
	 */
	public boolean readRichSequence(
			BufferedReader reader,
			SymbolTokenization symParser,
			RichSeqIOListener rsiol,
			Namespace ns
	)	throws
	IllegalSymbolException,
	IOException,
	ParseException {

		String line = reader.readLine();
		if (line == null) {
			throw new IOException("Premature stream end");
		}
		while(line.length() == 0) {
			line = reader.readLine();
			if (line == null) {
				throw new IOException("Premature stream end");
			}
		}
		if (!line.startsWith(">")) {
			throw new IOException("Stream does not appear to contain FASTA formatted data: " + line);
		}

		rsiol.startSequence();

		processHeader(line,rsiol,ns);

		StringBuffer seq = new StringBuffer();
		boolean hasMoreSeq = true;
		while (hasMoreSeq) {
			reader.mark(500);
			line = reader.readLine();
			if (line!=null) {
				line = line.trim();
				if (line.length() > 0 && line.charAt(0)=='>') {
					reader.reset();
					hasMoreSeq = false;
				} else {
					seq.append(line);
				}
			} else {
				hasMoreSeq = false;
			}
		}
		if (!this.getElideSymbols()) {
			try {
				SymbolList sl = new SimpleSymbolList(symParser,
						seq.toString().replaceAll("\\s+","").replaceAll("[\\.|~]","-"));
				rsiol.addSymbols(symParser.getAlphabet(),
						(Symbol[])(sl.toList().toArray(new Symbol[0])),
						0, sl.length());
			} catch (Exception e) {
				// do not know name and gi any longer, replace them with empty string.
				// why does the rsiol only have setter methods, but not getter???
				String message = ParseException.newMessage(this.getClass(), "", "", "problem parsing symbols", seq.toString());
				throw new ParseException(e, message);
			}
		}

		rsiol.endSequence();

		return line!=null;
	}

	/** Parse the Header information from the Fasta Description line
	 * 
	 * @param line
	 * @param rsiol
	 * @param ns
	 * @throws IOException
	 * @throws ParseException
	 */
	public void processHeader(String line,RichSeqIOListener rsiol,Namespace ns) 
	throws IOException, ParseException {
		Matcher m = hp.matcher(line);
		if (!m.matches()) {
			throw new IOException("Stream does not appear to contain FASTA formatted data: " + line);
		}

		String name = m.group(1);
		String desc = m.group(3);
		String gi = null;

		m = dp.matcher(name);
		if (m.matches()) {
			gi = m.group(2);
			String namespace = m.group(3);
			String accession = m.group(4);
			String verString = m.group(6);
			int version = verString==null?0:Integer.parseInt(verString);
			name = m.group(7);
			if (name==null) name=accession;

			rsiol.setAccession(accession);
			rsiol.setVersion(version);
			if (gi!=null) rsiol.setIdentifier(gi);
			if (ns==null) rsiol.setNamespace((Namespace)RichObjectFactory.getObject(SimpleNamespace.class,new Object[]{namespace}));
			else rsiol.setNamespace(ns);
		} else {
			rsiol.setAccession(name);
			rsiol.setNamespace((ns==null?RichObjectFactory.getDefaultNamespace():ns));
		}
		rsiol.setName(name);
		if (!this.getElideComments()) rsiol.setDescription(desc);

	}

	/**
	 * {@inheritDoc}
	 */
	public void	writeSequence(Sequence seq, PrintStream os) throws IOException {
		if (this.getPrintStream()==null) this.setPrintStream(os);
		this.writeSequence(seq, RichObjectFactory.getDefaultNamespace());
	}

	/**
	 * {@inheritDoc}
	 */
	public void writeSequence(Sequence seq, String format, PrintStream os) throws IOException {
		if (this.getPrintStream()==null) this.setPrintStream(os);
		if (!format.equals(this.getDefaultFormat())) throw new IllegalArgumentException("Unknown format: "+format);
		this.writeSequence(seq, RichObjectFactory.getDefaultNamespace());
	}


	/**
	 * {@inheritDoc}
	 * If namespace is null, then the sequence's own namespace is used.
	 */
	public void writeSequence(Sequence seq, Namespace ns) throws IOException {
		RichSequence rs;
		try {
			if (seq instanceof RichSequence) rs = (RichSequence)seq;
			else rs = RichSequence.Tools.enrich(seq);
		} catch (ChangeVetoException e) {
			IOException e2 = new IOException("Unable to enrich sequence");
			e2.initCause(e);
			throw e2;
		}

		StringBuilder sb = new StringBuilder();
		sb.append(">");

		String identifier = rs.getIdentifier();
		if (header.isShowIdentifier() && identifier!=null && !"".equals(identifier)) {
			sb.append("gi|");
			sb.append(identifier);
			sb.append("|");
		}
		if(header.isShowNamespace()){
			sb.append((ns==null?rs.getNamespace().getName():ns.getName()));
			sb.append("|");
		}
		if(header.isShowAccession()){
			sb.append(rs.getAccession());
			if(header.isShowVersion()){
				sb.append(".");
			}
		}
		if(header.isShowVersion()){
			sb.append(rs.getVersion());
			sb.append("|");
		}
		if(header.isShowName()){
			sb.append(rs.getName());
			sb.append(" ");
		}else{
			sb.append(" "); //in case the show the description there needs to be space
		}
		if(header.isShowDescription()){
			String desc = rs.getDescription();
			if (desc!=null && !"".equals(desc)) sb.append(desc.replaceAll("\\n"," "));
		}
		if(sb.charAt(sb.length() -1) == '|'){
			sb.deleteCharAt(sb.length() -1);
		}
		this.getPrintStream().print(sb.toString());
		this.getPrintStream().println();

		int length = rs.length();

		for (int pos = 1; pos <= length; pos += this.getLineWidth()) {
			int end = Math.min(pos + this.getLineWidth() - 1, length);
			this.getPrintStream().println(rs.subStr(pos, end));
		}
	}

	/**
	 * {@inheritDoc}
	 */
	public String getDefaultFormat() {
		return FASTA_FORMAT;
	}

	public FastaHeader getHeader() {
		return header;
	}

	public void setHeader(FastaHeader header) {
		this.header = header;
	}
}
