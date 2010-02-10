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
package org.biojavax.bio.phylo.io.nexus;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.Reader;
import java.io.Writer;
import java.util.Iterator;
import java.util.Stack;
import java.util.StringTokenizer;

import org.biojava.bio.seq.io.ParseException;

/**
 * Reads/writes Nexus files and fires events at a NexusFileListener object.
 * Blocks are parsed using NexusBlockParser objects provided at runtime. Each of
 * those objects should probably have a NexusBlockListener object associated
 * with them that receives events generated from the processed data in the
 * block.
 * 
 * @author Richard Holland
 * @author Tobias Thierer
 * @author Jim Balhoff
 * @since 1.6
 */
public class NexusFileFormat {

	/**
	 * New-line symbol.
	 */
	public static final String NEW_LINE = System.getProperty("line.separator");

	// Prevent instances.
	private NexusFileFormat() {
	}

	/**
	 * Parse a file and send events to the given listener.
	 * 
	 * @param listener
	 *            the listener that will receive events.
	 * @param inputFile
	 *            the file to parse.
	 * @throws IOException
	 *             if anything goes wrong with reading the file.
	 * @throws ParseException
	 *             if the file format is incorrect.
	 */
	public static void parseFile(final NexusFileListener listener,
			final File inputFile) throws IOException, ParseException {
		final FileReader fr = new FileReader(inputFile);
		try {
			NexusFileFormat.parseReader(listener, fr);
		} finally {
			fr.close();
		}
	}

	/**
	 * Parse a stream and send events to the given listener.
	 * 
	 * @param listener
	 *            the listener that will receive events.
	 * @param inputStream
	 *            the stream to parse.
	 * @throws IOException
	 *             if anything goes wrong with reading the stream.
	 * @throws ParseException
	 *             if the stream format is incorrect.
	 */
	public static void parseInputStream(final NexusFileListener listener,
			final InputStream inputStream) throws IOException, ParseException {
		NexusFileFormat.parseReader(listener,
				new InputStreamReader(inputStream));
	}

	/**
	 * Parse a reader and send events to the given listener.
	 * 
	 * @param listener
	 *            the listener that will receive events.
	 * @param inputReader
	 *            the file to parse.
	 * @throws IOException
	 *             if anything goes wrong with reading the reader.
	 * @throws ParseException
	 *             if the reader format is incorrect.
	 */
	public static void parseReader(final NexusFileListener listener,
			final Reader inputReader) throws IOException, ParseException {
		NexusFileFormat
				.parse(
						listener,
						inputReader instanceof BufferedReader ? (BufferedReader) inputReader
								: new BufferedReader(inputReader));
	}

	// Do the work!
	private static void parse(final NexusFileListener listener,
			final BufferedReader reader) throws IOException, ParseException {
		// What are our delims?
		String space = " ";
		String tab = "\t";
		String beginComment = "[";
		String endComment = "]";
		String singleQuote = "'";
		String underscore = "_";
		String endTokenGroup = ";";
		String openBracket = "(";
		String closeBracket = ")";
		String openBrace = "{";
		String closeBrace = "}";
		String newLine = "\n";
		String allDelims = space + tab + beginComment + endComment
				+ singleQuote + underscore + endTokenGroup + newLine
				+ openBracket + closeBracket + openBrace + closeBrace;

		// Reset status flags.
		int inComment = 0;
		boolean inSingleQuotes = false;
		boolean inDoubleQuotes = false;
		boolean singleQuoteOpened = false;
		TokenParser parser = new TokenParser(listener);

		// Read the file line-by-line.
		final Stack parsedTokBufferStack = new Stack();
		StringBuffer parsedTokBuffer = new StringBuffer();
		String line;
		while ((line = reader.readLine()) != null) {
			final StringTokenizer tokenizer = new StringTokenizer(
					line.replaceAll("\\r\\n|\\r", "\n")
					+ "\n", allDelims, true);
			while (tokenizer.hasMoreTokens()) {
				final String tok = tokenizer.nextToken();

				// Process token.
				if (allDelims.indexOf(tok) >= 0) {
					// Process double quotes by flipping inside quote
					// status and appending the quote to the end of
					// the current token buffer then skipping to the
					// next parsed token.
					if (singleQuoteOpened && singleQuote.equals(tok)) {
						inSingleQuotes = !inSingleQuotes;
						parsedTokBuffer.append(singleQuote);
					}
					// Stuff inside comments.
					else if (inComment > 0) {
						// Start or end quotes?
						if (singleQuote.equals(tok))
							inSingleQuotes = !inSingleQuotes;
						// Nested comment.
						else if (beginComment.equals(tok) && !inSingleQuotes
								&& !inDoubleQuotes) {
							// Flush any existing comment text.
							if (parsedTokBuffer.length() > 0) {
								listener
										.commentText(parsedTokBuffer.toString());
								parsedTokBuffer.setLength(0);
							}
							// Start the new comment.
							inComment++;
							listener.beginComment();
							parsedTokBufferStack.push(parsedTokBuffer);
							parsedTokBuffer = new StringBuffer();
						}
						// Closing comment, not inside quotes. This
						// fires the current token buffer contents
						// as plain text at the listener, then clears
						// the buffer.
						else if (endComment.equals(tok) && !inSingleQuotes
								&& !inDoubleQuotes) {
							inComment--;
							if (parsedTokBuffer.length() > 0)
								listener
										.commentText(parsedTokBuffer.toString());
							listener.endComment();
							parsedTokBuffer = (StringBuffer) parsedTokBufferStack
									.pop();
						}
						// All other tokens are appended to the comment
						// buffer.
						else
							parsedTokBuffer.append(tok);
					}
					// Delimiter inside quotes.
					else if (inSingleQuotes) {
						// Closing quote puts us outside quotes.
						if (singleQuote.equals(tok))
							inSingleQuotes = false;
						// All other delimiters copied verbatim.
						else
							parsedTokBuffer.append(tok);
					}
					// Delimiter outside quote or comment.
					else {
						// Begin comment.
						if (beginComment.equals(tok)) {
							// Start the new comment.
							inComment++;
							listener.beginComment();
							// Preserve any existing part-built tag.
							parsedTokBufferStack.push(parsedTokBuffer);
							parsedTokBuffer = new StringBuffer();
						}
						// Start quoted string.
						else if (singleQuote.equals(tok))
							inSingleQuotes = true;
						// Convert underscores to spaces.
						else if (underscore.equals(tok))
							parsedTokBuffer.append(space);
						// Brackets. Pass through as tokens if
						// the client wishes, or just append to buffer
						// if client doesn't care.
						else if (openBracket.equals(tok)
								|| closeBracket.equals(tok)
								|| openBrace.equals(tok)
								|| closeBrace.equals(tok)) {
							if (listener.wantsBracketsAndBraces()) {
								// Dump buffer so far.
								final String parsedTok = parsedTokBuffer
										.toString();
								parsedTokBuffer.setLength(0);
								parser.parseToken(parsedTok);
								// Parse bracket/brace itself.
								listener.parseToken(tok);
							} else
								parsedTokBuffer.append(tok);
						}
						// Use whitespace/semi-colon to indicate end
						// of current token.
						else if (space.equals(tok) || tab.equals(tok)
								|| endTokenGroup.equals(tok)
								|| newLine.equals(tok)) {
							// Don't bother checking token buffer contents if
							// the buffer is empty.
							if (parsedTokBuffer.length() > 0) {
								final String parsedTok = parsedTokBuffer
										.toString();
								parsedTokBuffer.setLength(0);
								parser.parseToken(parsedTok);
							}

							// If this was an end-line, let the listeners know.
							if (endTokenGroup.equals(tok))
								listener.endTokenGroup();
							// Otherwise pass all whitespace through as
							// additional tokens.
							else
								listener.parseToken(tok);
						}
					}
				}
				// Process all non-delimiter tokens.
				else
					// Add token to buffer so far.
					parsedTokBuffer.append(tok);

				// Update double quote status. The next token is a potential
				// double
				// quote if the previous token was NOT a quote but this one IS.
				singleQuoteOpened = !singleQuoteOpened
						&& singleQuote.equals(tok);
			}
		}

		// End the listener.
		listener.endFile();
	}

	private static class TokenParser {

		private boolean expectingHeader = true;

		private boolean expectingBeginTag = false;

		private boolean expectingBeginName = false;

		private boolean expectingBlockContents = false;

		private NexusFileListener listener;

		private TokenParser(final NexusFileListener listener) {
			this.listener = listener;
		}

		private void parseToken(final String parsedTok) throws ParseException {

			// Expecting header?
			if (this.expectingHeader && "#NEXUS".equalsIgnoreCase(parsedTok)) {
				this.expectingHeader = false;
				this.expectingBeginTag = true;
				this.listener.startFile();
			}

			// Expecting a BEGIN tag?
			else if (this.expectingBeginTag
					&& "BEGIN".equalsIgnoreCase(parsedTok)) {
				this.expectingBeginTag = false;
				this.expectingBeginName = true;
			}

			// Expecting a name for a BEGIN block?
			else if (this.expectingBeginName) {
				this.listener.startBlock(parsedTok);
				this.expectingBeginName = false;
				this.expectingBlockContents = true;
			}

			// Looking for block contents?
			else if (this.expectingBlockContents) {
				// End tag?
				if ("END".equalsIgnoreCase(parsedTok)) {
					this.listener.endBlock();
					this.expectingBlockContents = false;
					this.expectingBeginTag = true;
				}
				// Or just normal token.
				else
					this.listener.parseToken(parsedTok);
			}

			// All other situations.
			else
				throw new ParseException(
						"Parser in unknown state when parsing token \""
								+ parsedTok + "\"");
		}
	}

	/**
	 * Writes the given Nexus output to a file.
	 * 
	 * @param file
	 *            the file to write to.
	 * @param nexusFile
	 *            the Nexus output to write.
	 * @throws IOException
	 *             if there is a problem during writing.
	 */
	public static void writeFile(final File file, final NexusFile nexusFile)
			throws IOException {
		final FileWriter fw = new FileWriter(file);
		try {
			NexusFileFormat.writeWriter(fw, nexusFile);
		} finally {
			fw.close();
		}
	}

	/**
	 * Writes the given Nexus output to a stream.
	 * 
	 * @param os
	 *            the stream to write to.
	 * @param nexusFile
	 *            the Nexus output to write.
	 * @throws IOException
	 *             if there is a problem during writing.
	 */
	public static void writeStream(final OutputStream os,
			final NexusFile nexusFile) throws IOException {
		final OutputStreamWriter ow = new OutputStreamWriter(os);
		NexusFileFormat.writeWriter(ow, nexusFile);
	}

	/**
	 * Writes the given Nexus output to a writer.
	 * 
	 * @param writer
	 *            the writer to write to.
	 * @param nexusFile
	 *            the Nexus output to write.
	 * @throws IOException
	 *             if there is a problem during writing.
	 */
	public static void writeWriter(final Writer writer,
			final NexusFile nexusFile) throws IOException {
		writer.write("#NEXUS");
		writer.write(NexusFileFormat.NEW_LINE);
		for (final Iterator i = nexusFile.objectIterator(); i.hasNext();) {
			((NexusObject) i.next()).writeObject(writer);
			writer.write(NexusFileFormat.NEW_LINE);
		}
		writer.flush();
	}
}
