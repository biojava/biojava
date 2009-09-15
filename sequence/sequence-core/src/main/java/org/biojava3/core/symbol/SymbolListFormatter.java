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
package org.biojava3.core.symbol;

import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.Reader;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;

/**
 * The symbol list formatter can produce symbol lists from input streams and
 * strings, and can write them out to strings and output streams. It behaves in
 * one of three ways - either representing symbol lists as 1 character per
 * symbol, or by placing separators between the symbols, or by wrapping each
 * symbol in a start and end character. It provides for escaping the
 * separator/start/end characters in case they need to be appear inside the
 * individual symbol representations. The formatter requires an {@link Alphabet}
 * so that it knows how to decode or encode symbols from the stream using the
 * {@link SymbolEncoder} and {@link SymbolDecoder} instances from that alphabet.
 * 
 * @author Richard Holland
 * @since 3.0
 */
public class SymbolListFormatter implements Serializable {
	
	private static final long serialVersionUID = 1L;

	private Alphabet alpha;
	private String sep = null;
	private String start = null;
	private String end = null;
	private String escape = "\\";
	private char sepChar;
	private char startChar;
	private char endChar;
	private char escapeChar = '\\';

	/**
	 * Construct a formatter which will iterate over individual characters when
	 * inputting strings, and will output all symbols concatenated directly
	 * together when outputting strings.
	 * 
	 * @param alpha
	 *            the alphabet to provide the encoder and decoder for converting
	 *            between strings and symbols.
	 */
	public SymbolListFormatter(Alphabet alpha) {
		if (alpha == null) {
			throw new NullPointerException("Alphabet cannot be null.");
		}
		this.alpha = alpha;
	}

	/**
	 * Construct a formatter which will iterate over blocks of characters
	 * separated by the given separator character when inputting strings, and
	 * will output all symbols concatenated directly together when outputting
	 * strings.
	 * 
	 * @param alpha
	 *            the alphabet to provide the encoder and decoder for converting
	 *            between strings and symbols.
	 * @param sepChar
	 *            the character to use as a separator between symbols.
	 */
	public SymbolListFormatter(Alphabet alpha, char sepChar) {
		this(alpha);
		this.sep = "" + sepChar;
		this.sepChar = sepChar;
	}

	/**
	 * Construct a formatter which will iterate over blocks of characters
	 * wrapped by the given start/end characters when inputting strings, and
	 * will output all symbols wrapped with the start/end characters when
	 * outputting strings.
	 * 
	 * @param alpha
	 *            the alphabet to provide the encoder and decoder for converting
	 *            between strings and symbols.
	 * @param startChar
	 *            the character that denotes the start of a symbol.
	 * @param endChar
	 *            the character that denotes the end of a symbol.
	 */
	public SymbolListFormatter(Alphabet alpha, char startChar, char endChar) {
		this(alpha);
		this.start = "" + startChar;
		this.end = "" + endChar;
		this.startChar = startChar;
		this.endChar = endChar;
	}

	/**
	 * The escape character is prefixed to any occurrences of the start, end or
	 * separator characters within the strings returned by the alphabet's
	 * {@link SymbolEncoder}.
	 * 
	 * @param escapeChar
	 *            the escape character to use.
	 */
	public void setEscapeChar(char escapeChar) {
		this.escape = "" + escapeChar;
		this.escapeChar = escapeChar;
	}

	/**
	 * Converts the symbol into a string with all necessary start/end characters
	 * appended and escaping performed.
	 * 
	 * @param sym
	 *            the symbol.
	 * @return the string.
	 * @throws IOException
	 *             should never happen!
	 */
	public String formatSymbol(Symbol sym) throws IOException {
		StringBuilder sb = new StringBuilder();
		this.formatSymbol(sym, sb);
		return sb.toString();
	}

	/**
	 * Converts the symbol into a string with all necessary start/end characters
	 * appended and escaping performed.
	 * 
	 * @param sym
	 *            the symbol.
	 * @param a
	 *            the place to write the output to.
	 * @throws if the output failed.
	 */
	public void formatSymbol(Symbol sym, Appendable a) throws IOException {
		if (this.start != null) {
			a.append(this.start);
		}
		a.append(this.getEscapedSymbol(sym));
		if (this.end != null) {
			a.append(this.end);
		}
	}

	/**
	 * Converts the symbol into a string with all necessary start/end characters
	 * appended and escaping performed.
	 * 
	 * @param sym
	 *            the symbol.
	 * @param os
	 *            the place to write the output to.
	 * @throws if the output failed.
	 */
	public void formatSymbol(Symbol sym, OutputStream os) throws IOException {
		if (this.start != null) {
			os.write(this.start.getBytes());
		}
		os.write(this.getEscapedSymbol(sym).getBytes());
		if (this.end != null) {
			os.write(this.end.getBytes());
		}
	}

	/**
	 * Converts a symbol into an escaped string but without the start/end chars.
	 * 
	 * @param sym
	 *            the symbol.
	 * @return the escaped string.
	 */
	private String getEscapedSymbol(Symbol sym) {
		String str = this.alpha.getCoder().encodeSymbol(sym);
		str = str.replaceAll(this.escape, this.escape + this.escape);
		if (this.start != null) {
			str = str.replaceAll(this.start, this.escape + this.start);
		}
		if (this.end != null) {
			str = str.replaceAll(this.end, this.escape + this.end);
		}
		if (this.sep != null) {
			str = str.replaceAll(this.sep, this.escape + this.sep);
		}
		return str;
	}

	/**
	 * Converts a list of symbols into a string.
	 * 
	 * @param syms
	 *            the symbols.
	 * @return the string.
	 * @throws IOException
	 *             should never happen!
	 */
	public String formatSymbols(List<Symbol> syms) throws IOException {
		StringBuilder sb = new StringBuilder();
		this.formatSymbols(syms, sb);
		return sb.toString();
	}

	/**
	 * Converts a list of symbols into an output.
	 * 
	 * @param syms
	 *            the symbols.
	 * @param a
	 *            the output.
	 * @throws IOException
	 *             if the output fails.
	 */
	public void formatSymbols(List<Symbol> syms, Appendable a)
			throws IOException {
		boolean first = true;
		for (Symbol sym : syms) {
			if (!first && this.sep != null) {
				a.append(this.sep);
			}
			this.formatSymbol(sym, a);
			first = false;
		}
	}

	/**
	 * Converts a list of symbols into an output.
	 * 
	 * @param syms
	 *            the symbols.
	 * @param os
	 *            the output.
	 * @throws IOException
	 *             if the output fails.
	 */
	public void formatSymbols(List<Symbol> syms, OutputStream os)
			throws IOException {
		boolean first = true;
		for (Symbol sym : syms) {
			if (!first && this.sep != null) {
				os.write(this.sep.getBytes());
			}
			this.formatSymbol(sym, os);
			first = false;
		}
	}

	/**
	 * Attempts to read a symbol from the string.
	 * 
	 * @param str
	 *            a string.
	 * @return a symbol, or {@code null} if none could be read.
	 */
	public Symbol parseSymbol(CharSequence str) {
		Iterator<Symbol> iter = this.parseSymbols(str);
		if (iter.hasNext()) {
			return iter.next();
		}
		return null;
	}

	/**
	 * Attempts to read all symbols from the string.
	 * 
	 * @param str
	 *            a string.
	 * @return an iterator over all symbols.
	 */
	public Iterator<Symbol> parseSymbols(final CharSequence str) {
		return new SymbolIterator() {

			private int pos = 0;

			public boolean hasNextChar() {
				return this.pos < str.length();
			}

			public char nextChar() {
				return str.charAt(this.pos++);
			}
		};
	}

	/**
	 * Attempts to read all symbols from the input.
	 * 
	 * @param is
	 *            the input.
	 * @return an iterator over all symbols.
	 */
	public Iterator<Symbol> parseSymbols(final InputStream is) {
		return this.parseSymbols(new InputStreamReader(is));
	}

	/**
	 * Attempts to read all symbols from the input.
	 * 
	 * @param is
	 *            the input.
	 * @return an iterator over all symbols.
	 */
	public Iterator<Symbol> parseSymbols(final Reader r) {
		return new SymbolIterator() {

			private static final int BUFSIZE = 1024;
			private char[] buf = new char[BUFSIZE];
			private int idx = 0;
			private int size = 0;

			private void updateQueue() throws IOException {
				if (this.size == -1) {
					return;
				}
				this.size = r.read(buf, 0, BUFSIZE);
				this.idx = 0;
			}

			public boolean hasNextChar() {
				if (this.idx < this.size) {
					try {
						this.updateQueue();
					} catch (IOException ioe) {
						// Can't do anything inside the iterator
						// except pretend that it's just the end
						// of the stream.
					}
				}
				return this.idx < this.size;
			}

			public char nextChar() {
				return this.buf[this.idx++];
			}
		};
	}

	/**
	 * This internal class performs a read-ahead-one iteration over a source of
	 * characters. When created it will immediately parse the first symbol, then
	 * waits to be told to move on to the next one. This allows the
	 * {@link SymbolIterator#hasNext()} method to function in the way the user
	 * would expect.
	 * 
	 * @author Richard Holland
	 * @since 3.0
	 */
	private abstract class SymbolIterator implements Iterator<Symbol> {

		/**
		 * Called to see if there are any more characters waiting to be read.
		 * 
		 * @return {@code true} if there are.
		 */
		protected abstract boolean hasNextChar();

		/**
		 * Called to get the next character.
		 * 
		 * @return the next character.
		 */
		protected abstract char nextChar();

		private Symbol nextSym = null;

		public boolean hasNext() {
			if (this.nextSym == null) {
				this.nextSym = this.nextSymbol();
			}
			return this.nextSym != null;
		}

		/**
		 * Looks for the next symbol in the stream.
		 * 
		 * @return the symbol, or {@code null} if there isn't one.
		 */
		private Symbol nextSymbol() {
			if (SymbolListFormatter.this.start != null) {
				// Skip through until hit start symbol.
				boolean foundStart = false;
				boolean inEscape = false;
				while (this.hasNextChar() && !foundStart) {
					char ch = this.nextChar();
					// Take note of escape character states.
					if (inEscape) {
						inEscape = false;
					} else {
						if (ch == SymbolListFormatter.this.escapeChar) {
							inEscape = true;
						} else if (ch == SymbolListFormatter.this.startChar) {
							foundStart = true;
						}
					}
				}
			}
			// If hit end of stream, return null.
			if (!this.hasNextChar()) {
				return null;
			}
			// Read until hit sep or end.
			StringBuilder sb = new StringBuilder();
			boolean expectingEnd = SymbolListFormatter.this.end != null;
			boolean expectingSep = SymbolListFormatter.this.sep != null;
			boolean foundEnd = false;
			boolean foundSep = false;
			boolean inEscape = false;
			while (this.hasNextChar() && !(foundEnd || foundSep)) {
				char ch = this.nextChar();
				// Take note of escape character states.
				if (inEscape) {
					inEscape = false;
					sb.append(ch);
				} else {
					if ((expectingEnd || expectingSep)
							&& ch == SymbolListFormatter.this.escapeChar) {
						inEscape = true;
					} else {
						if (expectingEnd
								&& ch == SymbolListFormatter.this.endChar) {
							foundEnd = true;
						} else if (expectingSep
								&& ch == SymbolListFormatter.this.sepChar) {
							foundSep = true;
						} else {
							sb.append(ch);
							// No end/sep = single-char parsing.
							if (!expectingEnd && !expectingSep) {
								break;
							}
						}
					}
				}
			}
			// Take note of escape character states.
			// If expecting end but hit didn't find it, then treat as if found
			// nothing.
			if (expectingEnd && !foundEnd) {
				return null;
			}
			// Return symbol.
			return SymbolListFormatter.this.alpha.getCoder().decodeSymbol(
					sb.toString());
		}

		public Symbol next() {
			if (!this.hasNext()) {
				throw new NoSuchElementException();
			}
			Symbol currSym = this.nextSym;
			this.nextSym = this.nextSymbol();
			return currSym;
		}

		public void remove() {
			throw new UnsupportedOperationException(
					"Cannot remove symbols from input iterator.");
		}
	}

	/**
	 * Converts an iterator into a symbol list.
	 * 
	 * @param iter
	 *            the iterator.
	 * @return a symbol list.
	 */
	private List<Symbol> makeList(Iterator<Symbol> iter) {
		List<Symbol> list = new ArrayList<Symbol>();
		this.copyIntoList(iter, list);
		return list;
	}

	/**
	 * Copies symbols from an iterator into a symbol list.
	 * 
	 * @param iter
	 *            the iterator.
	 * @param list
	 *            the list to copy into.
	 */
	private void copyIntoList(Iterator<Symbol> iter, List<Symbol> list) {
		while (iter.hasNext()) {
			list.add(iter.next());
		}
	}

	/**
	 * Convert the string into a symbol list.
	 * 
	 * @param str
	 *            the string.
	 * @return the symbol list.
	 */
	public List<Symbol> parseSymbolList(CharSequence str) {
		return this.makeList(this.parseSymbols(str));
	}

	/**
	 * Convert the input into a symbol list.
	 * 
	 * @param is
	 *            the input.
	 * @return the symbol list.
	 */
	public List<Symbol> parseSymbolList(InputStream is) {
		return this.makeList(this.parseSymbols(is));
	}

	/**
	 * Convert the input into a symbol list.
	 * 
	 * @param r
	 *            the input.
	 * @return the symbol list.
	 */
	public List<Symbol> parseSymbolList(Reader r) {
		return this.makeList(this.parseSymbols(r));
	}

	/**
	 * Convert the string into a symbol list.
	 * 
	 * @param str
	 *            the string.
	 * @param symList
	 *            the list to append the symbols to.
	 */
	public void parseSymbolList(CharSequence str, List<Symbol> symList) {
		this.copyIntoList(this.parseSymbols(str), symList);
	}

	/**
	 * Convert the input into a symbol list.
	 * 
	 * @param is
	 *            the input.
	 * @param symList
	 *            the list to append the symbols to.
	 */
	public void parseSymbolList(InputStream is, List<Symbol> symList) {
		this.copyIntoList(this.parseSymbols(is), symList);
	}

	/**
	 * Convert the input into a symbol list.
	 * 
	 * @param r
	 *            the input.
	 * @param symList
	 *            the list to append the symbols to.
	 */
	public void parseSymbolList(Reader r, List<Symbol> symList) {
		this.copyIntoList(this.parseSymbols(r), symList);
	}
}
