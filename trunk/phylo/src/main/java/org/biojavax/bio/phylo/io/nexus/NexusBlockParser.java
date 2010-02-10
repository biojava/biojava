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

import org.biojava.bio.seq.io.ParseException;

/**
 * Parses Nexus blocks. Each instance should parse tokens into events that can
 * be fired at some kind of NexusBlockListener. An Abstract parser is provided
 * from which all implementations should derive.
 * 
 * @author Richard Holland
 * @author Tobias Thierer
 * @author Jim Balhoff
 * @since 1.6
 */
public interface NexusBlockParser {

	/**
	 * The name for an unknown block parser.
	 */
	public static final String UNKNOWN_BLOCK = "__UNKNOWN";

	/**
	 * Notifies the parser that a new block is starting.
	 * 
	 * @param blockName
	 *            the name of the block.
	 */
	public void startBlock(String blockName);

	/**
	 * Notifies the parser that a block is ending.
	 */
	public void endBlock();

	/**
	 * Notifies the parser of the next token. Comment tokens will already have
	 * been parsed out and sent separately to the text() method of the listener.
	 * Quoted strings will have been parsed and underscores converted. What this
	 * token contains is the full string, after removal of quotes if necessary.
	 * The token will never be only whitespace.
	 * 
	 * @param token
	 *            the token to parse.
	 * @throws ParseException
	 *             if the token is unparseable.
	 */
	public void parseToken(String token) throws ParseException;

	/**
	 * Opening a comment tag.
	 */
	public void beginComment();

	/**
	 * Closing a comment tag.
	 */
	public void endComment();

	/**
	 * Closing a line (semi-colon encountered). This indicates that anything
	 * received after it is on the next logical line of the block.
	 */
	public void endTokenGroup();

	/**
	 * Receiving free text inside a comment tag.
	 * 
	 * @param comment
	 *            the text of the comment.
	 */
	public void commentText(String comment) throws ParseException;

	/**
	 * Obtain the listener for this parser.
	 * 
	 * @return the listener.
	 */
	public NexusBlockListener getBlockListener();

	/**
	 * Does the listener want to know about brackets and braces as separate
	 * tokens?
	 * 
	 * @return <tt>true</tt> if it does.
	 */
	public boolean wantsBracketsAndBraces();

	/**
	 * All block parsers should derive from this abstract parser.
	 */
	public abstract class Abstract implements NexusBlockParser {
		private NexusBlockListener blockListener;

		private String blockName;

		public Abstract(final NexusBlockListener blockListener) {
			this.blockListener = blockListener;
			this.blockName = null;
		}

		public NexusBlockListener getBlockListener() {
			return this.blockListener;
		}

		public void startBlock(final String blockName) {
			this.resetStatus();
			this.blockName = blockName;
			this.blockListener.startBlock(blockName);
		}

		/**
		 * This function is called when the parser is reset before starting a
		 * new block.
		 */
		protected abstract void resetStatus();

		protected String getBlockName() {
			return this.blockName;
		}

		public void endBlock() {
			this.blockName = null;
			this.blockListener.endBlock();
		}

		public void beginComment() {
			this.blockListener.beginComment();
		}

		public void endComment() {
			this.blockListener.endComment();
		}

		public void endTokenGroup() {
			this.blockListener.endTokenGroup();
		}

		public void commentText(final String comment) throws ParseException {
			this.blockListener.commentText(comment);
		}

		public abstract void parseToken(final String token)
				throws ParseException;

		public boolean wantsBracketsAndBraces() {
			return false;
		}
	}
}
