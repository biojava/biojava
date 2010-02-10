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

import java.util.HashMap;
import java.util.Map;

import org.biojava.bio.seq.io.ParseException;

/**
 * Listens to events fired by the Nexus parser. Use these events to handle data
 * directly or construct objects.
 * 
 * @author Richard Holland
 * @author Tobias Thierer
 * @author Jim Balhoff
 * @since 1.6
 */
public interface NexusFileListener {

	/**
	 * About to start a new file.
	 */
	public void startFile();

	/**
	 * Finished reading a file.
	 */
	public void endFile();

	/**
	 * Opening a comment tag.
	 */
	public void beginComment();

	/**
	 * Receiving free text inside a comment tag.
	 * 
	 * @param comment
	 *            the text of the comment.
	 */
	public void commentText(String comment) throws ParseException;

	/**
	 * Closing a comment tag.
	 */
	public void endComment();

	/**
	 * Closing a line (semi-colon encountered). This indicates that anything
	 * received after it is on the next logical line of the file.
	 */
	public void endTokenGroup();

	/**
	 * Causes the default block parsers to be assigned. This is called by the
	 * constructor of the abstract implementation. If it is not called, then at
	 * least the unknown block parser must be set by other means.
	 */
	public void setDefaultBlockParsers();

	/**
	 * Sets the parser to use for a given block.
	 * 
	 * @param blockName
	 *            the name of the block.
	 * @param parser
	 *            the parser to use. Use <tt>null</tt> to unset an existing
	 *            one and use the default one for that block instead.
	 */
	public void setBlockParser(String blockName, NexusBlockParser parser);

	/**
	 * Gets the parser to use for a given block.
	 * 
	 * @param blockName
	 *            the name of the block. return parser the parser to use. Is
	 *            never null.
	 */
	public NexusBlockParser getBlockParser(String blockName);

	/**
	 * About to start a new block.
	 * 
	 * @param blockName
	 *            the name of the new block.
	 */
	public void startBlock(String blockName);

	/**
	 * Finished reading a block.
	 */
	public void endBlock();

	/**
	 * Encountered a token.
	 * 
	 * @param token
	 *            the token.
	 * @throws ParseException
	 *             if the token is invalid.
	 */
	public void parseToken(String token) throws ParseException;

	/**
	 * Does the listener want to know about brackets and braces as separate
	 * tokens?
	 * 
	 * @return <tt>true</tt> if it does.
	 */
	public boolean wantsBracketsAndBraces();

	/**
	 * Example abstract implementation which all others should extend.
	 */
	public abstract class Abstract implements NexusFileListener {

		private static final NexusBlockParser ignoreUnknownBlocks = new NexusBlockParser.Abstract(
				new NexusBlockListener() {
					public void beginComment() {
					}

					public void commentText(String comment)
							throws ParseException {
					}

					public void endBlock() {
					}

					public void endComment() {
					}

					public void endTokenGroup() {
					}

					public void startBlock(final String blockName) {
					}
				}) {
			public void resetStatus() {
			}

			public void parseToken(final String token) throws ParseException {
			}
		};

		private Map blockParsers = new HashMap();

		private NexusBlockParser blockParser;

		public Abstract() {
			this.setDefaultBlockParsers();
		}

		public void beginComment() {
			if (this.blockParser != null)
				this.blockParser.beginComment();
			else
				this.beginFileComment();
		}

		/**
		 * This method will get called when a comment is started on the file,
		 * and not any block within it.
		 */
		protected abstract void beginFileComment();

		public void commentText(String comment) throws ParseException {
			if (this.blockParser != null)
				this.blockParser.commentText(comment);
			else
				this.fileCommentText(comment);
		}

		/**
		 * This method will get called when comment text is found on the file,
		 * and not any block within it.
		 * 
		 * @param comment
		 *            the comment text.
		 */
		protected abstract void fileCommentText(String comment);

		public void endComment() {
			if (this.blockParser != null)
				this.blockParser.endComment();
			else
				this.endFileComment();
		}

		/**
		 * This method will get called when a comment is ended on the file, and
		 * not any block within it.
		 */
		protected abstract void endFileComment();

		public void endBlock() {
			this.blockParser.endBlock();
			this.blockEnded(this.blockParser);
			this.blockParser = null;
		}

		/**
		 * This method gets called when the block parser is expected to have
		 * finished parsing a block.
		 * 
		 * @param blockParser
		 *            the parser that has finished.
		 */
		protected abstract void blockEnded(NexusBlockParser blockParser);

		public boolean wantsBracketsAndBraces() {
			return this.blockParser != null
					&& this.blockParser.wantsBracketsAndBraces();
		}

		public void setDefaultBlockParsers() {
			this.setBlockParser(NexusBlockParser.UNKNOWN_BLOCK,
					NexusFileListener.Abstract.ignoreUnknownBlocks);
		}

		public NexusBlockParser getBlockParser(String blockName) {
			blockName = blockName.toUpperCase();
			return blockParser = this.blockParsers.containsKey(blockName) ? (NexusBlockParser) this.blockParsers
					.get(blockName)
					: (NexusBlockParser) this.blockParsers
							.get(NexusBlockParser.UNKNOWN_BLOCK);
		}

		public void endTokenGroup() {
			// Only blocks care about semi-colons.
			if (this.blockParser != null)
				this.blockParser.endTokenGroup();
		}

		public void parseToken(String token) throws ParseException {
			// Only blocks can parse tokens.
			if (this.blockParser != null)
				this.blockParser.parseToken(token);
		}

		public void setBlockParser(String blockName, NexusBlockParser parser) {
			this.blockParsers.put(blockName.toUpperCase(), parser);
		}

		public void startBlock(String blockName) {
			this.blockParser = this.getBlockParser(blockName);
			this.blockParser.startBlock(blockName);
		}

	}
}
