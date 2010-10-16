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
 * Listens to events from NexusBlockParser objects to create objects. This empty
 * interface needs to be implemented/extended for each block type supported.
 * 
 * @author Richard Holland
 * @author Tobias Thierer
 * @author Jim Balhoff
 * @since 1.6
 */
public interface NexusBlockListener {

	/**
	 * Notifies the parser that a new block is starting.
	 * 
	 * @param blockName
	 *            the name of the newly started block.
	 */
	public void startBlock(String blockName);

	/**
	 * Notifies the parser that a block is ending.
	 */
	public void endBlock();

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
}
