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

import org.biojavax.bio.phylo.io.nexus.TreesBlock.NewickTreeString;

/**
 * Builds Nexus taxa blocks.
 * 
 * @author Richard Holland
 * @author Tobias Thierer
 * @author Jim Balhoff
 * @since 1.6
 */
public class TreesBlockBuilder extends NexusBlockBuilder.Abstract implements
		TreesBlockListener {

	private TreesBlock block;

	protected void addComment(final NexusComment comment) {
		this.block.addComment(comment);
	}

	protected NexusBlock startBlockObject() {
		this.block = new TreesBlock();
		this.resetStatus();
		return this.block;
	}

	private void resetStatus() {
		// Nothing to do.
	}

	public void endBlock() {
		// Don't care.
	}

	public void endTokenGroup() {
		// Nothing to do.
	}

	public void addTranslation(String label, String taxa) {
		this.block.addTranslation(label, taxa);
	}

	public void addTree(String label, NewickTreeString tree) {
		this.block.addTree(label, tree);
	}

}
