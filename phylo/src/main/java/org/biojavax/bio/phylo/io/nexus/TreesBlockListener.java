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
 * Listens to events that represent Nexus trees blocks.
 * 
 * @author Richard Holland
 * @author Tobias Thierer
 * @author Jim Balhoff
 * @since 1.6
 */
public interface TreesBlockListener extends NexusBlockListener {

	/**
	 * Add a translation.
	 * 
	 * @param label
	 *            the label to add.
	 * @param taxa
	 *            the taxa name this label will represent.
	 */
	public void addTranslation(final String label, final String taxa);

	/**
	 * Adds a tree.
	 * 
	 * @param label
	 *            the label to give the tree.
	 * @param tree
	 *            the tree to add.
	 */
	public void addTree(final String label, final NewickTreeString tree);
}
