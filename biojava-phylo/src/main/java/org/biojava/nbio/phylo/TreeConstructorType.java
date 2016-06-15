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
package org.biojava.nbio.phylo;

/**
 * The TreeConstructorType specifies the aligorithm used to construct the tree
 * (clustering algorithm). Only support for the NJ algorithm (from forester) is
 * currently implemented.
 *
 * @author Aleix Lafita
 * @since 4.1.1
 *
 */
public enum TreeConstructorType {

	/** Neighbor Joining Algorithm */
	NJ,

	/** Unweighted Pair-Group Method with Arithmetic mean */
	UPGMA,

	/** What does this stand for? (Aleix: Nov 2015) */
	AV;
}
