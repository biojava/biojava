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
 * The TreeType specifies the optimization criteria used to generate the tree.
 *
 * @author Aleix Lafita
 * @since 4.1.1
 *
 */
public enum TreeType {

	/** Maximum Likelihood Tree */
	ML("ML-Tree"),

	/** Distance Tree */
	DISTANCE("Distance-Tree"),

	/** Parsimony Tree */
	PARSIMONY("Parsimony-Tree");

	/** Description name of the Tree Type */
	protected final String name;

	private TreeType(String name){
		this.name = name;
	}
}
