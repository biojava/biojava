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
/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.biojava.nbio.phylo;

/**
 * The TreeType specified the aligorithm used to construct the tree. Only suppo
 * for one algorithm is implemented, that is why this Enum is not used. It can
 * be used in the future as an interface to add more algorithm types.
 * 
 * @author willishf
 * @author Aleix Lafita
 * 
 */
public enum TreeType {

	/** Neighbor Joining Algorithm */
	NJ,

	/** Unweighted Pair Group Method with Arithmetic Mean */
	UPGMA,

	/** What does this stand for? */
	AV;
}
