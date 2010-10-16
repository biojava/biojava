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

/**
 * Represents Nexus data blocks.
 * 
 * @author Richard Holland
 * @author Tobias Thierer
 * @author Jim Balhoff
 * @since 1.6
 */
public class DataBlock extends CharactersBlock {

	/**
	 * A constant representing the name of Data blocks.
	 */
	public static final String DATA_BLOCK = "DATA";

	/**
	 * Delegates to NexusBlock.Abstract constructor using DataBlock.DATA_BLOCK
	 * as the name.
	 */
	public DataBlock() {
		super(DataBlock.DATA_BLOCK);
	}

}
