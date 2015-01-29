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
 * Created on Sep 16, 2009
 * Author: Andreas Prlic 
 *
 */

package org.biojava.bio.structure.align.client;

public class JobKillException extends Exception {

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	public JobKillException(String message){
		super(message);
	}
}
