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
 * Created on Mar 9, 2010
 * Author: Spencer Bliven 
 *
 */

package org.biojava.bio.structure.align.ce;


import org.biojava.bio.structure.align.ce.CeMain;

/** 
 * A wrapper for {@link CeMain} which sets default parameters to be appropriate for finding
 * circular permutations
 * 
 * @author Spencer Bliven.
 *
 */
public class CeCPMain extends CeMain {

	public static final String algorithmName = "jCE-CP";

	public static final String version = "1.0";

	public CeCPMain(){
		super();
		this.params.setCheckCircular(true);
		this.params.setMaxGapSize(0);
	}
	
	public String getAlgorithmName() {

		return algorithmName;
	}
	
	public String getVersion() {
		return version;
	}

}
