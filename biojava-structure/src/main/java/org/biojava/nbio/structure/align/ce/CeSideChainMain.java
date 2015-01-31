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
package org.biojava.nbio.structure.align.ce;

import org.biojava.nbio.structure.align.StructureAlignment;
import org.biojava.nbio.structure.align.util.ConfigurationException;


public class CeSideChainMain  extends CeMain implements StructureAlignment {

	public static final String algorithmName = "jCE-sidechain";

	/**
	 *  version history:
	 *  2.4 - Added more parameters to the command line, including -maxOptRMSD
	 *  2.3 - Initial version
	 */
	private static final String version = "2.3";

	public CeSideChainMain(){
		super();

		if ( params == null) {
			CeSideChainUserArgumentProcessor proc = new CeSideChainUserArgumentProcessor();
			params = (CeParameters) proc.getParameters();
		}
	}

	public static void main(String[] args) throws ConfigurationException {
		CeSideChainUserArgumentProcessor processor = new CeSideChainUserArgumentProcessor();
		processor.process(args);
	}

	@Override
	public String getAlgorithmName() {

		return algorithmName;
	}

	@Override
	public ConfigStrucAligParams getParameters() {

		return params;
	}

	@Override
	public void setParameters(ConfigStrucAligParams params){
		System.out.println("setting params : " + params);
		if (! (params instanceof CeParameters )){
			throw new IllegalArgumentException("provided parameter object is not of type CeParameter");
		}
		this.params = (CeParameters) params;
	}

	@Override
	public String getVersion() {
		return version;
	}

}
