/**
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
 * Created on Apr 24, 2012
 * Created by Andreas Prlic
 *
 * @since 3.0.2
 */
package org.biojava.nbio.structure.align.ce;

import org.biojava.nbio.structure.align.StructureAlignment;
import org.biojava.nbio.structure.align.ce.CECPParameters.DuplicationHint;

public class CeCPUserArgumentProcessor extends CeUserArgumentProcessor {

	protected class CeCPStartupParams extends CeStartupParams {
		protected DuplicationHint duplicationHint;
		protected Integer minCPLength;

		public CeCPStartupParams() {
			duplicationHint = DuplicationHint.SHORTER;
			minCPLength = CECPParameters.DEFAULT_MIN_CP_LENGTH;
			maxGapSize = 0;
		}

		public DuplicationHint getDuplicationHint() {
			return duplicationHint;
		}

		public void setDuplicationHint(DuplicationHint duplicationHint) {
			this.duplicationHint = duplicationHint;
		}

		public Integer getMinCPLength() {
			return minCPLength;
		}

		public void setMinCPLength(Integer minCPLength) {
			this.minCPLength = minCPLength;
		}

		@Override
		public String toString() {
			StringBuilder builder = new StringBuilder();
			builder.append("CeCPStartupParams [duplicationHint=")
					.append(duplicationHint).append(", minCPLength=")
					.append(minCPLength).append("]");
			return builder.toString();
		}
	}

	@Override
	protected StartupParameters getStartupParametersInstance() {
		return  new CeCPStartupParams();
	}
	@Override
	public StructureAlignment getAlgorithm() {
		return new CeCPMain();
	}

	@Override
	public Object getParameters() {
		CECPParameters aligParams = (CECPParameters) super.getParameters();
		CeCPStartupParams startParams = (CeCPStartupParams) params;

		if ( aligParams == null)
			aligParams = new CECPParameters();

		// Copy relevant parameters from the startup parameters
		aligParams.setDuplicationHint(startParams.getDuplicationHint());
		aligParams.setMinCPLength(startParams.getMinCPLength());
		return aligParams;
	}
}

