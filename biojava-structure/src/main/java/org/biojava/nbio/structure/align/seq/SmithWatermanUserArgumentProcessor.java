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
package org.biojava.nbio.structure.align.seq;

import org.biojava.nbio.structure.align.StructureAlignment;
import org.biojava.nbio.structure.align.ce.AbstractUserArgumentProcessor;
import org.biojava.nbio.structure.align.ce.StartupParameters;

public class SmithWatermanUserArgumentProcessor extends AbstractUserArgumentProcessor{


	protected static class SmithWatermanStartupParams extends StartupParameters {
		
		private short gapOpen;
		private short gapExtend;
		private double maxRmsd;
		private int minLen;

		public SmithWatermanStartupParams() {
			super();
		}

		public short getGapOpen() {
			return gapOpen;
		}

		public void setGapOpen(short gapOpen) {
			this.gapOpen = gapOpen;
		}

		public short getGapExtend() {
			return gapExtend;
		}

		public void setGapExtend(short gapExtend) {
			this.gapExtend = gapExtend;
		}

		
		public double getMaxRmsd() {
			return maxRmsd;
		}

		public void setMaxRmsd(double maxRmsd) {
			this.maxRmsd = maxRmsd;
		}

		public int getMinLen() {
			return minLen;
		}

		public void setMinLen(int minLen) {
			this.minLen = minLen;
		}

		@Override
		public String toString() {
			StringBuilder builder = new StringBuilder();
			builder.append("SmithWatermanStartupParams [gapOpen=")
			.append(gapOpen).append(", gapExtend=").append(gapExtend)
			.append("]").append(", maxRmsd=").append(maxRmsd)
			.append(", minLen=").append(minLen).append("]");
			return builder.toString();
		}
	}


	@Override
	public StructureAlignment getAlgorithm() {
		return new SmithWaterman3Daligner();
	}



	@Override
	public Object getParameters() {
		StructureAlignment alignment = getAlgorithm();

		SmithWaterman3DParameters p = (SmithWaterman3DParameters) alignment.getParameters();
		SmithWatermanStartupParams startup = (SmithWatermanStartupParams) params;

		if ( p == null)
			p = new SmithWaterman3DParameters();

		p.setGapExtend(startup.getGapExtend());
		p.setGapOpen(startup.getGapOpen());
		p.setMaxRmsd(startup.getMaxRmsd());
		p.setMinLen(startup.getMinLen());

		return p;
	}

	@Override
	public String getDbSearchLegend(){
		String legend = "# name1\tname2\tscore\tz-score\trmsd\tlen1\tlen2\tsim1\tsim2\t " ;
		return legend;
	}



	@Override
	protected StartupParameters getStartupParametersInstance() {
		return new SmithWatermanStartupParams();
	}

}
