package org.biojava.bio.structure.align.seq;


import org.biojava.bio.structure.align.StructureAlignment;
import org.biojava.bio.structure.align.ce.AbstractUserArgumentProcessor;
import org.biojava.bio.structure.align.ce.StartupParameters;



public class SmithWatermanUserArgumentProcessor extends AbstractUserArgumentProcessor{


	protected static class SmithWatermanStartupParams extends StartupParameters {
		private short gapOpen;
		private short gapExtend;

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

		@Override
		public String toString() {
			StringBuilder builder = new StringBuilder();
			builder.append("SmithWatermanStartupParams [gapOpen=")
			.append(gapOpen).append(", gapExtend=").append(gapExtend)
			.append("]");
			return builder.toString();
		}
	}


	public StructureAlignment getAlgorithm() {
		return new SmithWaterman3Daligner();
	}



	@Override
	public Object getParameters() {
		StructureAlignment alignment = getAlgorithm();

		SmithWaterman3DParameters p = (SmithWaterman3DParameters) alignment.getParameters();
		SmithWatermanStartupParams startup = (SmithWatermanStartupParams)params;

		if ( p == null)
			p = new SmithWaterman3DParameters();

		p.setGapExtend(startup.getGapExtend());
		p.setGapOpen(startup.getGapOpen());

		return p;
	}

	public String getDbSearchLegend(){
		String legend = "# name1\tname2\tscore\tz-score\trmsd\tlen1\tlen2\tsim1\tsim2\t " ;
		return legend;
	}



	@Override
	protected StartupParameters getStartupParametersInstance() {
		return new SmithWatermanStartupParams();
	}



}
