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
 * Created on Nov 2, 2009
 * Author: Andreas Prlic 
 *
 */

package org.biojava.bio.structure.align.ce;


import org.biojava.bio.structure.align.StructureAlignment;

/** process the arguments from command line
 * 
 * @author Andreas Prlic
 *
 */
public class CeUserArgumentProcessor extends AbstractUserArgumentProcessor {
	
	protected class CeStartupParams extends StartupParameters {
		//protected int maxGapSize; //TODO move from superclass
		protected int winSize;
		protected int scoringStrategy;
		//protected double maxOptRMSD; //TODO
		protected double gapOpen;
		protected double gapExtension;
		//protected boolean showAFPRanges; //TODO

		public CeStartupParams() {
			super();
			//maxGapSize = 30;
			winSize = 8;
			scoringStrategy = CeParameters.DEFAULT_SCORING_STRATEGY;
			//showAFPRanges = false;
			//maxOptRMSD = 99d;
			gapOpen = CeParameters.DEFAULT_GAP_OPEN;
			gapExtension = CeParameters.DEFAULT_GAP_EXTENSION;
		}

		public int getWinSize() {
			return winSize;
		}

		public void setWinSize(int winSize) {
			this.winSize = winSize;
		}

		public int getScoringStrategy() {
			return scoringStrategy;
		}

		public void setScoringStrategy(int scoringStrategy) {
			this.scoringStrategy = scoringStrategy;
		}

		public double getGapOpen() {
			return gapOpen;
		}

		public void setGapOpen(double gapOpen) {
			this.gapOpen = gapOpen;
		}

		public double getGapExtension() {
			return gapExtension;
		}

		public void setGapExtension(double gapExtension) {
			this.gapExtension = gapExtension;
		}

		@Override
		public String toString() {
			StringBuilder builder = new StringBuilder();
			builder.append("CeStartupParams [winSize=").append(winSize)
					.append(", scoringStrategy=").append(scoringStrategy)
					.append(", gapOpen=").append(gapOpen)
					.append(", gapExtension=").append(gapExtension)
					.append(", pdbFilePath=").append(pdbFilePath)
					.append(", cacheFilePath=").append(cacheFilePath)
					.append(", outFile=").append(outFile).append(", pdb1=")
					.append(pdb1).append(", pdb2=").append(pdb2)
					.append(", file1=").append(file1).append(", file2=")
					.append(file2).append(", showDBresult=")
					.append(showDBresult).append(", printXML=")
					.append(printXML).append(", printFatCat=")
					.append(printFatCat).append(", show3d=").append(show3d)
					.append(", autoFetch=").append(autoFetch)
					.append(", flexible=").append(flexible)
					.append(", pdbDirSplit=").append(pdbDirSplit)
					.append(", printCE=").append(printCE).append(", showMenu=")
					.append(showMenu).append(", showAFPRanges=")
					.append(showAFPRanges).append(", printPDB=")
					.append(printPDB).append(", maxOptRMSD=")
					.append(maxOptRMSD).append(", isDomainSplit=")
					.append(isDomainSplit).append(", alignPairs=")
					.append(alignPairs).append(", searchFile=")
					.append(searchFile).append(", saveOutputDir=")
					.append(saveOutputDir).append(", nrCPU=").append(nrCPU)
					.append(", maxGapSize=").append(maxGapSize).append("]");
			return builder.toString();
		}

	}
	
	@Override
	protected StartupParameters getStartupParametersInstance() {
		return new CeStartupParams();
	}

	@Override
	public StructureAlignment getAlgorithm() {
		return new CeMain();
	}


	@Override
	public Object getParameters() {
		
		StructureAlignment alignment = getAlgorithm();
		
		CeParameters aligParams = (CeParameters) alignment.getParameters();
		CeStartupParams startParams = (CeStartupParams) params;
		
		if ( aligParams == null)
			aligParams = new CECPParameters();
		
		// Copy relevant parameters from the startup parameters
		aligParams.setMaxGapSize(startParams.getMaxGapSize());
		aligParams.setWinSize(startParams.getWinSize());
		aligParams.setScoringStrategy(startParams.getScoringStrategy());
		aligParams.setMaxOptRMSD(startParams.getMaxOptRMSD());
		aligParams.setGapOpen(startParams.getGapOpen());
		aligParams.setGapExtension(startParams.getGapExtension());
		aligParams.setShowAFPRanges(startParams.isShowAFPRanges());
		return aligParams;
	}


	public String getDbSearchLegend(){
		//String legend = "# name1\tname2\tscore\tz-score\trmsd\tlen1\tlen2\tsim1\tsim2\t " ;
		//return legend;
		
		return "# name1\tname2\tscore\tz-score\trmsd\tlen1\tlen2\tcov1\tcov2\t%ID\tDescription\t " ;
		
		
	}
	

}
