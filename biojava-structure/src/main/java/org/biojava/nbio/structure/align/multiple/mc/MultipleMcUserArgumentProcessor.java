package org.biojava.nbio.structure.align.multiple.mc;

import java.util.ArrayList;
import java.util.List;

import org.biojava.nbio.structure.align.MultipleStructureAligner;
import org.biojava.nbio.structure.align.ce.CeCPMain;
import org.biojava.nbio.structure.align.ce.UserArgumentProcessor;
import org.biojava.nbio.structure.align.util.ConfigurationException;

/** 
 * Process the arguments from the command line for {@link MultipleAlignment}s. 
 * <p>
 * Big TODO. Take as template {@link AbstractUserArgumentProcessor} but accept more than
 * one structure.
 * 
 * @author Aleix Lafita
 *
 */
public class MultipleMcUserArgumentProcessor implements UserArgumentProcessor {

	private MultipleStartupParameters params ;
	public static final List<String> mandatoryArgs= new ArrayList<String>();

	protected MultipleMcUserArgumentProcessor(){
		params = getStartupParametersInstance();
	}
	
	protected static class MultipleMcStartupParams extends MultipleStartupParameters {
		
		//Parameters to expose to the GUI
		private int randomSeed;
		private int minBlockLen;
		private int minAlignedStructures;
		private double gapOpen;
		private double gapExtension;
		private int convergenceSteps;
		private String pairwiseAlgorithm;
		
		public MultipleMcStartupParams() {
			super();
			randomSeed = 0;
			minBlockLen = 15;
			minAlignedStructures = 0;
			gapOpen = 10.0;
			gapExtension = 5.0;
			convergenceSteps = 0;
			pairwiseAlgorithm = CeCPMain.algorithmName;
		}

		@Override
		public String toString() {
			return "MultipleMcStartupParams [randomSeed=" + randomSeed
					+ ", minBlockLen=" + minBlockLen
					+ ", minAlignedStructures=" + minAlignedStructures
					+ ", gapOpen=" + gapOpen + ", gapExtension=" + gapExtension
					+ ", convergenceSteps=" + convergenceSteps
					+ ", pairwiseAlgorithm=" + pairwiseAlgorithm + "]";
		}

		public int getRandomSeed() {
			return randomSeed;
		}

		public void setRandomSeed(int randomSeed) {
			this.randomSeed = randomSeed;
		}

		public int getMinBlockLen() {
			return minBlockLen;
		}

		public void setMinBlockLen(int minBlockLen) {
			this.minBlockLen = minBlockLen;
		}

		public int getMinAlignedStructures() {
			return minAlignedStructures;
		}

		public void setMinAlignedStructures(int minAlignedStructures) {
			this.minAlignedStructures = minAlignedStructures;
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

		public int getConvergenceSteps() {
			return convergenceSteps;
		}

		public void setConvergenceSteps(int convergenceSteps) {
			this.convergenceSteps = convergenceSteps;
		}

		public String getPairwiseAlgorithm() {
			return pairwiseAlgorithm;
		}

		public void setPairwiseAlgorithm(String pairwiseAlgorithm) {
			this.pairwiseAlgorithm = pairwiseAlgorithm;
		}
		
	}
	
	protected MultipleStartupParameters getStartupParametersInstance() {
		return new MultipleMcStartupParams();
	}
	
	public MultipleStructureAligner getAlgorithm() {
		return new MultipleMcMain();
	}
	
	public Object getParameters() {
		
		MultipleStructureAligner alignment = getAlgorithm();
		
		MultipleMcParameters aligParams = (MultipleMcParameters) alignment.getParameters();
		MultipleMcStartupParams startParams = (MultipleMcStartupParams) params;
		
		if ( aligParams == null)
			aligParams = new MultipleMcParameters();
		
		// Copy relevant parameters from the startup parameters
		aligParams.setRandomSeed(startParams.getRandomSeed());
		aligParams.setMinBlockLen(startParams.getMinBlockLen());
		aligParams.setMinAlignedStructures(startParams.getMinAlignedStructures());
		aligParams.setGapOpen(startParams.getGapOpen());
		aligParams.setGapExtension(startParams.getGapExtension());
		aligParams.setConvergenceSteps(startParams.getConvergenceSteps());
		aligParams.setPairwiseAlgorithm(startParams.getPairwiseAlgorithm());
		
		return aligParams;
	}

	@Override
	public void process(String[] argv) throws ConfigurationException {
		//Big TODO		
	}

	@Override
	public String printHelp() {
		//Big TODO
		return null;
	}
}
