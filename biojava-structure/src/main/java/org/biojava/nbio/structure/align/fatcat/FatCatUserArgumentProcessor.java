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

package org.biojava.nbio.structure.align.fatcat;


import org.biojava.nbio.structure.align.StructureAlignment;
import org.biojava.nbio.structure.align.ce.AbstractUserArgumentProcessor;
import org.biojava.nbio.structure.align.ce.StartupParameters;
import org.biojava.nbio.structure.align.fatcat.calc.FatCatParameters;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;


public class FatCatUserArgumentProcessor extends AbstractUserArgumentProcessor {
	Logger logger = LoggerFactory.getLogger(FatCatUserArgumentProcessor.class);

	protected class FatCatStartupParams extends StartupParameters {
		int fragLen;
		Double rmsdCut;
		double disCut;
		int maxTra;
		boolean flexible;

		public FatCatStartupParams() {
			// Defaults should match those in FatCatParameters.reset()
			fragLen = FatCatParameters.DEFAULT_FRAGLEN;
			rmsdCut = 3.0;
			disCut = 5.0;
			maxTra = 5;
			flexible = false;
		}

		public int getFragLen() {
			return fragLen;
		}
		public void setFragLen(int fragLen) {
			this.fragLen = fragLen;
		}
		public Double getRmsdCut() {
			return rmsdCut;
		}
		public void setRmsdCut(Double rmsdCut) {
			this.rmsdCut = rmsdCut;
		}
		public double getDisCut() {
			return disCut;
		}
		public void setDisCut(double disCut) {
			this.disCut = disCut;
		}
		public int getMaxTra() {
			return maxTra;
		}
		public void setMaxTra(int maxTra) {
			this.maxTra = maxTra;
		}
		public boolean isFlexible() {
			return flexible;
		}
		public void setFlexible(boolean flexible) {
			this.flexible = flexible;
		}
	}

	@Override
	protected StartupParameters getStartupParametersInstance() {
		return new FatCatStartupParams();
	}

	@Override
	public StructureAlignment getAlgorithm() {
		StructureAlignment algorithm = null;
		if ( params != null && ((FatCatStartupParams)params).isFlexible()) {
			logger.info("running flexible alignment");
			algorithm = new FatCatFlexible();
		}
		else {
			logger.info("running rigid alignment");
			algorithm = new FatCatRigid();
		}
		return algorithm;

	}

	@Override
	public Object getParameters() {
		StructureAlignment alignment = getAlgorithm();

		FatCatParameters aligParams = (FatCatParameters) alignment.getParameters();
		FatCatStartupParams startParams = (FatCatStartupParams) params;

		if ( aligParams == null)
			aligParams = new FatCatParameters();

		aligParams.setFragLen(startParams.getFragLen());
		aligParams.setRmsdCut(startParams.getRmsdCut());
		aligParams.setDisCut(startParams.getDisCut());
		aligParams.setMaxTra(startParams.getMaxTra());

		return aligParams;
	}

	@Override
	public String getDbSearchLegend(){

		return "# name1\tname2\tscore\tprobability\trmsd\tlen1\tlen2\tcov1\tcov2\t%ID\tDescription\t " ;

	}

}
