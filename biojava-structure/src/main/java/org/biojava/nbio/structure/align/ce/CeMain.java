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
 * Created on Sep 15, 2009
 * Author: Andreas Prlic
 *
 */

package org.biojava.nbio.structure.align.ce;


import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.AbstractStructureAlignment;
import org.biojava.nbio.structure.align.StructureAlignment;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.jama.Matrix;

/** 
 * The main class of the Java implementation of the Combinatorial Extension Algorithm (CE),
 * as has been originally developed by I. Shindyalov and P.Bourne (1998).
 * The original CE paper is available from here: <a href="http://peds.oxfordjournals.org/cgi/content/short/11/9/739">http://peds.oxfordjournals.org/cgi/content/short/11/9/739</a>
 *
 * For a demo of how to use this algorithm, visit the BioJava web site:
 * <a href="">CE usage example</a>.
 *
 * The BioJava CE version is based on CE version 2.3 (2003 or 2004).
 *
 * @author Andreas Prlic.
 *
 */
public class CeMain extends AbstractStructureAlignment implements StructureAlignment {

	public static final String algorithmName = "jCE";

	/**
	 *  version history:
	 *  1.2 - Added more parameters to the command line, including -maxOptRMSD
	 *  1.1 - Additional parameters
	 *  1.0 - Initial port from C code
	 */
	public static final String version = "1.2";

	protected CeParameters params;
	protected CECalculator calculator;
	private Atom[] ca2clone;

	public CeMain(){
		super();
		params = new CeParameters();
		calculator = new CECalculator(params);
	}


	/**
	 * Example Parameters:
	 *
	 * -pdbFilePath /tmp -autoFetch -printCE -pdb1 1cnv -pdb2 3cna
	 *
	 */
	public static void main(String[] args) throws Exception {
		CeUserArgumentProcessor processor = new CeUserArgumentProcessor(); //Responsible for creating a CeMain instance
		processor.process(args);
	}

	/**
	 * Align ca2 onto ca1.
	 */
	@Override
	public AFPChain align(Atom[] ca1, Atom[] ca2, Object param) throws StructureException{
		if ( ! (param instanceof CeParameters))
			throw new IllegalArgumentException("CE algorithm needs an object of call CeParameters as argument.");

		params = (CeParameters) param;

		// we don't want to rotate input atoms, do we?
		ca2clone = new Atom[ca2.length];

		int pos = 0;
		for (Atom a : ca2){
			Group g = (Group)a.getGroup().clone(); // works because each group has only a CA atom

			ca2clone[pos] = g.getAtom(a.getName());

			pos++;
		}

		calculator = new CECalculator(params);

		//Build alignment ca1 to ca2-ca2
		AFPChain afpChain = new AFPChain(algorithmName);
		afpChain = calculator.extractFragments(afpChain, ca1, ca2clone);
		calculator.traceFragmentMatrix( afpChain,ca1, ca2clone);
		calculator.nextStep( afpChain,ca1, ca2clone);

		afpChain.setAlgorithmName(getAlgorithmName());
		afpChain.setVersion(version);

		// Try to guess names

		if (ca1.length!=0 && ca1[0].getGroup().getChain()!=null && ca1[0].getGroup().getChain().getStructure()!=null)
			afpChain.setName1(ca1[0].getGroup().getChain().getStructure().getName());

		if (ca2.length!=0 && ca2[0].getGroup().getChain()!=null && ca2[0].getGroup().getChain().getStructure()!=null)
			afpChain.setName2(ca2[0].getGroup().getChain().getStructure().getName());

		if ( afpChain.getNrEQR() == 0)
		   return afpChain;

		// Set the distance matrix

		int winSize = params.getWinSize();
		int winSizeComb1 = (winSize-1)*(winSize-2)/2;
		double[][] m = calculator.initSumOfDistances(ca1.length, ca2.length, winSize, winSizeComb1, ca1, ca2clone);
		afpChain.setDistanceMatrix(new Matrix(m));
		afpChain.setSequentialAlignment(true);

		return afpChain;
	}




	@Override
	public AFPChain align(Atom[] ca1, Atom[] ca2) throws StructureException {

		if (params == null)
			params = new CeParameters();

		return align(ca1,ca2,params);
	}

	@Override
	public String getAlgorithmName() {

		return CeMain.algorithmName;
	}

	@Override
	public ConfigStrucAligParams getParameters() {

		return params;
	}

	@Override
	public void setParameters(ConfigStrucAligParams params){
		if (! (params instanceof CeParameters )){
			throw new IllegalArgumentException("provided parameter object is not of type CeParameter");
		}
		this.params = (CeParameters) params;
	}

	@Override
	public String getVersion() {
		return CeMain.version;
	}

	public CECalculator getCECalculator() {
		return calculator;
	}
}
