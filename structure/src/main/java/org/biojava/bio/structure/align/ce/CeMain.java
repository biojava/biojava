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

package org.biojava.bio.structure.align.ce;


import org.biojava.bio.structure.Atom;

import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.align.AbstractStructureAlignment;
import org.biojava.bio.structure.align.StructureAlignment;


import org.biojava.bio.structure.align.model.AFPChain;

/** The main class of the Java implementation of the Combinatorial Extension Algorithm (CE).
 * 
 * The original CE paper is available from here: <a href="http://peds.oxfordjournals.org/cgi/content/short/11/9/739">http://peds.oxfordjournals.org/cgi/content/short/11/9/739</a>
 * 
 * For a demo of how to use this algorithm, visit the BioJava web site:
 * <a href="">CE usage example</a>.
 * 
 * @author Andreas Prlic.
 *
 */
public class CeMain extends AbstractStructureAlignment implements StructureAlignment {

	public static final String algorithmName = "jCE";

	public static final String version = "1.0";

	private CeParameters params;

	Atom[] ca2clone;

	public CeMain(){
		super();
		params = new CeParameters();
	}
	

	public static void main(String[] args){

		CeMain ce = new CeMain();
		if (args.length  == 0 ) {			
			System.out.println(ce.printHelp());
			return;			
		}

		if ( args.length == 1){
			if (args[0].equalsIgnoreCase("-h") || args[0].equalsIgnoreCase("-help")|| args[0].equalsIgnoreCase("--help")){
				System.out.println(ce.printHelp());								
				return;
			}
			
		}

		CeUserArgumentProcessor processor = new CeUserArgumentProcessor();
		processor.process(args);

	}

	public AFPChain align(Atom[] ca1, Atom[] ca2, Object param) throws StructureException{

		if ( ! (param instanceof CeParameters))
			throw new IllegalArgumentException("CE algorithm needs an object of call CeParameters as argument.");

		params = (CeParameters) param;

		AFPChain afpChain = new AFPChain();
		afpChain.setCa1Length(ca1.length);
		afpChain.setCa2Length(ca2.length);

		// we don;t want to rotate input atoms, do we?
		ca2clone = new Atom[ca2.length];

		int pos = 0;
		for (Atom a : ca2){
			Group g = (Group)a.getParent().clone();

			ca2clone[pos] = g.getAtom("CA");

			pos++;
		}

		CECalculator calculator = new CECalculator();

		calculator.extractFragments(params, afpChain,ca1, ca2clone);
		calculator.traceFragmentMatrix(params, afpChain,ca1, ca2clone);
		calculator.nextStep(params, afpChain,ca1, ca2clone);

		afpChain.setAlgorithmName(algorithmName);
		afpChain.setVersion(version);

		return afpChain;

	}



	public AFPChain align(Atom[] ca1, Atom[] ca2) throws StructureException {

		if (params == null)
			params = new CeParameters();
		
		return align(ca1,ca2,params);
	}

	public String getAlgorithmName() {

		return algorithmName;
	}

	public ConfigStrucAligParams getParameters() {

		return params;
	}
	
	public void setParameters(ConfigStrucAligParams params){
		if (! (params instanceof CeParameters )){
			throw new IllegalArgumentException("provided parameter object is not of type CeParameter");
		}
		this.params = (CeParameters) params;
	}
	
	public String getVersion() {
		return version;
	}
}
