package org.biojava.bio.structure.align.ce;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.align.AbstractStructureAlignment;
import org.biojava.bio.structure.align.StructureAlignment;
import org.biojava.bio.structure.align.ce.sidechain.CeSideChainCalculator;
import org.biojava.bio.structure.gui.AlignmentGui;
import org.biojava.bio.structure.align.model.AFPChain;

public class CeSideChainMain  extends AbstractStructureAlignment implements StructureAlignment {

	public static final String algorithmName = "jCE-sidechain";

	private static final String version = "2.3";

	private CeParameters params;

	Atom[] ca2clone;
	

	public CeSideChainMain(){
		super();
		params = new CeParameters();
	}
	
	public static void main(String[] args){

		CeSideChainMain ce = new CeSideChainMain();
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

		CeSideChainUserArgumentProcessor processor = new CeSideChainUserArgumentProcessor();
		processor.process(args);

	}

	public AFPChain align(Atom[] ca1, Atom[] ca2, Object param) throws StructureException{

		if ( ! (param instanceof CeParameters))
			throw new IllegalArgumentException("CE algorithm needs an object of call CeParameters as argument.");

		params = (CeParameters) param;

		System.out.println("Using atoms : " + params.getAlignmentAtoms());
		
		AFPChain afpChain = new AFPChain();
		afpChain.setCa1Length(ca1.length);
		afpChain.setCa2Length(ca2.length);

		// we don;t want to rotate input atoms, do we?
		ca2clone = new Atom[ca2.length];

		
		// this implementaiton allows to use any set of Atom
		// i.e. CA, CB, N, C,O atoms
		// but it is ok to put in only CA atoms,
		// the other atoms will be fetched from the ca.getParent() group
		
		int pos = 0;
		for (Atom a : ca2){
			Group g = (Group)a.getParent().clone();

			ca2clone[pos] = g.getAtom("CA");

			pos++;
		}

	
		CeSideChainCalculator calculator = new CeSideChainCalculator();

		calculator.extractFragments(params, afpChain,ca1, ca2clone);
		calculator.traceFragmentMatrix(params, afpChain,ca1, ca2clone);
		calculator.nextStep(params, afpChain,ca1, ca2clone);

		afpChain.setAlgorithmName(algorithmName);
		afpChain.setVersion(version);

		return afpChain;

	}



	public AFPChain align(Atom[] ca1, Atom[] ca2) throws StructureException {
		CeSideChainUserArgumentProcessor proc = new CeSideChainUserArgumentProcessor();
		Object params = proc.getParameters();
				
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
		params = (CeParameters) params;
	}
	
	public String getVersion() {
		return version;
	}
	
}
