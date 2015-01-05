package org.biojava.bio.structure.align.ce;

import org.biojava.bio.structure.align.StructureAlignment;


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

	public static void main(String[] args){
		CeSideChainUserArgumentProcessor processor = new CeSideChainUserArgumentProcessor();
		processor.process(args);
	}

	//   public AFPChain align(Atom[] ca1, Atom[] ca2, Object param) throws StructureException{
	//
	//      if ( ! (param instanceof CeParameters))
	//         throw new IllegalArgumentException("CE algorithm needs an object of call CeParameters as argument.");
	//
	//      params = (CeParameters) param;
	//
	//      System.out.println("Running alignment with side chain directions");
	//
	//      AFPChain afpChain = new AFPChain();
	//      afpChain.setCa1Length(ca1.length);
	//      afpChain.setCa2Length(ca2.length);
	//
	//      // we don;t want to rotate input atoms, do we?
	//      ca2clone = new Atom[ca2.length];
	//
	//      // this implementaiton allows to use any set of Atom
	//      // i.e. CA, CB, N, C,O atoms
	//      // but it is ok to put in only CA atoms,
	//      // the other atoms will be fetched from the ca.getParent() group
	//
	//      int pos = 0;
	//      for (Atom a : ca2){
	//         Group g = (Group)a.getParent().clone();
	//
	//         ca2clone[pos] = g.getAtom(StructureTools.caAtomName);
	//
	//         pos++;
	//      }
	//
	//      CECalculator calculator = new CECalculator(params);
	//
	//      calculator.extractFragments( afpChain,ca1, ca2clone);
	//      calculator.traceFragmentMatrix( afpChain,ca1, ca2clone);
	//      calculator.nextStep( afpChain,ca1, ca2clone);
	//
	//      afpChain.setAlgorithmName(algorithmName);
	//      afpChain.setVersion(version);
	//
	//      return afpChain;
	//
	//   }
	//
	//
	//
	//   public AFPChain align(Atom[] ca1, Atom[] ca2) throws StructureException {
	//      if (params == null){
	//         //CeSideChainUserArgumentProcessor proc = new CeSideChainUserArgumentProcessor();
	//         //params = (CeParameters) proc.getParameters();
	//         System.out.println("using default CE Parameters...");
	//         params = new CeParameters();
	//      }
	//      System.out.println("alignment with show ranges? " + params.isShowAFPRanges());
	//      return align(ca1,ca2,params);
	//   }

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
