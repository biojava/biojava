package demo;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.biojava.bio.structure.AsaCalculator;
import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.Element;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.GroupType;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava3.structure.StructureIO;

public class DemoAsa {

	private static final boolean hetAtoms = false;
	
	public static void main(String[] args) throws IOException, StructureException {

		String pdbCode = args[0];
		int numThreads = Integer.parseInt(args[1]);
		
		demoAsa(pdbCode, numThreads);
	}
	
	private static void demoAsa(String pdbCode, int numThreads) throws IOException, StructureException {
		
		AtomCache cache = new AtomCache("/tmp",false);
		cache.setUseMmCif(true);
		
		StructureIO.setAtomCache(cache); 
		
		Structure structure = StructureIO.getStructure(pdbCode);
			
		long start = System.currentTimeMillis();
		
		AsaCalculator asaCalc = new AsaCalculator(structure, 
				AsaCalculator.DEFAULT_PROBE_SIZE, 
				1000, numThreads, hetAtoms);
		
		double[] asas = asaCalc.calculateAsa();
		
		long end = System.currentTimeMillis();

		
		double tot = 0;
		
		//int i = 0;
		
		Atom[] atoms = StructureTools.getAllNonHAtomArray(structure, hetAtoms);
		
		for (int i=0;i<atoms.length;i++) {
			tot+=asas[i];
		}
		
//		for (Chain chain:structure.getChains()) {
//			for (Group res:chain.getAtomGroups()) {
//				//List<Group> alllocs = new ArrayList<Group>();// = alt.getAltLocs();
//				//alllocs.add(alt);
//				//for(Group res: alllocs) {
//					if (!hetAtoms && res.getType().equals(GroupType.HETATM)) continue;
//
//					//System.out.printf("%1s\t%3s\t%s\n",chain.getChainID(),res.getResidueNumber().toString(),res.getPDBName());
//					double perResidue = 0;
//					for (Atom atom:res.getAtoms()) {
//						if (atom.getElement()!=Element.H) {
//							System.out.println("    "+atom.getName()+" "+atom.getAltLoc()+ " "+atom.getPDBserial());
//							perResidue += asas[i];
//							i++;
//						}
//
//					}
//					System.out.printf("%1s\t%3d\t%s\t%6.2f\n",chain.getChainID(),res.getResidueNumber().getSeqNum(),res.getPDBName(),perResidue);
//					tot+=perResidue;
//				//}
//			}
//		}
		
		
		System.out.printf("Total area: %9.2f\n",tot);
		System.out.printf("Time: %4.1fs\n",((end-start)/1000.0));
		
		
		System.out.println("Testing scaling: ");
		double[] runTimes = new double[numThreads];
		for (int nThreads=1;nThreads<=numThreads;nThreads++) {
			start = System.currentTimeMillis();
			
			asaCalc = new AsaCalculator(structure, 
					AsaCalculator.DEFAULT_PROBE_SIZE, 
					1000, numThreads, hetAtoms);
			
			asas = asaCalc.calculateAsa();
			
			end = System.currentTimeMillis();
			runTimes[nThreads-1] = (end-start)/1000.0;
			
		}
		for (int nThreads=1;nThreads<=numThreads;nThreads++) {
			System.out.printf(nThreads+" threads, time: %4.1fs -- x%2.1f\n",runTimes[nThreads-1],runTimes[0]/runTimes[nThreads-1]);
		}
		
		
		
	}

	
}
