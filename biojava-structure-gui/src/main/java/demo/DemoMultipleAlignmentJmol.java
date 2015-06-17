package demo;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.gui.MultipleAlignmentDisplay;
import org.biojava.nbio.structure.align.multiple.Block;
import org.biojava.nbio.structure.align.multiple.BlockImpl;
import org.biojava.nbio.structure.align.multiple.BlockSet;
import org.biojava.nbio.structure.align.multiple.BlockSetImpl;
import org.biojava.nbio.structure.align.multiple.MultipleAlignment;
import org.biojava.nbio.structure.align.multiple.MultipleAlignmentEnsemble;
import org.biojava.nbio.structure.align.multiple.MultipleAlignmentEnsembleImpl;
import org.biojava.nbio.structure.align.multiple.MultipleAlignmentImpl;
import org.biojava.nbio.structure.align.multiple.MultipleAlignmentScorer;
import org.biojava.nbio.structure.align.multiple.MultipleAlignmentWriter;
import org.biojava.nbio.structure.align.multiple.MultipleSuperimposer;
import org.biojava.nbio.structure.align.multiple.ReferenceSuperimposer;
import org.biojava.nbio.structure.align.util.AtomCache;

/**
 * Demo for visualizing the results of a Multiple Alignment, from a sample MultipleAlignment object.
 * 
 * @author Aleix Lafita
 * 
 */
public class DemoMultipleAlignmentJmol {

	public static void main(String[] args) throws IOException, StructureException {
		
		//Specify the structures to align
		//List<String> names = Arrays.asList("1tim.a", "1vzw", "1nsj", "3tha.a");	//TIM barrels
		List<String> names = Arrays.asList("1mbc", "1hlb", "1thb.a", "1ith.a");		//globins
		
		//Load the CA atoms of the structures
		AtomCache cache = new AtomCache();
		List<Atom[]> atomArrays = new ArrayList<Atom[]>();
		for (String name:names) atomArrays.add(cache.getAtoms(name));
		
		//Here the multiple structural alignment algorithm comes in place to generate the alignment object
		MultipleAlignment fakeMultAln = fakeMultipleAlignment("globins", atomArrays);
		fakeMultAln.getEnsemble().setStructureNames(names);
		
		//Generate a pairwise alignment and convert it to a MultipleAlignment
		//FatCat fatcat  = new FatCat();
		//AFPChain afpChain = fatcat.alignRigid(atomArrays.get(0),atomArrays.get(1));
		//MultipleAlignmentEnsemble ensemble = new MultipleAlignmentEnsembleImpl(afpChain, atomArrays.get(0),atomArrays.get(1));
		//MultipleAlignment pairwise = ensemble.getMultipleAlignments().get(0);
		
		//System.out.println(MultipleAlignmentWriter.toFASTA(fakeMultAln));
		MultipleAlignmentDisplay.display(fakeMultAln);
		//StructureAlignmentDisplay.display(pairwise);
		//For comparison display the original AFP
		//StructureAlignmentDisplay.display(afpChain,atomArrays.get(0),atomArrays.get(1));
	}
	
	/**
	 * Method that constructs a fake MultipleAlignment with two BlockSets, with two and one Blocks respectively. Used to
	 * test the correctness of the DataStructure. In the future it will be in the Test packages.
	 * @param family name of the protein family
	 * @param atomArrays
	 * @return MultipleAlignment
	 * @throws StructureException
	 * @throws StructureAlignmentException
	 */
	private static MultipleAlignment fakeMultipleAlignment(String family, List<Atom[]>atomArrays) throws StructureException {
		
		//Initialize the multiple alignment parent ensemble
		MultipleAlignmentEnsemble ensemble = new MultipleAlignmentEnsembleImpl();
		ensemble.setAtomArrays(atomArrays);
		ensemble.setAlgorithmName("fakeAlgorithm");
		ensemble.setVersion("1.0");
		MultipleAlignment fakeMultAln = new MultipleAlignmentImpl(ensemble);
		
		if (family == "globins"){
			
			BlockSet blockSet1 = new BlockSetImpl(fakeMultAln); //first BlockSet with 2 Blocks
			BlockSet blockSet2 = new BlockSetImpl(fakeMultAln); //second BlockSet with 1 Block
			
			Block block1 = new BlockImpl(blockSet1);
			Block block2 = new BlockImpl(blockSet1);
			
			Block block3 = new BlockImpl(blockSet2);
			
			//Alignment obtained from MUSTANG multiple alignment (just some of the residues, not the whole alignment)
			List<Integer> aligned11 = Arrays.asList(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21);
			List<Integer> aligned12 = Arrays.asList(29,30,31,32,33,34,35,36,38);
			List<Integer> aligned13 = Arrays.asList(123,124,125,126,127,128,129,130,131,132,133,134);
			
			List<Integer> aligned21 = Arrays.asList(10,11,12,13,null,15,16,17,null,19,20,21,22,23,24,25,null,27,28,29,30,null);
			List<Integer> aligned22 = Arrays.asList(39,40,41,42,43,44,45,46,48);
			List<Integer> aligned23 = Arrays.asList(133,134,135,136,137,138,139,140,141,142,143,144);
			
			List<Integer> aligned31 = Arrays.asList(0,1,2,3,null,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21);
			List<Integer> aligned32 = Arrays.asList(29,30,31,32,33,34,35,36,38);
			List<Integer> aligned33 = Arrays.asList(117,118,119,120,121,122,123,124,125,126,127,128);
			
			List<Integer> aligned41 = Arrays.asList(0,1,2,3,null,5,6,7,8,9,10,11,12,13,14,15,null,17,18,19,20,21);
			List<Integer> aligned42 = Arrays.asList(30,31,32,33,34,35,36,37,39);
			List<Integer> aligned43 = Arrays.asList(121,122,123,124,125,126,127,128,129,130,131,132);
			
			block1.getAlignRes().add(aligned11);
			block1.getAlignRes().add(aligned21);
			block1.getAlignRes().add(aligned31);
			block1.getAlignRes().add(aligned41);
			
			block2.getAlignRes().add(aligned12);
			block2.getAlignRes().add(aligned22);
			block2.getAlignRes().add(aligned32);
			block2.getAlignRes().add(aligned42);
			
			block3.getAlignRes().add(aligned13);
			block3.getAlignRes().add(aligned23);
			block3.getAlignRes().add(aligned33);
			block3.getAlignRes().add(aligned43);
			
			//Calculating all information in the alignment is as easy as that line, once the residue equivalencies are set
			MultipleSuperimposer imposer = new ReferenceSuperimposer();
			imposer.superimpose(fakeMultAln);
			MultipleAlignmentScorer.calculateScores(fakeMultAln);
		}
		return fakeMultAln;
	}

}
