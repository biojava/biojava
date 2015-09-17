package org.biojava.nbio.structure.symmetry.internal;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import javax.vecmath.Matrix4d;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.multiple.Block;
import org.biojava.nbio.structure.align.multiple.BlockImpl;
import org.biojava.nbio.structure.align.multiple.BlockSet;
import org.biojava.nbio.structure.align.multiple.BlockSetImpl;
import org.biojava.nbio.structure.align.multiple.MultipleAlignment;
import org.biojava.nbio.structure.align.multiple.MultipleAlignmentImpl;
import org.biojava.nbio.structure.align.multiple.util.MultipleAlignmentScorer;
import org.biojava.nbio.structure.symmetry.utils.SymmetryTools;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Iterative version of CeSymm that aims at identifying all symmetry axis 
 * (internal or quaternary) of a particular structure.
 * <p>
 * Works in the following way:
 * <ul><li>Run CeSymm on the original structure.
 * <li>Calculate the symmetric unit boundaries.
 * <li>Run CeSymm on one of the symmetric units to find further symmetries.
 * <li>Repeat the last two steps until no more significant results are found.
 * <li>Map back all residues in a multiple alignment of the subunits.
 * <li>Run a final optimization of all symmetric units correctly superimposed.
 * </ul></li>
 * 
 * @author Aleix Lafita
 * @since 4.2.0
 */
public class CeSymmIterative {
	
	private static Logger logger = 
			LoggerFactory.getLogger(CeSymmIterative.class);

	private CESymmParameters params;
	private MultipleAlignment msa;
	private Atom[] allAtoms;
	private String name;
	private SymmetryAxes axes;
	private List<List<Integer>> alignment;
	private List<MultipleAlignment> levels;

	/**
	 * For the iterative algorithm to work properly the refinement and 
	 * optimization options should be turned on, because the alignment
	 * has to be consistent at every recursive step.
	 * 
	 * @param params CeSymm parameters, make sure they are cloned
	 */
	public CeSymmIterative(CESymmParameters params) {
		//Disable the iteration mode, because otherwise we get infinite calls
		this.params = params;
		params.setMultipleAxes(false);

		msa = new MultipleAlignmentImpl();
		msa.getEnsemble().setAtomArrays(new ArrayList<Atom[]>());
		msa.getEnsemble().setAlgorithmName(CeSymm.algorithmName);
		msa.getEnsemble().setVersion(CeSymm.version);
		msa.getEnsemble().setStructureNames(new ArrayList<String>());

		BlockSet bs = new BlockSetImpl(msa);
		Block b = new BlockImpl(bs);
		b.setAlignRes(new ArrayList<List<Integer>>());

		alignment = new ArrayList<List<Integer>>();
		levels = new ArrayList<MultipleAlignment>();
		axes = new SymmetryAxes();
		name = null;
	}

	/**
	 * This method uses iteratively CeSymm to calculate all 
	 * symmetries in the input array of atoms and organize 
	 * them in a multiple alignment of the subunits.
	 * 
	 * @param atoms atoms 
	 * @return MultipleAlignment of the subunits
	 * 
	 * @throws StructureException  
	 */
	public MultipleAlignment execute(Atom[] atoms) 
			throws StructureException {

		allAtoms = atoms;
		for (Integer res=0; res<allAtoms.length; res++){
			alignment.add(new ArrayList<Integer>());
		}

		iterate(atoms);
		buildAlignment();
		recoverAxes();

		return msa;
	}

	/**
	 * This method runs iteratively CeSymm on the symmetric units
	 * until no more symmetries exist.
	 * 
	 * @param atoms Coordinates of the structure atoms
	 * @throws StructureException
	 */
	private void iterate(Atom[] atoms) throws StructureException {
		
		logger.debug("Starting new iteration...");
		
		if( atoms.length <= params.getWinSize() || 
				atoms.length <= params.getMinSubunitLength()) {
			logger.debug("Aborting iteration due to insufficient Atom "
					+ "array length: %d", atoms.length);
			return;
		}
		
		//Perform the CeSymm alignment
		CeSymm aligner = new CeSymm();
		MultipleAlignment align = aligner.analyze(atoms, params);
		if (name == null) 
			name = align.getEnsemble().getStructureNames().get(0);

		//End iterations if asymmetric
		if (!SymmetryTools.isRefined(align)) return;
		else if (align.getScore(MultipleAlignmentScorer.AVGTM_SCORE) < 
				params.getSymmetryThreshold() || 
				align.getCoreLength() < params.getMinSubunitLength()) {
			return;
		}

		levels.add(align);
		//If symmetric store the residue dependencies in graph
		Block b = align.getBlock(0);
		for (int pos=0; pos<b.length(); pos++){
			for (int su=0; su<b.size()-1; su++){
				Integer pos1 = b.getAlignRes().get(su).get(pos);
				Integer pos2 = b.getAlignRes().get(su+1).get(pos);
				//Add edge from lower to higher positions
				if (pos1 != null && pos2 != null){
					alignment.get(pos1).add(pos2);
				}
			}
		}

		//Generate the Atoms of one of the symmetric subunit
		Integer start = null;
		int it = 0;
		while (start == null){
			start = align.getBlocks().get(0).getAlignRes().get(0).get(it);
			it++;
		}
		Integer end = null;
		it = align.getBlocks().get(0).getAlignRes().get(0).size()-1;
		while (end == null){
			end = align.getBlocks().get(0).getAlignRes().get(0).get(it);
			it--;
		}
		//Iterate further
		Atom[] atomsR = Arrays.copyOfRange(atoms, start, end+1);
		
		iterate(atomsR);
	}

	private void buildAlignment() throws StructureException {

		List<List<Integer>> groups = new ArrayList<List<Integer>>();
		List<Integer> alreadySeen = new ArrayList<Integer>();
		int size = 0;

		//Calculate the connected groups of the alignment graph
		for (int i=0; i<alignment.size(); i++){
			if (!alreadySeen.contains(i)){
				List<Integer> group = new ArrayList<Integer>();
				List<Integer> residues = new ArrayList<Integer>();
				residues.add(i);

				while (residues.size() > 0){
					List<Integer> newResidues = new ArrayList<Integer>();
					for (Integer residue : residues){
						group.add(residue);
						alreadySeen.add(residue);
						List<Integer> children = 
								alignment.get(residue);
						newResidues.addAll(children);
					}
					residues = newResidues;
				}
				Collections.sort(group);
				groups.add(group);
				if (group.size() > size) size = group.size();
			}
		}

		Block b = msa.getBlock(0);
		//Construct the MultipleAlignment
		for (int su=0; su<size; su++) {
			msa.getEnsemble().getStructureNames().add(name);
			msa.getEnsemble().getAtomArrays().add(allAtoms);
			b.getAlignRes().add(new ArrayList<Integer>());

			for (List<Integer> group : groups){
				if (group.size() != size) continue;
				b.getAlignRes().get(su).add(group.get(su));
			}
		}
	}

	private void recoverAxes(){

		int size = msa.size();
		int parents = 1;

		for (int m=0; m<levels.size(); m++){

			MultipleAlignment align = levels.get(m);
			Matrix4d axis = align.getBlockSet(0).getTransformations().get(1);

			int subsize = align.size();
			parents *= subsize;
			size /= subsize;

			List<Integer> subunitTransform = new ArrayList<Integer>();
			for (int i=0; i<size*parents; i++){
				subunitTransform.add(0);
			}

			List<List<Integer>> superpose = new ArrayList<List<Integer>>();
			superpose.add(new ArrayList<Integer>());
			superpose.add(new ArrayList<Integer>());

			for (int su=0; su<subsize-1; su++){
				for (int s=0; s<size; s++){
					Integer subIndex1 = su*size+s;
					Integer subIndex2 = (su+1)*size+s;
					superpose.get(0).add(subIndex1);
					superpose.get(1).add(subIndex2);
				}
			}

			for (int p=0; p<parents; p++){
				for (int s=0; s<size; s++){
					subunitTransform.set(p*size+s, p%subsize);
				}
			}
			axes.addAxis(axis, superpose, subunitTransform, subsize);
		}
	}

	/**
	 * Return the symmetry axes.
	 * @return SymmetryAxes
	 */
	public SymmetryAxes getSymmetryAxes(){
		return axes;
	}
	
}
