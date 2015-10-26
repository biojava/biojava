package org.biojava.nbio.structure.symmetry.internal;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.biojava.nbio.structure.align.util.AlignmentTools;
import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.symmetry.utils.SymmetryTools;

/**
 * Refines an alignment with open symmetry, that is, 
 * with tanslational and rotational components in the 
 * symmetry operation, so that the aligned residues 
 * are consistent. 
 * <p>
 * This refiner also determines the order of symmetry 
 * (number of repeats), by the heuristic: the order 
 * that contains maximum number of residues aligned.
 * <p>
 * The goal is to find the internal repeats without 
 * ths constraint that the first and last subunit align,
 * so no CP has to be present in the CeSymm alignment.
 * <p>
 * The symmetry axis indicates also the translation of 
 * the subunits in this case. It does not work for all
 * close symmetry cases, because it assumes that no CP 
 * ocurred in the first alignment.
 * 
 * @author Aleix Lafita
 * @since 4.2.0
 * 
 */
public class OpenRefiner implements Refiner {

	@Override
	public AFPChain refine(List<AFPChain> afpAlignments, Atom[] atoms,
			int order) throws RefinerFailedException, StructureException {

		//The two vertices of the graph mean (previous, next)
		List<List<Integer>> graph = 
				SymmetryTools.buildSymmetryGraph(afpAlignments, atoms, true);
		
		AFPChain afpChain = afpAlignments.get(afpAlignments.size()-1);
		List<Integer> alreadySeen = new ArrayList<Integer>();

		//Calculate the connected groups of the alignment graph
		List<List<Integer>> groups = new ArrayList<List<Integer>>();
		for (int i=0; i<graph.size(); i++){
			if (!alreadySeen.contains(i)){
				List<Integer> group = new ArrayList<Integer>();
				int residue = i;
				
				while (residue != -1 && !alreadySeen.contains(residue)){
					group.add(residue);
					List<Integer> neigh = graph.get(residue);
					//Go to the next residue in sequence
					if (neigh.size() > 1) {
						if (neigh.get(1) > residue) {
							residue = neigh.get(1);
						} else residue = -1;
					}
					else if (neigh.size() > 0) {
						if (neigh.get(0) > residue) {
							residue = neigh.get(0);
						} else residue = -1;
					}
					else residue = -1; //residue does not have a next
				}
				Collections.sort(group);
				if (group.size() > 1){
					groups.add(group);
					alreadySeen.addAll(group);
				}
			}
		}

		//Determine the order
		if (order == 0) {
			//Calculate the most common group size
			List<Integer> sizes = new ArrayList<Integer>(atoms.length);
			for (int p=0; p<atoms.length; p++) sizes.add(0);
	
			for (int i=0; i<groups.size(); i++){
				int gorder = groups.get(i).size();
				sizes.set(gorder, sizes.get(gorder)+1);
			}
			int maxNr = 0; //the total number of residues aligned
			for (int s=2; s<sizes.size(); s++){
				if (sizes.get(s)*s > maxNr) {
					order = s;
					maxNr = sizes.get(s)*s;
				}
			}
		}

		//Now create the new AFP alignment from the selected groups
		List<List<Integer>> subunits = new ArrayList<List<Integer>>();
		for (List<Integer> g:groups) if (g.size() == order) subunits.add(g);

		//Delete all inconsistent groups
		List<Integer> deleteIndices = new ArrayList<Integer>();
		for (int i=1; i<subunits.size(); i++){
			for (int j=0; j<subunits.get(i).size()-1; j++){
				if (subunits.get(i).get(j) > subunits.get(0).get(j+1)){
					deleteIndices.add(i);
					break;
				}
			}
		}
		for (int i=deleteIndices.size()-1; i>=0; i--){
			int index = deleteIndices.get(i);
			subunits.remove(index);
		}

		//From the groups of higher order take the consistent residues also
		for (List<Integer> g:groups){
			if (g.size() > order){
				List<Integer> group = new ArrayList<Integer>();
				for (int pos=0; pos<g.size() && group.size() < order; pos++){
					boolean consistent = true;
					for (List<Integer> sub:subunits){
						if (sub.get(group.size()) > g.get(pos)) 
							consistent = false;
						if (group.size() < order-1) {
							if (sub.get(group.size()+1) < g.get(pos)) 
								consistent = false;
						}
					}
					if (consistent && group.size()<order) 
						group.add(g.get(pos));
				}
				if (group.size()==order) subunits.add(group);
			}
		}
		
		if (subunits.size() == 0){
			throw new RefinerFailedException("Empty alignment");
		}

		int[][][] optAln = new int[order][2][subunits.size()];
		for (int bk=0; bk<order; bk++){
			optAln[bk] = new int[2][];
			optAln[bk][0] = new int[subunits.size()];
			optAln[bk][1] = new int[subunits.size()];
			for (int su=0; su<subunits.size(); su++){
				optAln[bk][0][su] = subunits.get(su).get(bk);
				optAln[bk][1][su] = subunits.get(su).get((bk+1)%order);
			}
		}
		afpChain = AlignmentTools.replaceOptAln(optAln, afpChain, atoms, atoms);
		return afpChain;
	}
	
}
