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
 */

package org.biojava.bio.structure.quaternary;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.io.mmcif.model.PdbxStructAssembly;
import org.biojava.bio.structure.io.mmcif.model.PdbxStructAssemblyGen;
import org.biojava.bio.structure.io.mmcif.model.PdbxStructOperList;

/** 
 * Reconstructs the quaternary structure of a protein from an asymmetric unit
 * 
 * @author Peter Rose
 * @author Andreas Prlic
 *
 */
public class BiologicalAssemblyBuilder {
	private OperatorResolver operatorResolver;
	private List<PdbxStructAssemblyGen> psags;

	private List<BiologicalAssemblyTransformation> modelTransformations;

	public BiologicalAssemblyBuilder(){
		init();
	}

	public Structure rebuildQuaternaryStructure(Structure asymUnit, List<BiologicalAssemblyTransformation> transformations){
		// ensure that new chains are build in the same order as they appear in the asymmetric unit
        orderTransformationsByChainId(asymUnit, transformations);
        
		Structure s = asymUnit.clone();
		s.setChains(new ArrayList<Chain>());

		for (BiologicalAssemblyTransformation transformation : transformations){
			for (Chain c : asymUnit.getChains()){
                
				String intChainID = c.getInternalChainID();
				if (intChainID == null) {
					//System.err.println("no internal chain ID found, using " + c.getChainID() + " ( while looking for " + max.ndbChainId+")");
					intChainID = c.getChainID();
				}

				if (transformation.getChainId().equals(intChainID)){
					Chain chain = (Chain)c.clone();
		//			System.out.println("Applying transformation to: " + intChainID);
		//			System.out.println(m);
		//			System.out.println(Arrays.toString(v.getCoords()));
					for (Group g : chain.getAtomGroups()) {
	//					System.out.println("before");
						for (Atom a: g.getAtoms()) {
		//					System.out.println(a);
							transformation.transformPoint(a.getCoords());
//						System.out.println("after");
//						for (Atom a: g.getAtoms()) {
//							System.out.println(a);
						}
					}

					int modelNumber = Integer.parseInt(transformation.getId());	
	//				System.out.println("model: " + modelNumber);
					addChainAndModel(s, chain, modelNumber);
				}								
			}
		}

		s.setBiologicalAssembly(true);
		return s;
	}
	
	/**
	 * Orders model transformations by chain ids in the same order as in the asymmetric unit
	 * @param asymUnit
	 * @param transformations
	 */
	private void orderTransformationsByChainId(Structure asymUnit, List<BiologicalAssemblyTransformation> transformations) {
		final List<String> chainIds = getChainIds(asymUnit);
		Collections.sort(transformations, new Comparator<BiologicalAssemblyTransformation>() {
			public int compare(BiologicalAssemblyTransformation t1, BiologicalAssemblyTransformation t2) {
				// set sort order only if the two ids are identical
				if (t1.getId().equals(t2.getId())) {
					 return chainIds.indexOf(t1.getChainId()) - chainIds.indexOf(t2.getChainId());
				}
			    return 0;
		    }
		});
	}
	
	/**
	 * Returns a list of chain ids in the order they are specified in the ATOM
	 * records in the asymmetric unit
	 * @param asymUnit
	 * @return
	 */
	private List<String> getChainIds(Structure asymUnit) {
		List<String> chainIds = new ArrayList<String>();
		for ( Chain c : asymUnit.getChains()){      
			String intChainID = c.getInternalChainID();
			if ( intChainID == null) {
				//System.err.println("no internal chain ID found, using " + c.getChainID() + " ( while looking for " + max.ndbChainId+")");
				intChainID = c.getChainID();
			}
			chainIds.add(intChainID);
		}
		return chainIds;
	}
	
	private void addChainAndModel(Structure s, Chain newChain, int modelCount) {
		if (modelCount == 0) {
			s.addChain(newChain);
		} else if (modelCount > s.nrModels()) {
			List<Chain> newModel = new ArrayList<Chain>();
			newModel.add(newChain);
			s.addModel(newModel);
		} else {
			s.addChain(newChain, modelCount-1);
		}
	}

	/**
	 * Returns a list of transformation matrices for the generation of a macromolecular
	 * assembly for the specified assembly Id. 
	 * 
	 * @param assemblyId Id of the macromolecular assembly to be generated
	 * @return list of transformation matrices to generate macromolecular assembly
	 */
	public ArrayList<BiologicalAssemblyTransformation> getBioUnitTransformationList(PdbxStructAssembly psa, List<PdbxStructAssemblyGen> psags, List<PdbxStructOperList> operators) {
		//System.out.println("Rebuilding " + psa.getDetails() + " | " + psa.getOligomeric_details() + " | " + psa.getOligomeric_count());
		//System.out.println(psag);
		init();
		this.psags = psags;

		//psa.getId();
		
		for (PdbxStructOperList oper: operators){
			BiologicalAssemblyTransformation transform = new BiologicalAssemblyTransformation();
			transform.setId(oper.getId());
			transform.setRotationMatrix(oper.getMatrix());
			transform.setTranslation(oper.getVector());
//			transform.setTransformationMatrix(oper.getMatrix(), oper.getVector());
			modelTransformations.add(transform);
		}

		ArrayList<BiologicalAssemblyTransformation> transformations = getBioUnitTransformationsListUnaryOperators(psa.getId());
		transformations.addAll(getBioUnitTransformationsListBinaryOperators(psa.getId()));
		transformations.trimToSize();
		return transformations;
	}


	private ArrayList<BiologicalAssemblyTransformation> getBioUnitTransformationsListBinaryOperators(String assemblyId) {

		ArrayList<BiologicalAssemblyTransformation> transformations = new ArrayList<BiologicalAssemblyTransformation>();

		List<OrderedPair<String>> operators = operatorResolver.getBinaryOperators();
	

		for ( PdbxStructAssemblyGen psag : psags){
			if ( psag.getAssembly_id().equals(assemblyId)) {

				List<String>asymIds= Arrays.asList(psag.getAsym_id_list().split(","));

				operatorResolver.parseOperatorExpressionString(psag.getOper_expression());

				// apply binary operators to the specified chains
				// Example 1M4X: generates all products of transformation matrices (1-60)(61-88)
				for (String chainId : asymIds) {

					int modelNumber = 1;
					for (OrderedPair<String> operator : operators) {
						BiologicalAssemblyTransformation original1 = getModelTransformationMatrix(operator.getElement1());
						BiologicalAssemblyTransformation original2 = getModelTransformationMatrix(operator.getElement2());
			//			ModelTransformationMatrix transform = ModelTransformationMatrix.multiply4square_x_4square2(original1, original2);
						BiologicalAssemblyTransformation transform = BiologicalAssemblyTransformation.combine(original1, original2);
						transform.setChainId(chainId);
				//		transform.setId(original1.getId() + "x" + original2.getId());
						transform.setId(String.valueOf(modelNumber));
						transformations.add(transform);
						modelNumber++;
					}
				}
			}

		}

		return transformations;
	}

	private BiologicalAssemblyTransformation getModelTransformationMatrix(String operator) {
		for (BiologicalAssemblyTransformation transform: modelTransformations) {
			if (transform.getId().equals(operator)) {
				return transform;
			}
		}
		System.err.println("Could not find modelTransformationmatrix for " + operator);
		return new BiologicalAssemblyTransformation();
	}

	private ArrayList<BiologicalAssemblyTransformation> getBioUnitTransformationsListUnaryOperators(String assemblyId) {	
		ArrayList<BiologicalAssemblyTransformation> transformations = new ArrayList<BiologicalAssemblyTransformation>();

		
		for ( PdbxStructAssemblyGen psag : psags){
			if ( psag.getAssembly_id().equals(assemblyId)) {
		
				operatorResolver.parseOperatorExpressionString(psag.getOper_expression());
				List<String> operators = operatorResolver.getUnaryOperators();
				
				List<String>asymIds= Arrays.asList(psag.getAsym_id_list().split(","));
				
				// apply unary operators to the specified chains
				for (String chainId : asymIds) {
					for (String operator : operators) {
				
						BiologicalAssemblyTransformation original = getModelTransformationMatrix(operator);
						BiologicalAssemblyTransformation transform = new BiologicalAssemblyTransformation(original);
						transform.setChainId(chainId);
						transform.setId(operator);
						transformations.add(transform);
					}
				}
			}
		}

		return transformations;
	}
	
	private void init(){
		operatorResolver= new OperatorResolver();
		modelTransformations = new ArrayList<BiologicalAssemblyTransformation>(1);
	}
}
