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
import org.biojava.bio.structure.AtomImpl;
import org.biojava.bio.structure.Calc;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.io.mmcif.model.PdbxStructAssembly;
import org.biojava.bio.structure.io.mmcif.model.PdbxStructAssemblyGen;
import org.biojava.bio.structure.io.mmcif.model.PdbxStructOperList;
import org.biojava.bio.structure.jama.Matrix;

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

	private List<ModelTransformationMatrix> modelTransformations;

	public BiologicalAssemblyBuilder(){
		init();
	}

	public Structure rebuildQuaternaryStructure(Structure asymUnit, List<ModelTransformationMatrix> transformations){
		// ensure that new chains are build in the same order as they appear in the asymmetric unit
        orderTransformationsByChainId(asymUnit, transformations);
        
		Structure s = asymUnit.clone();
		List<Chain> transformedChains = new ArrayList<Chain>();
		s.setChains(transformedChains);

		for (ModelTransformationMatrix max : transformations){
			for (Chain c : asymUnit.getChains()){
                
				String intChainID = c.getInternalChainID();
				if ( intChainID == null) {
					//System.err.println("no internal chain ID found, using " + c.getChainID() + " ( while looking for " + max.ndbChainId+")");
					intChainID = c.getChainID();
				}

				if (max.getChainId().equals(intChainID)){
					Chain chain = (Chain)c.clone();
					Matrix m = max.getMatrix();
					double[] vector = max.getVector();
					Atom v = new AtomImpl();
					v.setCoords(vector);
					for (Group g : chain.getAtomGroups()) {
						Calc.rotate(g, m);
						Calc.shift(g, v);
					}

					int modelNumber = Integer.parseInt(max.getId());	
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
	private void orderTransformationsByChainId(Structure asymUnit, List<ModelTransformationMatrix> transformations) {
		final List<String> chainIds = getChainIds(asymUnit);
		Collections.sort(transformations, new Comparator<ModelTransformationMatrix>() {
			public int compare(ModelTransformationMatrix t1, ModelTransformationMatrix t2) {
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
	public ArrayList<ModelTransformationMatrix> getBioUnitTransformationList(PdbxStructAssembly psa, List<PdbxStructAssemblyGen> psags, List<PdbxStructOperList> operators) {
		//System.out.println("Rebuilding " + psa.getDetails() + " | " + psa.getOligomeric_details() + " | " + psa.getOligomeric_count());
		//System.out.println(psag);
		init();
		this.psags = psags;

		psa.getId();
		
		for (PdbxStructOperList oper: operators){
			ModelTransformationMatrix transform = new ModelTransformationMatrix();
			transform.setId(oper.getId());
			transform.setTransformationMatrix(oper.getMatrix(), oper.getVector());
			modelTransformations.add(transform);
		}

		ArrayList<ModelTransformationMatrix> transformations = getBioUnitTransformationsListUnaryOperators(psa.getId());
		transformations.addAll(getBioUnitTransformationsListBinaryOperators(psa.getId()));
		transformations.trimToSize();
		return transformations;
	}


	private ArrayList<ModelTransformationMatrix> getBioUnitTransformationsListBinaryOperators(String assemblyId) {

		ArrayList<ModelTransformationMatrix> transformations = new ArrayList<ModelTransformationMatrix>();

		List<OrderedPair<String>> operators = operatorResolver.getBinaryOperators();

		for ( PdbxStructAssemblyGen psag : psags){
			if ( psag.getAssembly_id().equals(assemblyId)) {

				List<String>asymIds= Arrays.asList(psag.getAsym_id_list().split(","));

				operatorResolver.parseOperatorExpressionString(psag.getOper_expression());

				// apply binary operators to the specified chains
				// Example 1M4X: generates all products of transformation matrices (1-60)(61-88)
				for (String chainId : asymIds) {

					for (OrderedPair<String> operator : operators) {
						ModelTransformationMatrix original1 = getModelTransformationMatrix(operator.getElement1());
						ModelTransformationMatrix original2 = getModelTransformationMatrix(operator.getElement2());
						ModelTransformationMatrix transform = ModelTransformationMatrix.multiply4square_x_4square2(original1, original2);
						transform.setChainId(chainId);
						transform.setId(original1.getId() + "x" + original2.getId());
						transformations.add(transform);
					}
				}
			}

		}

		return transformations;
	}

	private ModelTransformationMatrix getModelTransformationMatrix(String operator) {
		for (ModelTransformationMatrix transform: modelTransformations) {
			if (transform.getId().equals(operator)) {
				return transform;
			}
		}
		System.err.println("Could not find modelTransformationmatrix for " + operator);
		return new ModelTransformationMatrix();
	}

	private ArrayList<ModelTransformationMatrix> getBioUnitTransformationsListUnaryOperators(String assemblyId) {	
		ArrayList<ModelTransformationMatrix> transformations = new ArrayList<ModelTransformationMatrix>();

		
		for ( PdbxStructAssemblyGen psag : psags){
			if ( psag.getAssembly_id().equals(assemblyId)) {

				operatorResolver.parseOperatorExpressionString(psag.getOper_expression());
				List<String> operators = operatorResolver.getUnaryOperators();
				
				List<String>asymIds= Arrays.asList(psag.getAsym_id_list().split(","));
				
				// apply unary operators to the specified chains
				for (String chainId : asymIds) {
					for (String operator : operators) {
						//System.out.println("transforming " + chainId + " " + operator);
						ModelTransformationMatrix original = getModelTransformationMatrix(operator);
						ModelTransformationMatrix transform = new ModelTransformationMatrix(original);
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
		modelTransformations = new ArrayList<ModelTransformationMatrix>(1);
	}
}
