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

package org.biojava.nbio.structure.quaternary;

import org.biojava.nbio.structure.Calc;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.io.mmcif.model.PdbxStructAssembly;
import org.biojava.nbio.structure.io.mmcif.model.PdbxStructAssemblyGen;
import org.biojava.nbio.structure.io.mmcif.model.PdbxStructOperList;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import javax.vecmath.Matrix4d;
import java.util.*;

/**
 * Reconstructs the quaternary structure of a protein from an asymmetric unit
 *
 * @author Peter Rose
 * @author Andreas Prlic
 * @author Jose Duarte
 *
 */
public class BiologicalAssemblyBuilder {

	private static final Logger logger = LoggerFactory.getLogger(BiologicalAssemblyBuilder.class);

	private OperatorResolver operatorResolver;


	/**
	 * All matrix operators present in _pdbx_struct_oper_list.
	 * Identifiers (_pdbx_struct_oper_list.id) to matrix operators.
	 */
	private Map<String, Matrix4d> allTransformations;

	private List<String> modelIndex = new ArrayList<String>();

	public BiologicalAssemblyBuilder(){
		init();
	}

	/**
	 * Builds a Structure object containing the quaternary structure built from given asymUnit and transformations,
	 * by adding symmetry partners as new models.
	 * The output Structure will be different depending on the multiModel parameter:
	 * <li>
	 * the symmetry-expanded chains are added as new models, one per transformId. All original models but 
	 * the first one are discarded.
	 * </li>
	 * <li>
	 * as original with symmetry-expanded chains added with renamed chain ids and names (in the form 
	 * originalAsymId_transformId and originalAuthId_transformId)
	 * </li>
	 * @param asymUnit
	 * @param transformations
	 * @param useAsymIds if true use {@link Chain#getId()} to match the ids in the BiologicalAssemblyTransformation (needed if data read from mmCIF), 
	 * if false use {@link Chain#getName()} for the chain matching (needed if data read from PDB).
	 * @param multiModel if true the output Structure will be a multi-model one with one transformId per model, 
	 * if false the outputStructure will be as the original with added chains with renamed asymIds (in the form originalAsymId_transformId and originalAuthId_transformId). 
	 * @return
	 */
	public Structure rebuildQuaternaryStructure(Structure asymUnit, List<BiologicalAssemblyTransformation> transformations, boolean useAsymIds, boolean multiModel) {
		
		// ensure that new chains are build in the same order as they appear in the asymmetric unit
		orderTransformationsByChainId(asymUnit, transformations);

		Structure s = asymUnit.clone();
		

		// this resets all models (not only the first one): this is important for NMR (multi-model)
		// like that we can be sure we start with an empty structures and we add models or chains to it
		s.resetModels();

		for (BiologicalAssemblyTransformation transformation : transformations){

			List<Chain> chainsToTransform = new ArrayList<>();
			
			// note: for NMR structures (or any multi-model) we use the first model only and throw away the rest
			if (useAsymIds) {
				Chain c = asymUnit.getChain(transformation.getChainId());
				chainsToTransform.add(c);
			} else {
				Chain polyC = asymUnit.getPolyChainByPDB(transformation.getChainId());
				List<Chain> nonPolyCs = asymUnit.getNonPolyChainsByPDB(transformation.getChainId());
				Chain waterC = asymUnit.getWaterChainByPDB(transformation.getChainId());
				if (polyC!=null) 
					chainsToTransform.add(polyC);
				if (!nonPolyCs.isEmpty()) 
					chainsToTransform.addAll(nonPolyCs);
				if (waterC!=null) 
					chainsToTransform.add(waterC);
			}
			
			for (Chain c: chainsToTransform) {

				Chain chain = (Chain)c.clone();
				
				Calc.transform(chain, transformation.getTransformationMatrix());

				String transformId = transformation.getId();

				// note that the Structure.addChain/Structure.addModel methods set the parent reference to the new Structure
				
				// TODO set entities properly in the new structures! at the moment they are a mess... - JD 2016-05-19
				
				if (multiModel) 
					addChainMultiModel(s, chain, transformId);
				else 
					addChainFlattened(s, chain, transformId);

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
			@Override
			public int compare(BiologicalAssemblyTransformation t1, BiologicalAssemblyTransformation t2) {
				// set sort order only if the two ids are identical
				if (t1.getId().equals(t2.getId())) {
					 return chainIds.indexOf(t1.getChainId()) - chainIds.indexOf(t2.getChainId());
				} else {
					return t1.getId().compareTo(t2.getId());
				}
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
			String intChainID = c.getId();
			chainIds.add(intChainID);
		}
		return chainIds;
	}

	/**
	 * Adds a chain to the given structure to form a biological assembly,
	 * adding the symmetry expanded chains as new models per transformId.
	 * @param s
	 * @param newChain
	 * @param transformId
	 */
	private void addChainMultiModel(Structure s, Chain newChain, String transformId) {

		// multi-model bioassembly

		if ( modelIndex.size() == 0)
			modelIndex.add("PLACEHOLDER FOR ASYM UNIT");

		int modelCount = modelIndex.indexOf(transformId);
		if ( modelCount == -1)  {
			modelIndex.add(transformId);
			modelCount = modelIndex.indexOf(transformId);
		}

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
	 * Adds a chain to the given structure to form a biological assembly,
	 * adding the symmetry-expanded chains as new chains with renamed 
	 * chain ids and names (in the form originalAsymId_transformId and originalAuthId_transformId).
	 * @param s
	 * @param newChain
	 * @param transformId
	 */
	private void addChainFlattened(Structure s, Chain newChain, String transformId) {
		newChain.setId(newChain.getId()+"_"+transformId);
		newChain.setName(newChain.getName()+"_"+transformId);
		s.addChain(newChain);		
	}

	/**
	 * Returns a list of transformation matrices for the generation of a macromolecular
	 * assembly for the specified assembly Id.
	 *
	 * @param psa
	 * @param psags
	 * @param operators
	 * @return list of transformation matrices to generate macromolecular assembly
	 */
	public ArrayList<BiologicalAssemblyTransformation> getBioUnitTransformationList(PdbxStructAssembly psa, List<PdbxStructAssemblyGen> psags, List<PdbxStructOperList> operators) {
		init();

		// first we populate the list of all operators from pdbx_struct_oper_list so that we can then
		// get them from getBioUnitTransformationsListUnaryOperators() and getBioUnitTransformationsListBinaryOperators()
		for (PdbxStructOperList oper: operators){
			try {
				Matrix4d m = new Matrix4d();
				m.m00 = Double.parseDouble(oper.getMatrix11());
				m.m01 = Double.parseDouble(oper.getMatrix12());
				m.m02 = Double.parseDouble(oper.getMatrix13());

				m.m10 = Double.parseDouble(oper.getMatrix21());
				m.m11 = Double.parseDouble(oper.getMatrix22());
				m.m12 = Double.parseDouble(oper.getMatrix23());

				m.m20 = Double.parseDouble(oper.getMatrix31());
				m.m21 = Double.parseDouble(oper.getMatrix32());
				m.m22 = Double.parseDouble(oper.getMatrix33());

				m.m03 = Double.parseDouble(oper.getVector1());
				m.m13 = Double.parseDouble(oper.getVector2());
				m.m23 = Double.parseDouble(oper.getVector3());

				m.m30 = 0;
				m.m31 = 0;
				m.m32 = 0;
				m.m33 = 1;

				allTransformations.put(oper.getId(), m);

			} catch (NumberFormatException e) {
				logger.warn("Could not parse a matrix value from pdbx_struct_oper_list for id {}. The operator id will be ignored. Error: {}", oper.getId(), e.getMessage());
			}
		}

		ArrayList<BiologicalAssemblyTransformation> transformations = getBioUnitTransformationsListUnaryOperators(psa.getId(), psags);
		transformations.addAll(getBioUnitTransformationsListBinaryOperators(psa.getId(), psags));
		transformations.trimToSize();
		return transformations;
	}


	private ArrayList<BiologicalAssemblyTransformation> getBioUnitTransformationsListBinaryOperators(String assemblyId, List<PdbxStructAssemblyGen> psags) {

		ArrayList<BiologicalAssemblyTransformation> transformations = new ArrayList<BiologicalAssemblyTransformation>();

		List<OrderedPair<String>> operators = operatorResolver.getBinaryOperators();


		for ( PdbxStructAssemblyGen psag : psags){
			if ( psag.getAssembly_id().equals(assemblyId)) {

				List<String>asymIds= Arrays.asList(psag.getAsym_id_list().split(","));

				operatorResolver.parseOperatorExpressionString(psag.getOper_expression());

				// apply binary operators to the specified chains
				// Example 1M4X: generates all products of transformation matrices (1-60)(61-88)
				for (String chainId : asymIds) {

					for (OrderedPair<String> operator : operators) {
						Matrix4d original1 = allTransformations.get(operator.getElement1());
						Matrix4d original2 = allTransformations.get(operator.getElement2());
						if (original1 == null || original2 == null) {
							logger.warn("Could not find matrix operator for operator id {} or {}. Assembly id {} will not contain the composed operator.", operator.getElement1(), operator.getElement2(), assemblyId);
							continue;
						}
						Matrix4d composed = new Matrix4d(original1);
						composed.mul(original2);
						BiologicalAssemblyTransformation transform = new BiologicalAssemblyTransformation();
						transform.setChainId(chainId);
						transform.setId(operator.getElement1() + "x" + operator.getElement2());
						transform.setTransformationMatrix(composed);
						transformations.add(transform);
					}
				}
			}

		}

		return transformations;
	}

	private ArrayList<BiologicalAssemblyTransformation> getBioUnitTransformationsListUnaryOperators(String assemblyId, List<PdbxStructAssemblyGen> psags) {
		ArrayList<BiologicalAssemblyTransformation> transformations = new ArrayList<BiologicalAssemblyTransformation>();


		for ( PdbxStructAssemblyGen psag : psags){
			if ( psag.getAssembly_id().equals(assemblyId)) {

				operatorResolver.parseOperatorExpressionString(psag.getOper_expression());
				List<String> operators = operatorResolver.getUnaryOperators();

				List<String>asymIds= Arrays.asList(psag.getAsym_id_list().split(","));

				// apply unary operators to the specified chains
				for (String chainId : asymIds) {
					for (String operator : operators) {
						Matrix4d original = allTransformations.get(operator);
						if (original == null) {
							logger.warn("Could not find matrix operator for operator id {}. Assembly id {} will not contain the operator.", operator, assemblyId);
							continue;
						}
						BiologicalAssemblyTransformation transform = new BiologicalAssemblyTransformation();
						transform.setChainId(chainId);
						transform.setId(operator);
						transform.setTransformationMatrix(original);
						transformations.add(transform);
					}
				}
			}
		}

		return transformations;
	}

	private void init(){
		operatorResolver= new OperatorResolver();
		allTransformations = new HashMap<>();
	}
}
