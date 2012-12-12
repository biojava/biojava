package org.biojava.bio.structure.quaternary;

import java.util.ArrayList;
import java.util.Arrays;
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



/** Reconstructs the quaternary structure of a protein from an asymmetric unit
 * 
 * @author Peter Rose
 * @author Andreas Prlic
 *
 */
public class BiologicalAssemblyBuilder {

	OperatorResolver operatorResolver;

	String asymId ;

	PdbxStructAssembly psa;
	List<PdbxStructAssemblyGen> psags;
	//List<PdbxStructOperList> operators;

	List<ModelTransformationMatrix> modelTransformations;

	//List<String> asymIds;
	public BiologicalAssemblyBuilder(){
		init();
	}

	private void init(){
		operatorResolver= new OperatorResolver();
		modelTransformations = new ArrayList<ModelTransformationMatrix>(1);
		
	}




	public Structure rebuildQuaternaryStructure(Structure asymUnit, List<ModelTransformationMatrix> transformations){

		Structure s = asymUnit.clone();
		List<Chain> transformedChains = new ArrayList<Chain>();
		s.setChains(transformedChains);
		//System.out.print("rebuilding " + s.getPDBCode() + " ");
		//for (Chain c : s.getChains()) {
			//System.out.print(c.getChainID());
			//if ( c.getInternalChainID() != null) {
				//System.out.print("("+c.getInternalChainID()+")");
			//}

		//}

		//asymUnit.getPDBHeader().setBioUnitTranformations(transformations);


		//double[] tmpcoords = new double[3]; 
		for (ModelTransformationMatrix max : transformations){
			boolean foundChain = false;
			for ( Chain c : asymUnit.getChains()){

				String intChainID = c.getInternalChainID();
				if ( intChainID == null) {
					//System.err.println("no internal chain ID found, using " + c.getChainID() + " ( while looking for " + max.ndbChainId+")");
					intChainID = c.getChainID();
				}

				if ( max.ndbChainId.equals(intChainID)){

					/*System.out.println("transforming " + max.ndbChainId + " " + c.getChainID());
					System.out.println(max );
					System.out.println();
					 */
					foundChain = true;
					Chain newChain = (Chain)c.clone();
					Matrix m = max.getMatrix();
					//m = m.transpose();
					double[] vector = max.getVector();
					Atom v = new AtomImpl();
					v.setCoords(vector);
					for ( Group g :newChain.getAtomGroups()) {
						for ( Atom a : g.getAtoms()) {
							//max.transformPoint(a.getCoords(), tmpcoords);
							//a.setX(tmpcoords[0]);
							//a.setY(tmpcoords[1]);
							//a.setZ(tmpcoords[2]);
							Calc.rotate(a, m);
							Calc.shift(a, v);
						}
					}

					
					addCheckChainModel(s,newChain);

				}								
			}
			if (! foundChain){
				//System.err.println("could not transform chain: " + max.ndbChainId);
			}
		}

		//s.setChains(transformedChains);
		s.setBiologicalAssembly(true);
		return s;


	}


	private void addCheckChainModel(Structure s, Chain newChain) {
		for ( int i = 0 ; i < s.nrModels() ; i ++){
			List<Chain> model = s.getModel(i);
			boolean found = false;
			for ( Chain c : model){
				if ( c.getChainID().equals(newChain.getChainID())) {
					found = true;
					break;
				}
			}
			if ( ! found){
			
				model.add(newChain);
				return;
			}
			
		}
		
		// all existing models contain already a chain with the same ID. add a new model
		List<Chain> newModel = new ArrayList<Chain>();
		newModel.add(newChain);
		s.addModel(newModel);
		
		
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
		this.psa=psa;
		this.psags = psags;

		//this.operators = operators;

		asymId = psa.getId();
		
		for (PdbxStructOperList oper: operators){
			ModelTransformationMatrix transform = new ModelTransformationMatrix();
			transform.id = oper.getId();
			transform.setTransformationMatrix(oper.getMatrix(), oper.getVector());
			modelTransformations.add(transform);
		}

		///

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
						transform.ndbChainId = chainId;
						transform.id = original1.id + "x" + original2.id;
						transformations.add(transform);
					}
				}

			}

		}

		// apply binary operators to the specified chains
		// Example 1M4X: generates all products of transformation matrices (1-60)(61-88)

		//System.out.println("BioUnitTrasnformListBinaryOperators " + assemblyId + " " + transformations.size());
		return transformations;
	}

	private ModelTransformationMatrix getModelTransformationMatrix(String operator) {
		for (ModelTransformationMatrix transform: modelTransformations) {
			if (transform.id.equals(operator)) {
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
						transform.ndbChainId = chainId;
						transform.id = operator;
						transformations.add(transform);
					}
				}
			}
		}
		//System.out.println("BioUnitTrasnformListUnary " + assemblyId + " " + transformations.size());
		return transformations;
	}
}
