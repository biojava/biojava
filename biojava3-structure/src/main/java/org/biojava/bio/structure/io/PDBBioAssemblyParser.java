package org.biojava.bio.structure.io;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.biojava.bio.structure.jama.Matrix;
import org.biojava.bio.structure.quaternary.ModelTransformationMatrix;

public class PDBBioAssemblyParser {

	
	Integer currentBioMolecule = null;

	List<String> currentChainIDs;
	
	Matrix currentMatrix ;
	
	double[] shift;
	
	Map<Integer,List<ModelTransformationMatrix>> transformationMap = new HashMap<Integer, List<ModelTransformationMatrix>>();
	
	int currentIndex;
	
	List<ModelTransformationMatrix> transformations;
	int correct ;
	
	private boolean DEBUG = false;
	private boolean justCommitted = false;
	
	public PDBBioAssemblyParser() {
		
		currentChainIDs = new ArrayList<String>();
		
		currentIndex = 1;
		
		transformations = new ArrayList<ModelTransformationMatrix>();
		
		currentMatrix = Matrix.identity(3,3);
		currentBioMolecule = null;
		shift = new double[3];
		correct = 0;
		
	}
	
	/** See here for the current spec:
	 * 
	 * http://www.wwpdb.org/documentation/format33/remarks2.html
	 * 
	 * @param line
	 */
	public void pdb_REMARK_350_Handler(String line) {
				
		if (line.startsWith("REMARK 350 BIOMOLECULE:")) {
		
			String nr = line.substring(24,line.length()).trim();
			
			if ( currentBioMolecule != null){
				finalizeCurrentBioMolecule();
			}
			//System.out.println(line);
			currentBioMolecule = Integer.parseInt(nr);
			currentIndex = 1;
			currentChainIDs.clear();
	
		} else if ( line.startsWith("REMARK 350 APPLY THE FOLLOWING TO CHAINS:")) {
			//System.out.println("NEW SECTION " + currentChainIDs.size()+ " " + line);
			if ( currentChainIDs.size() > 0){
				addNewMatrix();
				currentChainIDs.clear();
				currentIndex++;
			}
			
			if ( currentChainIDs.size() > 0) {
				addNewMatrix();
				
				currentMatrix = Matrix.identity(3,3);
				
				shift = new double[3];
				currentChainIDs.clear();
				justCommitted = true;
			}
			addToCurrentChainList(line);
		
		} else if ( line.startsWith("REMARK 350                    AND CHAINS:")) {
		
			addToCurrentChainList(line);
		
		} else if ( line.startsWith("REMARK 350   BIOMT")) {
		
			readMatrix(line);
	
		}
		
	}

	private void readMatrix(String line) {
		
		//System.out.println(line);
		String pos = line.substring(18,19);
		int i = Integer.parseInt(pos);
		
		
		String index = line.substring(21,24+correct).trim();
		
		int id = Integer.parseInt(index);
		
		if ( id != currentIndex && (! justCommitted)) {
			if ( id > currentIndex + 1) {
				//  2JBP
				// todo: still need to fix:  2J4Z bioassembly 2  somehow...
				System.err.println("WARNING ID " + id + " > " + currentIndex + " " + currentChainIDs);
			} else {
				addNewMatrix();
			}
			currentIndex = id;
			
			currentMatrix = Matrix.identity(3,3);
			
			shift = new double[3];
			
		}
		
		
		String x = line.substring(24+correct,33+correct);
		
		
		String y = line.substring(34+correct,43+correct);
		
		
		String z = line.substring(44+correct,53+correct);
		
		
		String vec = line.substring(58+correct,line.length()).trim();
	
		if (DEBUG)
			System.out.println(id + " |" + i + "| |" + x + "| |" + y + "| |" + z + "| >" + vec+"<");
		currentMatrix.set(0,(i-1),Float.parseFloat(x));
		currentMatrix.set(1,(i-1),Float.parseFloat(y));
		currentMatrix.set(2,(i-1),Float.parseFloat(z));
		
		shift[i-1] = Float.parseFloat(vec);
		//System.out.println(shift[i-1]);
		//System.out.println(currentMatrix);
		justCommitted = false;
		
	}

	private void addNewMatrix() {
		//System.out.println("adding new matrix " + currentIndex + " for " + currentChainIDs);
		
		
		ModelTransformationMatrix max = new ModelTransformationMatrix();
		
		// this can happen and be valid e.g 1hv4
//		if ( currentBioMolecule > currentIndex) {
//			System.err.println("WARNING current molecule index " + currentBioMolecule +" > global index " + currentIndex + "! current chains: " + Arrays.toString(currentChainIDs.toArray()) );
//		}
		if ( DEBUG) {
			System.out.println("AddNewMatrix bio ass: " + currentBioMolecule + " index: " + currentIndex + " " +  currentChainIDs);
			System.out.println(currentMatrix);
			System.out.println(Arrays.toString(shift));
		}
		max.setMatrix(currentMatrix);
		max.setVector(shift);
		max.id = currentIndex+"";
		
			
		for ( String chainId : currentChainIDs) {
			ModelTransformationMatrix m = (ModelTransformationMatrix) max.clone();
			m.setNdbChainId(chainId);
			transformations.add(m);
			
			//System.out.println(m);
		}
		
		currentMatrix = Matrix.identity(3,3);
		shift = new double[3];
		
		if ( currentIndex == 999){
			//System.err.println("CORRECTION!");
			correct +=1;
		}
		
		justCommitted  = true;
		
		
		
	}

	private void addToCurrentChainList(String line) {
		
		
		
		String chainIds = line.substring(41,line.length()).trim().replaceAll(",","");
		
		String[] spl = chainIds.split(" ");
		for ( String chainId : spl){
			currentChainIDs.add(chainId);
		}
		
		
	}

	public void finalizeCurrentBioMolecule() {
		if ( DEBUG )
			System.out.println("finalizing biomolecule..." + currentBioMolecule + " (current index: " + currentIndex +")");
		// the last matrix has not been added at this stage...
		addNewMatrix();
		
		//System.out.println("transformations for current molec: " + transformations);	
		transformationMap.put(currentBioMolecule,transformations);
		
		transformations = new ArrayList<ModelTransformationMatrix>();
		//System.out.println("BioMolecule " + currentBioMolecule + " has chainIDs: " + currentChainIDs);
			
		currentChainIDs.clear();
	}

	public Map<Integer,List<ModelTransformationMatrix>> getTransformationMap() {
		//System.out.println(transformationMap);
		return transformationMap;
	}

	public void setTransformationMap(Map<Integer,List<ModelTransformationMatrix>> transformationMap) {
		this.transformationMap = transformationMap;
	}
	
	public int getNrBioAssemblies(){
		
		if ( currentBioMolecule == null)
			return 0;
		return currentBioMolecule;
	}

	
}
