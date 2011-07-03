package org.biojava3.aaproperties.xml;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.xml.bind.annotation.XmlRootElement;

@XmlRootElement(name="AminoAcidsTable", namespace="http://biojava.org")
public class AminoAcidCompositionTable {
	
	/**
	 * Contains the list of amino acid composition 
	 */
	private List<AminoAcidComposition> aminoacid;
	
	/**
	 * Stores the mapping of amino acid symbol to its molecular weight
	 */
	private Map<String, Double> aaSymbol2MolecularWeight;
	
	public AminoAcidCompositionTable(){}
	
	public AminoAcidCompositionTable(List<AminoAcidComposition> aaList){
		this.setAminoacid(aaList);
	}

	public List<AminoAcidComposition> getAminoacid() {
		return aminoacid;
	}

	public void setAminoacid(List<AminoAcidComposition> aminoacid) {
		this.aminoacid = aminoacid;
	}
	
	/**
	 * Computes and store the molecular weight of each amino acid by its symbol in aaSymbol2MolecularWeight
	 * 
	 * @param eTable
	 * 	Stores the mass of elements and isotopes
	 */
	public void computeMolecularWeight(ElementTable eTable){
		this.aaSymbol2MolecularWeight = new HashMap<String, Double>();
		for(AminoAcidComposition a:aminoacid){
			double total = 0.0;
			for(String eName:a.getElement2Count().keySet()){
				Element e = eTable.getElement(eName);
				e.getMass();
				total += eTable.getElement(eName).getMass() * a.getElementCount(eName);
			}
			for(String iName:a.getIsotope2Count().keySet()){
				total += eTable.getIsotope(iName).weight() * a.getIsotopeCount(iName);
			}
			this.aaSymbol2MolecularWeight.put(a.getSymbol(), total);
		}
	}
	
	/**
	 * @param aaSymbol
	 * 	Standard symbol of Amino Acid
	 * @return the molecular weight given its symbol
	 * @throws NullPointerException 
	 * 	thrown if AminoAcidCompositionTable.computeMolecularWeight(ElementTable) is not called before this method
	 */
	public double getMolecularWeight(String aaSymbol) throws NullPointerException{
		if(this.aaSymbol2MolecularWeight == null){
			throw new NullPointerException("Please call AminoAcidCompositionTable.computeMolecularWeight(ElementTable) before this method");
		}
		Double d = this.aaSymbol2MolecularWeight.get(aaSymbol);
		if(d == null)
			return 0;
		else
			return d; 
	}
}
