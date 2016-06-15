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
package org.biojava.nbio.aaproperties.xml;

import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;


@XmlRootElement(name="compoundtable", namespace="http://biojava.org")
@XmlAccessorType(XmlAccessType.NONE)
public class AminoAcidCompositionTable {

	/**
	 * Contains the list of amino acid composition
	 */
	@XmlElement(name = "compound", required = true)
	private List<AminoAcidComposition> aminoacid;

	/**
	 * Defines the amino acid compound set unique to this table
	 */
	private ModifiedAminoAcidCompoundSet modifiedAminoAcidCompoundSet;

	/**
	 * Stores the mapping of amino acid symbol to its molecular weight
	 */
	private Map<Character, Double> aaSymbol2MolecularWeight;

	public AminoAcidCompositionTable(){}

	public AminoAcidCompositionTable(List<AminoAcidComposition> aaList){
		this.setAminoacid(aaList);
	}

	public ModifiedAminoAcidCompoundSet getAminoAcidCompoundSet(){
		return this.modifiedAminoAcidCompoundSet;
	}

	public List<AminoAcidComposition> getAminoacid() {
		return aminoacid;
	}

	public void setAminoacid(List<AminoAcidComposition> aminoacid) {
		this.aminoacid = aminoacid;
	}

	public Set<Character> getSymbolSet(){
		return this.aaSymbol2MolecularWeight.keySet();
	}

	private void generatesAminoAcidCompoundSet(){
		this.modifiedAminoAcidCompoundSet = new ModifiedAminoAcidCompoundSet(this.aminoacid, this.aaSymbol2MolecularWeight);
	}

	/**
	 * Computes and store the molecular weight of each amino acid by its symbol in aaSymbol2MolecularWeight.
	 *
	 * @param eTable
	 * 		Stores the mass of elements and isotopes
	 */
	public void computeMolecularWeight(ElementTable eTable){
		this.aaSymbol2MolecularWeight = new HashMap<Character, Double>();
		for(AminoAcidComposition a:aminoacid){
			//Check to ensure that the symbol is of single character
			if(a.getSymbol().length() != 1){
				throw new Error(a.getSymbol() + " is not allowed. Symbols must be single character.\r\nPlease check AminoAcidComposition XML file");
			}
			//Check to ensure that the symbols are not repeated
			char c = a.getSymbol().charAt(0);
			if(this.aaSymbol2MolecularWeight.keySet().contains(c)){
				throw new Error("Symbol " + c + " is repeated.\r\n" +
						"Please check AminoAcidComposition XML file to ensure there are no repeated symbols. Note that this is case-insensitive.\r\n" +
						"This means that having 'A' and 'a' would be repeating.");
			}
			double total = 0.0;
			if(a.getElementList() != null){
				for(Name2Count element:a.getElementList()){
					element.getName();
					if(eTable.getElement(element.getName()) == null){
						throw new Error("Element " + element.getName() + " could not be found. " +
								"\r\nPlease ensure that its name is correct in AminoAcidComposition.xml and is defined in ElementMass.xml.");
					}
					eTable.getElement(element.getName()).getMass();
					total += eTable.getElement(element.getName()).getMass() * element.getCount();
				}
			}
			if(a.getIsotopeList() != null){
				for(Name2Count isotope:a.getIsotopeList()){
					isotope.getName();
					if(eTable.getIsotope(isotope.getName()) == null){
						throw new Error("Isotope " + isotope.getName() + " could not be found. " +
							"\r\nPlease ensure that its name is correct in AminoAcidComposition.xml and is defined in ElementMass.xml.");
					}
					eTable.getIsotope(isotope.getName()).getMass();
					total += eTable.getIsotope(isotope.getName()).getMass() * isotope.getCount();
				}
			}
			c = a.getSymbol().charAt(0);
			this.aaSymbol2MolecularWeight.put(c, total);
		}
		generatesAminoAcidCompoundSet();
	}

	/**
	 * @param aaSymbol
	 * 	Standard symbol of Amino Acid
	 * @return the molecular weight given its symbol
	 * @throws NullPointerException
	 * 	thrown if AminoAcidCompositionTable.computeMolecularWeight(ElementTable) is not called before this method
	 */
	public double getMolecularWeight(Character aaSymbol) throws NullPointerException{
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
