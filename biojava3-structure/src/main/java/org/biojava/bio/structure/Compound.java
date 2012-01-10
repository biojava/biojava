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
 * Created on 22.01.2007
 *
 */


package org.biojava.bio.structure;


import java.io.Serializable;
import java.lang.reflect.Method;
import java.util.ArrayList;
import java.util.List;

/**
 * An object to contain the info from the PDB header for a Molecule.
 *
 * Now PDB file format 3.2 aware - contains the new TAX_ID fields for the
 * organism studied and the expression system.
 *
 * @author Jules Jacobsen
 * @since 1.5
 */
/**
 * Created by IntelliJ IDEA.
 * User: jacobsen
 * Date: 22-Jan-2007
 * Time: 18:14:25
 */
public class Compound implements Cloneable, Serializable {

	/**
    *
    */
   private static final long serialVersionUID = 2991897825657586356L;
   private List<Chain> chainList = new ArrayList<Chain>();
	private List<String> chainId = null;
	private String refChainId = null;
	private String molId = "0";
	//String molId = null;
	private String molName = null;
	private String title = null;
	private List<String> synonyms = null;
	private List<String> ecNums = null;
	private String engineered = null;
	private String mutation = null;
	private String biologicalUnit = null;
	private String details = null;

	private String numRes = null;
	private String resNames = null;

	private String headerVars = null;

	private String synthetic = null;
	private String fragment = null;
	private String organismScientific = null;
    private String organismTaxId = null;
	private String organismCommon = null;
	private String strain = null;
	private String variant = null;
	private String cellLine = null;
	private String atcc = null;
	private String organ = null;
	private String tissue = null;
	private String cell = null;
	private String organelle = null;
	private String secretion = null;
	private String gene = null;
	private String cellularLocation = null;
	private String expressionSystem = null;
    private String expressionSystemTaxId = null;
	private String expressionSystemStrain = null;
	private String expressionSystemVariant = null;
	private String expressionSystemCellLine = null;
	private String expressionSystemAtccNumber = null;
	private String expressionSystemOrgan = null;
	private String expressionSystemTissue = null;
	private String expressionSystemCell = null;
	private String expressionSystemOrganelle = null;
	private String expressionSystemCellularLocation = null;
	private String expressionSystemVectorType = null;
	private String expressionSystemVector = null;
	private String expressionSystemPlasmid = null;
	private String expressionSystemGene = null;
	private String expressionSystemOtherDetails = null;

	private Long id;

	@SuppressWarnings("unchecked")
	public String toString(){
		StringBuffer buf = new StringBuffer();
		buf.append("Compound: " + molId + " " +molName + " ");
		/* disabled for the moment

    	 buf.append(" chains: " );
    	Iterator<Chain> iter = chainList.iterator();
    	while (iter.hasNext()){
    		Chain c = iter.next();
    		buf.append (c.getName() + " ");
    	}

		 */
		try {
			@SuppressWarnings("rawtypes")
			Class c = Class.forName("org.biojava.bio.structure.Compound");
			Method[] methods  = c.getMethods();

			for (int i = 0; i < methods.length; i++) {
				Method m = methods[i];

				String name = m.getName();
				if ( name.substring(0,3).equals("get")) {
					if (name.equals("getMolId"))
						continue;
					if ( name.equals("getMolName"))
						continue;

					Object o  = m.invoke(this, new Object[]{});
					if ( o instanceof String){
						if ( o != null)
							buf.append(name.substring(3, name.length())+": "+ o + " ");
					}
					if ( o instanceof List){
						if ( o != null)
                            buf.append(name.substring(3, name.length())).append(": ");

						List<Object>lst = (List<Object>)o;
						for (Object obj : lst){
							if ( obj instanceof Chain){
								continue;
							}
                            buf.append(obj).append(" ");
						}

					}
				}

			}

		} catch (Exception e){
			e.printStackTrace();
		}


		//if ( organismScientific != null)
		//	buf.append(" organism scientific: " + organismScientific);


		return buf.toString();
	}

	/** get the ID used by Hibernate
	 *
	 * @return the ID used by Hibernate
	 */
	public Long getId() {
		return id;
	}

	/** set the ID used by Hibernate
	 *
	 * @param id
	 */
	public void setId(Long id) {
		this.id = id;
	}

	/**
	 * Print some debug statements to System.out
	 *
	 *
	 */
	public void showHeader(){
		this.showCompound();
		this.showSource();
	}

	public void showCompound() {
		System.out.println("COMPOUND INFO:");
		if (this.molId != null) {
			System.out.println("Mol ID: " + this.molId);
		}
		if (this.chainId != null) {
			System.out.println("Chain: " + this.chainId);
			//this.refChainId = chainId
		}
		if (this.molName != null) {
			System.out.println("Mol Name: " + this.molName);
		}
		if (this.title != null) {
			System.out.println("Title: " + this.title);
		}
		if (this.synonyms != null) {
			for (String x : this.synonyms) {
				System.out.println("Synomym: " + x);
			}
		}
		if (this.ecNums != null) {
			for (String x : this.ecNums) {
				System.out.println("EC: " + x);
			}
		}
		if (this.fragment != null) {
			System.out.println("Fragment? " + this.fragment);
		}
		if (this.engineered != null) {
			System.out.println("Engineered? " + this.engineered);
		}
		if (this.mutation != null) {
			System.out.println("Mutation? " + this.mutation);
		}
		if (this.biologicalUnit != null) {
			System.out.println("Biological Unit: " + this.biologicalUnit);
		}
		if (this.details != null) {
			System.out.println("Details: " + this.details);
		}
		if (this.numRes != null) {
			System.out.println("No. Residues: " + this.numRes);
		}
		//System.out.println( "\n"

	}


	public void showSource() {
		System.out.println("SOURCE INFO:");
		if (this.synthetic != null) {
			System.out.println("Synthetic? " + this.synthetic);
		}
		if (this.fragment != null) {
			System.out.println("Fragment? " + this.fragment);
		}
		if (this.organismScientific != null) {
			System.out.println("Organism Scientific: " + this.organismScientific);
		}
        if (this.organismTaxId != null) {
			System.out.println("Organism Tax Id: " + this.organismTaxId);
		}
		if (this.organismCommon != null) {
			System.out.println("Organism Common: " + this.organismCommon);
		}
		if (this.strain != null) {
			System.out.println("Strain: " + this.strain);
		}
		if (this.variant != null) {
			System.out.println("Variant: " + this.variant);
		}
		if (this.cellLine != null) {
			System.out.println("Cell Line: " + this.cellLine);
		}
		if (this.atcc != null) {
			System.out.println("ATCC: " + this.atcc);
		}
		if (this.organ != null) {
			System.out.println("Organ: " + this.organ);
		}
		if (this.tissue != null) {
			System.out.println("Tissue: " + this.tissue);
		}
		if (this.cell != null) {
			System.out.println("Cell: " + this.cell);
		}
		if (this.organelle != null) {
			System.out.println("Organelle: " + this.organelle);
		}
		if (this.secretion != null) {
			System.out.println("Secretion: " + this.secretion);
		}
		if (this.gene != null) {
			System.out.println("Gene: " + this.gene);
		}
		if (this.cellularLocation != null) {
			System.out.println("Cellular Location: " + this.cellularLocation);
		}
		if (this.expressionSystem != null) {
			System.out.println("Expression System: " + this.expressionSystem);
		}
        if (this.expressionSystemTaxId != null) {
			System.out.println("Expression System Tax Id: " + this.expressionSystemTaxId);
		}
		if (this.expressionSystemStrain != null) {
			System.out.println("Expression System Strain: " + this.expressionSystemStrain);
		}
		if (this.expressionSystemVariant != null) {
			System.out.println("Expression System Variant: " + this.expressionSystemVariant);
		}
		if (this.expressionSystemCellLine != null) {
			System.out.println("Expression System Cell Line: " + this.expressionSystemCellLine);
		}
		if (this.expressionSystemAtccNumber != null) {
			System.out.println("Expression System ATCC Number: " + this.expressionSystemAtccNumber);
		}
		if (this.expressionSystemOrgan != null) {
			System.out.println("Expression System Organ: " + this.expressionSystemOrgan);
		}
		if (this.expressionSystemTissue != null) {
			System.out.println("Expression System Tissue: " + this.expressionSystemTissue);
		}
		if (this.expressionSystemCell != null) {
			System.out.println("Expression System Cell: " + this.expressionSystemCell);
		}
		if (this.expressionSystemOrganelle != null) {
			System.out.println("Expression System Organelle: " + this.expressionSystemOrganelle);
		}
		if (this.expressionSystemCellularLocation != null) {
			System.out.println("Expression System Cellular Location: " + this.expressionSystemCellularLocation);
		}
		if (this.expressionSystemVectorType != null) {
			System.out.println("Expression System Vector Type: " + this.expressionSystemVectorType);
		}
		if (this.expressionSystemVector != null) {
			System.out.println("Expression System Vector: " + this.expressionSystemVector);
		}
		if (this.expressionSystemPlasmid != null) {
			System.out.println("Expression System Plasmid: " + this.expressionSystemPlasmid);
		}
		if (this.expressionSystemGene != null) {
			System.out.println("Expression System Gene: " + this.expressionSystemGene);
		}
		if (this.expressionSystemOtherDetails != null) {
			System.out.println("Expression System Other Details: " + this.expressionSystemOtherDetails);
		}
	}

	/**
	 * Returns the chain id value.
	 * @return the list of ChainIDs that are described by this Compound
	 * @see #setChainId(List)
	 */
	public List<String> getChainId() {
		return chainId;
	}

	/**
	 * Sets the list of chain IDs.
	 * @param chainId  the list of ChainIDs that are described by this Compound
	 * @see #getChainId()
	 */
	public void setChainId(List<String> chainId) {
		this.chainId = chainId;
	}

	/**
	 * Returns the ref chain id value.
	 * @return the RefChainID
	 * @see #setRefChainId(String)
	 */
	public String getRefChainId() {
		return refChainId;
	}

	/**
	 * Returns the ref chain id value.
	 * @param refChainId the RefChainID
	 * @see #getRefChainId()
	 */
	public void setRefChainId(String refChainId) {
		this.refChainId = refChainId;
	}

	/**
	 * Returns the mol id value.
	 * @return the MolId value
	 * @see #setMolId(String)
	 */
	public String getMolId() {
		return molId;
	}

	/**
	 * Set the mol id value.
	 * @param molId the MolId value
	 * @see #getMolId()
	 */
	public void setMolId(String molId) {
		this.molId = molId;
	}

	public String getMolName() {
		return molName;
	}

	public void setMolName(String molName) {
		this.molName = molName;
	}

	public String getTitle() {
		return title;
	}

	public void setTitle(String title) {
		this.title = title;
	}

	public List<String> getSynonyms() {
		return synonyms;
	}

	public void setSynonyms(List<String> synonyms) {
		this.synonyms = synonyms;
	}

	public List<String> getEcNums() {
		return ecNums;
	}

	public void setEcNums(List<String> ecNums) {
		this.ecNums = ecNums;
	}

	public String getEngineered() {
		return engineered;
	}

	public void setEngineered(String engineered) {
		this.engineered = engineered;
	}

	public String getMutation() {
		return mutation;
	}

	public void setMutation(String mutation) {
		this.mutation = mutation;
	}

	public String getBiologicalUnit() {
		return biologicalUnit;
	}

	public void setBiologicalUnit(String biologicalUnit) {
		this.biologicalUnit = biologicalUnit;
	}

	public String getDetails() {
		return details;
	}

	public void setDetails(String details) {
		this.details = details;
	}

	public String getNumRes() {
		return numRes;
	}

	public void setNumRes(String numRes) {
		this.numRes = numRes;
	}

	public String getResNames() {
		return resNames;
	}

	public void setResNames(String resNames) {
		this.resNames = resNames;
	}

	public String getHeaderVars() {
		return headerVars;
	}

	public void setHeaderVars(String headerVars) {
		this.headerVars = headerVars;
	}

	public String getSynthetic() {
		return synthetic;
	}

	public void setSynthetic(String synthetic) {
		this.synthetic = synthetic;
	}

	public String getFragment() {
		return fragment;
	}

	public void setFragment(String fragment) {
		this.fragment = fragment;
	}

	public String getOrganismScientific() {
		return organismScientific;
	}

	public void setOrganismScientific(String organismScientific) {
		this.organismScientific = organismScientific;
	}

    public String getOrganismTaxId() {
        return organismTaxId;
    }

    public void setOrganismTaxId(String organismTaxId) {
        this.organismTaxId = organismTaxId;
    }

	public String getOrganismCommon() {
		return organismCommon;
	}

	public void setOrganismCommon(String organismCommon) {
		this.organismCommon = organismCommon;
	}

	public String getStrain() {
		return strain;
	}

	public void setStrain(String strain) {
		this.strain = strain;
	}

	public String getVariant() {
		return variant;
	}

	public void setVariant(String variant) {
		this.variant = variant;
	}

	public String getCellLine() {
		return cellLine;
	}

	public void setCellLine(String cellLine) {
		this.cellLine = cellLine;
	}

	public String getAtcc() {
		return atcc;
	}

	public void setAtcc(String atcc) {
		this.atcc = atcc;
	}

	public String getOrgan() {
		return organ;
	}

	public void setOrgan(String organ) {
		this.organ = organ;
	}

	public String getTissue() {
		return tissue;
	}

	public void setTissue(String tissue) {
		this.tissue = tissue;
	}

	public String getCell() {
		return cell;
	}

	public void setCell(String cell) {
		this.cell = cell;
	}

	public String getOrganelle() {
		return organelle;
	}

	public void setOrganelle(String organelle) {
		this.organelle = organelle;
	}

	public String getSecretion() {
		return secretion;
	}

	public void setSecretion(String secretion) {
		this.secretion = secretion;
	}

	public String getGene() {
		return gene;
	}

	public void setGene(String gene) {
		this.gene = gene;
	}

	public String getCellularLocation() {
		return cellularLocation;
	}

	public void setCellularLocation(String cellularLocation) {
		this.cellularLocation = cellularLocation;
	}

	public String getExpressionSystem() {
		return expressionSystem;
	}

    public String getExpressionSystemTaxId() {
        return expressionSystemTaxId;
    }

    public void setExpressionSystemTaxId(String expressionSystemTaxId) {
        this.expressionSystemTaxId = expressionSystemTaxId;
    }

	public void setExpressionSystem(String expressionSystem) {
		this.expressionSystem = expressionSystem;
	}

	public String getExpressionSystemStrain() {
		return expressionSystemStrain;
	}

	public void setExpressionSystemStrain(String expressionSystemStrain) {
		this.expressionSystemStrain = expressionSystemStrain;
	}

	public String getExpressionSystemVariant() {
		return expressionSystemVariant;
	}

	public void setExpressionSystemVariant(String expressionSystemVariant) {
		this.expressionSystemVariant = expressionSystemVariant;
	}

	public String getExpressionSystemCellLine() {
		return expressionSystemCellLine;
	}

	public void setExpressionSystemCellLine(String expressionSystemCellLine) {
		this.expressionSystemCellLine = expressionSystemCellLine;
	}

	public String getExpressionSystemAtccNumber() {
		return expressionSystemAtccNumber;
	}

	public void setExpressionSystemAtccNumber(String expressionSystemAtccNumber) {
		this.expressionSystemAtccNumber = expressionSystemAtccNumber;
	}

	public String getExpressionSystemOrgan() {
		return expressionSystemOrgan;
	}

	public void setExpressionSystemOrgan(String expressionSystemOrgan) {
		this.expressionSystemOrgan = expressionSystemOrgan;
	}

	public String getExpressionSystemTissue() {
		return expressionSystemTissue;
	}

	public void setExpressionSystemTissue(String expressionSystemTissue) {
		this.expressionSystemTissue = expressionSystemTissue;
	}

	public String getExpressionSystemCell() {
		return expressionSystemCell;
	}

	public void setExpressionSystemCell(String expressionSystemCell) {
		this.expressionSystemCell = expressionSystemCell;
	}

	public String getExpressionSystemOrganelle() {
		return expressionSystemOrganelle;
	}

	public void setExpressionSystemOrganelle(String expressionSystemOrganelle) {
		this.expressionSystemOrganelle = expressionSystemOrganelle;
	}

	public String getExpressionSystemCellularLocation() {
		return expressionSystemCellularLocation;
	}

	public void setExpressionSystemCellularLocation(String expressionSystemCellularLocation) {
		this.expressionSystemCellularLocation = expressionSystemCellularLocation;
	}

	public String getExpressionSystemVectorType() {
		return expressionSystemVectorType;
	}

	public void setExpressionSystemVectorType(String expressionSystemVectorType) {
		this.expressionSystemVectorType = expressionSystemVectorType;
	}

	public String getExpressionSystemVector() {
		return expressionSystemVector;
	}

	public void setExpressionSystemVector(String expressionSystemVector) {
		this.expressionSystemVector = expressionSystemVector;
	}

	public String getExpressionSystemPlasmid() {
		return expressionSystemPlasmid;
	}

	public void setExpressionSystemPlasmid(String expressionSystemPlasmid) {
		this.expressionSystemPlasmid = expressionSystemPlasmid;
	}

	public String getExpressionSystemGene() {
		return expressionSystemGene;
	}

	public void setExpressionSystemGene(String expressionSystemGene) {
		this.expressionSystemGene = expressionSystemGene;
	}

	public String getExpressionSystemOtherDetails() {
		return expressionSystemOtherDetails;
	}

	public void setExpressionSystemOtherDetails(String expressionSystemOtherDetails) {
		this.expressionSystemOtherDetails = expressionSystemOtherDetails;
	}

	public Compound clone() throws CloneNotSupportedException {
		Compound newMolId = (Compound) super.clone();
		return newMolId;
	}

	/** get the chains that are part of this Compound
	 *
	 * @return a List of Chain objects
	 */
	 public List<Chain> getChains(){
		return this.chainList;
	}

	public void addChain(Chain chain){
		this.chainList.add(chain);
	}

	public void setChains(List<Chain> chains){
		this.chainList = chains;
	}
}