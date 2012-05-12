package org.biojava3.aaproperties;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import javax.xml.bind.JAXBException;

import org.biojava3.aaproperties.xml.AminoAcidCompositionTable;
import org.biojava3.core.sequence.ProteinSequence;
import org.biojava3.core.sequence.compound.AminoAcidCompound;
import org.biojava3.core.sequence.compound.AminoAcidCompoundSet;

/**
 * This is an adaptor class which enable the ease of generating protein properties.
 * At least one adaptor method is written for each available properties provided in IPeptideProperties. 
 * 
 * @author kohchuanhock
 * @version 2011.08.22
 * @since 3.0.2
 * @see IPeptideProperties
 * @see PeptidePropertiesImpl
 */
public class PeptideProperties {
	/**
	 * Enumeration of 20 standard amino acid code
	 */
	public enum SingleLetterAACode { W, C, M, H, Y, F, Q, N, I, R, D, P, T, K, E, V, S, G, A, L}

	/**
	 * Contains the 20 standard AA code in a set
	 */
	public static Set<Character> standardAASet;

	/**
	 * To initialize the standardAASet
	 */
	static{
		standardAASet = new HashSet<Character>();
		for(SingleLetterAACode c:SingleLetterAACode.values()) standardAASet.add(c.toString().charAt(0));
	}

	/**
	 * An adaptor method to return the molecular weight of sequence. 
	 * The sequence argument must be a protein sequence consisting of only non-ambiguous characters.
	 * This method will sum the molecular weight of each amino acid in the
	 * sequence. Molecular weights are based on <a href="http://web.expasy.org/findmod/findmod_masses.html">here</a>.
	 * 
	 * @param sequence
	 * 	a protein sequence consisting of non-ambiguous characters only
	 * @return the total molecular weight of sequence + weight of water molecule
	 */
	public static final double getMolecularWeight(String sequence){
		sequence = Utils.checkSequence(sequence);
		ProteinSequence pSequence = new ProteinSequence(sequence);
		IPeptideProperties pp = new PeptidePropertiesImpl();
		return pp.getMolecularWeight(pSequence);
	}

	/**
	 * An adaptor method to return the molecular weight of sequence. 
	 * The sequence argument must be a protein sequence consisting of only non-ambiguous characters.
	 * This method will sum the molecular weight of each amino acid in the
	 * sequence. Molecular weights are based on the input xml file.
	 * 
	 * @param sequence
	 * 	a protein sequence consisting of non-ambiguous characters only
	 * @param elementMassFile 
	 * 	xml file that details the mass of each elements and isotopes
	 * @param aminoAcidCompositionFile 
	 * 	xml file that details the composition of amino acids
	 * @return the total molecular weight of sequence + weight of water molecule
	 * @throws FileNotFoundException 
	 * 	thrown if either elementMassFile or aminoAcidCompositionFile are not found
	 * @throws JAXBException 
	 * 	thrown if unable to properly parse either elementMassFile or aminoAcidCompositionFile
	 */
	public static final double getMolecularWeight(String sequence, File elementMassFile, File aminoAcidCompositionFile) 
	throws FileNotFoundException, JAXBException{
		sequence = Utils.checkSequence(sequence);
		ProteinSequence pSequence = new ProteinSequence(sequence);
		IPeptideProperties pp = new PeptidePropertiesImpl();
		return pp.getMolecularWeight(pSequence, elementMassFile, aminoAcidCompositionFile);
	}

	/**
	 * An adaptor method to return the molecular weight of sequence. The sequence argument must be a protein sequence consisting of only non-ambiguous characters.
	 * This method will sum the molecular weight of each amino acid in the
	 * sequence. Molecular weights are based on the input files. These input files must be XML using the defined schema.
	 * Note that it assumes that ElementMass.xml file can be found in default location.
	 * 
	 * @param sequence
	 * 	a protein sequence consisting of non-ambiguous characters only
	 * 	xml file that details the mass of each elements and isotopes
	 * @param aminoAcidCompositionFile
	 * 	xml file that details the composition of amino acids
	 * @return the total molecular weight of sequence + weight of water molecule
	 * @throws JAXBException
	 * 	thrown if unable to properly parse either elementMassFile or aminoAcidCompositionFile
	 * @throws FileNotFoundException
	 * 	thrown if either elementMassFile or aminoAcidCompositionFile are not found
	 */
	public static final double getMolecularWeight(String sequence, File aminoAcidCompositionFile) throws FileNotFoundException, JAXBException{
		sequence = Utils.checkSequence(sequence);
		ProteinSequence pSequence = new ProteinSequence(sequence);
		IPeptideProperties pp = new PeptidePropertiesImpl();
		return pp.getMolecularWeight(pSequence, aminoAcidCompositionFile);
	}

	/**
	 * An adaptor method would initialize amino acid composition table based on the input xml files and stores the table for usage in future calls to 
	 * IPeptideProperties.getMolecularWeightBasedOnXML(ProteinSequence, AminoAcidCompositionTable).
	 * Note that ElementMass.xml is assumed to be able to be seen in default location.
	 * 
	 * @param aminoAcidCompositionFile
	 * 	xml file that details the composition of amino acids
	 * @return the initialized amino acid composition table
	 * @throws JAXBException
	 * 	thrown if unable to properly parse either elementMassFile or aminoAcidCompositionFile
	 * @throws FileNotFoundException
	 * 	thrown if either elementMassFile or aminoAcidCompositionFile are not found
	 */
	public static final AminoAcidCompositionTable obtainAminoAcidCompositionTable(File aminoAcidCompositionFile) 
	throws JAXBException, FileNotFoundException{
		IPeptideProperties pp = new PeptidePropertiesImpl();
		return pp.obtainAminoAcidCompositionTable(aminoAcidCompositionFile);
	}

	/**
	 * An adaptor method would initialize amino acid composition table based on the input xml files and stores the table for usage in future calls to 
	 * IPeptideProperties.getMolecularWeightBasedOnXML(ProteinSequence, AminoAcidCompositionTable).
	 * 
	 * @param elementMassFile
	 * 	xml file that details the mass of each elements and isotopes
	 * @param aminoAcidCompositionFile
	 * 	xml file that details the composition of amino acids
	 * @return the initialized amino acid composition table
	 * @throws JAXBException
	 * 	thrown if unable to properly parse either elementMassFile or aminoAcidCompositionFile
	 * @throws FileNotFoundException
	 * 	thrown if either elementMassFile or aminoAcidCompositionFile are not found
	 */
	public static final AminoAcidCompositionTable obtainAminoAcidCompositionTable(File elementMassFile, File aminoAcidCompositionFile)
	throws JAXBException, FileNotFoundException{
		IPeptideProperties pp = new PeptidePropertiesImpl();
		return pp.obtainAminoAcidCompositionTable(elementMassFile, aminoAcidCompositionFile);
	}

	/**
	 * An adaptor method that returns the molecular weight of sequence. The sequence argument must be a protein sequence consisting of only non-ambiguous characters.
	 * This method will sum the molecular weight of each amino acid in the
	 * sequence. Molecular weights are based on the AminoAcidCompositionTable. 
	 * Those input files must be XML using the defined schema.
	 * 
	 * @param sequence
	 * 	a protein sequence consisting of non-ambiguous characters only
	 * @param aminoAcidCompositionTable
	 * 	a amino acid composition table obtained by calling IPeptideProperties.obtainAminoAcidCompositionTable
	 * @return the total molecular weight of sequence + weight of water molecule
	 * 	thrown if the method IPeptideProperties.setMolecularWeightXML(File, File) is not successfully called before calling this method.
	 */
	public static double getMolecularWeightBasedOnXML(String sequence, AminoAcidCompositionTable aminoAcidCompositionTable){
		sequence = Utils.checkSequence(sequence, aminoAcidCompositionTable.getSymbolSet());
		ProteinSequence pSequence = new ProteinSequence(sequence, aminoAcidCompositionTable.getAminoAcidCompoundSet());
		IPeptideProperties pp = new PeptidePropertiesImpl();
		return pp.getMolecularWeightBasedOnXML(pSequence, aminoAcidCompositionTable);
	}

	/**
	 * An adaptor method to returns the absorbance (optical density) of sequence. The sequence argument
	 * must be a protein sequence consisting of only non-ambiguous characters.
	 * The computation of absorbance (optical density) follows the
	 * documentation in <a href="http://web.expasy.org/protparam/protparam-doc.html">here</a>.
	 * 
	 * @param sequence
	 * 	a protein sequence consisting of non-ambiguous characters only
	 * @param assumeCysReduced
	 * 	true if Cys are assumed to be reduced and false if Cys are assumed to form cystines
	 * @return the absorbance (optical density) of sequence
	 */
	public static final double getAbsorbance(String sequence, boolean assumeCysReduced){
		sequence = Utils.checkSequence(sequence);
		ProteinSequence pSequence = new ProteinSequence(sequence);
		IPeptideProperties pp = new PeptidePropertiesImpl();
		return pp.getAbsorbance(pSequence, assumeCysReduced);
	}

	/**
	 * An adaptor method to return the extinction coefficient of sequence. The sequence argument
	 * must be a protein sequence consisting of only non-ambiguous characters.
	 * The extinction coefficient indicates how much light a protein absorbs at
	 * a certain wavelength. It is useful to have an estimation of this
	 * coefficient for following a protein which a spectrophotometer when
	 * purifying it. The computation of extinction coefficient follows the
	 * documentation in <a href="http://web.expasy.org/protparam/protparam-doc.html">here</a>.
	 * 
	 * @param sequence
	 *            a protein sequence consisting of non-ambiguous characters only
	 * @param assumeCysReduced
	 *            true if Cys are assumed to be reduced and false if Cys are
	 *            assumed to form cystines
	 * @return the extinction coefficient of sequence
	 */
	public static final double getExtinctionCoefficient(String sequence, boolean assumeCysReduced) {
		sequence = Utils.checkSequence(sequence);
		ProteinSequence pSequence = new ProteinSequence(sequence);
		IPeptideProperties pp = new PeptidePropertiesImpl();
		return pp.getExtinctionCoefficient(pSequence, assumeCysReduced);
	}

	/**
	 * An adaptor method to return the instability index of sequence. The sequence argument must be
	 * a protein sequence consisting of only non-ambiguous characters.
	 * The instability index provides an estimate of the stability of your
	 * protein in a test tube. The computation of instability index follows the
	 * documentation in <a href="http://web.expasy.org/protparam/protparam-doc.html">here</a>.
	 * 
	 * @param sequence
	 * 		a protein sequence consisting of non-ambiguous characters only
	 * @return the instability index of sequence
	 */
	public static final double getInstabilityIndex(String sequence) {
		sequence = Utils.checkSequence(sequence);
		ProteinSequence pSequence = new ProteinSequence(sequence);
		IPeptideProperties pp = new PeptidePropertiesImpl();
		return pp.getInstabilityIndex(pSequence);
	}

	/**
	 * An adaptor method to return the apliphatic index of sequence. The sequence argument must be a
	 * protein sequence consisting of only non-ambiguous characters.
	 * The aliphatic index of a protein is defined as the relative volume
	 * occupied by aliphatic side chains (alanine, valine, isoleucine, and
	 * leucine). It may be regarded as a positive factor for the increase of
	 * thermostability of globular proteins. The computation of aliphatic index
	 * follows the documentation in <a href="http://web.expasy.org/protparam/protparam-doc.html">here</a>.
	 * A protein whose instability index is smaller than 40 is predicted as stable, a value above 40 predicts that the protein may be unstable.
	 * 
	 * @param sequence
	 * 		a protein sequence consisting of non-ambiguous characters only
	 * @return the aliphatic index of sequence
	 */
	public static final double getApliphaticIndex(String sequence) {
		sequence = Utils.checkSequence(sequence);
		ProteinSequence pSequence = new ProteinSequence(sequence);
		IPeptideProperties pp = new PeptidePropertiesImpl();
		return pp.getApliphaticIndex(pSequence);
	}

	/**
	 * An adaptor method to return the average hydropathy value of sequence. The sequence argument
	 * must be a protein sequence consisting of only non-ambiguous characters.
	 * The average value for a sequence is calculated as the sum of hydropathy
	 * values of all the amino acids, divided by the number of residues in the
	 * sequence. Hydropathy values are based on (Kyte, J. and Doolittle, R.F.
	 * (1982) A simple method for displaying the hydropathic character of a
	 * protein. J. Mol. Biol. 157, 105-132).
	 * 
	 * @param sequence 
	 * 		a protein sequence consisting of non-ambiguous characters only
	 * @return the average hydropathy value of sequence
	 */
	public static final double getAvgHydropathy(String sequence) {
		sequence = Utils.checkSequence(sequence);
		ProteinSequence pSequence = new ProteinSequence(sequence);
		IPeptideProperties pp = new PeptidePropertiesImpl();
		return pp.getAvgHydropathy(pSequence);
	}

	/**
	 * An adaptor method to return the isoelectric point of sequence. The sequence argument must be
	 * a protein sequence consisting of only non-ambiguous characters.
	 * The isoelectric point is the pH at which the protein carries no net
	 * electrical charge. The isoelectric point will be computed based on
	 * approach stated in 
	 * <a href="http://www.innovagen.se/custom-peptide-synthesis/peptide-property-calculator/peptide-property-calculator-notes.asp#PI">here</a>
	 * 
	 * pKa values used will be either 
	 * those used by Expasy which referenced "Electrophoresis 1994, 15, 529-539" 
	 * OR
	 * A.Lehninger, Principles of Biochemistry, 4th Edition (2005), Chapter 3, page78, Table 3-1.
	 * 
	 * @param sequence
	 * 		a protein sequence consisting of non-ambiguous characters only
	 * @param useExpasyValues
	 * 		whether to use Expasy values (Default) or Innovagen values
	 * @return the isoelectric point of sequence
	 */
	public static final double getIsoelectricPoint(String sequence, boolean useExpasyValues) {
		sequence = Utils.checkSequence(sequence);
		ProteinSequence pSequence = new ProteinSequence(sequence);
		IPeptideProperties pp = new PeptidePropertiesImpl();
		return pp.getIsoelectricPoint(pSequence, useExpasyValues);
	}
	
	public static final double getIsoelectricPoint(String sequence){
		return getIsoelectricPoint(sequence, true);
	}

	/**
	 * An adaptor method to return the net charge of sequence at pH 7. The sequence argument must be
	 * a protein sequence consisting of only non-ambiguous characters.
	 * The net charge will be computed using the approach stated in 
	 * <a href="http://www.innovagen.se/custom-peptide-synthesis/peptide-property-calculator/peptide-property-calculator-notes.asp#PI">here</a>
	 * 
	 * pKa values used will be either 
	 * those used by Expasy which referenced "Electrophoresis 1994, 15, 529-539" 
	 * OR
	 * A.Lehninger, Principles of Biochemistry, 4th Edition (2005), Chapter 3, page78, Table 3-1.
	 * 
	 * @param sequence
	 * 		a protein sequence consisting of non-ambiguous characters only
	 * @param useExpasyValues
	 * 		whether to use Expasy values (Default) or Innovagen values
	 * @param pHPoint
	 * 		the pH value to use for computation of the net charge. Default at 7.
	 * @return the net charge of sequence at given pHPoint
	 */
	public static final double getNetCharge(String sequence, boolean useExpasyValues, double pHPoint){
		sequence = Utils.checkSequence(sequence);
		ProteinSequence pSequence = new ProteinSequence(sequence);
		IPeptideProperties pp = new PeptidePropertiesImpl();
		return pp.getNetCharge(pSequence, useExpasyValues, pHPoint);
	}
	
	public static final double getNetCharge(String sequence, boolean useExpasyValues) {
		return getNetCharge(sequence, useExpasyValues, 7.0);
	}
	
	public static final double getNetCharge(String sequence){
		return getNetCharge(sequence, true);
	}

	/**
	 * An adaptor method to return the composition of specified amino acid in the sequence. The
	 * sequence argument must be a protein sequence consisting of only
	 * non-ambiguous characters. The aminoAcidCode must be a non-ambiguous
	 * character.
	 * The composition of an amino acid is the total number of its occurrence,
	 * divided by the total length of the sequence.
	 * 
	 * @param sequence
	 *            a protein sequence consisting of non-ambiguous characters only
	 * @param aminoAcidCode
	 *            the code of the amino acid to compute
	 * @return the composition of specified amino acid in the sequence
	 * @see SingleLetterAACode
	 */
	public static final double getEnrichment(String sequence, SingleLetterAACode aminoAcidCode) {
		return getEnrichment(sequence, aminoAcidCode.toString());
	}

	/**
	 * An adaptor method to return the composition of specified amino acid in the sequence. The
	 * sequence argument must be a protein sequence consisting of only
	 * non-ambiguous characters. The aminoAcidCode must be a non-ambiguous
	 * character.
	 * The composition of an amino acid is the total number of its occurrence,
	 * divided by the total length of the sequence.
	 * 
	 * @param sequence
	 * 		a protein sequence consisting of non-ambiguous characters only
	 * @param aminoAcidCode
	 * 		the code of the amino acid to compute
	 * @return the composition of specified amino acid in the sequence
	 */
	public static final double getEnrichment(String sequence, char aminoAcidCode){
		return getEnrichment(sequence, aminoAcidCode + "");
	}

	/**
	 * An adaptor method to return the composition of specified amino acid in the sequence. The
	 * sequence argument must be a protein sequence consisting of only
	 * non-ambiguous characters. The aminoAcidCode must be a non-ambiguous
	 * character.
	 * The composition of an amino acid is the total number of its occurrence,
	 * divided by the total length of the sequence.
	 * 
	 * @param sequence
	 *            a protein sequence consisting of non-ambiguous characters only
	 * @param aminoAcidCode
	 *            the code of the amino acid to compute
	 * @return the composition of specified amino acid in the sequence
	 */
	public static final double getEnrichment(String sequence, String aminoAcidCode){
		sequence = Utils.checkSequence(sequence);
		ProteinSequence pSequence = new ProteinSequence(sequence);
		IPeptideProperties pp = new PeptidePropertiesImpl();
		AminoAcidCompoundSet aaSet = new AminoAcidCompoundSet();
		return pp.getEnrichment(pSequence, aaSet.getCompoundForString(aminoAcidCode));
	}

	/**
	 * An adaptor method to return the composition of the 20 standard amino acid in the sequence.
	 * The sequence argument must be a protein sequence consisting of only
	 * non-ambiguous characters.
	 * The composition of an amino acid is the total number of its occurrence,
	 * divided by the total length of the sequence.
	 * 
	 * @param sequence
	 *            a protein sequence consisting of non-ambiguous characters only
	 * @return the composition of the 20 standard amino acid in the sequence
	 * @see AminoAcidCompound
	 */
	public static final Map<AminoAcidCompound, Double> getAAComposition(String sequence) {
		sequence = Utils.checkSequence(sequence);
		ProteinSequence pSequence = new ProteinSequence(sequence);
		IPeptideProperties pp = new PeptidePropertiesImpl();
		return pp.getAAComposition(pSequence);
	}

	/**
	 * An adaptor method to return the composition of the 20 standard amino acid in the sequence.
	 * The sequence argument must be a protein sequence consisting of only
	 * non-ambiguous characters.
	 * The composition of an amino acid is the total number of its occurrence,
	 * divided by the total length of the sequence.
	 * 
	 * @param sequence
	 * 		a protein sequence consisting of non-ambiguous characters only
	 * @return the composition of the 20 standard amino acid in the sequence
	 */
	public static final Map<String, Double> getAACompositionString(String sequence){
		Map<AminoAcidCompound, Double> aa2Composition = getAAComposition(sequence);
		Map<String, Double> aaString2Composition = new HashMap<String, Double>();
		for(AminoAcidCompound aaCompound:aa2Composition.keySet()){
			aaString2Composition.put(aaCompound.getShortName(), aa2Composition.get(aaCompound));
		}
		return aaString2Composition;
	}

	/**
	 * An adaptor method to return the composition of the 20 standard amino acid in the sequence.
	 * The sequence argument must be a protein sequence consisting of only
	 * non-ambiguous characters.
	 * The composition of an amino acid is the total number of its occurrence,
	 * divided by the total length of the sequence.
	 * 
	 * @param sequence
	 * 		a protein sequence consisting of non-ambiguous characters only
	 * @return the composition of the 20 standard amino acid in the sequence
	 */
	public static final Map<Character, Double> getAACompositionChar(String sequence){
		Map<AminoAcidCompound, Double> aa2Composition = getAAComposition(sequence);
		Map<Character, Double> aaChar2Composition = new HashMap<Character, Double>();
		for(AminoAcidCompound aaCompound:aa2Composition.keySet()){
			aaChar2Composition.put(aaCompound.getShortName().charAt(0), aa2Composition.get(aaCompound));
		}
		return aaChar2Composition;
	}
}	
