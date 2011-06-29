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
 * Created on 2011.05.09 by kohchuanhock
 *
 */
package org.biojava3.aaproperties;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Map;

import javax.xml.bind.JAXBException;

import org.biojava3.aaproperties.xml.AminoAcidCompositionTable;
import org.biojava3.core.sequence.ProteinSequence;
import org.biojava3.core.sequence.compound.AminoAcidCompound;

/**
 * TODO
 * 4) Write tutorials for BioJava - Cookbook BioJava Wiki -> ProteinStructure
 * 5) Make the xml file more restrictive.
 * 6) Simplify the MW.xml and give a isotope example.
 * 
 * DONE
 * 9) Next week to be on wednesday
 * 
 * Question
 * 1) Where should the method auto locate the elementMass file
 * 
 * 
 * 
 * 
 * An interface to generate some basic physico-chemical properties of protein sequences.<br/>
 * The following properties could be generated:
 * <p/>
 * Molecular weight<br/>
 * Absorbance<br/>
 * Extinction coefficient<br/>
 * Instability index<br/>
 * Apliphatic index<br/>
 * Average hydropathy value<br/>
 * Isoelectric point<br/>
 * Net charge at pH 7<br/>
 * Composition of specified amino acid<br/>
 * Composition of the 20 standard amino acid<br/>
 * @author kohchuanhock
 * @version 2011.05.09
 * @see PeptideProperties
 */
public interface IPeptideProperties{
	/**
	 * Returns the molecular weight of sequence. The sequence argument must be a protein sequence consisting of only non-ambiguous characters.
	 * This method will sum the molecular weight of each amino acid in the
	 * sequence. Molecular weights are based on <a href="http://au.expasy.org/tools/findmod/findmod_masses.html#AA">here</a>.
	 * 
	 * @param sequence
	 * 	a protein sequence consisting of non-ambiguous characters only
	 * @return the total molecular weight of sequence + weight of water molecule
	 * @see ProteinSequence
	 */
	public double getMolecularWeight(ProteinSequence sequence);
	
	/**
	 * Returns the molecular weight of sequence. The sequence argument must be a protein sequence consisting of only non-ambiguous characters.
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
	public double getMolecularWeight(ProteinSequence sequence, File aminoAcidCompositionFile) throws JAXBException, FileNotFoundException;
	
	/**
	 * Returns the molecular weight of sequence. The sequence argument must be a protein sequence consisting of only non-ambiguous characters.
	 * This method will sum the molecular weight of each amino acid in the
	 * sequence. Molecular weights are based on the input files. These input files must be XML using the defined schema. 
	 * 
	 * @param sequence
	 * 	a protein sequence consisting of non-ambiguous characters only
	 * @param elementMassFile
	 * 	xml file that details the mass of each elements and isotopes
	 * @param aminoAcidCompositionFile
	 * 	xml file that details the composition of amino acids
	 * @return the total molecular weight of sequence + weight of water molecule
	 * @throws JAXBException
	 * 	thrown if unable to properly parse either elementMassFile or aminoAcidCompositionFile
	 * @throws FileNotFoundException
	 * 	thrown if either elementMassFile or aminoAcidCompositionFile are not found
	 */
	public double getMolecularWeight(ProteinSequence sequence, File elementMassFile, File aminoAcidCompositionFile) throws JAXBException, FileNotFoundException;
	
	/**
	 * Returns the molecular weight of sequence. The sequence argument must be a protein sequence consisting of only non-ambiguous characters.
	 * This method will sum the molecular weight of each amino acid in the
	 * sequence. Molecular weights are based on the AminoAcidCompositionTable. 
	 * Those input files must be XML using the defined schema.
	 * 
	 * @param sequence
	 * 	a protein sequence consisting of non-ambiguous characters only
	 * @param aminoAcidCompositionTable
	 * 	a amino acid composition table obtained by calling IPeptideProperties.obtainAminoAcidCompositionTable
	 * @return the total molecular weight of sequence + weight of water molecule
	 * @throws Exception
	 * 	thrown if the method IPeptideProperties.setMolecularWeightXML(File, File) is not successfully called before calling this method.
	 */
	public double getMolecularWeightBasedOnXML(ProteinSequence sequence, AminoAcidCompositionTable aminoAcidCompositionTable) throws Exception;
	
	/**
	 * This method would initialize amino acid composition table based on the input xml files and stores the table for usage in future calls to 
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
	public AminoAcidCompositionTable obtainAminoAcidCompositionTable(File elementMassFile, File aminoAcidCompositionFile) throws JAXBException, FileNotFoundException;

	/**
	 * Returns the extinction coefficient of sequence. The sequence argument
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
	 * @see ProteinSequence
	 */
	public double getExtinctionCoefficient(ProteinSequence sequence, boolean assumeCysReduced);

	/**
	 * Returns the absorbance (optical density) of sequence. The sequence argument
	 * must be a protein sequence consisting of only non-ambiguous characters.
	 * The computation of absorbance (optical density) follows the
	 * documentation in <a href="http://web.expasy.org/tools/protparam-doc.html">here</a>.
	 * 
	 * @param sequence
	 *            a protein sequence consisting of non-ambiguous characters only
	 * @param assumeCysReduced
	 *            true if Cys are assumed to be reduced and false if Cys are
	 *            assumed to form cystines
	 * @return the absorbance (optical density) of sequence
	 * @see ProteinSequence
	 */
	public double getAbsorbance(ProteinSequence sequence, boolean assumeCysReduced);
	
	/**
	 * Returns the instability index of sequence. The sequence argument must be
	 * a protein sequence consisting of only non-ambiguous characters.
	 * The instability index provides an estimate of the stability of your
	 * protein in a test tube. The computation of instability index follows the
	 * documentation in <a href="http://web.expasy.org/tools/protparam-doc.html">here</a>.
	 * 
	 * @param sequence
	 *            a protein sequence consisting of non-ambiguous characters only
	 * @return the instability index of sequence
	 * @see ProteinSequence
	 */
	public double getInstabilityIndex(ProteinSequence sequence);

	/**
	 * Returns the apliphatic index of sequence. The sequence argument must be a
	 * protein sequence consisting of only non-ambiguous characters.
	 * The aliphatic index of a protein is defined as the relative volume
	 * occupied by aliphatic side chains (alanine, valine, isoleucine, and
	 * leucine). It may be regarded as a positive factor for the increase of
	 * thermostability of globular proteins. The computation of aliphatic index
	 * follows the documentation in <a href="http://web.expasy.org/tools/protparam-doc.html">here</a>.
	 * A protein whose instability index is smaller than 40 is predicted as stable, a value above 40 predicts that the protein may be unstable.
	 * 
	 * @param sequence
	 *            a protein sequence consisting of non-ambiguous characters only
	 * @return the aliphatic index of sequence
	 * @see ProteinSequence
	 */
	public double getApliphaticIndex(ProteinSequence sequence);

	/**
	 * Returns the average hydropathy value of sequence. The sequence argument
	 * must be a protein sequence consisting of only non-ambiguous characters.
	 * The average value for a sequence is calculated as the sum of hydropathy
	 * values of all the amino acids, divided by the number of residues in the
	 * sequence. Hydropathy values are based on (Kyte, J. and Doolittle, R.F.
	 * (1982) A simple method for displaying the hydropathic character of a
	 * protein. J. Mol. Biol. 157, 105-132).
	 * 
	 * @param sequence
	 *            a protein sequence consisting of non-ambiguous characters only
	 * @return the average hydropathy value of sequence
	 * @see ProteinSequence
	 */
	public double getAvgHydropathy(ProteinSequence sequence);

	/**
	 * Returns the isoelectric point of sequence. The sequence argument must be
	 * a protein sequence consisting of only non-ambiguous characters.
	 * The isoelectric point is the pH at which the protein carries no net
	 * electrical charge. The isoelectric point will be computed based on
	 * approach stated in 
	 * <a href="http://www.innovagen.se/custom-peptide-synthesis/peptide-property-calculator/peptide-property-calculator-notes.asp#PI">here</a>
	 * 
	 * @param sequence
	 *            a protein sequence consisting of non-ambiguous characters only
	 * @return the isoelectric point of sequence
	 * @see ProteinSequence
	 */
	public double getIsoelectricPoint(ProteinSequence sequence);

	/**
	 * Returns the net charge of sequence at pH 7. The sequence argument must be
	 * a protein sequence consisting of only non-ambiguous characters.
	 * The net charge will be computed using the approach stated in 
	 * <a href="http://www.innovagen.se/custom-peptide-synthesis/peptide-property-calculator/peptide-property-calculator-notes.asp#NetCharge>here</a>
	 * 
	 * @param sequence
	 *            a protein sequence consisting of non-ambiguous characters only
	 * @return the net charge of sequence at pH 7
	 * @see ProteinSequence
	 */
	public double getNetCharge(ProteinSequence sequence);

	/**
	 * Returns the composition of specified amino acid in the sequence. The
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
	 * @see ProteinSequence 
	 * @see AminoAcidCompound
	 */
	public double getEnrichment(ProteinSequence sequence, AminoAcidCompound aminoAcidCode);

	/**
	 * Returns the composition of the 20 standard amino acid in the sequence.
	 * The sequence argument must be a protein sequence consisting of only
	 * non-ambiguous characters.
	 * The composition of an amino acid is the total number of its occurrence,
	 * divided by the total length of the sequence.
	 * 
	 * @param sequence
	 *            a protein sequence consisting of non-ambiguous characters only
	 * @return the composition of the 20 standard amino acid in the sequence
	 * @see ProteinSequence 
	 * @see AminoAcidCompound
	 */
	public Map<AminoAcidCompound, Double> getAAComposition(ProteinSequence sequence);
}
