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
 * Created on July 26, 2010
 * Author: Mark Chapman
 */

package org.biojava3.alignment;

import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.Map;

import org.biojava3.alignment.aaindex.AAindexFactory;
import org.biojava3.alignment.template.SubstitutionMatrix;
import org.biojava3.core.sequence.compound.AmbiguityDNACompoundSet;
import org.biojava3.core.sequence.compound.AminoAcidCompound;
import org.biojava3.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava3.core.sequence.compound.NucleotideCompound;


/**
 * Static utility to access substitution matrices that come bundled with BioJava.  All matrices were downloaded from
 * ftp://ftp.ncbi.nih.gov/blast/matrices/
 *
 * @author Mark Chapman
 */
public class SubstitutionMatrixHelper {

    private static Map<String, SubstitutionMatrix<AminoAcidCompound>> aminoAcidMatrices =
            new HashMap<String, SubstitutionMatrix<AminoAcidCompound>>();
    private static Map<String, SubstitutionMatrix<NucleotideCompound>> nucleotideMatrices =
            new HashMap<String, SubstitutionMatrix<NucleotideCompound>>();

    // prevents instantiation
    private SubstitutionMatrixHelper() { }

    
    /** Returns any matrix from the AAINDEX database file
     * 
     * @param matrixName
     * @return a {@link SubstitutionMatrix}
     */
    public static SubstitutionMatrix<AminoAcidCompound> getMatrixFromAAINDEX(String matrixName){
    	
    	return AAindexFactory.getAAIndexProvider().getMatrix(matrixName);
    	
    }
    
    
    /**
     * Returns Blosum 100 matrix by Henikoff & Henikoff
     * @return Blosum 100 matrix
     */
    public static SubstitutionMatrix<AminoAcidCompound> getBlosum100() {
        return getAminoAcidCompoundSubstitutionMatrix("blosum100");
    }

    /**
     * Returns Blosum 30 matrix by Henikoff & Henikoff
     * @return Blosum 30 matrix
     */
    public static SubstitutionMatrix<AminoAcidCompound> getBlosum30() {
        return getAminoAcidCompoundSubstitutionMatrix("blosum30");
    }

    /**
     * Returns Blosum 35 matrix by Henikoff & Henikoff
     * @return Blosum 35 matrix
     */
    public static SubstitutionMatrix<AminoAcidCompound> getBlosum35() {
        return getAminoAcidCompoundSubstitutionMatrix("blosum35");
    }

    /**
     * Returns Blosum 40 matrix by Henikoff & Henikoff
     * @return Blosum 40 matrix
     */
    public static SubstitutionMatrix<AminoAcidCompound> getBlosum40() {
        return getAminoAcidCompoundSubstitutionMatrix("blosum40");
    }

    /**
     * Returns Blosum 45 matrix by Henikoff & Henikoff
     * @return Blosum 45 matrix
     */
    public static SubstitutionMatrix<AminoAcidCompound> getBlosum45() {
        return getAminoAcidCompoundSubstitutionMatrix("blosum45");
    }

    /**
     * Returns Blosum 50 matrix by Henikoff & Henikoff
     * @return Blosum 50 matrix
     */
    public static SubstitutionMatrix<AminoAcidCompound> getBlosum50() {
        return getAminoAcidCompoundSubstitutionMatrix("blosum50");
    }

    /**
     * Returns Blosum 55 matrix by Henikoff & Henikoff
     * @return Blosum 55 matrix
     */
    public static SubstitutionMatrix<AminoAcidCompound> getBlosum55() {
        return getAminoAcidCompoundSubstitutionMatrix("blosum55");
    }

    /**
     * Returns Blosum 60 matrix by Henikoff & Henikoff
     * @return Blosum 60 matrix
     */
    public static SubstitutionMatrix<AminoAcidCompound> getBlosum60() {
        return getAminoAcidCompoundSubstitutionMatrix("blosum60");
    }

    /**
     * Returns Blosum 62 matrix by Henikoff & Henikoff
     * @return Blosum 62 matrix
     */
    public static SubstitutionMatrix<AminoAcidCompound> getBlosum62() {
        return getAminoAcidCompoundSubstitutionMatrix("blosum62");
    }

    /**
     * Returns Blosum 65 matrix by Henikoff & Henikoff
     * @return Blosum 65 matrix
     */
    public static SubstitutionMatrix<AminoAcidCompound> getBlosum65() {
        return getAminoAcidCompoundSubstitutionMatrix("blosum65");
    }

    /**
     * Returns Blosum 70 matrix by Henikoff & Henikoff
     * @return Blosum 70 matrix
     */
    public static SubstitutionMatrix<AminoAcidCompound> getBlosum70() {
        return getAminoAcidCompoundSubstitutionMatrix("blosum70");
    }

    /**
     * Returns Blosum 75 matrix by Henikoff & Henikoff
     * @return Blosum 75 matrix
     */
    public static SubstitutionMatrix<AminoAcidCompound> getBlosum75() {
        return getAminoAcidCompoundSubstitutionMatrix("blosum75");
    }

    /**
     * Returns Blosum 80 matrix by Henikoff & Henikoff
     * @return Blosum 80 matrix
     */
    public static SubstitutionMatrix<AminoAcidCompound> getBlosum80() {
        return getAminoAcidCompoundSubstitutionMatrix("blosum80");
    }

    /**
     * Returns Blosum 85 matrix by Henikoff & Henikoff
     * @return Blosum 85 matrix
     */
    public static SubstitutionMatrix<AminoAcidCompound> getBlosum85() {
        return getAminoAcidCompoundSubstitutionMatrix("blosum85");
    }

    /**
     * Returns Blosum 90 matrix by Henikoff & Henikoff
     * @return Blosum 90 matrix
     */
    public static SubstitutionMatrix<AminoAcidCompound> getBlosum90() {
        return getAminoAcidCompoundSubstitutionMatrix("blosum90");
    }

    /**
     * Returns PAM 250 matrix by Gonnet, Cohen & Benner
     * @return Gonnet 250 matrix
     */
    public static SubstitutionMatrix<AminoAcidCompound> getGonnet250() {
        return getAminoAcidCompoundSubstitutionMatrix("gonnet250");
    }

    /**
     * Returns Nuc 4.2 matrix by Lowe
     * @return Nuc 4.2 matrix
     */
    public static SubstitutionMatrix<NucleotideCompound> getNuc4_2() {
        return getNucleotideCompoundSubstitutionMatrix("nuc-4_2");
    }

    /**
     * Returns Nuc 4.4 matrix by Lowe
     * @return Nuc 4.4 matrix
     */
    public static SubstitutionMatrix<NucleotideCompound> getNuc4_4() {
        return getNucleotideCompoundSubstitutionMatrix("nuc-4_4");
    }

    /**
     * Returns PAM 250 matrix by Dayhoff
     * @return PAM 250 matrix
     */
    public static SubstitutionMatrix<AminoAcidCompound> getPAM250() {
        return getAminoAcidCompoundSubstitutionMatrix("pam250");
    }

    // helper methods

    // reads in an amino acid substitution matrix, if necessary
    private static SubstitutionMatrix<AminoAcidCompound> getAminoAcidCompoundSubstitutionMatrix(String file) {
        if (!aminoAcidMatrices.containsKey(file)) {
            aminoAcidMatrices.put(file, new SimpleSubstitutionMatrix<AminoAcidCompound>(
                    AminoAcidCompoundSet.getAminoAcidCompoundSet(), getReader(file), file));
        }
        return aminoAcidMatrices.get(file);
    }

    // reads in a nucleotide substitution matrix, if necessary
    private static SubstitutionMatrix<NucleotideCompound> getNucleotideCompoundSubstitutionMatrix(String file) {
        if (!nucleotideMatrices.containsKey(file)) {
            nucleotideMatrices.put(file, new SimpleSubstitutionMatrix<NucleotideCompound>(
                    AmbiguityDNACompoundSet.getDNACompoundSet(), getReader(file), file));
        }
        return nucleotideMatrices.get(file);
    }

    // reads in a substitution matrix from a resource file
    private static InputStreamReader getReader(String file) {
        return new InputStreamReader(SubstitutionMatrixHelper.class.getResourceAsStream(String.format("/%s.txt",
                file)));
    }

}
