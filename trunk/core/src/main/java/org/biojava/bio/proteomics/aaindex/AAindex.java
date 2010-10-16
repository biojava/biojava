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

/*
 *    AAindex.java
 */
package org.biojava.bio.proteomics.aaindex;

import java.util.HashMap;
import java.util.Map;

import org.biojava.bio.seq.ProteinTools;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.SimpleSymbolPropertyTable;

/**
 * Symbol property table based on the Amino Acid Index Database. Each 
 * <code>AAindex</code> object represents a single entry of an AAindex1 file.
 * Each entry contains twenty numeric values for the twenty amino acids, e.g.
 * describing the hydrophobicity of an amino acid. To get this value for a
 * certain amino acid call the 
 * {@link org.biojava.bio.symbol.SymbolPropertyTable#getDoubleValue(org.biojava.bio.symbol.Symbol)}
 * method with the appropriate symbol, e.g. 
 * <code>aaindex.getDoubleValue(ProteinTools.gln())</code>.
 * The remaining data fields, i.e. object properties, are fully described in the
 * <a href="http://www.genome.ad.jp/dbget-bin/show_man?aaindex">AAindex manual
 * </a>.
 * <p><b>References:</b></p>
 * <p><a href="http://www.genome.ad.jp/dbget/aaindex.html">AAindex web 
 * site</a>.</p>
 * <p>Kawashima, S. and Kanehisa, M.; AAindex: amino acid index database. 
 * Nucleic Acids Res. 28, 374 (2000).</p>
 * <p>Tomii, K. and Kanehisa, M.;  Analysis of amino acid indices and mutation 
 * matrices for sequence comparison and structure prediction of proteins. 
 * Protein Eng. 9, 27-36 (1996).</p>
 * <p>Nakai, K., Kidera, A., and Kanehisa, M.;  Cluster analysis of amino acid 
 * indices for prediction of protein structure and function.  
 * Protein Eng. 2, 93-100 (1988)</p>
 * @author <a href="mailto:Martin.Szugat@GMX.net">Martin Szugat</a>
 * @version $Revision$
 */
public class AAindex extends SimpleSymbolPropertyTable {

    /* PUBLIC CONSTANTS */
    
    /**
     * The alphabet of the symbol property table, that is 
     * {@linkplain ProteinTools#getAlphabet() protein}.
     */
    public static final Alphabet PROTEIN_ALPHABET = ProteinTools.getAlphabet(); 

    /* PRIVATE FIELDS */

    /**
     * The accession number of the AAindex entry.
     */
    private String accessionNumber = null;

    /**
     * The description of the AAindex entry.
     */
    private String description = null;

    /**
     * Literature database entry numbers for the AAindex entry.
     */
    private String[] litdbEntryNumbers = null;

    /**
     * The authors of the article first explaining this AAindex entry.
     */
    private String articleAuthors = null;

    /**
     * The title of the article.
     */
    private String articleTitle = null;

    /**
     * The reference to the journal which published the article.
     */
    private String journalReference = null;

    /**
     * User commments.
     */
    private String comment = null;

    /**
     * A map of similar AAindex entries with correlation coefficients.
     */
    private Map similarEntries = new HashMap();
    
    /* PUBLIC CONSTRUCTOR */

    /**
     * Initializes the AAindex symbol property table.
     * @param accessionNumber the AAindex accession number (same as the table
     * name)
     * @throws NullPointerException if <code>accessionNumber</code> is
     * <code>null</code>.
     */
    public AAindex(String accessionNumber) throws NullPointerException {
        super(PROTEIN_ALPHABET, accessionNumber);
        if (accessionNumber == null) {
            throw new NullPointerException("accessionNumber is null.");
        }
        this.accessionNumber = accessionNumber;
    }
    
    /* PUBLIC PROPERTIES */

    /**
     * Gets the accession number of the AAindex entry.
     * @return the accession number (same as 
     * {@link org.biojava.bio.symbol.SymbolPropertyTable#getName()}
     */
    public String accessionNumber() {
        return accessionNumber;
    }

    /**
     * Gets the names of the authors which first published an article about the
     * AAindex entry.
     * @return a list of names. May be <code>null</code>.
     */
    public String getArticleAuthors() {
        return articleAuthors;
    }

    /**
     * Sets the names of the authors which first published an article about the
     * AAindex entry.
     * @param articleAuthors May be <code>null</code>.
     */
    public void setArticleAuthors(String articleAuthors) {
        this.articleAuthors = articleAuthors;
    }

    /**
     * Gets the user comment for the AAindex entry.
     * @return free text. May be <code>null</code>.
     */
    public String getComment() {
        return comment;
    }

    /**
     * Sets the user comment for the AAindex entry.
     * @param comment free text. May be <code>null</code>.
     */
    public void setComment(String comment) {
        this.comment = comment;
    }

    /**
     * Gets the title of the article which describes the AAindex entry.
     * @return the article title. May be <code>null</code>.
     */
    public String getArticleTitle() {
        return articleTitle;
    }

    /**
     * Sets the title of the article which describes the AAindex entry.
     * @param articleTitle the article title. May be <code>null</code>.
     */
    public void setArticleTitle(String articleTitle) {
        this.articleTitle = articleTitle;
    }

    /**
     * Gets the description for the AAindex entry.
     * @return a human readable description. May be <code>null</code>.
     */
    public String getDescription() {
        return description;
    }

    /**
     * Sets the description for the AAindex entry.
     * @param description a human readable description. 
     * May be <code>null</code>.
     */
    public void setDescription(String description) {
        this.description = description;
    }

    /**
     * Gets a reference to the journal which published the article about the
     * AAindex entry.
     * @return the journal reference. May be <code>null</code>.
     */
    public String getJournalReference() {
        return journalReference;
    }

    /**
     * Sets a reference to the journal which published the article about the
     * AAindex entry.
     * @param journalReference the journal reference. May be <code>null</code>.
     */
    public void setJournalReference(String journalReference) {
        this.journalReference = journalReference;
    }

    /**
     * Gets the list of literature database identifiers for the AAindex entry.
     * @return a list of identifiers. May be <code>null</code>.
     */
    public String[] getLITDBEntryNumbers() {
        return (litdbEntryNumbers == null ? null 
                : (String[]) litdbEntryNumbers.clone());
    }

    /**
     * Sets the list of literature database identifiers for the AAindex entry.
     * @param litdbEntryNumbers a list of identifiers
     */
    public void setLITDBEntryNumbers(String[] litdbEntryNumbers) {
        if (litdbEntryNumbers == null) {
            this.litdbEntryNumbers = null;
        } else {
            this.litdbEntryNumbers = (String[]) litdbEntryNumbers.clone();
        }
    }
    
    /**
     * Returns a map with the names of similar AAindex entries and its 
     * correlation coefficients. 
     * @return a map which keys are the names of the similar AAindex entries and
     * which values are the corresponding correlation coefficients
     */
    public Map similarEntries() {
        return similarEntries;
    }

}
