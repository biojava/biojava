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

package org.biojava.bio.symbol;

/** A genetic code translation table representing a translation table in the 
 * DDBJ/EMBL/GenBank Feature Table (appendix V).
 *
 * @author gwaldon
 * @since 1.5
 */
//PENDING a GeneticCodeTable interface that includes start codons.
public class SimpleGeneticCodeTable extends SimpleManyToOneTranslationTable {
    
    private int table_num;
    private String description;
    
    /** Creates a new instance of SimpleGeneticCodeTable */
    public SimpleGeneticCodeTable(FiniteAlphabet source, FiniteAlphabet target) {
        super(source,target);
    }
    
    public void setTableNumber(int num) {
        table_num = num;
    }
   
    /**
     * @return the value for the feature qualifier table_num  
     * found in the DDBJ/EMBL/GenBank Feature Table. The associated
     * feature key is CDS.
     */
    public int getTableNumber() {
        return table_num;
    }
               
    public void setDescription(String description) {
        this.description = description;
    }
    
    /**
     * @return A string descripting this table, normally the one found in
     * the DDBJ/EMBL/GenBank Feature Table.
     */
    public String getDescription() {
        if(description==null)
            return("");
        return description;
    }
}
