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

package org.biojavax.bio.taxa.io;

import java.io.BufferedReader;
import java.io.IOException;

import org.biojava.bio.seq.io.ParseException;
import org.biojavax.bio.taxa.NCBITaxon;

/**
 * Implementors are able to load taxonomy files and generate sets of NCBITaxon objects
 * that represent them. Taxon objects should be generated using RichObjectFactory.
 * @author Richard Holland
 * @since 1.5
 */
public interface NCBITaxonomyLoader {
    
    /**
     * Reads the next entry from the nodes.dmp file and returns the corresponding 
     * NCBITaxon object.
     * @param nodes something that reads the nodes.dmp file
     * @return the next NCBITaxon object in the file, or null if the file has ended.
     */
    public NCBITaxon readNode(BufferedReader nodes) throws IOException, ParseException;
    
    /**
     * Reads the next entry from the names.dmp file and returns the corresponding 
     * NCBITaxon object with the name added in already. Note that this does not clear
     * out existing names from the taxon, it only adds them. Use the code snippet
     * below if you want to clear out the names first:
     * <code>
     * for (Iterator i = taxon.getNameClasses().iterator(); i.hasNext(); ) {
     *      taxon.getNames((String)i.next()).clear();
     * }
     * </code>
     * @param names something that reads the names.dmp file
     * @return the NCBITaxon object corresponding to the next entry in the file, 
     * or null if the file has ended.
     */
    public NCBITaxon readName(BufferedReader names) throws IOException, ParseException;
}
