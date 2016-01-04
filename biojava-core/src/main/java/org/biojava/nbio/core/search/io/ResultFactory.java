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
package org.biojava.nbio.core.search.io;

import java.io.File;
import java.io.IOException;
import java.text.ParseException;
import java.util.List;
import org.biojava.nbio.core.sequence.template.Sequence;

/**
 * Designed by Paolo Pavan.
 * You may want to find my contacts on Github and LinkedIn for code info 
 * or discuss major changes.
 * https://github.com/paolopavan
 * 
 * @author Paolo Pavan
 */

public interface ResultFactory {
    /**
     * returns a list of file extensions associated to this ResultFactory
     * 
     * @return 
     */
    List<String> getFileExtensions();
    void setFile(File f);
    /**
     * Launch the parsing and get back a list of Result objects representing the
     * search result in the specified file.
     * 
     * @param maxEScore
     * @return
     * @throws Exception 
     */
    List<Result> createObjects(double maxEScore) throws IOException, ParseException;
    /**
     * The factory that implements this method will be able to save the Search results
     * to a file in the same format that it is able to read.
     * 
     * @param results
     * @throws Exception 
     */
    void storeObjects(List<Result> results) throws IOException, ParseException;
    
    /**
     * Specify the collection of sequences objects used as queries in the Search run. 
     * They will be associated back to the query during the construction of the Result object.
     * @param sequences 
     */
    void setQueryReferences(List<Sequence> sequences);
    /**
     * Specify the collection of sequences objects used as database in the Search run. 
     * They will be associated back to the Hit during the construction of the Hit object.
     * @param sequences 
     */
    void setDatabaseReferences(List<Sequence> sequences);
}
