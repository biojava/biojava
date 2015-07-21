/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.biojava.nbio.core.search.io;

import java.io.File;
import java.util.List;

/**
 *
 * @author pavanpa
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
    List<Result> createObjects(double maxEScore) throws Exception;
    /**
     * The factory that implements this method will be able to save the Search results
     * to a file in the same format that it is able to read.
     * 
     * @param results
     * @throws Exception 
     */
    void storeObjects(List<Result> results) throws Exception;
    
    /**
     * Specify the collection of sequences objects used as queries in the Search run. 
     * They will be associated back to the query during the construction of the Result object.
     * @param sequences 
     */
    void setQueryReferences(List sequences);
    /**
     * Specify the collection of sequences objects used as database in the Search run. 
     * They will be associated back to the Hit during the construction of the Hit object.
     * @param sequences 
     */
    void setDatabaseReferences(List sequences);
}
