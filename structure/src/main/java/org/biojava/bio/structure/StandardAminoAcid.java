/*
 *                  BioJava development code
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
 * Created on Dec 21, 2005
 *
 */
package org.biojava.bio.structure;

import java.io.InputStream;
import java.util.HashMap;
import java.util.Map;


import org.biojava.bio.structure.io.PDBFileParser;


/** A class that provides a set of standard amino acids.
 * 
 * 
 * 
 * @author Andreas Prlic
 * @author Tamas Horvath provided the standard amino acids.
 *   
 *
 */
public final class StandardAminoAcid {
    
    static private Map<String,AminoAcid> aminoAcids;
    
    /** can not be instanciated
     * 
     */
    private StandardAminoAcid() {
        super();
        
    }
    
    /**
     * <p>
     * Initialize the static StandardAminoAcid resource.
     * </p>
     *
     * <p>
     * This parses the resource
     * <code>org/biojava/bio/structure/standardaminos.pdb</code>
     * and builds a basic set of amino acids.
     *</p>
     * @author Tamas Horvath provided the standard amino acids 
     */
    static {
        aminoAcids = new HashMap<String,AminoAcid>();
    
        try {
            InputStream fileStream = StandardAminoAcid.class.getClassLoader().getResourceAsStream(
                    "org/biojava/bio/structure/standardaminos.pdb"
            );
            if (fileStream == null) {
                throw new Exception("Couldn't locate standardaminos.pdb.  This probably means that your biojava.jar file is corrupt or incorrectly built.");
            }
            
            PDBFileParser parser = new PDBFileParser();
            Structure s = parser.parsePDBFile(fileStream);

       
            GroupIterator iter = new GroupIterator(s);
            while (iter.hasNext()){
                Group g = (Group) iter.next();
               
                if ( g instanceof AminoAcid){
                    AminoAcid aa = (AminoAcid)g;
       
                    aminoAcids.put(aa.getPDBName(),aa);
                    aminoAcids.put(aa.getAminoType().toString(),aa);
                    
                }
            }
            
        } catch (Exception t) {
            throw new RuntimeException( "Unable to initialize standard aminoacids", t);
        }
    }
    
    /** get a standard amino acid.
     * 
     * @param name the 3- or 1-letter representation of the amino acid.
     * @return the amino acids, or null if the name can not be matched
     */
    public static AminoAcid getAminoAcid(String name){
        
        return (AminoAcid) aminoAcids.get(name);
    }
    
}
