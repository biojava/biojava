

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
 * Created on 26.04.2004
 * @author Andreas Prlic
 *
 */
package org.biojava.bio.structure.io;

import java.io.IOException;

import org.biojava.bio.structure.Structure;


/**
 * Defines the interface how to access Structure (- PDB file) input,output readers, writers.
 *
 * @author Andreas Prlic
 * @version %I% %G%
 * @deprecated Do no use any longer, StructureIO is now a static utility calss on the org.biojava3.structure level..
 */

public interface StructureIO {
    
    /** Get a structure by providing a PDB code.
     * 
     * @param pdbId  a String specifying the id value (PDB code)
     * @return a Structure object
     * @throws IOException ...
     */

    public Structure getStructureById(String pdbId) throws IOException;


    
}
  
