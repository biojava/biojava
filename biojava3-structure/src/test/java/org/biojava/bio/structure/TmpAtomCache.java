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
 * Created on May 11, 2010
 * Author: Andreas Prlic 
 *
 */

package org.biojava.bio.structure;

import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.io.FileParsingParameters;

/**
 * @deprecated The zero-argument constructor to AtomCache mostly replaces this.
 *  The singleton pattern is not really needed--just make a new AtomCache
 *  instance and set the updateRemediatedFiles parameter if needed.
 *  (deprecated by Spencer Bliven)
 */
@Deprecated
public class TmpAtomCache
{
   
   public static AtomCache cache = new AtomCache();
   static {
      FileParsingParameters params = new FileParsingParameters();
      params.setLoadChemCompInfo(false); //default value
      params.setUpdateRemediatedFiles(true); //non-default value
      cache.setFileParsingParams(params);
      cache.setAutoFetch(true);
   }
}
