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

package org.biojava3.protmod.structure;

import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.io.FileParsingParameters;

public class TmpAtomCache
{
   static String tmpDir = System.getProperty("java.io.tmpdir");
   public static AtomCache cache = new AtomCache(tmpDir, true);
   static {
	   FileParsingParameters params = new FileParsingParameters();
	   params.setLoadChemCompInfo(true);
	   params.setAlignSeqRes(true);
	   params.setParseSecStruc(false);
	   params.setUpdateRemediatedFiles(true);
	   cache.setFileParsingParams(params);
	   cache.setAutoFetch(true);
   }
}
