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
 * Created on Jun 17, 2010
 * Author: ap3 
 *
 */

package org.biojava.bio.structure;

import java.io.Serializable;

/** Everything that is needed to uniquely describe a residue position
 * 
 * @author Andreas Prlic
 *
 */
public class PDBResidueNumber implements Serializable
{
   
   /**
    * 
    */
   private static final long serialVersionUID = 1773011704758536083L;
   String chainId;
   String insCode;
   Integer residueNumber;
   
   public String getChainId()
   {
      return chainId;
   }
   public void setChainId(String chainId)
   {
      this.chainId = chainId;
   }
   public String getInsCode()
   {
      return insCode;
   }
   public void setInsCode(String insCode)
   {
      this.insCode = insCode;
   }
   public Integer getResidueNumber()
   {
      return residueNumber;
   }
   public void setResidueNumber(Integer residueNumber)
   {
      this.residueNumber = residueNumber;
   }

}
