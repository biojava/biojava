/*
 *                    PDB web development code
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
 *
 * Created on Jul 25, 2009
 * Created by Andreas Prlic
 *
 */

package org.biojava.bio.structure.io.mmcif.model;

public class AuditAuthor
{
   String name;
   String pdbx_ordinal;
   public String getName()
   {
      return name;
   }
   public void setName(String name)
   {
      this.name = name;
   }
   public String getPdbx_ordinal()
   {
      return pdbx_ordinal;
   }
   public void setPdbx_ordinal(String pdbx_ordinal)
   {
      this.pdbx_ordinal = pdbx_ordinal;
   }


}
