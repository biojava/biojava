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
 * Created on Feb 9, 2006
 *
 */
package org.biojava.dasobert.das2;

import org.biojava.dasobert.dasregistry.DasSource;

public interface Das2Source
extends DasSource {

   public Das2Capability[] getDas2Capabilities();
   public void setDas2Capabilities(Das2Capability[] capabilities);
   
   /** test if this is a DAS1 source represented as a DAS2 source
    *  if true - this source can be converted into a DAS1 source by using
    *  DasSourceConverter.toDas1(Das2Source);
    *
    * @return true if the DasSource has DAS1 capabilties
    */
   public boolean hasDas1Capabilities();
}
