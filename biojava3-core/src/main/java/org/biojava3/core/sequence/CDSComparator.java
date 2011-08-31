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
 * Created on 3/1/2010
 * @author Scooter Willis <willishf at gmail dot com>
 */

package org.biojava3.core.sequence;

import java.util.Comparator;



    public class CDSComparator implements Comparator<CDSSequence>{

/**
 * Used to sort two CDSSequences where Negative Strand makes it tough
 * @param o1
 * @param o2
 * @return val
 */
        public int compare(CDSSequence o1, CDSSequence o2) {
            if(o1.getStrand() != o2.getStrand()){
                return o1.getBioBegin() - o2.getBioBegin();
            }
            if(o1.getStrand() == Strand.NEGATIVE){
                return -1 * (o1.getBioBegin() - o2.getBioBegin());
            }

            return o1.getBioBegin() - o2.getBioBegin();
        }

    }
