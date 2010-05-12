/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.biojava3.core.sequence;

import java.util.Comparator;


/**
 * 
 * @author Scooter Willis <willishf at gmail dot com>
 */
    public class CDSComparator implements Comparator<CDSSequence>{

        @Override
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
