/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.biojava3.core.sequence;

import java.util.Comparator;
import org.biojava3.core.sequence.TranscriptSequence.Sense;

/**
 *
 * @author Scooter Willis <willishf at gmail dot com>
 */
    public class ExonComparator implements Comparator<ExonSequence>{

        @Override
        public int compare(ExonSequence o1, ExonSequence o2) {
            if(o1.sense != o2.sense){
                return o1.getBegin() - o2.getBegin();
            }
            if(o1.sense == Sense.NEGATIVE){
                return -1 * (o1.getBegin() - o2.getBegin());
            }

            return o1.getBegin() - o2.getBegin();
        }

    }
