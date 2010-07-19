/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.biojava3.core.sequence;

import java.util.Comparator;
import org.biojava3.core.sequence.template.AbstractSequence;


/**
 * 
 * @author Scooter Willis <willishf at gmail dot com>
 */
    public class SequenceComparator implements Comparator<AbstractSequence>{

        @Override
        public int compare(AbstractSequence o1, AbstractSequence o2) {
          
            return o1.getBioBegin() - o2.getBioBegin();
        }

    }
